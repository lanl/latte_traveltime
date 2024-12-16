#
#(c) 2024. Triad National Security, LLC. All rights reserved.
#
# This program was produced under U.S. Government contract 89233218CNA000001 
# for Los Alamos National Laboratory (LANL), which is operated by 
# Triad National Security, LLC for the U.S. Department of Energy/National Nuclear 
# Security Administration. All rights in the program are reserved by 
# Triad National Security, LLC, and the U.S. Department of Energy/National 
# Nuclear Security Administration. The Government is granted for itself and 
# others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide 
# license in this material to reproduce, prepare. derivative works, 
# distribute copies to the public, perform publicly and display publicly, 
# and to permit others to do so.
#
# Author:
#    Kai Gao, kaigao@lanl.gov
#

import argparse
import torch
import lightning as pl
from utility import *
from model2 import *

#==============================================================================
parser = argparse.ArgumentParser(description='MTL inference-refinement options')
parser.add_argument('--nodes', '-nodes', type=int, default=1, help='number of nodes')
parser.add_argument('--gpus_per_node', '-gpus_per_node', type=int, default=1, help='number of gpus per node')
parser.add_argument('--n1', '-n1', type=int, default=256, help='number of sampling points in x1')
parser.add_argument('--n2', '-n2', type=int, default=256, help='number of sampling points in x2')
parser.add_argument('--input', '-in', type=str, default=None, help='test model using test checkpoint')
parser.add_argument('--model_infer', '-model_infer', type=str, default=None, help='inference model')
parser.add_argument('--model_refine', '-model_refine', type=str, default=None, help='refinement model')
parser.add_argument('--output', '-out', type=str, default=None, help='output name')
parser.add_argument('--niter', '-niter', type=int, default=3, help='number of refinement cycles')
opts = parser.parse_args()

assert opts.n1 >= 1
assert opts.n2 >= 1

if torch.cuda.is_available() and opts.gpus_per_node >= 1:
    device = torch.device('cuda')
    print(date_time(), ' >> Using GPU')
else:
    device = torch.device('cpu')
    print(date_time(), ' >> Using CPU')

torch.set_float32_matmul_precision('high')

#==============================================================================
class mtlnet_infer(pl.LightningModule):
    def __init__(self, lr: float = 1.0e-4):

        super(mtlnet_infer, self).__init__()
        
        self.lr = lr
        self.in_channels = 1
        
        # encoder
        self.l1 = 16
        self.l2 = 32
        self.l3 = 64
       
        self.encoder1_meq = resu1(self.in_channels, self.l1)
        self.encoder2_meq = resu2(self.l1, self.l2)
        self.encoder3 = resu3(self.l2, self.l3)
        
        self.decoder = mtl_decoder(self.l1, self.l2, self.l3, out_ch=self.l1, last_kernel_size=3)
        self.subdecoder_fault_semantic = mtl_subdecoder(in_ch=self.l1,
                                                        out_ch=1,
                                                        bn=False,
                                                        mid_activation='relu',
                                                        activation='sigmoid')
        self.subdecoder_fault_dip = mtl_subdecoder(in_ch=self.l1,
                                                   out_ch=1,
                                                   bn=True,
                                                   mid_activation='relu',
                                                   activation='sigmoid')

    def forward(self, x):

        out_encoder1_meq = self.encoder1_meq(x['meq'])
        out_encoder2_meq = self.encoder2_meq(maxpool(out_encoder1_meq, 2))
        out_encoder3 = self.encoder3(maxpool(out_encoder2_meq, 2))
        
        # decoders
        out = {}
        out_fault = self.decoder(x['meq'], out_encoder1_meq, out_encoder2_meq, out_encoder3)
        out_fault_semantic = self.subdecoder_fault_semantic(out_fault)
        out_fault_dip = self.subdecoder_fault_dip(out_fault) * out_fault_semantic
        out['fsem'] = out_fault_semantic
        out['fdip'] = out_fault_dip

        return out

#==============================================================================
class mtlnet_refine(pl.LightningModule):
    def __init__(self, lr: float = 1.0e-4):

        super(mtlnet_refine, self).__init__()
        
        self.lr = lr
        self.in_ch = 1
        
        # encoder
        self.l1 = 16
        self.l2 = 32
        self.l3 = 64
       
        self.input1 = nn.Sequential(conv(self.in_ch, 4), conv(4, 4))
        self.input2 = nn.Sequential(conv(self.in_ch, 4), conv(4, 4))
        self.input3 = nn.Sequential(conv(self.in_ch, 4), conv(4, 4))
        
        self.encoder1 = resu1(4 * 3, self.l1)
        self.encoder2 = resu2(self.l1, self.l2)
        self.encoder3 = resu3(self.l2, self.l3)
        
        self.decoder = mtl_decoder(self.l1, self.l2, self.l3, out_ch=self.l1, last_kernel_size=3)
        self.subdecoder_fault_semantic = mtl_subdecoder(in_ch=self.l1,
                                                        out_ch=1,
                                                        bn=False,
                                                        mid_activation='relu',
                                                        activation='sigmoid')
        self.subdecoder_fault_dip = mtl_subdecoder(in_ch=self.l1,
                                                   out_ch=1,
                                                   bn=True,
                                                   mid_activation='relu',
                                                   activation='sigmoid')

    def forward(self, x):

        # encoders
        i1 = self.input1(x['meq'])
        i2 = self.input2(x['fsem'])
        i3 = self.input3(x['fdip'])

        out_encoder1 = self.encoder1(torch.cat((i1, i2, i3), dim=1))
        out_encoder2 = self.encoder2(maxpool(out_encoder1, 2))
        out_encoder3 = self.encoder3(maxpool(out_encoder2, 2))
        
        # decoders
        out = {}
        out_fault = self.decoder(x['meq'], out_encoder1, out_encoder2, out_encoder3)
        out_fault_semantic = self.subdecoder_fault_semantic(out_fault)
        out_fault_dip = self.subdecoder_fault_dip(out_fault) * out_fault_semantic
        out['fsem'] = out_fault_semantic
        out['fdip'] = out_fault_dip

        return out
    
#==============================================================================
if __name__ == '__main__':

    # Inference phase
    n1 = opts.n1
    n2 = opts.n2

    img = {}
    img['meq'] = read_array(opts.input, (1, 1, n1, n2))

    net = mtlnet_infer()
    net.load_state_dict(torch.load(opts.model_infer, map_location=device, weights_only=True)['state_dict'])
    net.to(device)

    print(date_time(), " >> Pretrained model loaded")

    with torch.no_grad():
        img['meq'] = img['meq'].to(device)
        predict = net(img)
        
    print(date_time(), ' >> Inference finished')

    # Refinement phase
    img['fsem'] = predict['fsem']
    img['fdip'] = predict['fdip']
    
    # Load trained model
    net = mtlnet_refine()
    net.load_state_dict(torch.load(opts.model_refine, map_location=device, weights_only=True)['state_dict'])
    net.to(device)

    print(date_time(), " >> Pretrained model loaded")

    for iter in range(opts.niter):
        with torch.no_grad():
            predict = net(img)
                            
            # Next iteration                
            print(date_time(), ' >> Iteration ' + str(iter + 1) + ' finished')
            img['fsem'] = predict['fsem']
            img['fdip'] = predict['fdip']
                
    # Save results
    write_array(get_numpy(predict['fsem']), opts.output + '.fsem')
    write_array(get_numpy(predict['fdip']), opts.output + '.fdip')

    print(date_time(), ' >> Refinement finished')