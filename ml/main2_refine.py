#
# (c) 2024. Triad National Security, LLC. All rights reserved.
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

import os
import sys
import warnings
import argparse
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset

import lightning as pl
from lightning.pytorch.callbacks import ModelCheckpoint
from lightning.pytorch.loggers import TensorBoardLogger
from lightning.pytorch.callbacks import LearningRateMonitor
from torchmetrics.functional.image import structural_similarity_index_measure as ssim
from torchmetrics.functional.image import multiscale_structural_similarity_index_measure as mssim
import torchmetrics.functional as mf

warnings.filterwarnings("ignore")

from utility import *
from model2 import *

#==============================================================================
parser = argparse.ArgumentParser(description='MTL-Net options')
parser.add_argument('--ntrain', type=int, default=1000, help='training size')
parser.add_argument('--nvalid', type=int, default=100, help='validation size')
parser.add_argument('--batch_train', type=int, default=1, help='training batch size')
parser.add_argument('--batch_valid', type=int, default=1, help='validation batch size')
parser.add_argument('--epochs', type=int, default=100, help='max number of epochs')
parser.add_argument('--lr', type=float, default=0.5e-4, help='learning rate')
parser.add_argument('--threads', type=int, default=8, help='number of threads for data loader')
parser.add_argument('--dir_output', type=str, default='./result', help='directory')
parser.add_argument('--dir_data_train', type=str, default='./dataset2/data_train', help='directory')
parser.add_argument('--dir_target_train', type=str, default='./dataset2/target_train', help='directory')
parser.add_argument('--dir_data_valid', type=str, default='./dataset2/data_valid', help='directory')
parser.add_argument('--dir_target_valid', type=str, default='./dataset2/target_valid', help='directory')
parser.add_argument('--resume', type=str, default=None, help='restart training from resume checkopoint')
parser.add_argument('--nodes', type=int, default=1, help='number of nodes')
parser.add_argument('--gpus_per_node', type=int, default=4, help='number of gpus per node')
parser.add_argument('--seed', type=int, default=12345, help='random seed for initialization')
parser.add_argument('--check', type=str, default=None, help='test model using test checkpoint')
parser.add_argument('--n1', '-n1', type=int, default=256, help='number of sampling points in x1')
parser.add_argument('--n2', '-n2', type=int, default=256, help='number of sampling points in x2')
parser.add_argument('--input', '-in', type=str, default=None, help='test model using test checkpoint')
parser.add_argument('--model', '-model', type=str, default=None, help='test model using test checkpoint')
parser.add_argument('--output', '-out', type=str, default=None, help='test model using test checkpoint')
parser.add_argument('--niter', '-niter', type=int, default=1, help='number of iterations of refinement')
parser.add_argument('--save_iter', '-save_iter', type=str2bool, default=None, help='save iteration results')

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
class BasicDataset(Dataset):
    def __init__(self, dir_data, dir_target, data_ids, dim=(opts.n1, opts.n2)):
        self.dir_data = dir_data
        self.dir_target = dir_target
        self.ids = data_ids
        self.n1 = dim[0]
        self.n2 = dim[1]

    def __len__(self):
        return len(self.ids)

    def __getitem__(self, i):

        idx = str(self.ids[i])
        
        data = {}
        data['meq'] = read_array(self.dir_data + '/' + idx + '_meq.bin', (1, self.n1, self.n2))
        data['fsem'] = read_array(self.dir_data + '/' + idx + '_fsem.bin', (1, self.n1, self.n2))
        data['fdip'] = read_array(self.dir_data + '/' + idx + '_fdip.bin', (1, self.n1, self.n2))

        target = {}
        target['fsem'] = read_array(self.dir_target + '/' + idx + '_fsem.bin', (1, self.n1, self.n2))
        target['fdip'] = read_array(self.dir_target + '/' + idx + '_fdip.bin', (1, self.n1, self.n2))

        return data, target


def custom_loss(y_pred, y_true):

    # # fault semantic
    mp = y_pred['fsem']
    mt = y_true['fsem']
    loss_fsem = 1.0 - (2.0 * torch.sum(mp * mt) + 1.0) / (torch.sum(mp + mt) + 1.0)

    # fault dip
    dp = y_pred['fdip']
    dt = y_true['fdip']
    loss_fdip = F.l1_loss(dp, dt) * 10
    
    loss = loss_fsem + loss_fdip
    
    return loss, loss_fsem, loss_fdip
    
def custom_accuracy(y_pred, y_true):
    
    # accuracy
    accuracy = mf.classification.binary_accuracy(y_pred, y_true)

    # ssim
    s = ssim(y_pred, y_true, data_range=1.0)
    
    # precision and recall
    precision = mf.classification.binary_precision(y_pred, y_true)
    recall = mf.classification.binary_recall(y_pred, y_true)

    return accuracy, precision, recall, s

#==============================================================================
class mtlnet(pl.LightningModule):
    def __init__(self, lr: float = 1.0e-4):

        super(mtlnet, self).__init__()
        
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

    def training_step(self, batch):

        x, y_true = batch
        y_pred = self.forward(x)

        loss, loss_fault_semantic, loss_fault_dip = custom_loss(y_pred, y_true)
        self.log("train_loss", loss, on_step=False, on_epoch=True, prog_bar=True)
        self.log("train_loss_fsem", loss_fault_semantic, on_step=False, on_epoch=True, prog_bar=False)
        self.log("train_loss_fdip", loss_fault_dip, on_step=False, on_epoch=True, prog_bar=False)

        accuracy, precision, recall, s = custom_accuracy(y_pred['fsem'], y_true['fsem'])
        self.log("train_accuracy", accuracy, on_step=False, on_epoch=True, prog_bar=False)
        self.log("train_precision", precision, on_step=False, on_epoch=True, prog_bar=False)
        self.log("train_recall", recall, on_step=False, on_epoch=True, prog_bar=False)
        self.log("train_ssim", s, on_step=False, on_epoch=True, prog_bar=False)
    
        return loss

    def validation_step(self, batch):

        x, y_true = batch
        y_pred = self.forward(x)

        loss, loss_fault_semantic, loss_fault_dip = custom_loss(y_pred, y_true)
        self.log("valid_loss", loss, on_step=False, on_epoch=True, prog_bar=True)
        self.log("valid_loss_fsem", loss_fault_semantic, on_step=False, on_epoch=True, prog_bar=False)
        self.log("valid_loss_fdip", loss_fault_dip, on_step=False, on_epoch=True, prog_bar=False)

        accuracy, precision, recall, s = custom_accuracy(y_pred['fsem'], y_true['fsem'])
        self.log("valid_accuracy", accuracy, on_step=False, on_epoch=True, prog_bar=False)
        self.log("valid_precision", precision, on_step=False, on_epoch=True, prog_bar=False)
        self.log("valid_recall", recall, on_step=False, on_epoch=True, prog_bar=False)
        self.log("valid_ssim", s, on_step=False, on_epoch=True, prog_bar=False)
        
        return loss

    def configure_optimizers(self):

        optimizer = torch.optim.Adam(self.parameters(), lr=self.lr)
        scheduler = {
            'scheduler': torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer),
            'name': 'lr',
            'monitor': 'valid_loss'
            }
        return {'optimizer': optimizer, 'lr_scheduler': scheduler}

#==============================================================================
if __name__ == '__main__':

    if opts.input is None and opts.check is None:
        ## Training phase

        logger = TensorBoardLogger(opts.dir_output)
        checkpoint_callback = ModelCheckpoint(monitor='valid_loss',
                                              dirpath=opts.dir_output,
                                              filename='{epoch:d}',
                                              mode='min',
                                              save_top_k=opts.epochs,
                                              save_last=True)

        lr_monitor = LearningRateMonitor(logging_interval='step')
        
        params = {
            'max_epochs': opts.epochs,
            'default_root_dir': opts.dir_output,
            'logger': logger,
            'callbacks': [checkpoint_callback, lr_monitor]
        }
        
        if torch.cuda.is_available() and opts.gpus_per_node >= 1:
            params['devices'] = opts.gpus_per_node
            params['num_nodes'] = opts.nodes
            params['accelerator'] = 'gpu'
            params['strategy'] = 'ddp'
            
        trainer = pl.Trainer(**params)

        t = BasicDataset(opts.dir_data_train,
                         opts.dir_target_train,
                         data_ids=np.arange(0, opts.ntrain),
                         dim=(opts.n1, opts.n2))
        train_loader = DataLoader(t, batch_size=opts.batch_train, num_workers=opts.threads, shuffle=True)

        v = BasicDataset(opts.dir_data_valid,
                         opts.dir_target_valid,
                         data_ids=np.arange(0, opts.nvalid),
                         dim=(opts.n1, opts.n2))
        valid_loader = DataLoader(v, batch_size=opts.batch_valid, num_workers=opts.threads)

        set_random_seed(opts.seed)
        params = {'lr': opts.lr}
        net = mtlnet(**params)

        if opts.resume:
            trainer.fit(net, train_loader, valid_loader, ckpt_path=opts.resume)
        else:
            trainer.fit(net, train_loader, valid_loader)

        print(date_time(), ' >> Training finished')

    if opts.input is None and opts.check is not None:
        ## Validation phase
        
        v = BasicDataset(opts.dir_data_valid,
                         opts.dir_target_valid,
                         data_ids=np.arange(0, opts.nvalid),
                         dim=(opts.n1, opts.n2))
        valid_loader = DataLoader(v, batch_size=opts.batch_valid, num_workers=opts.threads)

        net = mtlnet()
        net.load_state_dict(torch.load(opts.check, map_location=device, weights_only=True)['state_dict'])
        net.to(device)
        l = 1

        with tqdm(total=len(v), desc='', unit='image') as pbar:

            for (input, target) in valid_loader:

                with torch.no_grad():
                    input['meq'] = input['meq'].to(device)
                    input['fsem'] = input['fsem'].to(device)
                    input['fdip'] = input['fdip'].to(device)
                    predict = net(input)

                for i in range(0, input['meq'].shape[0]):

                    ir = (l - 1) * opts.batch_valid + i

                    fig, ax = plt.subplots(2, 3)
                    
                    mq = get_numpy(input['meq'][i])
                    dip = get_numpy(input['fdip'][i])
                    fsem = get_numpy(target['fsem'][i])
                    fdip = get_numpy(target['fdip'][i])
                    ax[0][0].imshow(mq, cmap='jet'), ax[0][0].set_title('meq')
                    ax[0][1].imshow(fsem, cmap='viridis'), ax[0][1].set_title('fsem')
                    ax[0][2].imshow(fdip, cmap='jet'), ax[0][2].set_title('fdip')
                    
                    fsem = get_numpy(predict['fsem'][i])
                    fdip = get_numpy(predict['fdip'][i])
                    ax[1][0].imshow(dip, cmap='jet'), ax[1][1].set_title('fdip gt')
                    ax[1][1].imshow(fsem, cmap='viridis'), ax[1][1].set_title('fsem')
                    ax[1][2].imshow(fdip, cmap='jet'), ax[1][2].set_title('fdip')

                    plt.show()

                pbar.update(input['meq'].shape[0])
                l = l + 1

        print(date_time(), ' >> Validation finished')

    if opts.input is not None:
        ## Inference phase

        # Read image
        n1 = opts.n1
        n2 = opts.n2

        img = {}
        img['meq'] = read_array(opts.input, (1, 1, n1, n2))
        img['fsem'] = read_array(opts.input + '.fsem', (1, 1, n1, n2))
        img['fdip'] = read_array(opts.input + '.fdip', (1, 1, n1, n2))

        # Load trained model
        net = mtlnet()
        net.load_state_dict(torch.load(opts.model, map_location=device, weights_only=True)['state_dict'])
        net.to(device)

        print(date_time(), " >> Pretrained model loaded")

        for iter in range(opts.niter):
            with torch.no_grad():
                img['meq'] = img['meq'].to(device)
                img['fsem'] = img['fsem'].to(device)
                img['fdip'] = img['fdip'].to(device)
                predict = net(img)
                                
                # Next iteration                
                if opts.niter > 1:
                    print(date_time(), ' >> Iteration ' + str(iter + 1) + ' finished')

                    img['fsem'] = predict['fsem']
                    img['fdip'] = predict['fdip']
                    
                    if opts.save_iter:
                        write_array(get_numpy(predict['fsem']), opts.output + '.fsem.iter' + str(iter))
                        write_array(get_numpy(predict['fdip']), opts.output + '.fdip.iter' + str(iter))

        if opts.save_iter is False:
            write_array(get_numpy(predict['fsem']), opts.output + '.fsem')
            write_array(get_numpy(predict['fdip']), opts.output + '.fdip')

        print(date_time(), ' >> Refinement finished')
