#
# © 2024. Triad National Security, LLC. All rights reserved.
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

import numpy as np
import os
import torch
from datetime import datetime


## Convert string to bool
def str2bool(v):

    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'y', 'true', 't', 'on', '1'):
        return True
    elif v.lower() in ('no', 'n', 'false', 'f', 'off', '0'):
        return False
    else:
        print(' Error: Argument must be one of yes/no, y/n, true/false, t/f, on/off, 1/0. ')
        exit()


## Read a raw binary array and convert to PyTorch tensor
def read_array(filename, shape, dtype=np.float32, totorch=True):

    x = np.fromfile(filename, count=np.prod(shape), dtype=dtype)
    x = np.reshape(x, shape[::-1])
    x = np.transpose(x)

    if totorch:
        x = torch.from_numpy(x).type(torch.FloatTensor)

    return x


## Write a binary array
def write_array(x, filename, dtype=np.float32):

    x = np.asarray(x, dtype=dtype)
    x = np.transpose(x)
    x.tofile(filename)


## Forward integer indices
def forward_range(start, n, step=1):

    r = np.zeros(n, dtype=np.int32)
    for k in range(n):
        r[k] = start + k * step

    return r


##  Set random seeds for training
def set_random_seed(seed=12345):
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    return None


## Get numpy array and transfer to CPU
def get_numpy(w):

    return w.squeeze().data.cpu().numpy()


## Get time stamp
def date_time():

    now = datetime.now()
    return now.strftime(" %Y/%m/%d %H:%M:%S ")
