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


import math
import numpy as np

def regtick(start, end, interval, mtick):
    
    tick = np.arange(start, end + 0.1 * abs(interval), interval)
    minor_tick_interval = interval / (mtick + 1.0)
    minor_tick = np.arange(start, end + 0.1 * np.abs(minor_tick_interval), minor_tick_interval)
    
    return tick, minor_tick

def regspace(start, stop, step=1):
    dir = 1 if (step > 0) else -1
    return np.arange(start, stop + dir, step)

def next_power_of_2(x):
    
    if x == 0:
        return 1
    else:
        return 2**math.ceil(math.log2(x))
    
def rescale(w, range=(0, 1)):
    
    if np.max(w) != np.min(w):
        
        dr = range[1] - range[0]
        w = w - np.min(w)
        w = w/np.max(w)
        w = w*dr + range[0]
    
    else:
        w = np.ones_like(w)*range[0]
    
    return w
    
## forward integer indices
def forward_range(start, n, step=1):
    
    r = np.zeros(n, dtype=np.int32)
    for k in range(n):
        r[k] = start + k*step
    
    return r

## backward integer indices
def backward_range(end, n, step=1):
    
    r = np.zeros(n, dtype=np.int32)
    for k in range(n):
        r[k] = end - k*step

    return r        

## strict start-end indices
def strict_range(start, end, step=1):
    
    if end < start:
        s = -np.abs(step)
    else:
        s = np.abs(step)
    r = np.zeros(np.int32(np.floor((end - start)/s)) + 1, dtype=np.int32)
    for k in range(np.size(r)):
        r[k] = start + k*s
        
    return r
    
## date and time
from datetime import datetime
def date_time():

	now = datetime.now()
	return now.strftime(" %Y/%m/%d %H:%M:%S ")
	
## convert number to string
def num2str(x, f=''):
    
    if type(x) is int and f == '':
        ff = "{:d}"
    elif type(x) is float and f == 'd':
        ff = "{:.0f}"
    else:
        ff = "{:" + f + "}"
    
    return ff.format(x)
