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
import torch


def count_nonempty_lines(filename):

    file = open(filename, "r")
    line_count = 0
    for line in file:
        if line != "\n":
            line_count += 1
    file.close()

    return line_count


def read_array(filename, shape, dtype=np.float32, ascii=False, totorch=False):

    if ascii:
        w = np.zeros(shape, dtype=dtype)
        f = open(filename, 'r')
        l = 1
        while l <= shape[0]:
            w[l - 1, :] = f.readline().strip().split()
            l = l + 1
        f.close()

    else:
        if type(shape) is tuple:
            w = np.fromfile(filename, count=np.prod(shape), dtype=dtype)
            w = np.reshape(w, shape[::-1])
            w = np.transpose(w)
        else:
            w = np.fromfile(filename, count=shape, dtype=dtype)

    if totorch:
        w = torch.from_numpy(w).type(torch.FloatTensor)

    return w


def write_array(w, filename, dtype=np.float32, ascii=False):

    w = np.asarray(w, dtype=dtype)
    w = np.transpose(w)

    if ascii:
        np.savetxt(filename, w.transpose(), fmt='%e', delimiter=' ')
    else:
        w.tofile(filename)
