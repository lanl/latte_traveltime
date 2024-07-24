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

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np

def upsample_like(src, tar, mode='trilinear'):

    if mode == 'trilinear':
        src = F.interpolate(src, size=tar.shape[2:], mode='trilinear', align_corners=True)
    else:
        src = F.interpolate(src, size=tar.shape[2:], mode='nearest')
    return src


def pad_like(source, target):

    pZ = target.size()[2] - source.size()[2]
    pY = target.size()[3] - source.size()[3]
    pX = target.size()[4] - source.size()[4]

    return F.pad(source,
                 [pX // 2, pX - pX // 2, pY // 2, pY - pY // 2, pZ // 2, pZ - pZ // 2])


def maxpool(l, kernel_size=2):

    return F.max_pool3d(l, kernel_size=kernel_size)


def upsample(l, scale_factor=2, mode='trilinear'):

    if mode == 'trilinear':
        return F.upsample(l, scale_factor=scale_factor, mode='trilinear', align_corners=True)
    else:
        return F.upsample(l, scale_factor=scale_factor, mode='nearest')


def taper(x, width, apply):

    n1, n2, n3 = x.shape

    if apply[0]:
        for i in range(0, width[0]):
            x[i, :, :] = x[i, :, :] * i * 1.0 / width[0]
    if apply[1]:
        for i in range(0, width[1]):
            x[n1 - 1 - i, :, :] = x[n1 - 1 - i, :, :] * i * 1.0 / width[1]

    if apply[2]:
        for i in range(0, width[2]):
            x[:, i, :] = x[:, i, :] * i * 1.0 / width[2]
    if apply[3]:
        for i in range(0, width[3]):
            x[:, n2 - 1 - i, :] = x[:, n2 - 1 - i, :] * i * 1.0 / width[3]

    if apply[4]:
        for i in range(0, width[4]):
            x[:, :, i] = x[:, :, i] * i * 1.0 / width[4]
    if apply[5]:
        for i in range(0, width[5]):
            x[:, :, n3 - 1 - i] = x[:, :, n3 - 1 - i] * i * 1.0 / width[5]

    return x

class conv(nn.Module):
    def __init__(self,
                 in_channels,
                 out_channels,
                 kernel_size=3,
                 dilation=1,
                 activation='relu',
                 bn=True,
                 bn_type='instance'):

        super(conv, self).__init__()

        layers = []
        layers.append(nn.Conv3d(in_channels, out_channels, kernel_size=kernel_size, dilation=dilation, padding='same'))

        if bn:
            if bn_type == 'batch':
                layers.append(nn.BatchNorm3d(out_channels))
            elif bn_type == 'instance':
                layers.append(nn.InstanceNorm3d(out_channels))
            else:
                layers.append(nn.InstanceNorm3d(out_channels))

        if activation == 'relu':
            layers.append(nn.ReLU(inplace=True))
        elif activation == 'leaky_relu':
            layers.append(nn.LeakyReLU(inplace=True))
        elif activation == 'sigmoid':
            layers.append(nn.Sigmoid())
        elif activation == 'softmax':
            layers.append(nn.Softmax())
        elif activation == 'tanh':
            layers.append(nn.Tanh())

        self.c = nn.Sequential(*layers)
        
    def forward(self, x):

        return self.c(x)


class resu1(nn.Module):
    def __init__(self, in_channels=1, out_channels=1):

        super(resu1, self).__init__()

        n = out_channels

        self.conv0 = conv(in_channels, n)
        self.conv1 = conv(n, n)
        self.conv2 = conv(n, 2 * n)
        self.conv3 = conv(2 * n, 4 * n)
        self.convd = conv(4 * n, 16 * n, dilation=2)
        self.up3 = conv(4 * n + 16 * n, 4 * n)
        self.up2 = conv(2 * n + 4 * n, 2 * n)
        self.up1 = conv(2 * n + n, n)
        self.up0 = conv(n, out_channels)

        self.pool = nn.MaxPool3d(kernel_size=2)
        self.up = nn.Upsample(scale_factor=2, mode='trilinear')

    def forward(self, x):

        conv0 = self.conv0(x)

        conv1 = self.conv1(conv0)
        pool1 = self.pool(conv1)

        conv2 = self.conv2(pool1)
        pool2 = self.pool(conv2)

        conv3 = self.conv3(pool2)
        pool3 = self.pool(conv3)

        d = self.convd(pool3)

        up3 = pad_like(self.up(d), conv3)
        up3 = self.up3(torch.cat((up3, conv3), dim=1))

        up2 = pad_like(self.up(up3), conv2)
        up2 = self.up2(torch.cat((up2, conv2), dim=1))

        up1 = pad_like(self.up(up2), conv1)
        up1 = self.up1(torch.cat((up1, conv1), dim=1))

        up = pad_like(self.up0(up1), conv0)
        up = up + conv0

        return up


class resu2(nn.Module):
    def __init__(self, in_channels=1, out_channels=1):

        super(resu2, self).__init__()

        n = out_channels

        self.conv0 = conv(in_channels, n)
        self.conv1 = conv(n, n)
        self.conv2 = conv(n, 2 * n)
        self.convd = conv(2 * n, 8 * n, dilation=2)
        self.up2 = conv(2 * n + 8 * n, 2 * n)
        self.up1 = conv(2 * n + n, n)
        self.up0 = conv(n, out_channels)

        self.pool = nn.MaxPool3d(kernel_size=2)
        self.up = nn.Upsample(scale_factor=2, mode='trilinear')

    def forward(self, x):

        conv0 = self.conv0(x)

        conv1 = self.conv1(conv0)
        pool1 = self.pool(conv1)

        conv2 = self.conv2(pool1)
        pool2 = self.pool(conv2)

        d = self.convd(pool2)

        up2 = pad_like(self.up(d), conv2)
        up2 = self.up2(torch.cat((up2, conv2), dim=1))

        up1 = pad_like(self.up(up2), conv1)
        up1 = self.up1(torch.cat((up1, conv1), dim=1))

        up = pad_like(self.up0(up1), conv0)
        up = up + conv0

        return up


class resu3(nn.Module):
    def __init__(self, in_channels=1, out_channels=1):

        super(resu3, self).__init__()

        n = out_channels

        self.conv0 = conv(in_channels, n)
        self.conv1 = conv(n, n)
        self.conv2 = conv(n, 2 * n, dilation=2)
        self.conv3 = conv(2 * n, 2 * n, dilation=4)
        self.convd = conv(2 * n, 4 * n, dilation=8)
        self.up3 = conv(2 * n + 4 * n, 2 * n, dilation=4)
        self.up2 = conv(2 * n + 2 * n, 2 * n, dilation=2)
        self.up1 = conv(2 * n + n, n)
        self.up0 = conv(n, out_channels)

    def forward(self, x):

        conv0 = self.conv0(x)
        conv1 = self.conv1(conv0)
        conv2 = self.conv2(conv1)
        conv3 = self.conv3(conv2)

        d = self.convd(conv3)

        up3 = self.up3(torch.cat((d, conv3), dim=1))
        up2 = self.up2(torch.cat((up3, conv2), dim=1))
        up1 = self.up1(torch.cat((up2, conv1), dim=1))
        up = self.up0(up1)
        up = up + conv0

        return up


## Decoder ##
class mtl_decoder(nn.Module):
    def __init__(self, l1, l2, l3, out_ch=1, bn=False, out_activation=None, last_kernel_size=1):

        super(mtl_decoder, self).__init__()

        self.out_ch = out_ch
        self.out_activation = out_activation
        self.last_kernel_size = last_kernel_size
        self.l1 = l1
        self.l2 = l2
        self.l3 = l3
        self.bn = bn

        self.decoder2 = resu2(self.l2 + self.l3, self.l2)
        self.decoder1 = resu1(self.l1 + self.l2, self.l1)

        layers = []
        layers.append(nn.Conv3d(self.l1, self.l1, 3, padding='same'))
        layers.append(nn.Conv3d(self.l1, self.l1, 3, padding='same'))
        layers.append(nn.Conv3d(self.l1, self.out_ch, self.last_kernel_size, padding='same'))

        if self.bn:
            layers.append(nn.InstanceNorm3d(self.out_ch))

        if self.out_activation == 'relu':
            layers.append(nn.ReLU(inplace=True))
        elif self.out_activation == 'softmax':
            layers.append(nn.Softmax())
        elif self.out_activation == 'sigmoid':
            layers.append(nn.Sigmoid())

        self.out = nn.Sequential(*layers)

    def forward(self, x, out_encoder1, out_encoder2, l):

        # pytorch uses chw
        up3 = upsample_like(l, out_encoder2)
        up3 = self.decoder2(torch.cat((up3, out_encoder2), dim=1))

        up2 = upsample_like(up3, out_encoder1)
        up2 = self.decoder1(torch.cat((up2, out_encoder1), dim=1))

        #side output
        out = pad_like(self.out(up2), x)

        return out


class mtl_subdecoder(nn.Module):
    def __init__(self, in_ch, out_ch=1, bn=False, mid_activation=None, activation=None):

        super(mtl_subdecoder, self).__init__()

        # self.layer = l
        self.in_ch = in_ch
        self.out_ch = out_ch
        self.activation = activation
        self.bn = bn
        self.mid_activation = mid_activation

        self.conv = nn.Sequential(conv(self.in_ch, 2 * self.in_ch, activation=self.mid_activation, bn=self.bn),
                                  conv(2 * self.in_ch, 2 * self.in_ch, activation=self.mid_activation, bn=self.bn),
                                  conv(2 * self.in_ch, 2 * self.in_ch, activation=self.mid_activation, bn=self.bn),
                                  conv(2 * self.in_ch, self.in_ch, activation=self.mid_activation, bn=self.bn))
        self.out = nn.Sequential(conv(self.in_ch, 1, activation=self.activation, bn=False))

    def forward(self, x):

        out = self.out(x + self.conv(x))

        return out
