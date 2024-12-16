import os
import warnings
import sys
import argparse
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from python.libpy_io import *
from python.libpy_filedir import *
from python.libpy_utility import *
import subprocess

# fonts
set_font()

plt.rcParams.update({'font.size': 14})  # Replace 14 with your desired font size


tags = ['ad', 'dd']

ns = 40
nr = 401

data = np.empty(0)
data_ad = np.empty(0)
data_dd = np.empty(0)
for i in regspace(1, ns, 1):
    data = np.append(data, read_array('./data/shot_' + num2str(i) + '_traveltime_p.bin', nr))
    data_ad = np.append(data_ad, read_array('./test_ad/iteration_100/synthetic/shot_' + num2str(i) + '_traveltime_p.bin', nr))
    data_dd = np.append(data_dd, read_array('./test_dd/iteration_100/synthetic/shot_' + num2str(i) + '_traveltime_p.bin', nr))
    
data_ad = data_ad - data
data_dd = data_dd - data

d = np.zeros((nr*ns, 2))
d[:, 0] = data_ad
d[:, 1] = data_dd

for i in regspace(0, 1, 1):
    
    f, ax = plt.subplots(1, 1, figsize=(4, 4), squeeze=True)
    cmin = -0.03
    cmax = 0.03

    data = d[:, i]
    mu, sigma = norm.fit(data)
    print(mu, sigma)

    vmin, vmax = np.amin(data), np.amax(data)
    lim = np.amax([np.abs(vmin), vmax]) * 1.1
    x, bins, p = ax.hist(
        data,
        regspace(-0.03, 0.03, 0.0025),
        density=True,
        facecolor="royalblue",
        alpha=1,
        ec="black",
        align="mid",
    )
    for item in p:
        item.set_height(item.get_height() / sum(x))
    ax.set_xlabel("Traveltime Misfit (s)", fontsize=14)
    ax.set_axisbelow(True)
    ax.set_aspect('auto')
    ax.grid(color="gray", linestyle="dashed", linewidth=1)
    ax.set_xlim([cmin, cmax])
    ax.set_ylim([0, 0.5])
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)
    label1 = "{:.3f}$".format(vmin)
    label2 = "{:.3f}$".format(vmax)
    ax.set_title(
        " Error Range = [" + label1 + ", " + label2 + "]\n" + "$\\mu$ = {:.3f}".format(mu) + ", " +
        "$\\sigma$ = {:.3f}".format(sigma),
        fontsize=14,
    )

    y = np.linspace(cmin, cmax, 100)
    p = norm.pdf(y, mu, sigma)
    ax.plot(y, p / np.sum(x), "r--", linewidth=2)
    ax.set_ylabel("Normalized Frequency Count", fontsize=14)

    plt.savefig('./fatt/misfit_stats_' + tags[i] + ".pdf", dpi=300, bbox_inches="tight", pad_inches=0.01)
    plt.close()
