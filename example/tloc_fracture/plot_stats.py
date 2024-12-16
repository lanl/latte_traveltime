import os
import warnings
import sys
import argparse
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from tqdm import tqdm
from pathlib import Path
warnings.filterwarnings("ignore")
sys.path.append(str(Path.home()) + "/src/python")
from libpy_utility import *
from libpy_filedir import *
from libpy_io import *
from libpy_getpar import *
import matplotlib as mplt
from matplotlib import rcParams
from scipy.stats import norm

basefamily = "sans-serif"
basefont = "Arial"
fontset = "custom"
rcParams["font.family"] = basefamily
rcParams["font." + basefamily] = basefont
mplt.rcParams["mathtext.fontset"] = fontset
mplt.rcParams["mathtext.rm"] = basefont
mplt.rcParams["mathtext.sf"] = basefont
mplt.rcParams["mathtext.it"] = basefont + ":italic"
mplt.rcParams["mathtext.bf"] = basefont + ":bold"

tags = ['loc_dd_acoustic', 'loc_dd_acoustic_reg']

ns = 40
nr = 1200

data = np.empty(0)
data_ad = np.empty(0)
data_dd = np.empty(0)

st0 = read_array('./model/st0.bin', nr)

for i in regspace(1, ns, 1):
    data = np.append(data, read_array('./test_loc_dd_acoustic/record_processed/shot_' + num2str(i) + '_traveltime_p.bin', nr) - st0)
    data_ad = np.append(data_ad, read_array('./test_loc_dd_acoustic/iteration_50/synthetic/shot_' + num2str(i) + '_traveltime_p.bin', nr))
    data_dd = np.append(data_dd, read_array('./test_loc_dd_acoustic_reg/iteration_50/synthetic/shot_' + num2str(i) + '_traveltime_p.bin', nr))
    
data_ad = data_ad - data
data_dd = data_dd - data

d = np.zeros((nr*ns, 2))
d[:, 0] = data_ad
d[:, 1] = data_dd

for i in regspace(0, 1, 1):
    
    f, ax = plt.subplots(1, 1, figsize=(5, 5), squeeze=True)
    cmin = -0.1
    cmax = 0.1

    data = d[:, i]
    mu, sigma = norm.fit(data)
    print(mu, sigma)

    vmin, vmax = np.amin(data), np.amax(data)
    lim = np.amax([np.abs(vmin), vmax]) * 1.1
    x, bins, p = ax.hist(
        data,
        10,
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
        " Error Range = [" + label1 + ", " + label2 + "]\n" + "$\mu$ = {:.3f}".format(mu) + ", " +
        "$\sigma$ = {:.3f}".format(sigma),
        fontsize=14,
    )

    y = np.linspace(cmin, cmax, 100)
    p = norm.pdf(y, mu, sigma)
    ax.plot(y, p / np.sum(x), "r--", linewidth=2)
    ax.set_ylabel("Normalized Count", fontsize=14)

    plt.savefig('./tloc_fracture/misfit_stats_' + tags[i] + ".pdf", dpi=300, bbox_inches="tight", pad_inches=0.01)
    plt.show()
