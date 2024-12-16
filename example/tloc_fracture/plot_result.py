
import matplotlib as mplt
import matplotlib.pyplot as plt
import numpy
from matplotlib.ticker import (MultipleLocator, FixedLocator, AutoMinorLocator)
from python.libpy_io import *
from python.libpy_filedir import *
from python.libpy_utility import *
import subprocess

# fonts
set_font()
plt.rcParams.update({'font.size': 14})  # Replace 14 with your desired font size
rcParams["xtick.top"] = True
rcParams["xtick.labeltop"] = True
rcParams["xtick.bottom"] = False
rcParams["xtick.labelbottom"] = False


# params
xrange = [900, 4100]
yrange = [1000, 3100]
nr = 1200

maxiter = [50, 50]
label = ['loc_dd_acoustic', 'loc_dd_acoustic_reg']
t0s = [False, False]
mms = [1.0e-4, 1.0e-4]

for itest in range(2):

    iter = maxiter[itest]
    iters = [5, 10, iter]
    dir = 'test_' + label[itest]
    output = 'tloc_fracture/' + label[itest]
    has_t0 = t0s[itest]
    has_s = False
    mm = mms[itest]

    subprocess.Popen("mkdir -p tloc_fracture", shell=True)
    subprocess.Popen("rm -rf " + output + "_invt_iter_*.pdf", shell=True)

    x_obs = read_array('./model/sx.bin', nr)
    z_obs = read_array('./model/sz.bin', nr)
    t_obs = read_array('./model/st0.bin', nr)
    p_obs = read_array(dir + '/record_processed/shot_2_traveltime_p.bin', nr) - t_obs
    if has_s:
        s_obs = read_array(dir + '/record_processed/shot_2_traveltime_s.bin', nr) - t_obs
    x_init = read_array(dir + '/iteration_0/model/sx.bin', nr)
    z_init = read_array(dir + '/iteration_0/model/sz.bin', nr)
    if has_t0:
        t_init = read_array(dir + '/iteration_0/model/st0.bin', nr)
    else:
        t_init = np.zeros(nr)
    p_init = read_array(dir + '/iteration_0/synthetic/shot_2_traveltime_p.bin', nr) - t_init
    if has_s:
        s_init = read_array(dir + '/iteration_0/synthetic/shot_2_traveltime_s.bin', nr) - t_init

    x_final = read_array(dir + '/iteration_' + num2str(iter) + '/model/updated_sx.bin', nr)
    z_final = read_array(dir + '/iteration_' + num2str(iter) + '/model/updated_sz.bin', nr)
    if has_t0:
        t_final = read_array(dir + '/iteration_' + num2str(iter) + '/model/updated_st0.bin', nr)
    else:
        t_final = np.zeros(nr)
    p_final = read_array(dir + '/iteration_' + num2str(iter) + '/synthetic/shot_2_traveltime_p.bin', nr) - t_final
    if has_s:
        s_final = read_array(dir + '/iteration_' + num2str(iter) + '/synthetic/shot_2_traveltime_s.bin', nr) - t_final

    dist_init = np.sqrt((x_init - x_obs)**2 + (z_init - z_obs)**2)
    dmin = np.min(dist_init)
    dmax = np.max(dist_init)
    print(dmin, dmax)

    dist_final = np.sqrt((x_final - x_obs)**2 + (z_final - z_obs)**2)
    dmin = np.min(dist_final)
    dmax = np.max(dist_final)
    print(dmin, dmax)

    norm = mplt.colors.Normalize(vmin=0, vmax=3600)

    #=================================================================== initial
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    ax.scatter(regspace(50.0, 4000.0, 100.0), np.zeros(40), 50, 'yellow', 'v', edgecolors='k')
    for i in range(len(x_init)):
        ax.arrow(x_obs[i], z_obs[i], x_init[i] - x_obs[i], z_init[i] - z_obs[i], color='gray', alpha=0.75, linewidth=0.5)
    ax.scatter(x_obs, z_obs, 10, 'b', 'o', label='Ground Truth')
    ax.scatter(x_init, z_init, 10, 'r', 'o', label='Initial Guess')

    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    ax.invert_yaxis()
    ax.set_xlabel('Horizontal Position (m)')
    ax.xaxis.set_label_position('top')
    ax.set_ylabel('Depth (m)')
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.yaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_minor_locator(MultipleLocator(100))
    ax.legend(ncols=2)
    # ax.grid(True) #, linestyle='dashed')
    plt.tight_layout()
    plt.savefig(output + '_init.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)

    #=================================================================== final
    for l in iters:

        x_final = read_array(dir + '/iteration_' + num2str(l) + '/model/updated_sx.bin', nr)
        z_final = read_array(dir + '/iteration_' + num2str(l) + '/model/updated_sz.bin', nr)
            
        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
        ax.scatter(regspace(50.0, 4000.0, 100.0), np.zeros(40), 50, 'yellow', 'v', edgecolors='k')
        for i in range(len(x_final)):
            ax.arrow(x_obs[i], z_obs[i], x_final[i] - x_obs[i], z_final[i] - z_obs[i], color='gray', alpha=0.75, linewidth=0.5)
        ax.scatter(x_obs, z_obs, 10, 'b', 'o', label='Ground Truth')
        ax.scatter(x_final, z_final, 10, 'r', 'o', label='Inverted')

        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
        ax.invert_yaxis()
        ax.set_xlabel('Horizontal Position (m)')
        ax.xaxis.set_label_position('top')
        ax.set_ylabel('Depth (m)')
        ax.xaxis.set_major_locator(MultipleLocator(1000))
        ax.xaxis.set_minor_locator(MultipleLocator(100))
        ax.yaxis.set_major_locator(MultipleLocator(1000))
        ax.yaxis.set_minor_locator(MultipleLocator(100))
        ax.legend(ncols=2)
        # ax.grid(True) #, linestyle='dashed')
        plt.tight_layout()
        plt.savefig(output + '_invt_iter_' + num2str(l) + '.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)

    #=================================================================== convergence
    if has_t0:
        fig, ax = plt.subplots(1, 1, figsize=(7, 4))

        ax.plot(regspace(0.0, nr - 1.0, 1.0), t_obs, 'b', label='Ground Truth')
        ax.plot(regspace(0.0, nr - 1.0, 1.0), t_init, 'lime', label='Initial Guess')
        ax.plot(regspace(0.0, nr - 1.0, 1.0), t_final, 'r', label='Inverted')
        ax.scatter(regspace(0.0, nr - 1.0, 1.0), t_final - t_obs, 5, 'gray', label='Final Misfit')

        # ax.scatter(regspace(0.0, nr - 1.0, 1.0), t_obs, 20, 'b', marker='o', label='Ground Truth')
        # ax.scatter(regspace(0.0, nr - 1.0, 1.0), t_init, 20, 'lime', marker='*', label='Initial Guess')
        # ax.scatter(regspace(0.0, nr - 1.0, 1.0), t_final, 20, 'r', marker='*', label='Inverted')

        ax.set_xlim([0, nr - 1.0])
        ax.set_ylim([-5, 11])
        ax.xaxis.set_major_locator(MultipleLocator(50))
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.yaxis.set_major_locator(MultipleLocator(3))
        ax.yaxis.set_minor_locator(MultipleLocator(1))
        ax.set_xlabel('Source Index')
        ax.xaxis.set_label_position('top')
        ax.set_ylabel('Origin Time (s)')
        ax.legend(ncols=2)
        # ax.grid(True) #, linestyle='dashed')

        # plt.show()
        plt.tight_layout()
        plt.savefig(output + '_st0.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)

    #=================================================================== traveltime
    fig, ax = plt.subplots(1, 1, figsize=(7, 4))

    ax.plot(regspace(0.0, nr - 1.0, 1.0), p_obs, 'b', linewidth=1, label='Ground Truth')
    ax.plot(regspace(0.0, nr - 1.0, 1.0), p_init, 'lime', linewidth=1, label='Initial Model')
    ax.plot(regspace(0.0, nr - 1.0, 1.0), p_final, 'r', linewidth=1, label='Inverted Model')
    ax.scatter(regspace(0.0, nr - 1.0, 1.0), p_final - p_obs, 5, 'gray', label='Final Misfit')


    ax.set_xlim([0, nr - 1.0])
    ax.set_ylim([-1.5, 3])
    ax.xaxis.set_major_locator(MultipleLocator(200))
    ax.xaxis.set_minor_locator(MultipleLocator(20))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.invert_yaxis()
    ax.set_xlabel('Source Index')
    ax.xaxis.set_label_position('top')
    ax.set_ylabel('Traveltime (s)')
    ax.legend(ncols=2, loc='lower left')
    ax.invert_yaxis()
    # ax.grid(True) #, linestyle='dashed')

    # plt.show()
    plt.tight_layout()
    plt.savefig(output + '_time_p.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)

    if has_s:

        fig, ax = plt.subplots(1, 1, figsize=(7, 4))

        ax.plot(regspace(0.0, nr - 1.0, 1.0), s_obs, 'b', linewidth=1, label='Ground Truth')
        ax.plot(regspace(0.0, nr - 1.0, 1.0), s_init, 'lime', linewidth=1, label='Initial Model')
        ax.plot(regspace(0.0, nr - 1.0, 1.0), s_final, 'r', linewidth=1, label='Inverted Model')
        ax.scatter(regspace(0.0, nr - 1.0, 1.0), s_final - s_obs, 5, 'gray', label='Final Misfit')

        ax.set_xlim([0, nr - 1.0])
        ax.set_ylim([-2.2, 4.5])
        ax.xaxis.set_major_locator(MultipleLocator(25))
        ax.xaxis.set_minor_locator(MultipleLocator(5))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.invert_yaxis()
        ax.set_xlabel('Source Index')
        ax.xaxis.set_label_position('top')
        ax.set_ylabel('Traveltime (s)')
        ax.legend(ncols=2, loc='lower left')
        ax.invert_yaxis()
        # ax.grid(True) #, linestyle='dashed')

        # plt.show()
        plt.tight_layout()
        plt.savefig(output + '_time_s.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)



    #=================================================================== convergence
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    c = read_array(dir + '/data_misfit.txt', (iter + 1, 3), ascii=True)
    ax.plot(c[:, 0], c[:, 2], 'b', linewidth=2)

    ax.set_xlim([0, 40])
    ax.set_ylim([mm, 1])
    ax.set_yscale('log')
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(MultipleLocator(2))
    ax.yaxis.set_major_locator(FixedLocator(locs = np.logspace(np.log10(mm), 0, np.int16(np.abs(np.log10(mm))) + 1)))
    # ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.set_xlabel('Number of Iterations')
    ax.xaxis.set_label_position('top')
    ax.set_ylabel('Normalized Data Misfit')
    ax.grid(True) #, linestyle='dashed')

    # plt.show()
    plt.tight_layout()
    plt.savefig(output + '_misfit.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
