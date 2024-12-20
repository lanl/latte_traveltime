
import matplotlib.pyplot as plt
import numpy
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
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
xrange = [-100, 3100]
yrange = [-100, 2100]
zrange = [-100, 1100]
nr = 1200


label = ['tloc', 'tloc_smooth', 'tloc_reg']

for itest in range(len(label)):

    iter = 50
    iters = [iter]  #[5, 10, iter]
    dir = 'test_' + label[itest]
    output = 'tloc_fault/' + label[itest]
    has_t0 = False
    has_s = True

    subprocess.Popen("mkdir -p tloc_fault", shell=True)
    subprocess.Popen("rm -rf " + output + "_invt_iter*.pdf", shell=True)

    x_obs = read_array('./model/sx.bin', nr)
    y_obs = read_array('./model/sy.bin', nr)
    z_obs = read_array('./model/sz.bin', nr)
    t_obs = read_array('./model/st0.bin', nr)

    p_obs = read_array(dir + '/record_processed/shot_2_traveltime_p.bin', nr) - t_obs
    if has_s:
        s_obs = read_array(dir + '/record_processed/shot_2_traveltime_s.bin', nr) - t_obs
        
    x_init = read_array(dir + '/iteration_0/model/sx.bin', nr)
    y_init = read_array(dir + '/iteration_0/model/sy.bin', nr)
    z_init = read_array(dir + '/iteration_0/model/sz.bin', nr)

    if has_t0:
        t_init = read_array(dir + '/iteration_0/model/st0.bin', nr)
    else:
        t_init = np.zeros(nr)
        
    p_init = read_array(dir + '/iteration_0/synthetic/shot_2_traveltime_p.bin', nr) - t_init
    if has_s:
        s_init = read_array(dir + '/iteration_0/synthetic/shot_2_traveltime_s.bin', nr) - t_init

    x_final = read_array(dir + '/iteration_' + num2str(iter) + '/model/updated_sx.bin', nr)
    y_final = read_array(dir + '/iteration_' + num2str(iter) + '/model/updated_sy.bin', nr)
    z_final = read_array(dir + '/iteration_' + num2str(iter) + '/model/updated_sz.bin', nr)
    if has_t0:
        t_final = read_array(dir + '/iteration_' + num2str(iter) + '/model/updated_st0.bin', nr)
    else:
        t_final = np.zeros(nr)
        
    p_final = read_array(dir + '/iteration_' + num2str(iter) + '/synthetic/shot_2_traveltime_p.bin', nr) - t_final
    if has_s:
        s_final = read_array(dir + '/iteration_' + num2str(iter) + '/synthetic/shot_2_traveltime_s.bin', nr) - t_final

    dist_init = np.sqrt((x_init - x_obs)**2 + (y_init - y_obs)**2 + (z_init - z_obs)**2)
    dmin = np.min(dist_init)
    dmax = np.max(dist_init)
    print(dmin, dmax)

    dist_final = np.sqrt((x_final - x_obs)**2 + (y_final - y_obs)**2 + (z_final - z_obs)**2)
    dmin = np.min(dist_final)
    dmax = np.max(dist_final)
    print(dmin, dmax)

    norm = mplt.colors.Normalize(vmin=0, vmax=3600)

    #=================================================================== gt
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 4))
    x, y = np.meshgrid(regspace(100.0, 3000.0, 200.0), regspace(100.0, 2000.0, 200.0))
    ax.scatter(x, y, 10, 'r', 'v', edgecolors=None)
    # for i in range(len(x_init)):
    #     ax.arrow(x_obs[i], y_obs[i], x_init[i] - x_obs[i], y_init[i] - y_obs[i], color='gray', alpha=0.75, linewidth=0.5)
    sc = ax.scatter(x_obs, y_obs, 10, z_obs, 'o', label='Ground Truth', vmin=0, vmax=1000)
    # ax.scatter(x_init, y_init, 10, z_init, '*', label='Initial Guess')

    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    ax.invert_yaxis()
    ax.set_xlabel('X (m)')
    ax.xaxis.set_label_position('top')
    ax.set_ylabel('Y (m)')
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.yaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_minor_locator(MultipleLocator(100))
    ax.legend(ncols=1, loc='lower right', fontsize=12)
    # ax.grid(True) #, linestyle='dashed')
    cb = plt.colorbar(sc, pad=0.01)
    cb.set_label('Depth (m)')
    cb.ax.yaxis.set_minor_locator(MultipleLocator(40))
    plt.tight_layout()
    plt.savefig(output + '_gt.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
    plt.close()

    #=================================================================== init
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 4))
    x, y = np.meshgrid(regspace(100.0, 3000.0, 200.0), regspace(100.0, 2000.0, 200.0))
    ax.scatter(x, y, 10, 'r', 'v', edgecolors=None)
    sc = ax.scatter(x_init, y_init, 10, z_init, 'o', label='Initial Guess', vmin=0, vmax=1000)

    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    ax.invert_yaxis()
    ax.set_xlabel('X (m)')
    ax.xaxis.set_label_position('top')
    ax.set_ylabel('Y (m)')
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.yaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_minor_locator(MultipleLocator(100))
    ax.legend(ncols=1, loc='lower right', fontsize=12)
    cb = plt.colorbar(sc, pad=0.01)
    cb.set_label('Depth (m)')
    cb.ax.yaxis.set_minor_locator(MultipleLocator(40))
    plt.tight_layout()
    plt.savefig(output + '_init.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
    plt.close()
    
    
    #=================================================================== final
    for l in iters:

        x_final = read_array(dir + '/iteration_' + num2str(l) + '/model/updated_sx.bin', nr)
        y_final = read_array(dir + '/iteration_' + num2str(l) + '/model/updated_sy.bin', nr)
        z_final = read_array(dir + '/iteration_' + num2str(l) + '/model/updated_sz.bin', nr)

        fig, ax = plt.subplots(1, 1, figsize=(5.5, 4))
        # sc = ax.scatter(x_final, y_final, 10, np.abs(z_final - z_obs), 'o', label='Inverted', vmin=0, vmax=100, cmap='magma_r')
        
        for i in range(len(x_init)):
            ax.arrow(x_obs[i], y_obs[i], x_final[i] - x_obs[i], y_final[i] - y_obs[i], color='gray', alpha=0.75, linewidth=0.5)
        sc = ax.scatter(x_obs, y_obs, 5, 'b', 'o', label='Ground Truth')
        ax.scatter(x_final, y_final, 5, 'r', 'o', label='Inverted')

        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
        ax.invert_yaxis()
        ax.set_xlabel('X (m)')
        ax.xaxis.set_label_position('top')
        ax.set_ylabel('Y (m)')
        ax.xaxis.set_major_locator(MultipleLocator(1000))
        ax.xaxis.set_minor_locator(MultipleLocator(100))
        ax.yaxis.set_major_locator(MultipleLocator(1000))
        ax.yaxis.set_minor_locator(MultipleLocator(100))
        ax.legend(ncols=1, loc='lower right', fontsize=12)
        plt.tight_layout()
        plt.savefig(output + '_invt_iter_' + num2str(l) + '_xy.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
                
        fig, ax = plt.subplots(1, 1, figsize=(6.5, 4))
        for i in range(len(x_init)):
            ax.arrow(x_obs[i], z_obs[i], x_final[i] - x_obs[i], z_final[i] - z_obs[i], color='gray', alpha=0.75, linewidth=0.5)
        sc = ax.scatter(x_obs, z_obs, 5, 'b', 'o', label='Ground Truth')
        ax.scatter(x_final, z_final, 5, 'r', 'o', label='Inverted')

        ax.set_xlim(xrange)
        ax.set_ylim(zrange)
        ax.invert_yaxis()
        ax.set_xlabel('X (m)')
        ax.xaxis.set_label_position('top')
        ax.set_ylabel('Depth (m)')
        ax.xaxis.set_major_locator(MultipleLocator(1000))
        ax.xaxis.set_minor_locator(MultipleLocator(100))
        ax.yaxis.set_major_locator(MultipleLocator(250))
        ax.yaxis.set_minor_locator(MultipleLocator(50))
        ax.legend(ncols=1, loc='lower right', fontsize=12)
        plt.tight_layout()
        plt.savefig(output + '_invt_iter_' + num2str(l) + '_xz.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
        
        fig, ax = plt.subplots(1, 1, figsize=(6.5, 4))
        for i in range(len(x_init)):
            ax.arrow(y_obs[i], z_obs[i], y_final[i] - y_obs[i], z_final[i] - z_obs[i], color='gray', alpha=0.75, linewidth=0.5)
        sc = ax.scatter(y_obs, z_obs, 5, 'b', 'o', label='Ground Truth')
        ax.scatter(y_final, z_final, 5, 'r', 'o', label='Inverted')
        
        ax.set_xlim(yrange)
        ax.set_ylim(zrange)
        ax.invert_yaxis()
        ax.set_xlabel('Y (m)')
        ax.xaxis.set_label_position('top')
        ax.set_ylabel('Depth (m)')
        ax.xaxis.set_major_locator(MultipleLocator(1000))
        ax.xaxis.set_minor_locator(MultipleLocator(100))
        ax.yaxis.set_major_locator(MultipleLocator(250))
        ax.yaxis.set_minor_locator(MultipleLocator(50))
        ax.legend(ncols=1, loc='lower right', fontsize=12)
        plt.tight_layout()
        plt.savefig(output + '_invt_iter_' + num2str(l) + '_yz.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
        plt.close()
        
    #=================================================================== final
    for l in iters:

        x_final = read_array(dir + '/iteration_' + num2str(l) + '/model/updated_sx.bin', nr)
        y_final = read_array(dir + '/iteration_' + num2str(l) + '/model/updated_sy.bin', nr)
        z_final = read_array(dir + '/iteration_' + num2str(l) + '/model/updated_sz.bin', nr)

        fig, ax = plt.subplots(1, 1, figsize=(6.5, 4))
        sc = ax.scatter(x_final, y_final, 10, np.abs(z_final - z_obs), 'o', label='Inverted', vmin=0, vmax=100, cmap='magma_r')
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
        ax.invert_yaxis()
        ax.set_xlabel('X (m)')
        ax.xaxis.set_label_position('top')
        ax.set_ylabel('Y (m)')
        ax.xaxis.set_major_locator(MultipleLocator(1000))
        ax.xaxis.set_minor_locator(MultipleLocator(100))
        ax.yaxis.set_major_locator(MultipleLocator(1000))
        ax.yaxis.set_minor_locator(MultipleLocator(100))
        ax.legend(ncols=1, loc='lower right', fontsize=12)
        cb = plt.colorbar(sc, pad=0.01)
        cb.set_label('Depth Error (m)')
        cb.ax.yaxis.set_minor_locator(MultipleLocator(10))
        plt.tight_layout()
        plt.savefig(output + '_invt_iter_' + num2str(l) + '_z_error.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
        
        fig, ax = plt.subplots(1, 1, figsize=(6.5, 4))
        sc = ax.scatter(x_final, y_final, 10, np.sqrt((x_final - x_obs)**2 + (y_final - y_obs)**2), 'o', label='Inverted', vmin=0, vmax=100, cmap='magma_r')
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
        ax.invert_yaxis()
        ax.set_xlabel('X (m)')
        ax.xaxis.set_label_position('top')
        ax.set_ylabel('Y (m)')
        ax.xaxis.set_major_locator(MultipleLocator(1000))
        ax.xaxis.set_minor_locator(MultipleLocator(100))
        ax.yaxis.set_major_locator(MultipleLocator(1000))
        ax.yaxis.set_minor_locator(MultipleLocator(100))
        ax.legend(ncols=1, loc='lower right', fontsize=12)
        cb = plt.colorbar(sc, pad=0.01)
        cb.set_label('Horizontal Error (m)')
        cb.ax.yaxis.set_minor_locator(MultipleLocator(10))
        plt.tight_layout()
        plt.savefig(output + '_invt_iter_' + num2str(l) + '_xy_error.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
        
        fig, ax = plt.subplots(1, 1, figsize=(6.5, 4))
        sc = ax.scatter(x_final, y_final, 10, np.sqrt((x_final - x_obs)**2 + (y_final - y_obs)**2 + (z_final - z_obs)**2), 'o', label='Inverted', vmin=0, vmax=200, cmap='magma_r')
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
        ax.invert_yaxis()
        ax.set_xlabel('X (m)')
        ax.xaxis.set_label_position('top')
        ax.set_ylabel('Y (m)')
        ax.xaxis.set_major_locator(MultipleLocator(1000))
        ax.xaxis.set_minor_locator(MultipleLocator(100))
        ax.yaxis.set_major_locator(MultipleLocator(1000))
        ax.yaxis.set_minor_locator(MultipleLocator(100))
        ax.legend(ncols=1, loc='lower right', fontsize=12)
        cb = plt.colorbar(sc, pad=0.01)
        cb.set_label('Absolute Error (m)')
        cb.ax.yaxis.set_minor_locator(MultipleLocator(10))
        plt.tight_layout()
        plt.savefig(output + '_invt_iter_' + num2str(l) + '_xyz_error.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
        plt.close()

    #=================================================================== convergence
    if has_t0:
        fig, ax = plt.subplots(1, 1, figsize=(8, 3))

        ax.plot(regspace(0.0, nr - 1.0, 1.0), t_obs, 'b', linewidth=0.5, label='Ground Truth')
        ax.plot(regspace(0.0, nr - 1.0, 1.0), t_init, 'lime', linewidth=0.5, label='Initial Guess')
        ax.plot(regspace(0.0, nr - 1.0, 1.0), t_final, 'r', linewidth=0.5, label='Inverted')

        # ax.scatter(regspace(0.0, nr - 1.0, 1.0), t_obs, 20, 'b', marker='o', label='Ground Truth')
        # ax.scatter(regspace(0.0, nr - 1.0, 1.0), t_init, 20, 'lime', marker='*', label='Initial Guess')
        # ax.scatter(regspace(0.0, nr - 1.0, 1.0), t_final, 20, 'r', marker='*', label='Inverted')

        ax.set_xlim([0, nr - 1.0])
        ax.set_ylim([-4, 11])
        ax.xaxis.set_major_locator(MultipleLocator(50))
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.yaxis.set_major_locator(MultipleLocator(2))
        ax.yaxis.set_minor_locator(MultipleLocator(1))
        ax.set_xlabel('Source Index')
        ax.xaxis.set_label_position('top')
        ax.set_ylabel('Origin Time (s)')
        ax.legend(ncols=3)
        # ax.grid(True) #, linestyle='dashed')

        # plt.show()
        plt.tight_layout()
        plt.savefig(output + '_st0.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
        plt.close()

    #=================================================================== convergence
    fig, ax = plt.subplots(1, 1, figsize=(9, 3))

    ax.plot(regspace(0.0, nr - 1.0, 1.0), p_obs, 'b', linewidth=1, label='Ground Truth')
    ax.plot(regspace(0.0, nr - 1.0, 1.0), p_init, 'lime', linewidth=1, label='Initial Model')
    ax.plot(regspace(0.0, nr - 1.0, 1.0), p_final, 'r', linewidth=1, label='Inverted Model')
    ax.plot(regspace(0.0, nr - 1.0, 1.0), p_final - p_obs, 'gray', linewidth=1, label='Final Misfit')

    ax.set_xlim([-1, nr - 1.0])
    ax.set_ylim([-1, 2.5])
    ax.xaxis.set_major_locator(MultipleLocator(200))
    ax.xaxis.set_minor_locator(MultipleLocator(20))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.invert_yaxis()
    ax.set_xlabel('Source Index')
    ax.xaxis.set_label_position('top')
    ax.set_ylabel('Traveltime (s)')
    ax.legend(ncols=4, loc='lower left', fontsize=12)
    ax.invert_yaxis()
    # ax.grid(True) #, linestyle='dashed')

    # plt.show()
    plt.tight_layout()
    plt.savefig(output + '_time_p.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)

    if has_s:

        fig, ax = plt.subplots(1, 1, figsize=(9, 3))

        ax.plot(regspace(0.0, nr - 1.0, 1.0), s_obs, 'b', linewidth=1, label='Ground Truth')
        ax.plot(regspace(0.0, nr - 1.0, 1.0), s_init, 'lime', linewidth=1, label='Initial Model')
        ax.plot(regspace(0.0, nr - 1.0, 1.0), s_final, 'r', linewidth=1, label='Inverted Model')
        ax.plot(regspace(0.0, nr - 1.0, 1.0), s_final - s_obs, 'r', linewidth=1, label='Final Misfit')

        ax.set_xlim([0, nr - 1.0])
        ax.set_ylim([-1, 4])
        ax.xaxis.set_major_locator(MultipleLocator(25))
        ax.xaxis.set_minor_locator(MultipleLocator(5))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.invert_yaxis()
        ax.set_xlabel('Source Index')
        ax.xaxis.set_label_position('top')
        ax.set_ylabel('Traveltime (s)')
        ax.legend(ncols=4, loc='lower left', fontsize=12)
        ax.invert_yaxis()
        # ax.grid(True) #, linestyle='dashed')

        # plt.show()
        plt.tight_layout()
        plt.savefig(output + '_time_s.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
        plt.close()



    #=================================================================== convergence
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    c = read_array(dir + '/data_misfit.txt', (iter + 1, 3), ascii=True)
    ax.plot(c[:, 0], c[:, 2], 'b', linewidth=2)

    ax.set_xlim([0, iter])
    ax.set_ylim([0, 1])
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_xlabel('Number of Iterations')
    ax.set_ylabel('Normalized Data Misfit')
    ax.grid(True) #, linestyle='dashed')

    # plt.show()
    plt.tight_layout()
    plt.savefig(output + '_misfit.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
    plt.close()
    
    #=================================================================== convergence
    x_final = read_array(dir + '/iteration_' + num2str(l) + '/model/updated_sx.bin', nr)
    y_final = read_array(dir + '/iteration_' + num2str(l) + '/model/updated_sy.bin', nr)
    z_final = read_array(dir + '/iteration_' + num2str(l) + '/model/updated_sz.bin', nr)
    data = np.sqrt((x_final - x_obs)**2 + (y_final - y_obs)**2 + (z_final - z_obs)**2)
    
    # from scipy.stats import f
    # d1, d2, loc, scale = f.fit(data)
    # xx = np.linspace(np.min(data), np.max(data), 100)
    # pdf = f.pdf(xx, d1, d2, scale=scale)

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    vmin, vmax = np.amin(data), np.amax(data)
    lim = np.amax([np.abs(vmin), vmax]) * 1.1
    x, bins, p = ax.hist(
        data,
        bins=regspace(0.0, 200.0, 10.0),
        density=True,
        facecolor="royalblue",
        alpha=1,
        ec="black",
        align="mid",
    )
    for item in p:
        item.set_height(item.get_height() / sum(x))
      
    axx = ax.twinx()  
    axx.hist(
        data,
        bins=regspace(0.0, 200.0, 10.0),
        density=True,
        cumulative=True,
        histtype='step',
        ec='r',
        linewidth=2
    )
    axx.set_ylim([0, 1])
    axx.yaxis.set_major_locator(MultipleLocator(0.2))
    axx.yaxis.set_minor_locator(MultipleLocator(0.05))
    axx.set_ylabel("Normalized Cumlative Count", fontsize=14, color='r')
    axx.tick_params(axis='y', labelcolor='r')
        
    # ax.plot(xx, pdf*10, linestyle='--', color='r')
    ax.set_xlabel("Location Error (m)", fontsize=14)
    ax.set_axisbelow(True)
    ax.set_aspect('auto')
    ax.grid(color="gray", linestyle="dashed", linewidth=1)
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.02))
    ax.set_xlim([0, 200])
    ax.set_ylim([0, 0.3])
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.xaxis.set_ticks_position('bottom')  # Move ticks to the top
    ax.xaxis.set_label_position('bottom')  # Move labels to the top
    ax.tick_params(axis='x', top=False, labeltop=False)
    label1 = "{:.3f}$".format(vmin)
    label2 = "{:.3f}$".format(vmax)
    ax.set_ylabel("Normalized Frequency Count", fontsize=14)
    plt.savefig(output + '_stats.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
    plt.close()

    # exit(0)
    
