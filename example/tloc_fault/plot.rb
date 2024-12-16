

dir_noreg = "./test_tloc"
dir_smooth = "./test_tloc_smooth"
dir_reg = "./test_tloc_reg"


system "x_showgraph -in=data_no_st0/shot_2_traveltime_p.bin.clean,data_no_st0/shot_2_traveltime_p.bin.noisy,data_no_st0/shot_2_traveltime_p.bin.noise -n1=150,150,150 -reverse2=y -x2beg=-0.1 -x2end=2 -tick2d=0.5 -tick2beg=-0.5 -mtick2=4 -tick1d=20 -mtick1=9 -linecolor=b,r,gray -plotlabelloc=lower_left -size1=5 -size2=3 -label1='Trace Number' -label2='Traveltime (s)' -plotlabel='Clean':'Noisy':'Added Noise' -out=tloc_fault/data1.pdf &"


system "x_showgraph -in=data_no_st0/shot_1112_traveltime_p.bin.clean,data_no_st0/shot_1112_traveltime_p.bin.noisy,data_no_st0/shot_1112_traveltime_p.bin.noise -n1=150,150,150 -reverse2=y -x2beg=-0.1 -x2end=2 -tick2d=0.5 -tick2beg=-0.5 -mtick2=4 -tick1d=20 -mtick1=9 -linecolor=b,r,gray -plotlabelloc=lower_right -size1=5 -size2=3 -label1='Trace Number' -label2='Traveltime (s)' -plotlabel='Clean':'Noisy':'Added Noise' -out=tloc_fault/data2.pdf &"


fsize = 22

opts = "-n1=101 -d1=0.01 -n2=201 -d2=0.01 -d3=0.01 -label1='Depth (km)' -label2='Y (km)' -label3='X (km)' -tick1d=0.25 -tick2d=0.5 -tick3d=1.0 -mtick1=4 -mtick2=4 -mtick3=9 -color=rainbowcmyk -legend=y -size1=3 -size2=4 -size2=5 -slice1=0.15 -slice2=1 -slice3=1.2 -lticksize=#{fsize} -unitsize=#{fsize} -tick1size=#{fsize} -tick2size=#{fsize} -tick3size=#{fsize} -lmtick=4 -label1size=#{fsize} -label2size=#{fsize} -label3size=#{fsize} "

pmin = 800
pmax = 3000
smin = 350
smax = 2100

# true
system "x_showslice -tr=tloc_fault/gt_fault_sxyz_1.jpg -in=model/vp.bin #{opts} -cmin=#{pmin} -cmax=#{pmax} -unit='P-wave Velocity (m/s)' -ld=500 -lmtick=4 -out=tloc_fault/gt_vp.pdf &"
system "x_showslice -tr=tloc_fault/gt_fault_sxyz_2.jpg -in=model/vp_init.bin #{opts} -cmin=#{pmin} -cmax=#{pmax} -unit='P-wave Velocity (m/s)' -ld=500 -lmtick=4 -out=tloc_fault/init_vp.pdf &"

system "x_showslice -in=#{dir_noreg}/iteration_50/model/updated_vp.bin #{opts} -cmin=#{pmin} -cmax=#{pmax} -unit='P-wave Velocity (m/s)' -ld=500 -lmtick=4 -out=tloc_fault/invt_vp_tomo_dd.pdf &"
system "x_showslice -in=#{dir_noreg}/iteration_50/model/updated_vs.bin #{opts} -cmin=#{smin} -cmax=#{smax} -unit='S-wave Velocity (m/s)' -ld=500 -lmtick=4 -out=tloc_fault/invt_vs_tomo_dd.pdf &"

system "x_showslice -in=#{dir_smooth}/iteration_50/model/updated_vp.bin #{opts} -cmin=#{pmin} -cmax=#{pmax} -unit='P-wave Velocity (m/s)' -ld=500 -lmtick=4 -out=tloc_fault/invt_vp_tomo_dd_smooth.pdf &"
system "x_showslice -in=#{dir_smooth}/iteration_50/model/updated_vs.bin #{opts} -cmin=#{smin} -cmax=#{smax} -unit='S-wave Velocity (m/s)' -ld=500 -lmtick=4 -out=tloc_fault/invt_vs_tomo_dd_smooth.pdf &"

system "x_showslice -in=#{dir_reg}/iteration_50/model/updated_vp.bin #{opts} -cmin=#{pmin} -cmax=#{pmax} -unit='P-wave Velocity (m/s)' -ld=500 -lmtick=4 -out=tloc_fault/invt_vp_tomo_dd_reg.pdf &"
system "x_showslice -in=#{dir_reg}/iteration_50/model/updated_vs.bin #{opts} -cmin=#{smin} -cmax=#{smax} -unit='S-wave Velocity (m/s)' -ld=500 -lmtick=4 -out=tloc_fault/invt_vs_tomo_dd_reg.pdf &"


abort

pmin = -300
pmax = 300
smin = -150
smax = 150

for v in ["vp", "vs"]

    system "x_diff model/#{v}.bin model/#{v}_init.bin >diff_#{v}_gt.bin"
    system "x_diff #{dir_noreg}/iteration_50/model/updated_#{v}.bin model/#{v}_init.bin >diff_#{v}_tloc.bin"
    system "x_diff #{dir_smooth}/iteration_50/model/updated_#{v}.bin model/#{v}_init.bin >diff_#{v}_tloc_smooth.bin"
    system "x_diff #{dir_reg}/iteration_50/model/updated_#{v}.bin model/#{v}_init.bin >diff_#{v}_tloc_reg.bin"

end

#system "rm -rf ./tloc_fault/*.pdf "

for i in [11, 35, 61, 88]
    for name in ["gt", "tloc", "tloc_reg", "tloc_smooth"]

        system "x_select <diff_vp_#{name}.bin n1=101 n2=201 s1=#{i},#{i} >diff_vp_#{name}_slice#{i}.bin "
        system "x_showmatrix -in=diff_vp_#{name}_slice#{i}.bin -n1=201 -d1=0.01 -d2=0.01 -label1='Y (km)' -label2='X (km)' -tick1d=1 -mtick1=9 -tick2d=1 -mtick2=9 -color=rainbowcmyk -legend=y -unit='P-wave Velocity Perturbation (m/s)' -cmin=#{pmin} -cmax=#{pmax} -out=tloc_fault/diff_vp_#{name}_slice#{i}.pdf &"

        system "x_select <diff_vs_#{name}.bin n1=101 n2=201 s1=#{i},#{i} >diff_vs_#{name}_slice#{i}.bin "
        system "x_showmatrix -in=diff_vs_#{name}_slice#{i}.bin -n1=201 -d1=0.01 -d2=0.01 -label1='Y (km)' -label2='X (km)' -tick1d=1 -mtick1=9 -tick2d=1 -mtick2=9 -color=rainbowcmyk -legend=y -unit='S-wave Velocity Perturbation (m/s)' -cmin=#{smin} -cmax=#{smax} -out=tloc_fault/diff_vs_#{name}_slice#{i}.pdf &"

    end
end
