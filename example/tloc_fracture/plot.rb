opts = "-n1=301 -d1=0.01 -d2=0.01 -label1='Depth (km)' -label2='Horizontal Position (km)' -size1=3 -size2=5 -tick1d=1 -tick2d=1 -mtick1=9 -mtick2=9 -color=rainbowcmyk -legend=y -curve=model/sxz.txt,model/rxz.txt -curvesize=10,10 -curveselect=2,1 -curvefacecolor=b,lime -plotlabel='Source':'Receiver' -curveorder=11,11 -curvestyle=scatter*,scatterv -curveedgecolor=none,none "

pmin = 1000
pmax = 3000

# true
system "x_showmatrix -in=model/vp.bin #{opts} -cmin=#{pmin} -cmax=#{pmax} -unit='P-wave Velocity (m/s)' -ld=500 -lmtick=4 -out=tloc_fracture/gt_vp.pdf &"
system "x_showmatrix -in=model/vp_init.bin #{opts} -cmin=#{pmin} -cmax=#{pmax} -unit='P-wave Velocity (m/s)' -ld=500 -lmtick=4 -out=tloc_fracture/init_vp.pdf &"


opts = opts + " -x1beg=1 -x2beg=1 -x2end=4 -curvefacecolor=r,lime "

dir = 'loc_dd_acoustic_reg'
iter = 5

system "x_showmatrix -in=test_#{dir}/iteration_#{iter}/model/source_image.bin #{opts} -cmin=0 -cmax=1 -unit='Source Probability' -ld=0.5 -lmtick=4 -out=tloc_fracture/#{dir}_source_iter_#{iter}.pdf &"

system "x_showmatrix -in=test_#{dir}/iteration_#{iter}/model/source_image.bin.fdip #{opts} -cmin=0 -cmax=1 -unit='Normalized Fault Dip' -ld=0.5 -lmtick=4 -out=tloc_fracture/#{dir}_fdip_iter_#{iter}.pdf &"

iter = 50

system "x_showmatrix -in=test_#{dir}/iteration_#{iter}/model/source_image.bin #{opts} -cmin=0 -cmax=1 -unit='Source Probability' -ld=0.5 -lmtick=4 -out=tloc_fracture/#{dir}_source_iter_#{iter}.pdf &"

system "x_showmatrix -in=test_#{dir}/iteration_#{iter}/model/source_image.bin.fdip #{opts} -cmin=0 -cmax=1 -unit='Normalized Fault Dip' -ld=0.5 -lmtick=4 -out=tloc_fracture/#{dir}_fdip_iter_#{iter}.pdf &"


system "x_showgraph -in=data_no_st0/shot_2_traveltime_p.bin.clean,data_no_st0/shot_2_traveltime_p.bin.noisy,data_no_st0/shot_2_traveltime_p.bin.noise -n1=80,80,80 -reverse2=y -x2beg=-0.1 -x2end=3 -tick2d=0.5 -tick2beg=-0.5 -mtick2=4 -tick1d=20 -mtick1=9 -linecolor=b,r,gray -plotlabelloc=lower_left -size1=4 -size2=3 -label1='Trace Number' -label2='Traveltime (s)' -plotlabel='Clean':'Noisy':'Added Noise' -out=tloc_fracture/data1.pdf &"


system "x_showgraph -in=data_no_st0/shot_12_traveltime_p.bin.clean,data_no_st0/shot_12_traveltime_p.bin.noisy,data_no_st0/shot_12_traveltime_p.bin.noise -n1=80,80,80 -reverse2=y -x2beg=-0.1 -x2end=3 -tick2d=0.5 -tick2beg=-0.5 -mtick2=4 -tick1d=20 -mtick1=9 -linecolor=b,r,gray -plotlabelloc=lower_right -size1=4 -size2=3 -label1='Trace Number' -label2='Traveltime (s)' -plotlabel='Clean':'Noisy':'Added Noise' -out=tloc_fracture/data2.pdf &"
