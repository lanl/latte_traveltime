
opts = "-n1=301 -d1=0.01 -d2=0.01 -label1='Depth (km)' -label2='Horizontal Position (km)' -size1=3 -size2=4 -tick1d=1 -tick2d=1 -mtick1=9 -mtick2=9 -color=rainbowcmyk -legend=y "

pmin = 1800
pmax = 2200
smin = 1050
smax = 1250

# true
system "x_showmatrix -in=model/vp.bin #{opts} -cmin=#{pmin} -cmax=#{pmax} -unit='P-wave Velocity (m/s)' -ld=100 -lmtick=9 -out=tloc/gt_vp.pdf &"

system "x_showmatrix -in=model/vs.bin #{opts} -cmin=#{smin} -cmax=#{smax} -unit='S-wave Velocity (m/s)' -ld=50 -lmtick=4 -out=tloc/gt_vs.pdf &"
system "x_showmatrix -in=model/vp_init.bin #{opts} -cmin=#{pmin} -cmax=#{pmax} -unit='P-wave Velocity (m/s)' -ld=100 -lmtick=9 -out=tloc/init_vp.pdf &"
system "x_showmatrix -in=model/vs_init.bin #{opts} -cmin=#{smin} -cmax=#{smax} -unit='S-wave Velocity (m/s)' -ld=50 -lmtick=9 -out=tloc/init_vs.pdf &"

dir = 'tomo_dd_elastic_reg'
iter = 100

system "x_showmatrix -in=test_#{dir}/iteration_#{iter}/model/updated_vp.bin #{opts} -cmin=#{pmin} -cmax=#{pmax} -unit='P-wave Velocity (m/s)' -ld=100 -lmtick=9 -out=tloc/#{dir}_vp.pdf &"
system "x_showmatrix -in=test_#{dir}/iteration_#{iter}/model/updated_vs.bin #{opts} -cmin=#{smin} -cmax=#{smax} -unit='S-wave Velocity (m/s)' -ld=50 -lmtick=4 -out=tloc/#{dir}_vs.pdf &"


dir = 'tomo_dd_elastic'
iter = 100

system "x_showmatrix -in=test_#{dir}/iteration_#{iter}/model/updated_vp.bin #{opts} -cmin=#{pmin} -cmax=#{pmax} -unit='P-wave Velocity (m/s)' -ld=100 -lmtick=9 -out=tloc/#{dir}_vp.pdf &"
system "x_showmatrix -in=test_#{dir}/iteration_#{iter}/model/updated_vs.bin #{opts} -cmin=#{smin} -cmax=#{smax} -unit='S-wave Velocity (m/s)' -ld=50 -lmtick=4 -out=tloc/#{dir}_vs.pdf &"
