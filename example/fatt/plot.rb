
n1 = 51
n2 = 401
s = 3
outdir = "./fatt"

# model and reflectors
opts = " -n1=#{n1} -d1=0.01 -d2=0.01 -label1='Depth (km)' -label2='Horizontal Position (km)' -size1=2 -size2=6.5 -tick1d=0.1 -mtick1=1 -tick2d=1 -mtick2=9 -legend=y -lmtick=4 -color=jet "

system "mkdir -p #{outdir}"

system "x_showmatrix -in=model/vp_init.bin #{opts} -cmin=500 -cmax=2250 -ld=500 -ltickbeg=500 -color=rainbowcmyk -unit='P-wave Velocity (m/s)' -out=#{outdir}/vp_init.pdf &"

system "x_showmatrix -in=model/vp.bin #{opts} -cmin=500 -cmax=2250 -ld=500 -ltickbeg=500 -color=rainbowcmyk -unit='P-wave Velocity (m/s)' -out=#{outdir}/vp.pdf &"

system "x_showmatrix -in=test_ad/iteration_50/model/updated_vp.bin #{opts} -cmin=500 -cmax=2250 -ld=500 -ltickbeg=500 -color=rainbowcmyk -unit='P-wave Velocity (m/s)' -out=#{outdir}/vp_ad.pdf &"

system "x_showmatrix -in=test_dd/iteration_50/model/updated_vp.bin #{opts} -cmin=500 -cmax=2250 -ld=500 -ltickbeg=500 -color=rainbowcmyk -unit='P-wave Velocity (m/s)' -out=#{outdir}/vp_dd.pdf &"

system "x_showgraph -in=test_ad/data_misfit.txt,test_dd/data_misfit.txt -ftype=ascii -ptype=2 -n1=101,101 -select=1,3 -size1=4 -size2=3 -norm2=log -linewidth=2,2 -x2beg=0.0005 -x2end=1 -x1beg=0 -x1end=100 -tick1d=20 -mtick1=9 -mtick2=9 -linecolor=b,r -plotlabel='Absolute Differnece (AD)':'Double Difference (DD)' -label1='Iteration Number' -label2='Normalized Data Misfit' -out=#{outdir}/misfit.pdf & "

opts = "-n1=401,401  -size1=4 -size2=3 -linewidth=1.5,1.5 -linecolor=b,r -label1='Trace Number' -tick1d=100 -mtick1=9 -label2='Traveltime (s)' -reverse2=y -x2beg=0 -x2end=3.0 -tick2d=0.5 -mtick2=4"
system "x_showgraph -in=data/shot_3_traveltime_p.bin,test_ad/iteration_0/synthetic/shot_3_traveltime_p.bin #{opts} -plotlabel='Ground Truth':'Initial Model' -out=#{outdir}/time_init.pdf & "
system "x_showgraph -in=data/shot_3_traveltime_p.bin,test_ad/iteration_50/synthetic/shot_3_traveltime_p.bin #{opts} -plotlabel='Ground Truth':'Inverted Model - AD' -out=#{outdir}/time_fatt_ad.pdf & "
system "x_showgraph -in=data/shot_3_traveltime_p.bin,test_dd/iteration_50/synthetic/shot_3_traveltime_p.bin #{opts} -plotlabel='Ground Truth':'Inverted Model - DD' -out=#{outdir}/time_fatt_dd.pdf & "
