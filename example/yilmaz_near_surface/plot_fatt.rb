
iter1 = 90
iter2 = 49

opts = "-n1=21 -label1='Depth (m)' -label2='Horizontal Position (m)' -legend=1 -legendloc=bottom -unit='P-wave Velocity (m/s)' -cmin=350 -cmax=4500 -interp=Gaussian -size1=2.5 -size2=6 -lmtick=9 -ld=1000 -tick1d=5 -mtick1=4 -tick2d=10 -mtick2=9 "
system "mkdir -p result"
system "x_showmatrix -in=vp.bin #{opts} -out=./result/vp_oz.pdf &"
system "x_showmatrix -in=vpgrad.bin #{opts} -out=./result/vp_init.pdf &"
system "x_showmatrix -in=fatt/iteration_#{iter1}/model/updated_vp.bin #{opts} -out=./result/vp_fatt.pdf &"
system "x_showmatrix -in=fatt_ozvp/iteration_#{iter2}/model/updated_vp.bin #{opts} -out=./result/vp_fatt_ozvp.pdf &"

opts = " -n1=48,48 -x2beg=0 -x2end=0.05 -tick2d=0.01 -mtick2=9 -label2='Time (s)' -label1='Trace Number' -tick1d=5 -mtick1=4 -linewidth=2,2 -linecolor=b,r -reverse2=1 -size1=4 -size2=2.5 -plotlabelloc=lower_right " #-marker=o,x -markersize=5,5 -linestyle=solid,solid
for i in 1..49
	system "x_showgraph -in=pick/shot_#{i}_traveltime_p.bin,fatt_ozvp/iteration_0/synthetic/shot_#{i}_traveltime_p.bin -plotlabel='Observed':'Oz Model' #{opts} -out=result/time_oz_#{i}.pdf "
	system "x_showgraph -in=pick/shot_#{i}_traveltime_p.bin,fatt/iteration_#{iter1}/synthetic/shot_#{i}_traveltime_p.bin -plotlabel='Observed':'FATT' #{opts} -out=result/time_fatt_#{i}.pdf "
end
