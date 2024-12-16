
# To keep codes clean, I have commented out the nearest-point methods, and the 
# codes just use interpolation schemes
# Therefore to generate the files in the below with names _neareast, you 
# have to recompile the codes. 

system "x_diff data_gradient_nearest/shot_1_traveltime_p.bin tref_gradient.bin >diff_nearest_gradient.bin "
system "x_diff data_gradient_interp/shot_1_traveltime_p.bin tref_gradient.bin >diff_interp_gradient.bin "

n1 = 500

a = [0]
b = [499]
e2 = [0.25]
dd2 = [0.1]

for i in 1..1

	opts = " -n1=#{n1},#{n1},#{n1} -x1beg=#{a[i - 1]} -x1end=#{b[i - 1]} -linewidth=1,0.5,0.5 -tick1d=100 -mtick1=9 -label1='Receiver Index' -label2='Traveltime (s)' -size1=6 -size2=2 -x2beg=0 "

	system "x_showgraph -in=tref_gradient.bin #{opts} -n1=#{n1} -linecolor=k -x2end=#{e2[i - 1]} -tick2d=#{dd2[i - 1]} -mtick2=4 -out=compare_interp_#{i}_gradient.pdf &"

	opts = " -n1=#{n1},#{n1} -linewidth=1,1 -x1beg=#{a[i - 1]} -x1end=#{b[i - 1]} -tick1d=100 -mtick1=9 -label1='Receiver Index' -label2='Traveltime Error (s)' -size1=6 -size2=2 "

	system "x_showgraph -in=diff_nearest_gradient.bin,diff_interp_gradient.bin #{opts} -linecolor=b,r -plotlabel='Nearest Grid Point':'Interpolation' -plotlabelloc=lower_center -plotlabelcol=2 -plotlabelsize=12 -x2beg=-0.006 -x2end=0.006 -tick2d=0.003 -mtick2=4 -out=compare_interp_error_#{i}_gradient.pdf &"

end

system "x_showmatrix -in=model/v_gradient.bin -n1=101 -size1=3.5 -size2=3.5 -d1=10 -d2=10 -label1='Depth (m)' -label2='Horizontal Position (m)' -curve=recr_gradient.txt,srcr_gradient.txt -curvestyle=scatterv,scatter* -curvesize=10,100 -curvefacecolor=k,r -curveedgecolor=none,k -curveselect=2,1 -mtick1=4 -mtick2=4 -out=recr_gradient.pdf -legend=y -unit='Velocity (m/s)' -ltickbeg=2500 -ld=250 -color=rainbowcmyk -lmtick=4 &"


system "x_diff data3_gradient_nearest/shot_1_traveltime_p.bin tref3_gradient.bin >diff3_nearest_gradient.bin "
system "x_diff data3_gradient_interp/shot_1_traveltime_p.bin tref3_gradient.bin >diff3_interp_gradient.bin "

n1 = 500

a = [0]
b = [499]
e2 = [0.3]
dd2 = [0.1]

for i in 1..1

	opts = " -n1=#{n1},#{n1},#{n1} -x1beg=#{a[i - 1]} -x1end=#{b[i - 1]} -linewidth=1,0.5,0.5 -tick1d=100 -mtick1=9 -label1='Receiver Index' -label2='Traveltime (s)' -size1=6 -size2=2 -x2beg=0 "

	system "x_showgraph -in=tref3_gradient.bin #{opts} -n1=#{n1} -linecolor=k -x2end=#{e2[i - 1]} -tick2d=#{dd2[i - 1]} -mtick2=4 -out=compare3_interp_#{i}_gradient.pdf &"

	opts = " -n1=#{n1},#{n1} -linewidth=1,1 -x1beg=#{a[i - 1]} -x1end=#{b[i - 1]} -tick1d=100 -mtick1=9 -label1='Receiver Index' -label2='Traveltime Error (s)' -size1=6 -size2=2 "

	system "x_showgraph -in=diff3_nearest_gradient.bin,diff3_interp_gradient.bin #{opts} -linecolor=b,r -plotlabel='Nearest Grid Point':'Interpolation' -plotlabelloc=lower_center -plotlabelcol=2 -plotlabelsize=12 -x2beg=-0.006 -x2end=0.006 -tick2d=0.003 -mtick2=4 -out=compare3_interp_error_#{i}_gradient.pdf &"

end

system "x_showcolorbar -cmin=2.2930000E+03 -cmax=4.0730000E+03 -unit='Velocity (m/s)' -ltickbeg=2000 -lmtick=4 -lloc=bottom -ld=500 -out=v3bar.png -lwidth=3.3  &"




system "x_diff data_homo_nearest/shot_1_traveltime_p.bin tref_homo.bin >diff_nearest_homo.bin "
system "x_diff data_homo_interp/shot_1_traveltime_p.bin tref_homo.bin >diff_interp_homo.bin "

n1 = 1000

a = [0, 600]
b = [599, 999]
e2 = [0.4, 0.7]
dd2 = [0.1, 0.2]

for i in 1..2

	opts = " -n1=#{n1},#{n1},#{n1} -x1beg=#{a[i - 1]} -x1end=#{b[i - 1]} -linewidth=1,0.5,0.5 -tick1d=100 -mtick1=9 -label1='Receiver Index' -label2='Traveltime (s)' -size1=6 -size2=2 -x2beg=0 "

	system "x_showgraph -in=tref_homo.bin #{opts} -n1=#{n1} -linecolor=k -x2end=#{e2[i - 1]} -tick2d=#{dd2[i - 1]} -mtick2=4 -out=compare_interp_#{i}_homo.pdf &"

	opts = " -n1=#{n1},#{n1} -linewidth=1,1 -x1beg=#{a[i - 1]} -x1end=#{b[i - 1]} -tick1d=100 -mtick1=9 -label1='Receiver Index' -label2='Traveltime Error (s)' -size1=6 -size2=2 "

	system "x_showgraph -in=diff_nearest_homo.bin,diff_interp_homo.bin #{opts} -linecolor=b,r -plotlabel='Nearest Grid Point':'Interpolation' -plotlabelloc=lower_center -plotlabelcol=2 -plotlabelsize=12 -x2beg=-0.005 -x2end=0.01 -tick2d=0.005 -mtick2=4 -out=compare_interp_error_#{i}_homo.pdf &"

end

system "x_showmatrix -in=model/v_homo.bin -n1=101 -size1=3.5 -size2=3.5 -d1=10 -d2=10 -label1='Depth (m)' -label2='Horizontal Position (m)' -curve=recr_homo.txt,srcr_homo.txt -curvestyle=scatterv,scatter* -curvesize=10,100 -curvefacecolor=yellow,r -curveedgecolor=none,none -curveselect=2,1 -mtick1=4 -mtick2=4 -out=recr_homo.pdf &"

