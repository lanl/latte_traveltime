
n1 = 201
n2 = 401
s = 2
outdir = "./eikonal"

# model and reflectors
opts = " -n1=#{n1} -d1=0.01 -d2=0.01 -label1='Depth (km)' -label2='Horizontal Position (km)' -size1=2.5 -size2=5 -tick1d=1 -mtick1=9 -tick2d=1 -mtick2=9 -legend=y -lmtick=4 -color=ncar "

system "mkdir -p #{outdir}"

system "x_showmatrix -in=model/vp.bin #{opts} -cmin=700 -cmax=3000 -ld=500 -ltickbeg=500 -color=rainbowcmyk -unit='P-wave Velocity (m/s)' -out=#{outdir}/vp.pdf &"
system "x_showmatrix -in=model/vs.bin #{opts} -cmin=400 -cmax=1700 -ld=500 -ltickbeg=500 -color=rainbowcmyk -unit='S-wave Velocity (m/s)' -out=#{outdir}/vs.pdf &"
system "x_showmatrix -in=model/refl.bin #{opts} -ld=1 -lmtick=0 -ncolor=3 -unit='Reflector Index' -out=#{outdir}/refl.pdf &"


# data
opts = " -n1=#{n2},#{n2},#{n2},#{n2},#{n2} -linecolor=k,b,r,b,r -reverse2=y -linestyle=solid,solid,solid,dashed,dashed -linewidth=2,2,2,2,2 -label1='Trace Number' -label2='Traveltime (s)' -size1=4 -size2=3 -tick1d=100 -mtick1=9 -tickbottom=n -ticktop=y -label1loc=top -label2pad=6 -label1pad=6 "

system "x_select <data_elastic_refl_p/shot_#{s}_traveltime_p.bin n1=#{n2} s2=1,3 >pp.bin "
system "x_select <data_elastic_refl_p/shot_#{s}_traveltime_s.bin n1=#{n2} s2=2,3 >ps.bin "
system "cat pp.bin ps.bin >p.bin "
system "x_showgraph -in=p.bin #{opts} -plotlabel='P':'PP1':'PP2':'PS1':'PS2' -x2beg=0 -x2end=4 -tick2d=1 -mtick2=9 -out=#{outdir}/data_p.pdf &"

system "x_select <data_elastic_refl_s/shot_#{s}_traveltime_s.bin n1=#{n2} s2=1,3 >ss.bin "
system "x_select <data_elastic_refl_s/shot_#{s}_traveltime_p.bin n1=#{n2} s2=2,3 >sp.bin "
system "cat ss.bin sp.bin >s.bin "
system "x_showgraph -in=s.bin #{opts} -plotlabel='S':'SS1':'SS2':'SP1':'SP2' -x2beg=0 -x2end=5 -tick2d=1 -mtick2=9 -out=#{outdir}/data_s.pdf &"


# p incident
opts = opts + " -unit='Time (s)' -ld=0.5 -cmin=0 -cmax=3 -lmtick=4 -contourlevel=0.2 -clabelsize=0 -contourfill=y -color=rainbowcmyk "
system "x_select <./snapshot_elastic_refl_p/shot_#{s}_traveltime_p.bin n1=#{n1} s2=1,#{n2} >p_p.bin; \
    x_showcontour -in=p_p.bin #{opts} -out=#{outdir}/elastic_p_p.pdf & "

opts = opts + " -unit='Time (s)' -ld=0.5 -cmin=0 -cmax=3 -lmtick=4 "
system "x_select <./snapshot_elastic_refl_p/shot_#{s}_traveltime_p.bin n1=#{n1} s2=#{n2 + 1},#{2*n2} >p_pp1.bin; \
    x_showcontour -in=p_pp1.bin #{opts} -out=#{outdir}/elastic_p_pp1.pdf & "
system "x_select <./snapshot_elastic_refl_p/shot_#{s}_traveltime_p.bin n1=#{n1} s2=#{2*n2 + 1},#{3 * n2} >p_pp2.bin; \
    x_showcontour -in=p_pp2.bin #{opts} -out=#{outdir}/elastic_p_pp2.pdf & "

opts = opts + " -unit='Time (s)' -ld=0.5 -cmin=0 -cmax=4 -lmtick=4 "
system "x_select <./snapshot_elastic_refl_p/shot_#{s}_traveltime_s.bin n1=#{n1} s2=#{n2 + 1},#{2*n2} >p_ps1.bin; \
    x_showcontour -in=p_ps1.bin #{opts} -out=#{outdir}/elastic_p_ps1.pdf & "
system "x_select <./snapshot_elastic_refl_p/shot_#{s}_traveltime_s.bin n1=#{n1} s2=#{2*n2 + 1},#{3 * n2} >p_ps2.bin; \
    x_showcontour -in=p_ps2.bin #{opts} -out=#{outdir}/elastic_p_ps2.pdf & "

# s incident
opts = opts + " -unit='Time (s)' -ld=0.5 -cmin=0 -cmax=5 -lmtick=4 "
system "x_select <./snapshot_elastic_refl_s/shot_#{s}_traveltime_s.bin n1=#{n1} s2=1,#{n2} >s_s.bin; \
    x_showcontour -in=s_s.bin #{opts} -out=#{outdir}/elastic_s_s.pdf & "

opts = opts + " -unit='Time (s)' -ld=0.5 -cmin=0 -cmax=5 -lmtick=4 "
system "x_select <./snapshot_elastic_refl_s/shot_#{s}_traveltime_s.bin n1=#{n1} s2=#{n2 + 1},#{2*n2} >s_ss1.bin; \
    x_showcontour -in=s_ss1.bin #{opts} -out=#{outdir}/elastic_s_ss1.pdf & "
system "x_select <./snapshot_elastic_refl_s/shot_#{s}_traveltime_s.bin n1=#{n1} s2=#{2*n2 + 1},#{3 * n2} >s_ss2.bin; \
    x_showcontour -in=s_ss2.bin #{opts} -out=#{outdir}/elastic_s_ss2.pdf & "

opts = opts + " -unit='Time (s)' -ld=0.5 -cmin=0 -cmax=4 -lmtick=4 "
system "x_select <./snapshot_elastic_refl_s/shot_#{s}_traveltime_p.bin n1=#{n1} s2=#{n2 + 1},#{2*n2} >s_sp1.bin; \
    x_showcontour -in=s_sp1.bin #{opts} -out=#{outdir}/elastic_s_sp1.pdf & "
system "x_select <./snapshot_elastic_refl_s/shot_#{s}_traveltime_p.bin n1=#{n1} s2=#{2*n2 + 1},#{3 * n2} >s_sp2.bin; \
    x_showcontour -in=s_sp2.bin #{opts} -out=#{outdir}/elastic_s_sp2.pdf & "
