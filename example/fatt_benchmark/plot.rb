
n1 = 201
n2 = 201
outdir = "./"

opts = " -n1=#{n1} -d1=0.01 -d2=0.01 -label1='Depth (km)' -label2='Horizontal Position (km)' -size1=3.5 -size2=3.5 -tick1d=0.5 -mtick1=4 -tick2d=0.5 -mtick2=4 -backlegend=y -unit='Adjoint-State Gradient' -lmtick=4 -curve=srcr.txt,recr.txt -curveselect=3,1 -curvestyle=scatter*,scatterv -curvefacecolor=k,k -curvesize=80,30 "

system "x_showcontour -background=test_low/iteration_1/model/grad_vp.bin #{opts} -backcolor=bwr -contours=1500 -backcmin=-0.1 -backcmax=0.1 -clabelsize=0 -ld=0.05 -in=model/vp_low.bin -text='$v = 1000$ m/s' -textloc=0.3,1 -contourwidth=2 -out=#{outdir}/grad_low.pdf &"

system "x_showcontour -background=test_high/iteration_1/model/grad_vp.bin #{opts} -backcolor=bwr -contours=2000 -backcmin=-0.2 -backcmax=0.2 -clabelsize=0 -ld=0.05 -in=model/vp_high.bin -text='$v = 3000$ m/s' -textloc=0.3,1 -contourwidth=2 -out=#{outdir}/grad_high.pdf &"
