set term postscript eps enhanced color blacktext "Helvetica" 20


filename1='track_truncated_map_madx.txt'

output1='madx_map_track_x.eps'

set output output1
set ylabel "PX (1 {/Symbol \264 10^{-3})"
set xlabel 'X (1 {/Symbol \264} 10^{-2}m)'
set title ''

set grid
#set logscale y
#set xrange [1:10]
set key left top
plot filename1 using ($2*1e2):($3*1e3) every 10 with points pointtype 1 pointsize 1 t ''



output2='madx_map_track_y.eps'

set output output2
set ylabel "PY (1 {/Symbol \264 10^{-3})"
set xlabel 'Y (1 {/Symbol \264} 10^{-2}m)'
set title ''

set grid
#set logscale y
#set xrange [1:10]
set key left top
plot filename1 using ($4*1e2):($5*1e3) every 10 with points pointtype 1 pointsize 1 t ''




output3='madx_map_track_s.eps'

set output output3
set ylabel "PT (1 {/Symbol \264 10^{-3})"
set xlabel 'T (1 {/Symbol \264} 10^{-4}m)'
set title ''

set grid
#set logscale y
#set xrange [1:10]
set key left top
plot filename1 using ($6*1e4):($7*1e3) every 10 with points pointtype 1 pointsize 1 t ''



filename1='track_symplectic_type_2_madx.txt'

output1='madx_gf2_track_x.eps'

set output output1
set ylabel "PX (1 {/Symbol \264 10^{-3})"
set xlabel 'X (1 {/Symbol \264} 10^{-2}m)'
set title ''

set grid
#set logscale y
#set xrange [1:10]
set key left top
plot filename1 using ($2*1e2):($3*1e3) every 10 with points pointtype 1 pointsize 1 t ''



output2='madx_gf2_track_y.eps'

set output output2
set ylabel "PY (1 {/Symbol \264 10^{-3})"
set xlabel 'Y (1 {/Symbol \264} 10^{-2}m)'
set title ''

set grid
#set logscale y
#set xrange [1:10]
set key left top
plot filename1 using ($4*1e2):($5*1e3) every 10 with points pointtype 1 pointsize 1 t ''




output3='madx_gf2_track_s.eps'

set output output3
set ylabel "PT (1 {/Symbol \264 10^{-3})"
set xlabel 'T (1 {/Symbol \264} 10^{-4}m)'
set title ''

set grid
#set logscale y
#set xrange [1:10]
set key left top
plot filename1 using ($6*1e4):($7*1e3) every 10 with points pointtype 1 pointsize 1 t ''




filename1='track_truncated_map_cosy.txt'

output1='cosy_map_track_x.eps'

set output output1
set ylabel "PX (1 {/Symbol \264 10^{-3})"
set xlabel 'X (1 {/Symbol \264} 10^{-2}m)'
set title ''

set grid
#set logscale y
#set xrange [1:10]
set key left top
plot filename1 using ($2*1e2):($3*1e3) every 10 with points pointtype 1 pointsize 1 t ''



output2='cosy_map_track_y.eps'

set output output2
set ylabel "PY (1 {/Symbol \264 10^{-3})"
set xlabel 'Y (1 {/Symbol \264} 10^{-2}m)'
set title ''

set grid
#set logscale y
#set xrange [1:10]
set key left top
plot filename1 using ($4*1e2):($5*1e3) every 10 with points pointtype 1 pointsize 1 t ''




filename1='track_symplectic_type_1_cosy.txt'

output1='cosy_gf1_track_x.eps'

set output output1
set ylabel "PX (1 {/Symbol \264 10^{-3})"
set xlabel 'X (1 {/Symbol \264} 10^{-2}m)'
set title ''

set grid
#set logscale y
#set xrange [1:10]
set key left top
plot filename1 using ($2*1e2):($3*1e3) every 10 with points pointtype 1 pointsize 1 t ''



output2='cosy_gf1_track_y.eps'

set output output2
set ylabel "PY (1 {/Symbol \264 10^{-3})"
set xlabel 'Y (1 {/Symbol \264} 10^{-2}m)'
set title ''

set grid
#set logscale y
#set xrange [1:10]
set key left top
plot filename1 using ($4*1e2):($5*1e3) every 10 with points pointtype 1 pointsize 1 t ''


filename1='exp03_track_symplectic_type_1_cosy.txt'

output1='ep03_cosy_gf1_track_x.eps'

set output output1
set ylabel "PX (1 {/Symbol \264 10^{-3})"
set xlabel 'X (1 {/Symbol \264} 10^{-2}m)'
set title ''

set grid
#set logscale y
#set xrange [1:10]
set key left top
plot filename1 using ($2*1e2):($3*1e3) every 10 with points pointtype 1 pointsize 1 t ''



output2='ep03_cosy_gf1_track_y.eps'

set output output2
set ylabel "PY (1 {/Symbol \264 10^{-3})"
set xlabel 'Y (1 {/Symbol \264} 10^{-2}m)'
set title ''

set grid
#set logscale y
#set xrange [1:10]
set key left top
plot filename1 using ($4*1e2):($5*1e3) every 10 with points pointtype 1 pointsize 1 t ''


