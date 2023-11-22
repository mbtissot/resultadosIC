if (ARG2 eq 'Fe'){
    set xlabel 'u_{/Symbol \136}' font ",22";
    set ylabel 'u_{/Symbol \174\174}' font ",22";
    set xrange [-12:12]
    set yrange [-12:12]
    set zrange [1e-15:1]
    set cntrparam levels discrete 1e-2,1e-4,1e-6,1e-8
} else {
    set xlabel 'q_{/Symbol \136}' font ",22";
    set ylabel 'q_{/Symbol \174\174}' font ",22";
    set xrange [-0.6:0.6]
    set yrange [-0.6:0.6]
    set cntrparam levels discrete 1e-5,2.5e-5,5e-5,7.5e-5,1e-4
}

set term postscript eps color blacktext "Helvetica" 14
set output ARG3.'.eps'
set contour
set logscale z
set view 55,65
set ztics font ",8"
set format z "%g"
unset key
sp ARG1 u 1:2:3 w lines lt rgb "#1735FF"
pause mouse close
