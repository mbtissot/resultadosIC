if (ARG2 eq 'Fe'){
    set xlabel 'u_{/Symbol \136}';
    set ylabel 'u_{/Symbol \174\174}';
} else {
    set xlabel 'q_{/Symbol \136}';
    set ylabel 'q_{/Symbol \174\174}';
}
set contour
set title ARG1
sp ARG1 u 1:2:3 w lines
pause mouse close
