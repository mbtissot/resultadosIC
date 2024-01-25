if (ARG2 eq 'Fe'){
    set title 'Fe'
    set xlabel 'u_{/Symbol \136}' font ",14";
    set ylabel 'u_{/Symbol \174\174}' font ",14";
    set xrange [-12:12]
    set yrange [-12:12]
    set zrange [1e-15:1]
    set contour
    set cntrparam levels discrete 1e-2,1e-4,1e-6,1e-8
}else{
    if(ARG2 eq 'IL'){
        set title 'IL'
        set xlabel 'q_{/Symbol \136}' font ",14";
        set ylabel 'q_{/Symbol \174\174}' font ",14";
        set xrange [-0.6:0.6]
        set yrange [-0.6:0.6]
        set contour
        set cntrparam levels discrete 1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2
    }else{
        if(ARG2 eq 'IS'){
            set title 'IS'
            set xlabel 'q_{/Symbol \136}' font ",14";
            set ylabel 'q_{/Symbol \174\174}' font ",14";
            set xrange [-0.6:0.6]
            set yrange [-0.6:0.6]
            set contour
            set cntrparam levels discrete 1e-5,2.5e-5,5e-5,7.5e-5,1e-4
        } else {
            set title 'IT'
            set xlabel 'q_{/Symbol \136}' font ",14";
            set ylabel 'q_{/Symbol \174\174}' font ",14";
            set xrange [-0.6:0.6]
            set yrange [-0.6:0.6]
            set contour
            set cntrparam levels discrete 1e-6,9e-7,8e-7,7e-7,6e-7,5e-7,4e-7,3e-7,2e-7,1e-7
        }
    }
}

set logscale z
set view 55,65
set ztics font ",8"
set format z "%g"
#unset key
sp ARG1 u 1:2:3 w lines lt rgb "#1735FF"
pause mouse close
