if (ARG2 eq 'Fe'){
    set title 'Fe'
    set xlabel 'u_{/Symbol \136}' font ",14";
    set ylabel 'u_{/Symbol \174\174}' font ",14";
}else{
    if(ARG2 eq 'IL'){
        set title 'IL'
        set xlabel 'q_{/Symbol \136}' font ",14";
        set ylabel 'q_{/Symbol \174\174}' font ",14";
    }else{
        if(ARG2 eq 'IS'){
            set title 'IS'
            set xlabel 'q_{/Symbol \136}' font ",14";
            set ylabel 'q_{/Symbol \174\174}' font ",14";
        } else {
            set title 'IT'
            set xlabel 'q_{/Symbol \136}' font ",14";
            set ylabel 'q_{/Symbol \174\174}' font ",14";
        }
    }
}

set logscale y
set format y "%g"
#unset key
p ARG1 u 1:2 w lines lt rgb "#1735FF"
pause mouse close
