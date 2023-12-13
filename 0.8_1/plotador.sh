#!/bin/bash
filename=$1
var=$2
saveAs=$3
#echo $var
gnuplot -persist -c /home/matheus/ShellScripts/script2.gp $filename $var $saveAs
