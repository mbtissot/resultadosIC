#!/bin/bash

regexFe="Fe+([0-9]?)+\.+wt$"
regexIL="IL+([0-9]?)+\.+wt$"
regexIS="IS+([0-9]?)+\.+wt$"
regexIT="IT+([0-9]?)+\.+wt$"

filename=$1
filename_noext=${filename::-3}
var="ILISIT"


if [[ $filename =~ $regexFe ]] 
then
#    echo "Deu Fe"
    var="Fe"
fi

if [[ $filename =~ $regexIL ]] 
then
#    echo "Deu IL"
    var="IL"
fi

if [[ $filename =~ $regexIS ]] 
then
#    echo "Deu IS"
    var="IS"
fi

if [[ $filename =~ $regexIT ]] 
then
#    echo "Deu IT"
    var="IT"
fi

gnuplot -persist -c plot.gp $filename $var $filename_noext
