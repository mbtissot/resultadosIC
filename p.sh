#!/bin/bash

regexFe="^Fe+([0-9]?)+\.+wt$"
regexIL="^IL+([0-9]?)+\.+wt$"
regexIS="^IS+([0-9]?)+\.+wt$"
regexIT="^IT+([0-9]?)+\.+wt$"
regexIIT="^IIT+([0-9]?)+\.+wt$"

filename=$1
filename_bom=$(basename $1)
filename_noext=${filename_bom::-3}
var="ILISIT"


if [[ $filename_bom =~ $regexFe ]] 
then
#    echo "Deu Fe"
    var="Fe"
fi

if [[ $filename_bom =~ $regexIL ]] 
then
#    echo "Deu IL"
    var="IL"
fi

if [[ $filename_bom =~ $regexIS ]] 
then
#    echo "Deu IS"
    var="IS"
fi

if [[ $filename_bom =~ $regexIT ]] 
then
#    echo "Deu IT"
    var="IT"
fi

if [[ $filename_bom =~ $regexIIT ]]
then
    echo $filename
    gnuplot -persist -c plot2.gp $filename $var $filename_noext
    return 0
fi

echo $filename


gnuplot -persist -c plot.gp $filename $var $filename_noext
