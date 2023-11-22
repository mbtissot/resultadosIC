#!/bin/bash

regexFe="Fe+([0-9]?)+\.+wt$"
regexIL="IL+([0-9]?)+\.+wt$"
regexIS="IS+([0-9]?)+\.+wt$"
regexIT="IT+([0-9]?)+\.+wt$"

filename=$(basename $1)
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

if [[ $filename =~ $regexFe ]] || [[ $filename =~ $regexIL ]] || [[ $filename =~ $regexIS ]] || [[ $filename =~ $regexIT ]]
then
    #echo "Chamei plotador"
	~/ShellScripts/plotador.sh $filename $var $filename_noext
else
	xed $filename
fi
