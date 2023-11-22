#! /usr/bin/bash

args=($@)
num=$#

usage () {
	echo "		USAGE: rmall [-optional]";
	echo "		Where [-optional] is -ism";
	echo "		Meaning that you want to remove all, except Ini, Start and Macro.";
}

if [ $num != 0 ]
then
	echo "Tem um ou mais argumentos"
else
	echo "Tem 0 argumentos"
#	usage;
fi

for (( var=0; var<$num; var++ ))
do
	echo "$var: ${args[var]}"
done
#Ini.wt
#Ini_backup.wt
#macro.sh
#Start*

#rm Fe* IL* IS* IT* status.wt phys* math* comm* sub_prog.mod Warning_qsimp Tempo.txt Ewave.wt a.out Ratios.wt

# if [ ${args[0]} == "abc" ]
# then
# 	echo "Deu";
# else
# 	echo "N deu";
# fi