#! /usr/bin/bash


echo "Tempos desta run: " > Tempo.txt

dt=$(date '+%d/%m/%Y %H:%M:%S');

echo "$dt - Copiando Ini para Ini_backup e criando Files"
cp Ini.wt Ini_backup.wt
mkdir Files

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando compilação"
echo "$dt - Começando compilação" >> Tempo.txt
gfortran wt_lst_aniso.f90

#---------------------------------------------
