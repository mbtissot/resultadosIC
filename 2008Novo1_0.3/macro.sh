#! /usr/bin/bash

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 8 execução"
echo "$dt - Começando 8 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe4000.wt
cp IL.wt IL4000.wt
cp IS.wt IS4000.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe4000.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 9 execução"
echo "$dt - Começando 9 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe4500.wt
cp IL.wt IL4500.wt
cp IS.wt IS4500.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe4500.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 10 execução"
echo "$dt - Começando 10 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe5000.wt
cp IL.wt IL5000.wt
cp IS.wt IS5000.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe5000.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

#paplay /usr/share/sounds/LinuxMint/stereo/desktop-login.ogg
sleep 600
shutdown now