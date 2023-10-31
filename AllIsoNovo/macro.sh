#! /usr/bin/bash


#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 18 execução"
echo "$dt - Começando 18 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1800.wt
cp IL.wt IL1800.wt
cp IS.wt IS1800.wt
cp IT.wt IT1800.wt
cp DDzzduz.wt DDzzduz1800.wt
cp Out.wt Out1800.wt
cp Ratios.wt Ratios1800.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe1800.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 19 execução"
echo "$dt - Começando 19 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1900.wt
cp IL.wt IL1900.wt
cp IS.wt IS1900.wt
cp IT.wt IT1900.wt
cp DDzzduz.wt DDzzduz1900.wt
cp Out.wt Out1900.wt
cp Ratios.wt Ratios1900.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe1900.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 20 execução"
echo "$dt - Começando 20 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe2000.wt
cp IL.wt IL2000.wt
cp IS.wt IS2000.wt
cp IT.wt IT2000.wt
cp DDzzduz.wt DDzzduz2000.wt
cp Out.wt Out2000.wt
cp Ratios.wt Ratios2000.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe2000.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi