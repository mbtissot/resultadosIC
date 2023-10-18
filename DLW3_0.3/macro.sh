#! /usr/bin/bash


#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 17 execução"
echo "$dt - Começando 17 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1700.wt
cp IL.wt IL1700.wt
cp IS.wt IS1700.wt
cp IT.wt IT1700.wt
cp DDzzduz.wt DDzzduz1700.wt
cp Out.wt Out1700.wt
cp Ratios.wt Ratios1700.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe1700.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

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

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 21 execução"
echo "$dt - Começando 21 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe2100.wt
cp IL.wt IL2100.wt
cp IS.wt IS2100.wt
cp IT.wt IT2100.wt
cp DDzzduz.wt DDzzduz2100.wt
cp Out.wt Out2100.wt
cp Ratios.wt Ratios2100.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe2100.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 22 execução"
echo "$dt - Começando 22 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe2200.wt
cp IL.wt IL2200.wt
cp IS.wt IS2200.wt
cp IT.wt IT2200.wt
cp DDzzduz.wt DDzzduz2200.wt
cp Out.wt Out2200.wt
cp Ratios.wt Ratios2200.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe2200.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 23 execução"
echo "$dt - Começando 23 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe2300.wt
cp IL.wt IL2300.wt
cp IS.wt IS2300.wt
cp IT.wt IT2300.wt
cp DDzzduz.wt DDzzduz2300.wt
cp Out.wt Out2300.wt
cp Ratios.wt Ratios2300.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe2300.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 24 execução"
echo "$dt - Começando 24 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe2400.wt
cp IL.wt IL2400.wt
cp IS.wt IS2400.wt
cp IT.wt IT2400.wt
cp DDzzduz.wt DDzzduz2400.wt
cp Out.wt Out2400.wt
cp Ratios.wt Ratios2400.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe2400.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 25 execução"
echo "$dt - Começando 25 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe2500.wt
cp IL.wt IL2500.wt
cp IS.wt IS2500.wt
cp IT.wt IT2500.wt
cp DDzzduz.wt DDzzduz2500.wt
cp Out.wt Out2500.wt
cp Ratios.wt Ratios2500.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe2500.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------
sleep 1200
shutdown now