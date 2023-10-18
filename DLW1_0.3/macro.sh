#! /usr/bin/bash


echo "Tempos desta run: " > Tempo.txt

dt=$(date '+%d/%m/%Y %H:%M:%S');

echo "$dt - Copiando Ini para Ini_backup"
cp Ini.wt Ini_backup.wt

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando compilação"
echo "$dt - Começando compilação" >> Tempo.txt
gfortran wt_lst_aniso.f90

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 1 execução"
echo "$dt - Começando 1 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0100.wt
cp IL.wt IL0100.wt
cp IS.wt IS0100.wt
cp IT.wt IT0100.wt
cp DDzzduz.wt DDzzduz0100.wt
cp Out.wt Out0100.wt
cp Ratios.wt Ratios0100.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe0100.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 2 execução"
echo "$dt - Começando 2 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0200.wt
cp IL.wt IL0200.wt
cp IS.wt IS0200.wt
cp IT.wt IT0200.wt
cp DDzzduz.wt DDzzduz0200.wt
cp Out.wt Out0200.wt
cp Ratios.wt Ratios0200.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe0200.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 3 execução"
echo "$dt - Começando 3 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0300.wt
cp IL.wt IL0300.wt
cp IS.wt IS0300.wt
cp IT.wt IT0300.wt
cp DDzzduz.wt DDzzduz0300.wt
cp Out.wt Out0300.wt
cp Ratios.wt Ratios0300.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe0300.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 4 execução"
echo "$dt - Começando 4 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0400.wt
cp IL.wt IL0400.wt
cp IS.wt IS0400.wt
cp IT.wt IT0400.wt
cp DDzzduz.wt DDzzduz0400.wt
cp Out.wt Out0400.wt
cp Ratios.wt Ratios0400.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe0400.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 5 execução"
echo "$dt - Começando 5 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0500.wt
cp IL.wt IL0500.wt
cp IS.wt IS0500.wt
cp IT.wt IT0500.wt
cp DDzzduz.wt DDzzduz0500.wt
cp Out.wt Out0500.wt
cp Ratios.wt Ratios0500.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe0500.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 6 execução"
echo "$dt - Começando 6 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0600.wt
cp IL.wt IL0600.wt
cp IS.wt IS0600.wt
cp IT.wt IT0600.wt
cp DDzzduz.wt DDzzduz0600.wt
cp Out.wt Out0600.wt
cp Ratios.wt Ratios0600.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe0600.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 7 execução"
echo "$dt - Começando 7 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0700.wt
cp IL.wt IL0700.wt
cp IS.wt IS0700.wt
cp IT.wt IT0700.wt
cp DDzzduz.wt DDzzduz0700.wt
cp Out.wt Out0700.wt
cp Ratios.wt Ratios0700.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe0700.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 8 execução"
echo "$dt - Começando 8 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0800.wt
cp IL.wt IL0800.wt
cp IS.wt IS0800.wt
cp IT.wt IT0800.wt
cp DDzzduz.wt DDzzduz0800.wt
cp Out.wt Out0800.wt
cp Ratios.wt Ratios0800.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe0800.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 9 execução"
echo "$dt - Começando 9 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0900.wt
cp IL.wt IL0900.wt
cp IS.wt IS0900.wt
cp IT.wt IT0900.wt
cp DDzzduz.wt DDzzduz0900.wt
cp Out.wt Out0900.wt
cp Ratios.wt Ratios0900.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe0900.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 10 execução"
echo "$dt - Começando 10 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1000.wt
cp IL.wt IL1000.wt
cp IS.wt IS1000.wt
cp IT.wt IT1000.wt
cp DDzzduz.wt DDzzduz1000.wt
cp Out.wt Out1000.wt
cp Ratios.wt Ratios1000.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe1000.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 11 execução"
echo "$dt - Começando 11 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1100.wt
cp IL.wt IL1100.wt
cp IS.wt IS1100.wt
cp IT.wt IT1100.wt
cp DDzzduz.wt DDzzduz1100.wt
cp Out.wt Out1100.wt
cp Ratios.wt Ratios1100.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe1100.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 12 execução"
echo "$dt - Começando 12 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1200.wt
cp IL.wt IL1200.wt
cp IS.wt IS1200.wt
cp IT.wt IT1200.wt
cp DDzzduz.wt DDzzduz1200.wt
cp Out.wt Out1200.wt
cp Ratios.wt Ratios1200.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe1200.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 13 execução"
echo "$dt - Começando 13 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1300.wt
cp IL.wt IL1300.wt
cp IS.wt IS1300.wt
cp IT.wt IT1300.wt
cp DDzzduz.wt DDzzduz1300.wt
cp Out.wt Out1300.wt
cp Ratios.wt Ratios1300.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe1300.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 14 execução"
echo "$dt - Começando 14 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1400.wt
cp IL.wt IL1400.wt
cp IS.wt IS1400.wt
cp IT.wt IT1400.wt
cp DDzzduz.wt DDzzduz1400.wt
cp Out.wt Out1400.wt
cp Ratios.wt Ratios1400.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe1400.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 15 execução"
echo "$dt - Começando 15 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1500.wt
cp IL.wt IL1500.wt
cp IS.wt IS1500.wt
cp IT.wt IT1500.wt
cp DDzzduz.wt DDzzduz1500.wt
cp Out.wt Out1500.wt
cp Ratios.wt Ratios1500.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe1500.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 16 execução"
echo "$dt - Começando 16 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1600.wt
cp IL.wt IL1600.wt
cp IS.wt IS1600.wt
cp IT.wt IT1600.wt
cp DDzzduz.wt DDzzduz1600.wt
cp Out.wt Out1600.wt
cp Ratios.wt Ratios1600.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe1600.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

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

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 26 execução"
echo "$dt - Começando 26 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe2600.wt
cp IL.wt IL2600.wt
cp IS.wt IS2600.wt
cp IT.wt IT2600.wt
cp DDzzduz.wt DDzzduz2600.wt
cp Out.wt Out2600.wt
cp Ratios.wt Ratios2600.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe2600.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 27 execução"
echo "$dt - Começando 27 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe2700.wt
cp IL.wt IL2700.wt
cp IS.wt IS2700.wt
cp IT.wt IT2700.wt
cp DDzzduz.wt DDzzduz2700.wt
cp Out.wt Out2700.wt
cp Ratios.wt Ratios2700.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe2700.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 28 execução"
echo "$dt - Começando 28 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe2800.wt
cp IL.wt IL2800.wt
cp IS.wt IS2800.wt
cp IT.wt IT2800.wt
cp DDzzduz.wt DDzzduz2800.wt
cp Out.wt Out2800.wt
cp Ratios.wt Ratios2800.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe2800.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 29 execução"
echo "$dt - Começando 29 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe2900.wt
cp IL.wt IL2900.wt
cp IS.wt IS2900.wt
cp IT.wt IT2900.wt
cp DDzzduz.wt DDzzduz2900.wt
cp Out.wt Out2900.wt
cp Ratios.wt Ratios2900.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe2900.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 30 execução"
echo "$dt - Começando 30 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe3000.wt
cp IL.wt IL3000.wt
cp IS.wt IS3000.wt
cp IT.wt IT3000.wt
cp DDzzduz.wt DDzzduz3000.wt
cp Out.wt Out3000.wt
cp Ratios.wt Ratios3000.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe3000.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 31 execução"
echo "$dt - Começando 31 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe3100.wt
cp IL.wt IL3100.wt
cp IS.wt IS3100.wt
cp IT.wt IT3100.wt
cp DDzzduz.wt DDzzduz3100.wt
cp Out.wt Out3100.wt
cp Ratios.wt Ratios3100.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe3100.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 32 execução"
echo "$dt - Começando 32 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe3200.wt
cp IL.wt IL3200.wt
cp IS.wt IS3200.wt
cp IT.wt IT3200.wt
cp DDzzduz.wt DDzzduz3200.wt
cp Out.wt Out3200.wt
cp Ratios.wt Ratios3200.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe3200.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 33 execução"
echo "$dt - Começando 33 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe3300.wt
cp IL.wt IL3300.wt
cp IS.wt IS3300.wt
cp IT.wt IT3300.wt
cp DDzzduz.wt DDzzduz3300.wt
cp Out.wt Out3300.wt
cp Ratios.wt Ratios3300.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe3300.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 34 execução"
echo "$dt - Começando 34 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe3400.wt
cp IL.wt IL3400.wt
cp IS.wt IS3400.wt
cp IT.wt IT3400.wt
cp DDzzduz.wt DDzzduz3400.wt
cp Out.wt Out3400.wt
cp Ratios.wt Ratios3400.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe3400.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 35 execução"
echo "$dt - Começando 35 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe3500.wt
cp IL.wt IL3500.wt
cp IS.wt IS3500.wt
cp IT.wt IT3500.wt
cp DDzzduz.wt DDzzduz3500.wt
cp Out.wt Out3500.wt
cp Ratios.wt Ratios3500.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe3500.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 36 execução"
echo "$dt - Começando 36 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe3600.wt
cp IL.wt IL3600.wt
cp IS.wt IS3600.wt
cp IT.wt IT3600.wt
cp DDzzduz.wt DDzzduz3600.wt
cp Out.wt Out3600.wt
cp Ratios.wt Ratios3600.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe3600.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 37 execução"
echo "$dt - Começando 37 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe3700.wt
cp IL.wt IL3700.wt
cp IS.wt IS3700.wt
cp IT.wt IT3700.wt
cp DDzzduz.wt DDzzduz3700.wt
cp Out.wt Out3700.wt
cp Ratios.wt Ratios3700.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe3700.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 38 execução"
echo "$dt - Começando 38 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe3800.wt
cp IL.wt IL3800.wt
cp IS.wt IS3800.wt
cp IT.wt IT3800.wt
cp DDzzduz.wt DDzzduz3800.wt
cp Out.wt Out3800.wt
cp Ratios.wt Ratios3800.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe3800.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 39 execução"
echo "$dt - Começando 39 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe3900.wt
cp IL.wt IL3900.wt
cp IS.wt IS3900.wt
cp IT.wt IT3900.wt
cp DDzzduz.wt DDzzduz3900.wt
cp Out.wt Out3900.wt
cp Ratios.wt Ratios3900.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe3900.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 40 execução"
echo "$dt - Começando 40 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe4000.wt
cp IL.wt IL4000.wt
cp IS.wt IS4000.wt
cp IT.wt IT4000.wt
cp DDzzduz.wt DDzzduz4000.wt
cp Out.wt Out4000.wt
cp Ratios.wt Ratios4000.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe4000.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 41 execução"
echo "$dt - Começando 41 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe4100.wt
cp IL.wt IL4100.wt
cp IS.wt IS4100.wt
cp IT.wt IT4100.wt
cp DDzzduz.wt DDzzduz4100.wt
cp Out.wt Out4100.wt
cp Ratios.wt Ratios4100.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe4100.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 42 execução"
echo "$dt - Começando 42 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe4200.wt
cp IL.wt IL4200.wt
cp IS.wt IS4200.wt
cp IT.wt IT4200.wt
cp DDzzduz.wt DDzzduz4200.wt
cp Out.wt Out4200.wt
cp Ratios.wt Ratios4200.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe4200.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 43 execução"
echo "$dt - Começando 43 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe4300.wt
cp IL.wt IL4300.wt
cp IS.wt IS4300.wt
cp IT.wt IT4300.wt
cp DDzzduz.wt DDzzduz4300.wt
cp Out.wt Out4300.wt
cp Ratios.wt Ratios4300.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe4300.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 44 execução"
echo "$dt - Começando 44 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe4400.wt
cp IL.wt IL4400.wt
cp IS.wt IS4400.wt
cp IT.wt IT4400.wt
cp DDzzduz.wt DDzzduz4400.wt
cp Out.wt Out4400.wt
cp Ratios.wt Ratios4400.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe4400.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 45 execução"
echo "$dt - Começando 45 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe4500.wt
cp IL.wt IL4500.wt
cp IS.wt IS4500.wt
cp IT.wt IT4500.wt
cp DDzzduz.wt DDzzduz4500.wt
cp Out.wt Out4500.wt
cp Ratios.wt Ratios4500.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe4500.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 46 execução"
echo "$dt - Começando 46 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe4600.wt
cp IL.wt IL4600.wt
cp IS.wt IS4600.wt
cp IT.wt IT4600.wt
cp DDzzduz.wt DDzzduz4600.wt
cp Out.wt Out4600.wt
cp Ratios.wt Ratios4600.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe4600.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 47 execução"
echo "$dt - Começando 47 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe4700.wt
cp IL.wt IL4700.wt
cp IS.wt IS4700.wt
cp IT.wt IT4700.wt
cp DDzzduz.wt DDzzduz4700.wt
cp Out.wt Out4700.wt
cp Ratios.wt Ratios4700.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe4700.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 48 execução"
echo "$dt - Começando 48 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe4800.wt
cp IL.wt IL4800.wt
cp IS.wt IS4800.wt
cp IT.wt IT4800.wt
cp DDzzduz.wt DDzzduz4800.wt
cp Out.wt Out4800.wt
cp Ratios.wt Ratios4800.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe4800.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 49 execução"
echo "$dt - Começando 49 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe4900.wt
cp IL.wt IL4900.wt
cp IS.wt IS4900.wt
cp IT.wt IT4900.wt
cp DDzzduz.wt DDzzduz4900.wt
cp Out.wt Out4900.wt
cp Ratios.wt Ratios4900.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe4900.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 50 execução"
echo "$dt - Começando 50 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe5000.wt
cp IL.wt IL5000.wt
cp IS.wt IS5000.wt
cp IT.wt IT5000.wt
cp DDzzduz.wt DDzzduz5000.wt
cp Out.wt Out5000.wt
cp Ratios.wt Ratios5000.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe5000.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

paplay /usr/share/sounds/LinuxMint/stereo/desktop-login.ogg
#shutdown now