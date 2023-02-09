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
mv Out.wt Ini.wt

if grep -q "NaN" Fe2500.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

paplay /usr/share/sounds/LinuxMint/stereo/desktop-login.ogg
#shutdown now