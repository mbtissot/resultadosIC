#! /usr/bin/bash


echo "Tempos desta run: " > Tempo.txt

dt=$(date '+%d/%m/%Y %H:%M:%S');

echo "$dt - Copiando Ini para Ini_backup"
cp Ini.wt Ini_backup.wt

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando compilação"
echo "$dt - Começando compilação" >> Tempo.txt
gfortran wt_2d.f90

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 1 execução"
echo "$dt - Começando 1 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0500.wt
cp IL.wt IL0500.wt
cp IS.wt IS0500.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe0500.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 2 execução"
echo "$dt - Começando 2 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1000.wt
cp IL.wt IL1000.wt
cp IS.wt IS1000.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe1000.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 3 execução"
echo "$dt - Começando 3 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1500.wt
cp IL.wt IL1500.wt
cp IS.wt IS1500.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe1500.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 4 execução"
echo "$dt - Começando 4 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe2000.wt
cp IL.wt IL2000.wt
cp IS.wt IS2000.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe2000.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 5 execução"
echo "$dt - Começando 5 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe2500.wt
cp IL.wt IL2500.wt
cp IS.wt IS2500.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe2500.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 6 execução"
echo "$dt - Começando 6 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe3000.wt
cp IL.wt IL3000.wt
cp IS.wt IS3000.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe3000.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 7 execução"
echo "$dt - Começando 7 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe3500.wt
cp IL.wt IL3500.wt
cp IS.wt IS3500.wt
mv Out.wt Ini.wt

if grep -q "NaN" Fe3500.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

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

paplay /usr/share/sounds/LinuxMint/stereo/desktop-login.ogg
#shutdown now