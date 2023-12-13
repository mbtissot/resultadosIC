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

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 1 execução"
echo "$dt - Começando 1 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0020.wt
cp IL.wt IL0020.wt
cp IS.wt IS0020.wt
cp IT.wt IT0020.wt
cp DDzzduz.wt DDzzduz0020.wt
cp Out.wt Out0020.wt
cp Ratios.wt Ratios0020.wt
cp energies.wt energies0020.wt
cp Ewave.wt Ewave0020.wt
mv DDzzduz0020.wt Out0020.wt Ratios0020.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0020.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 2 execução"
echo "$dt - Começando 2 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0040.wt
cp IL.wt IL0040.wt
cp IS.wt IS0040.wt
cp IT.wt IT0040.wt
cp DDzzduz.wt DDzzduz0040.wt
cp Out.wt Out0040.wt
cp Ratios.wt Ratios0040.wt
cp energies.wt energies0040.wt
cp Ewave.wt Ewave0040.wt
mv DDzzduz0040.wt Out0040.wt Ratios0040.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0040.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 3 execução"
echo "$dt - Começando 3 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0060.wt
cp IL.wt IL0060.wt
cp IS.wt IS0060.wt
cp IT.wt IT0060.wt
cp DDzzduz.wt DDzzduz0060.wt
cp Out.wt Out0060.wt
cp Ratios.wt Ratios0060.wt
cp energies.wt energies0060.wt
cp Ewave.wt Ewave0060.wt
mv DDzzduz0060.wt Out0060.wt Ratios0060.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0060.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 4 execução"
echo "$dt - Começando 4 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0080.wt
cp IL.wt IL0080.wt
cp IS.wt IS0080.wt
cp IT.wt IT0080.wt
cp DDzzduz.wt DDzzduz0080.wt
cp Out.wt Out0080.wt
cp Ratios.wt Ratios0080.wt
cp energies.wt energies0080.wt
cp Ewave.wt Ewave0080.wt
mv DDzzduz0080.wt Out0080.wt Ratios0080.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0080.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 5 execução"
echo "$dt - Começando 5 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0100.wt
cp IL.wt IL0100.wt
cp IS.wt IS0100.wt
cp IT.wt IT0100.wt
cp DDzzduz.wt DDzzduz0100.wt
cp Out.wt Out0100.wt
cp Ratios.wt Ratios0100.wt
cp energies.wt energies0100.wt
cp Ewave.wt Ewave0100.wt
mv DDzzduz0100.wt Out0100.wt Ratios0100.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0100.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 6 execução"
echo "$dt - Começando 6 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0120.wt
cp IL.wt IL0120.wt
cp IS.wt IS0120.wt
cp IT.wt IT0120.wt
cp DDzzduz.wt DDzzduz0120.wt
cp Out.wt Out0120.wt
cp Ratios.wt Ratios0120.wt
cp energies.wt energies0120.wt
cp Ewave.wt Ewave0120.wt
mv DDzzduz0120.wt Out0120.wt Ratios0120.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0120.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 7 execução"
echo "$dt - Começando 7 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0140.wt
cp IL.wt IL0140.wt
cp IS.wt IS0140.wt
cp IT.wt IT0140.wt
cp DDzzduz.wt DDzzduz0140.wt
cp Out.wt Out0140.wt
cp Ratios.wt Ratios0140.wt
cp energies.wt energies0140.wt
cp Ewave.wt Ewave0140.wt
mv DDzzduz0140.wt Out0140.wt Ratios0140.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0140.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 8 execução"
echo "$dt - Começando 8 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0160.wt
cp IL.wt IL0160.wt
cp IS.wt IS0160.wt
cp IT.wt IT0160.wt
cp DDzzduz.wt DDzzduz0160.wt
cp Out.wt Out0160.wt
cp Ratios.wt Ratios0160.wt
cp energies.wt energies0160.wt
cp Ewave.wt Ewave0160.wt
mv DDzzduz0160.wt Out0160.wt Ratios0160.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0160.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 9 execução"
echo "$dt - Começando 9 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0180.wt
cp IL.wt IL0180.wt
cp IS.wt IS0180.wt
cp IT.wt IT0180.wt
cp DDzzduz.wt DDzzduz0180.wt
cp Out.wt Out0180.wt
cp Ratios.wt Ratios0180.wt
cp energies.wt energies0180.wt
cp Ewave.wt Ewave0180.wt
mv DDzzduz0180.wt Out0180.wt Ratios0180.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0180.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 10 execução"
echo "$dt - Começando 10 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0200.wt
cp IL.wt IL0200.wt
cp IS.wt IS0200.wt
cp IT.wt IT0200.wt
cp DDzzduz.wt DDzzduz0200.wt
cp Out.wt Out0200.wt
cp Ratios.wt Ratios0200.wt
cp energies.wt energies0200.wt
cp Ewave.wt Ewave0200.wt
mv DDzzduz0200.wt Out0200.wt Ratios0200.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0200.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 11 execução"
echo "$dt - Começando 11 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0220.wt
cp IL.wt IL0220.wt
cp IS.wt IS0220.wt
cp IT.wt IT0220.wt
cp DDzzduz.wt DDzzduz0220.wt
cp Out.wt Out0220.wt
cp Ratios.wt Ratios0220.wt
cp energies.wt energies0220.wt
cp Ewave.wt Ewave0220.wt
mv DDzzduz0220.wt Out0220.wt Ratios0220.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0220.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 12 execução"
echo "$dt - Começando 12 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0240.wt
cp IL.wt IL0240.wt
cp IS.wt IS0240.wt
cp IT.wt IT0240.wt
cp DDzzduz.wt DDzzduz0240.wt
cp Out.wt Out0240.wt
cp Ratios.wt Ratios0240.wt
cp energies.wt energies0240.wt
cp Ewave.wt Ewave0240.wt
mv DDzzduz0240.wt Out0240.wt Ratios0240.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0240.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 13 execução"
echo "$dt - Começando 13 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0260.wt
cp IL.wt IL0260.wt
cp IS.wt IS0260.wt
cp IT.wt IT0260.wt
cp DDzzduz.wt DDzzduz0260.wt
cp Out.wt Out0260.wt
cp Ratios.wt Ratios0260.wt
cp energies.wt energies0260.wt
cp Ewave.wt Ewave0260.wt
mv DDzzduz0260.wt Out0260.wt Ratios0260.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0260.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 14 execução"
echo "$dt - Começando 14 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0280.wt
cp IL.wt IL0280.wt
cp IS.wt IS0280.wt
cp IT.wt IT0280.wt
cp DDzzduz.wt DDzzduz0280.wt
cp Out.wt Out0280.wt
cp Ratios.wt Ratios0280.wt
cp energies.wt energies0280.wt
cp Ewave.wt Ewave0280.wt
mv DDzzduz0280.wt Out0280.wt Ratios0280.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0280.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 15 execução"
echo "$dt - Começando 15 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0300.wt
cp IL.wt IL0300.wt
cp IS.wt IS0300.wt
cp IT.wt IT0300.wt
cp DDzzduz.wt DDzzduz0300.wt
cp Out.wt Out0300.wt
cp Ratios.wt Ratios0300.wt
cp energies.wt energies0300.wt
cp Ewave.wt Ewave0300.wt
mv DDzzduz0300.wt Out0300.wt Ratios0300.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0300.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 16 execução"
echo "$dt - Começando 16 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0320.wt
cp IL.wt IL0320.wt
cp IS.wt IS0320.wt
cp IT.wt IT0320.wt
cp DDzzduz.wt DDzzduz0320.wt
cp Out.wt Out0320.wt
cp Ratios.wt Ratios0320.wt
cp energies.wt energies0320.wt
cp Ewave.wt Ewave0320.wt
mv DDzzduz0320.wt Out0320.wt Ratios0320.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0320.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 17 execução"
echo "$dt - Começando 17 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0340.wt
cp IL.wt IL0340.wt
cp IS.wt IS0340.wt
cp IT.wt IT0340.wt
cp DDzzduz.wt DDzzduz0340.wt
cp Out.wt Out0340.wt
cp Ratios.wt Ratios0340.wt
cp energies.wt energies0340.wt
cp Ewave.wt Ewave0340.wt
mv DDzzduz0340.wt Out0340.wt Ratios0340.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0340.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 18 execução"
echo "$dt - Começando 18 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0360.wt
cp IL.wt IL0360.wt
cp IS.wt IS0360.wt
cp IT.wt IT0360.wt
cp DDzzduz.wt DDzzduz0360.wt
cp Out.wt Out0360.wt
cp Ratios.wt Ratios0360.wt
cp energies.wt energies0360.wt
cp Ewave.wt Ewave0360.wt
mv DDzzduz0360.wt Out0360.wt Ratios0360.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0360.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 19 execução"
echo "$dt - Começando 19 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0380.wt
cp IL.wt IL0380.wt
cp IS.wt IS0380.wt
cp IT.wt IT0380.wt
cp DDzzduz.wt DDzzduz0380.wt
cp Out.wt Out0380.wt
cp Ratios.wt Ratios0380.wt
cp energies.wt energies0380.wt
cp Ewave.wt Ewave0380.wt
mv DDzzduz0380.wt Out0380.wt Ratios0380.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0380.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 20 execução"
echo "$dt - Começando 20 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0400.wt
cp IL.wt IL0400.wt
cp IS.wt IS0400.wt
cp IT.wt IT0400.wt
cp DDzzduz.wt DDzzduz0400.wt
cp Out.wt Out0400.wt
cp Ratios.wt Ratios0400.wt
cp energies.wt energies0400.wt
cp Ewave.wt Ewave0400.wt
mv DDzzduz0400.wt Out0400.wt Ratios0400.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0400.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 21 execução"
echo "$dt - Começando 21 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0420.wt
cp IL.wt IL0420.wt
cp IS.wt IS0420.wt
cp IT.wt IT0420.wt
cp DDzzduz.wt DDzzduz0420.wt
cp Out.wt Out0420.wt
cp Ratios.wt Ratios0420.wt
cp energies.wt energies0420.wt
cp Ewave.wt Ewave0420.wt
mv DDzzduz0420.wt Out0420.wt Ratios0420.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0420.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 22 execução"
echo "$dt - Começando 22 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0440.wt
cp IL.wt IL0440.wt
cp IS.wt IS0440.wt
cp IT.wt IT0440.wt
cp DDzzduz.wt DDzzduz0440.wt
cp Out.wt Out0440.wt
cp Ratios.wt Ratios0440.wt
cp energies.wt energies0440.wt
cp Ewave.wt Ewave0440.wt
mv DDzzduz0440.wt Out0440.wt Ratios0440.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0440.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 23 execução"
echo "$dt - Começando 23 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0460.wt
cp IL.wt IL0460.wt
cp IS.wt IS0460.wt
cp IT.wt IT0460.wt
cp DDzzduz.wt DDzzduz0460.wt
cp Out.wt Out0460.wt
cp Ratios.wt Ratios0460.wt
cp energies.wt energies0460.wt
cp Ewave.wt Ewave0460.wt
mv DDzzduz0460.wt Out0460.wt Ratios0460.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0460.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 24 execução"
echo "$dt - Começando 24 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0480.wt
cp IL.wt IL0480.wt
cp IS.wt IS0480.wt
cp IT.wt IT0480.wt
cp DDzzduz.wt DDzzduz0480.wt
cp Out.wt Out0480.wt
cp Ratios.wt Ratios0480.wt
cp energies.wt energies0480.wt
cp Ewave.wt Ewave0480.wt
mv DDzzduz0480.wt Out0480.wt Ratios0480.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0480.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 25 execução"
echo "$dt - Começando 25 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0500.wt
cp IL.wt IL0500.wt
cp IS.wt IS0500.wt
cp IT.wt IT0500.wt
cp DDzzduz.wt DDzzduz0500.wt
cp Out.wt Out0500.wt
cp Ratios.wt Ratios0500.wt
cp energies.wt energies0500.wt
cp Ewave.wt Ewave0500.wt
mv DDzzduz0500.wt Out0500.wt Ratios0500.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0500.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 26 execução"
echo "$dt - Começando 26 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0520.wt
cp IL.wt IL0520.wt
cp IS.wt IS0520.wt
cp IT.wt IT0520.wt
cp DDzzduz.wt DDzzduz0520.wt
cp Out.wt Out0520.wt
cp Ratios.wt Ratios0520.wt
cp energies.wt energies0520.wt
cp Ewave.wt Ewave0520.wt
mv DDzzduz0520.wt Out0520.wt Ratios0520.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0520.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 27 execução"
echo "$dt - Começando 27 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0540.wt
cp IL.wt IL0540.wt
cp IS.wt IS0540.wt
cp IT.wt IT0540.wt
cp DDzzduz.wt DDzzduz0540.wt
cp Out.wt Out0540.wt
cp Ratios.wt Ratios0540.wt
cp energies.wt energies0540.wt
cp Ewave.wt Ewave0540.wt
mv DDzzduz0540.wt Out0540.wt Ratios0540.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0540.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 28 execução"
echo "$dt - Começando 28 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0560.wt
cp IL.wt IL0560.wt
cp IS.wt IS0560.wt
cp IT.wt IT0560.wt
cp DDzzduz.wt DDzzduz0560.wt
cp Out.wt Out0560.wt
cp Ratios.wt Ratios0560.wt
cp energies.wt energies0560.wt
cp Ewave.wt Ewave0560.wt
mv DDzzduz0560.wt Out0560.wt Ratios0560.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0560.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 29 execução"
echo "$dt - Começando 29 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0580.wt
cp IL.wt IL0580.wt
cp IS.wt IS0580.wt
cp IT.wt IT0580.wt
cp DDzzduz.wt DDzzduz0580.wt
cp Out.wt Out0580.wt
cp Ratios.wt Ratios0580.wt
cp energies.wt energies0580.wt
cp Ewave.wt Ewave0580.wt
mv DDzzduz0580.wt Out0580.wt Ratios0580.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0580.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 30 execução"
echo "$dt - Começando 30 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0600.wt
cp IL.wt IL0600.wt
cp IS.wt IS0600.wt
cp IT.wt IT0600.wt
cp DDzzduz.wt DDzzduz0600.wt
cp Out.wt Out0600.wt
cp Ratios.wt Ratios0600.wt
cp energies.wt energies0600.wt
cp Ewave.wt Ewave0600.wt
mv DDzzduz0600.wt Out0600.wt Ratios0600.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0600.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 31 execução"
echo "$dt - Começando 31 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0620.wt
cp IL.wt IL0620.wt
cp IS.wt IS0620.wt
cp IT.wt IT0620.wt
cp DDzzduz.wt DDzzduz0620.wt
cp Out.wt Out0620.wt
cp Ratios.wt Ratios0620.wt
cp energies.wt energies0620.wt
cp Ewave.wt Ewave0620.wt
mv DDzzduz0620.wt Out0620.wt Ratios0620.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0620.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 32 execução"
echo "$dt - Começando 32 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0640.wt
cp IL.wt IL0640.wt
cp IS.wt IS0640.wt
cp IT.wt IT0640.wt
cp DDzzduz.wt DDzzduz0640.wt
cp Out.wt Out0640.wt
cp Ratios.wt Ratios0640.wt
cp energies.wt energies0640.wt
cp Ewave.wt Ewave0640.wt
mv DDzzduz0640.wt Out0640.wt Ratios0640.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0640.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 33 execução"
echo "$dt - Começando 33 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0660.wt
cp IL.wt IL0660.wt
cp IS.wt IS0660.wt
cp IT.wt IT0660.wt
cp DDzzduz.wt DDzzduz0660.wt
cp Out.wt Out0660.wt
cp Ratios.wt Ratios0660.wt
cp energies.wt energies0660.wt
cp Ewave.wt Ewave0660.wt
mv DDzzduz0660.wt Out0660.wt Ratios0660.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0660.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 34 execução"
echo "$dt - Começando 34 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0680.wt
cp IL.wt IL0680.wt
cp IS.wt IS0680.wt
cp IT.wt IT0680.wt
cp DDzzduz.wt DDzzduz0680.wt
cp Out.wt Out0680.wt
cp Ratios.wt Ratios0680.wt
cp energies.wt energies0680.wt
cp Ewave.wt Ewave0680.wt
mv DDzzduz0680.wt Out0680.wt Ratios0680.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0680.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 35 execução"
echo "$dt - Começando 35 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0700.wt
cp IL.wt IL0700.wt
cp IS.wt IS0700.wt
cp IT.wt IT0700.wt
cp DDzzduz.wt DDzzduz0700.wt
cp Out.wt Out0700.wt
cp Ratios.wt Ratios0700.wt
cp energies.wt energies0700.wt
cp Ewave.wt Ewave0700.wt
mv DDzzduz0700.wt Out0700.wt Ratios0700.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0700.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 36 execução"
echo "$dt - Começando 36 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0720.wt
cp IL.wt IL0720.wt
cp IS.wt IS0720.wt
cp IT.wt IT0720.wt
cp DDzzduz.wt DDzzduz0720.wt
cp Out.wt Out0720.wt
cp Ratios.wt Ratios0720.wt
cp energies.wt energies0720.wt
cp Ewave.wt Ewave0720.wt
mv DDzzduz0720.wt Out0720.wt Ratios0720.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0720.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 37 execução"
echo "$dt - Começando 37 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0740.wt
cp IL.wt IL0740.wt
cp IS.wt IS0740.wt
cp IT.wt IT0740.wt
cp DDzzduz.wt DDzzduz0740.wt
cp Out.wt Out0740.wt
cp Ratios.wt Ratios0740.wt
cp energies.wt energies0740.wt
cp Ewave.wt Ewave0740.wt
mv DDzzduz0740.wt Out0740.wt Ratios0740.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0740.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 38 execução"
echo "$dt - Começando 38 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0760.wt
cp IL.wt IL0760.wt
cp IS.wt IS0760.wt
cp IT.wt IT0760.wt
cp DDzzduz.wt DDzzduz0760.wt
cp Out.wt Out0760.wt
cp Ratios.wt Ratios0760.wt
cp energies.wt energies0760.wt
cp Ewave.wt Ewave0760.wt
mv DDzzduz0760.wt Out0760.wt Ratios0760.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0760.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 39 execução"
echo "$dt - Começando 39 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0780.wt
cp IL.wt IL0780.wt
cp IS.wt IS0780.wt
cp IT.wt IT0780.wt
cp DDzzduz.wt DDzzduz0780.wt
cp Out.wt Out0780.wt
cp Ratios.wt Ratios0780.wt
cp energies.wt energies0780.wt
cp Ewave.wt Ewave0780.wt
mv DDzzduz0780.wt Out0780.wt Ratios0780.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0780.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 40 execução"
echo "$dt - Começando 40 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0800.wt
cp IL.wt IL0800.wt
cp IS.wt IS0800.wt
cp IT.wt IT0800.wt
cp DDzzduz.wt DDzzduz0800.wt
cp Out.wt Out0800.wt
cp Ratios.wt Ratios0800.wt
cp energies.wt energies0800.wt
cp Ewave.wt Ewave0800.wt
mv DDzzduz0800.wt Out0800.wt Ratios0800.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0800.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 41 execução"
echo "$dt - Começando 41 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0820.wt
cp IL.wt IL0820.wt
cp IS.wt IS0820.wt
cp IT.wt IT0820.wt
cp DDzzduz.wt DDzzduz0820.wt
cp Out.wt Out0820.wt
cp Ratios.wt Ratios0820.wt
cp energies.wt energies0820.wt
cp Ewave.wt Ewave0820.wt
mv DDzzduz0820.wt Out0820.wt Ratios0820.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0820.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 42 execução"
echo "$dt - Começando 42 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0840.wt
cp IL.wt IL0840.wt
cp IS.wt IS0840.wt
cp IT.wt IT0840.wt
cp DDzzduz.wt DDzzduz0840.wt
cp Out.wt Out0840.wt
cp Ratios.wt Ratios0840.wt
cp energies.wt energies0840.wt
cp Ewave.wt Ewave0840.wt
mv DDzzduz0840.wt Out0840.wt Ratios0840.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0840.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 43 execução"
echo "$dt - Começando 43 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0860.wt
cp IL.wt IL0860.wt
cp IS.wt IS0860.wt
cp IT.wt IT0860.wt
cp DDzzduz.wt DDzzduz0860.wt
cp Out.wt Out0860.wt
cp Ratios.wt Ratios0860.wt
cp energies.wt energies0860.wt
cp Ewave.wt Ewave0860.wt
mv DDzzduz0860.wt Out0860.wt Ratios0860.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0860.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 44 execução"
echo "$dt - Começando 44 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0880.wt
cp IL.wt IL0880.wt
cp IS.wt IS0880.wt
cp IT.wt IT0880.wt
cp DDzzduz.wt DDzzduz0880.wt
cp Out.wt Out0880.wt
cp Ratios.wt Ratios0880.wt
cp energies.wt energies0880.wt
cp Ewave.wt Ewave0880.wt
mv DDzzduz0880.wt Out0880.wt Ratios0880.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0880.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 45 execução"
echo "$dt - Começando 45 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0900.wt
cp IL.wt IL0900.wt
cp IS.wt IS0900.wt
cp IT.wt IT0900.wt
cp DDzzduz.wt DDzzduz0900.wt
cp Out.wt Out0900.wt
cp Ratios.wt Ratios0900.wt
cp energies.wt energies0900.wt
cp Ewave.wt Ewave0900.wt
mv DDzzduz0900.wt Out0900.wt Ratios0900.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0900.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 46 execução"
echo "$dt - Começando 46 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0920.wt
cp IL.wt IL0920.wt
cp IS.wt IS0920.wt
cp IT.wt IT0920.wt
cp DDzzduz.wt DDzzduz0920.wt
cp Out.wt Out0920.wt
cp Ratios.wt Ratios0920.wt
cp energies.wt energies0920.wt
cp Ewave.wt Ewave0920.wt
mv DDzzduz0920.wt Out0920.wt Ratios0920.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0920.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 47 execução"
echo "$dt - Começando 47 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0940.wt
cp IL.wt IL0940.wt
cp IS.wt IS0940.wt
cp IT.wt IT0940.wt
cp DDzzduz.wt DDzzduz0940.wt
cp Out.wt Out0940.wt
cp Ratios.wt Ratios0940.wt
cp energies.wt energies0940.wt
cp Ewave.wt Ewave0940.wt
mv DDzzduz0940.wt Out0940.wt Ratios0940.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0940.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 48 execução"
echo "$dt - Começando 48 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0960.wt
cp IL.wt IL0960.wt
cp IS.wt IS0960.wt
cp IT.wt IT0960.wt
cp DDzzduz.wt DDzzduz0960.wt
cp Out.wt Out0960.wt
cp Ratios.wt Ratios0960.wt
cp energies.wt energies0960.wt
cp Ewave.wt Ewave0960.wt
mv DDzzduz0960.wt Out0960.wt Ratios0960.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0960.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 49 execução"
echo "$dt - Começando 49 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe0980.wt
cp IL.wt IL0980.wt
cp IS.wt IS0980.wt
cp IT.wt IT0980.wt
cp DDzzduz.wt DDzzduz0980.wt
cp Out.wt Out0980.wt
cp Ratios.wt Ratios0980.wt
cp energies.wt energies0980.wt
cp Ewave.wt Ewave0980.wt
mv DDzzduz0980.wt Out0980.wt Ratios0980.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe0980.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 50 execução"
echo "$dt - Começando 50 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1000.wt
cp IL.wt IL1000.wt
cp IS.wt IS1000.wt
cp IT.wt IT1000.wt
cp DDzzduz.wt DDzzduz1000.wt
cp Out.wt Out1000.wt
cp Ratios.wt Ratios1000.wt
cp energies.wt energies1000.wt
cp Ewave.wt Ewave1000.wt
mv DDzzduz1000.wt Out1000.wt Ratios1000.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1000.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 51 execução"
echo "$dt - Começando 51 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1020.wt
cp IL.wt IL1020.wt
cp IS.wt IS1020.wt
cp IT.wt IT1020.wt
cp DDzzduz.wt DDzzduz1020.wt
cp Out.wt Out1020.wt
cp Ratios.wt Ratios1020.wt
cp energies.wt energies1020.wt
cp Ewave.wt Ewave1020.wt
mv DDzzduz1020.wt Out1020.wt Ratios1020.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1020.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 52 execução"
echo "$dt - Começando 52 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1040.wt
cp IL.wt IL1040.wt
cp IS.wt IS1040.wt
cp IT.wt IT1040.wt
cp DDzzduz.wt DDzzduz1040.wt
cp Out.wt Out1040.wt
cp Ratios.wt Ratios1040.wt
cp energies.wt energies1040.wt
cp Ewave.wt Ewave1040.wt
mv DDzzduz1040.wt Out1040.wt Ratios1040.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1040.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 53 execução"
echo "$dt - Começando 53 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1060.wt
cp IL.wt IL1060.wt
cp IS.wt IS1060.wt
cp IT.wt IT1060.wt
cp DDzzduz.wt DDzzduz1060.wt
cp Out.wt Out1060.wt
cp Ratios.wt Ratios1060.wt
cp energies.wt energies1060.wt
cp Ewave.wt Ewave1060.wt
mv DDzzduz1060.wt Out1060.wt Ratios1060.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1060.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 54 execução"
echo "$dt - Começando 54 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1080.wt
cp IL.wt IL1080.wt
cp IS.wt IS1080.wt
cp IT.wt IT1080.wt
cp DDzzduz.wt DDzzduz1080.wt
cp Out.wt Out1080.wt
cp Ratios.wt Ratios1080.wt
cp energies.wt energies1080.wt
cp Ewave.wt Ewave1080.wt
mv DDzzduz1080.wt Out1080.wt Ratios1080.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1080.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 55 execução"
echo "$dt - Começando 55 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1100.wt
cp IL.wt IL1100.wt
cp IS.wt IS1100.wt
cp IT.wt IT1100.wt
cp DDzzduz.wt DDzzduz1100.wt
cp Out.wt Out1100.wt
cp Ratios.wt Ratios1100.wt
cp energies.wt energies1100.wt
cp Ewave.wt Ewave1100.wt
mv DDzzduz1100.wt Out1100.wt Ratios1100.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1100.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 56 execução"
echo "$dt - Começando 56 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1120.wt
cp IL.wt IL1120.wt
cp IS.wt IS1120.wt
cp IT.wt IT1120.wt
cp DDzzduz.wt DDzzduz1120.wt
cp Out.wt Out1120.wt
cp Ratios.wt Ratios1120.wt
cp energies.wt energies1120.wt
cp Ewave.wt Ewave1120.wt
mv DDzzduz1120.wt Out1120.wt Ratios1120.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1120.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 57 execução"
echo "$dt - Começando 57 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1140.wt
cp IL.wt IL1140.wt
cp IS.wt IS1140.wt
cp IT.wt IT1140.wt
cp DDzzduz.wt DDzzduz1140.wt
cp Out.wt Out1140.wt
cp Ratios.wt Ratios1140.wt
cp energies.wt energies1140.wt
cp Ewave.wt Ewave1140.wt
mv DDzzduz1140.wt Out1140.wt Ratios1140.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1140.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 58 execução"
echo "$dt - Começando 58 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1160.wt
cp IL.wt IL1160.wt
cp IS.wt IS1160.wt
cp IT.wt IT1160.wt
cp DDzzduz.wt DDzzduz1160.wt
cp Out.wt Out1160.wt
cp Ratios.wt Ratios1160.wt
cp energies.wt energies1160.wt
cp Ewave.wt Ewave1160.wt
mv DDzzduz1160.wt Out1160.wt Ratios1160.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1160.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 59 execução"
echo "$dt - Começando 59 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1180.wt
cp IL.wt IL1180.wt
cp IS.wt IS1180.wt
cp IT.wt IT1180.wt
cp DDzzduz.wt DDzzduz1180.wt
cp Out.wt Out1180.wt
cp Ratios.wt Ratios1180.wt
cp energies.wt energies1180.wt
cp Ewave.wt Ewave1180.wt
mv DDzzduz1180.wt Out1180.wt Ratios1180.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1180.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 60 execução"
echo "$dt - Começando 60 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1200.wt
cp IL.wt IL1200.wt
cp IS.wt IS1200.wt
cp IT.wt IT1200.wt
cp DDzzduz.wt DDzzduz1200.wt
cp Out.wt Out1200.wt
cp Ratios.wt Ratios1200.wt
cp energies.wt energies1200.wt
cp Ewave.wt Ewave1200.wt
mv DDzzduz1200.wt Out1200.wt Ratios1200.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1200.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 61 execução"
echo "$dt - Começando 61 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1220.wt
cp IL.wt IL1220.wt
cp IS.wt IS1220.wt
cp IT.wt IT1220.wt
cp DDzzduz.wt DDzzduz1220.wt
cp Out.wt Out1220.wt
cp Ratios.wt Ratios1220.wt
cp energies.wt energies1220.wt
cp Ewave.wt Ewave1220.wt
mv DDzzduz1220.wt Out1220.wt Ratios1220.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1220.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 62 execução"
echo "$dt - Começando 62 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1240.wt
cp IL.wt IL1240.wt
cp IS.wt IS1240.wt
cp IT.wt IT1240.wt
cp DDzzduz.wt DDzzduz1240.wt
cp Out.wt Out1240.wt
cp Ratios.wt Ratios1240.wt
cp energies.wt energies1240.wt
cp Ewave.wt Ewave1240.wt
mv DDzzduz1240.wt Out1240.wt Ratios1240.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1240.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 63 execução"
echo "$dt - Começando 63 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1260.wt
cp IL.wt IL1260.wt
cp IS.wt IS1260.wt
cp IT.wt IT1260.wt
cp DDzzduz.wt DDzzduz1260.wt
cp Out.wt Out1260.wt
cp Ratios.wt Ratios1260.wt
cp energies.wt energies1260.wt
cp Ewave.wt Ewave1260.wt
mv DDzzduz1260.wt Out1260.wt Ratios1260.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1260.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 64 execução"
echo "$dt - Começando 64 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1280.wt
cp IL.wt IL1280.wt
cp IS.wt IS1280.wt
cp IT.wt IT1280.wt
cp DDzzduz.wt DDzzduz1280.wt
cp Out.wt Out1280.wt
cp Ratios.wt Ratios1280.wt
cp energies.wt energies1280.wt
cp Ewave.wt Ewave1280.wt
mv DDzzduz1280.wt Out1280.wt Ratios1280.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1280.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 65 execução"
echo "$dt - Começando 65 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1300.wt
cp IL.wt IL1300.wt
cp IS.wt IS1300.wt
cp IT.wt IT1300.wt
cp DDzzduz.wt DDzzduz1300.wt
cp Out.wt Out1300.wt
cp Ratios.wt Ratios1300.wt
cp energies.wt energies1300.wt
cp Ewave.wt Ewave1300.wt
mv DDzzduz1300.wt Out1300.wt Ratios1300.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1300.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 66 execução"
echo "$dt - Começando 66 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1320.wt
cp IL.wt IL1320.wt
cp IS.wt IS1320.wt
cp IT.wt IT1320.wt
cp DDzzduz.wt DDzzduz1320.wt
cp Out.wt Out1320.wt
cp Ratios.wt Ratios1320.wt
cp energies.wt energies1320.wt
cp Ewave.wt Ewave1320.wt
mv DDzzduz1320.wt Out1320.wt Ratios1320.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1320.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 67 execução"
echo "$dt - Começando 67 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1340.wt
cp IL.wt IL1340.wt
cp IS.wt IS1340.wt
cp IT.wt IT1340.wt
cp DDzzduz.wt DDzzduz1340.wt
cp Out.wt Out1340.wt
cp Ratios.wt Ratios1340.wt
cp energies.wt energies1340.wt
cp Ewave.wt Ewave1340.wt
mv DDzzduz1340.wt Out1340.wt Ratios1340.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1340.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 68 execução"
echo "$dt - Começando 68 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1360.wt
cp IL.wt IL1360.wt
cp IS.wt IS1360.wt
cp IT.wt IT1360.wt
cp DDzzduz.wt DDzzduz1360.wt
cp Out.wt Out1360.wt
cp Ratios.wt Ratios1360.wt
cp energies.wt energies1360.wt
cp Ewave.wt Ewave1360.wt
mv DDzzduz1360.wt Out1360.wt Ratios1360.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1360.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 69 execução"
echo "$dt - Começando 69 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1380.wt
cp IL.wt IL1380.wt
cp IS.wt IS1380.wt
cp IT.wt IT1380.wt
cp DDzzduz.wt DDzzduz1380.wt
cp Out.wt Out1380.wt
cp Ratios.wt Ratios1380.wt
cp energies.wt energies1380.wt
cp Ewave.wt Ewave1380.wt
mv DDzzduz1380.wt Out1380.wt Ratios1380.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1380.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 70 execução"
echo "$dt - Começando 70 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1400.wt
cp IL.wt IL1400.wt
cp IS.wt IS1400.wt
cp IT.wt IT1400.wt
cp DDzzduz.wt DDzzduz1400.wt
cp Out.wt Out1400.wt
cp Ratios.wt Ratios1400.wt
cp energies.wt energies1400.wt
cp Ewave.wt Ewave1400.wt
mv DDzzduz1400.wt Out1400.wt Ratios1400.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1400.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 71 execução"
echo "$dt - Começando 71 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1420.wt
cp IL.wt IL1420.wt
cp IS.wt IS1420.wt
cp IT.wt IT1420.wt
cp DDzzduz.wt DDzzduz1420.wt
cp Out.wt Out1420.wt
cp Ratios.wt Ratios1420.wt
cp energies.wt energies1420.wt
cp Ewave.wt Ewave1420.wt
mv DDzzduz1420.wt Out1420.wt Ratios1420.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1420.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 72 execução"
echo "$dt - Começando 72 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1440.wt
cp IL.wt IL1440.wt
cp IS.wt IS1440.wt
cp IT.wt IT1440.wt
cp DDzzduz.wt DDzzduz1440.wt
cp Out.wt Out1440.wt
cp Ratios.wt Ratios1440.wt
cp energies.wt energies1440.wt
cp Ewave.wt Ewave1440.wt
mv DDzzduz1440.wt Out1440.wt Ratios1440.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1440.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 73 execução"
echo "$dt - Começando 73 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1460.wt
cp IL.wt IL1460.wt
cp IS.wt IS1460.wt
cp IT.wt IT1460.wt
cp DDzzduz.wt DDzzduz1460.wt
cp Out.wt Out1460.wt
cp Ratios.wt Ratios1460.wt
cp energies.wt energies1460.wt
cp Ewave.wt Ewave1460.wt
mv DDzzduz1460.wt Out1460.wt Ratios1460.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1460.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 74 execução"
echo "$dt - Começando 74 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1480.wt
cp IL.wt IL1480.wt
cp IS.wt IS1480.wt
cp IT.wt IT1480.wt
cp DDzzduz.wt DDzzduz1480.wt
cp Out.wt Out1480.wt
cp Ratios.wt Ratios1480.wt
cp energies.wt energies1480.wt
cp Ewave.wt Ewave1480.wt
mv DDzzduz1480.wt Out1480.wt Ratios1480.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1480.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 75 execução"
echo "$dt - Começando 75 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1500.wt
cp IL.wt IL1500.wt
cp IS.wt IS1500.wt
cp IT.wt IT1500.wt
cp DDzzduz.wt DDzzduz1500.wt
cp Out.wt Out1500.wt
cp Ratios.wt Ratios1500.wt
cp energies.wt energies1500.wt
cp Ewave.wt Ewave1500.wt
mv DDzzduz1500.wt Out1500.wt Ratios1500.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1500.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 76 execução"
echo "$dt - Começando 76 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1520.wt
cp IL.wt IL1520.wt
cp IS.wt IS1520.wt
cp IT.wt IT1520.wt
cp DDzzduz.wt DDzzduz1520.wt
cp Out.wt Out1520.wt
cp Ratios.wt Ratios1520.wt
cp energies.wt energies1520.wt
cp Ewave.wt Ewave1520.wt
mv DDzzduz1520.wt Out1520.wt Ratios1520.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1520.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 77 execução"
echo "$dt - Começando 77 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1540.wt
cp IL.wt IL1540.wt
cp IS.wt IS1540.wt
cp IT.wt IT1540.wt
cp DDzzduz.wt DDzzduz1540.wt
cp Out.wt Out1540.wt
cp Ratios.wt Ratios1540.wt
cp energies.wt energies1540.wt
cp Ewave.wt Ewave1540.wt
mv DDzzduz1540.wt Out1540.wt Ratios1540.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1540.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 78 execução"
echo "$dt - Começando 78 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1560.wt
cp IL.wt IL1560.wt
cp IS.wt IS1560.wt
cp IT.wt IT1560.wt
cp DDzzduz.wt DDzzduz1560.wt
cp Out.wt Out1560.wt
cp Ratios.wt Ratios1560.wt
cp energies.wt energies1560.wt
cp Ewave.wt Ewave1560.wt
mv DDzzduz1560.wt Out1560.wt Ratios1560.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1560.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 79 execução"
echo "$dt - Começando 79 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1580.wt
cp IL.wt IL1580.wt
cp IS.wt IS1580.wt
cp IT.wt IT1580.wt
cp DDzzduz.wt DDzzduz1580.wt
cp Out.wt Out1580.wt
cp Ratios.wt Ratios1580.wt
cp energies.wt energies1580.wt
cp Ewave.wt Ewave1580.wt
mv DDzzduz1580.wt Out1580.wt Ratios1580.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1580.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 80 execução"
echo "$dt - Começando 80 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1600.wt
cp IL.wt IL1600.wt
cp IS.wt IS1600.wt
cp IT.wt IT1600.wt
cp DDzzduz.wt DDzzduz1600.wt
cp Out.wt Out1600.wt
cp Ratios.wt Ratios1600.wt
cp energies.wt energies1600.wt
cp Ewave.wt Ewave1600.wt
mv DDzzduz1600.wt Out1600.wt Ratios1600.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1600.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 81 execução"
echo "$dt - Começando 81 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1620.wt
cp IL.wt IL1620.wt
cp IS.wt IS1620.wt
cp IT.wt IT1620.wt
cp DDzzduz.wt DDzzduz1620.wt
cp Out.wt Out1620.wt
cp Ratios.wt Ratios1620.wt
cp energies.wt energies1620.wt
cp Ewave.wt Ewave1620.wt
mv DDzzduz1620.wt Out1620.wt Ratios1620.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1620.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 82 execução"
echo "$dt - Começando 82 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1640.wt
cp IL.wt IL1640.wt
cp IS.wt IS1640.wt
cp IT.wt IT1640.wt
cp DDzzduz.wt DDzzduz1640.wt
cp Out.wt Out1640.wt
cp Ratios.wt Ratios1640.wt
cp energies.wt energies1640.wt
cp Ewave.wt Ewave1640.wt
mv DDzzduz1640.wt Out1640.wt Ratios1640.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1640.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 83 execução"
echo "$dt - Começando 83 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1660.wt
cp IL.wt IL1660.wt
cp IS.wt IS1660.wt
cp IT.wt IT1660.wt
cp DDzzduz.wt DDzzduz1660.wt
cp Out.wt Out1660.wt
cp Ratios.wt Ratios1660.wt
cp energies.wt energies1660.wt
cp Ewave.wt Ewave1660.wt
mv DDzzduz1660.wt Out1660.wt Ratios1660.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1660.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 84 execução"
echo "$dt - Começando 84 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1680.wt
cp IL.wt IL1680.wt
cp IS.wt IS1680.wt
cp IT.wt IT1680.wt
cp DDzzduz.wt DDzzduz1680.wt
cp Out.wt Out1680.wt
cp Ratios.wt Ratios1680.wt
cp energies.wt energies1680.wt
cp Ewave.wt Ewave1680.wt
mv DDzzduz1680.wt Out1680.wt Ratios1680.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1680.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 85 execução"
echo "$dt - Começando 85 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1700.wt
cp IL.wt IL1700.wt
cp IS.wt IS1700.wt
cp IT.wt IT1700.wt
cp DDzzduz.wt DDzzduz1700.wt
cp Out.wt Out1700.wt
cp Ratios.wt Ratios1700.wt
cp energies.wt energies1700.wt
cp Ewave.wt Ewave1700.wt
mv DDzzduz1700.wt Out1700.wt Ratios1700.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1700.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 86 execução"
echo "$dt - Começando 86 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1720.wt
cp IL.wt IL1720.wt
cp IS.wt IS1720.wt
cp IT.wt IT1720.wt
cp DDzzduz.wt DDzzduz1720.wt
cp Out.wt Out1720.wt
cp Ratios.wt Ratios1720.wt
cp energies.wt energies1720.wt
cp Ewave.wt Ewave1720.wt
mv DDzzduz1720.wt Out1720.wt Ratios1720.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1720.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 87 execução"
echo "$dt - Começando 87 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1740.wt
cp IL.wt IL1740.wt
cp IS.wt IS1740.wt
cp IT.wt IT1740.wt
cp DDzzduz.wt DDzzduz1740.wt
cp Out.wt Out1740.wt
cp Ratios.wt Ratios1740.wt
cp energies.wt energies1740.wt
cp Ewave.wt Ewave1740.wt
mv DDzzduz1740.wt Out1740.wt Ratios1740.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1740.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 88 execução"
echo "$dt - Começando 88 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1760.wt
cp IL.wt IL1760.wt
cp IS.wt IS1760.wt
cp IT.wt IT1760.wt
cp DDzzduz.wt DDzzduz1760.wt
cp Out.wt Out1760.wt
cp Ratios.wt Ratios1760.wt
cp energies.wt energies1760.wt
cp Ewave.wt Ewave1760.wt
mv DDzzduz1760.wt Out1760.wt Ratios1760.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1760.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 89 execução"
echo "$dt - Começando 89 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1780.wt
cp IL.wt IL1780.wt
cp IS.wt IS1780.wt
cp IT.wt IT1780.wt
cp DDzzduz.wt DDzzduz1780.wt
cp Out.wt Out1780.wt
cp Ratios.wt Ratios1780.wt
cp energies.wt energies1780.wt
cp Ewave.wt Ewave1780.wt
mv DDzzduz1780.wt Out1780.wt Ratios1780.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1780.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 90 execução"
echo "$dt - Começando 90 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1800.wt
cp IL.wt IL1800.wt
cp IS.wt IS1800.wt
cp IT.wt IT1800.wt
cp DDzzduz.wt DDzzduz1800.wt
cp Out.wt Out1800.wt
cp Ratios.wt Ratios1800.wt
cp energies.wt energies1800.wt
cp Ewave.wt Ewave1800.wt
mv DDzzduz1800.wt Out1800.wt Ratios1800.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1800.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 91 execução"
echo "$dt - Começando 91 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1820.wt
cp IL.wt IL1820.wt
cp IS.wt IS1820.wt
cp IT.wt IT1820.wt
cp DDzzduz.wt DDzzduz1820.wt
cp Out.wt Out1820.wt
cp Ratios.wt Ratios1820.wt
cp energies.wt energies1820.wt
cp Ewave.wt Ewave1820.wt
mv DDzzduz1820.wt Out1820.wt Ratios1820.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1820.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 92 execução"
echo "$dt - Começando 92 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1840.wt
cp IL.wt IL1840.wt
cp IS.wt IS1840.wt
cp IT.wt IT1840.wt
cp DDzzduz.wt DDzzduz1840.wt
cp Out.wt Out1840.wt
cp Ratios.wt Ratios1840.wt
cp energies.wt energies1840.wt
cp Ewave.wt Ewave1840.wt
mv DDzzduz1840.wt Out1840.wt Ratios1840.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1840.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 93 execução"
echo "$dt - Começando 93 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1860.wt
cp IL.wt IL1860.wt
cp IS.wt IS1860.wt
cp IT.wt IT1860.wt
cp DDzzduz.wt DDzzduz1860.wt
cp Out.wt Out1860.wt
cp Ratios.wt Ratios1860.wt
cp energies.wt energies1860.wt
cp Ewave.wt Ewave1860.wt
mv DDzzduz1860.wt Out1860.wt Ratios1860.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1860.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 94 execução"
echo "$dt - Começando 94 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1880.wt
cp IL.wt IL1880.wt
cp IS.wt IS1880.wt
cp IT.wt IT1880.wt
cp DDzzduz.wt DDzzduz1880.wt
cp Out.wt Out1880.wt
cp Ratios.wt Ratios1880.wt
cp energies.wt energies1880.wt
cp Ewave.wt Ewave1880.wt
mv DDzzduz1880.wt Out1880.wt Ratios1880.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1880.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 95 execução"
echo "$dt - Começando 95 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1900.wt
cp IL.wt IL1900.wt
cp IS.wt IS1900.wt
cp IT.wt IT1900.wt
cp DDzzduz.wt DDzzduz1900.wt
cp Out.wt Out1900.wt
cp Ratios.wt Ratios1900.wt
cp energies.wt energies1900.wt
cp Ewave.wt Ewave1900.wt
mv DDzzduz1900.wt Out1900.wt Ratios1900.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1900.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 96 execução"
echo "$dt - Começando 96 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1920.wt
cp IL.wt IL1920.wt
cp IS.wt IS1920.wt
cp IT.wt IT1920.wt
cp DDzzduz.wt DDzzduz1920.wt
cp Out.wt Out1920.wt
cp Ratios.wt Ratios1920.wt
cp energies.wt energies1920.wt
cp Ewave.wt Ewave1920.wt
mv DDzzduz1920.wt Out1920.wt Ratios1920.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1920.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 97 execução"
echo "$dt - Começando 97 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1940.wt
cp IL.wt IL1940.wt
cp IS.wt IS1940.wt
cp IT.wt IT1940.wt
cp DDzzduz.wt DDzzduz1940.wt
cp Out.wt Out1940.wt
cp Ratios.wt Ratios1940.wt
cp energies.wt energies1940.wt
cp Ewave.wt Ewave1940.wt
mv DDzzduz1940.wt Out1940.wt Ratios1940.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1940.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 98 execução"
echo "$dt - Começando 98 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1960.wt
cp IL.wt IL1960.wt
cp IS.wt IS1960.wt
cp IT.wt IT1960.wt
cp DDzzduz.wt DDzzduz1960.wt
cp Out.wt Out1960.wt
cp Ratios.wt Ratios1960.wt
cp energies.wt energies1960.wt
cp Ewave.wt Ewave1960.wt
mv DDzzduz1960.wt Out1960.wt Ratios1960.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1960.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 99 execução"
echo "$dt - Começando 99 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe1980.wt
cp IL.wt IL1980.wt
cp IS.wt IS1980.wt
cp IT.wt IT1980.wt
cp DDzzduz.wt DDzzduz1980.wt
cp Out.wt Out1980.wt
cp Ratios.wt Ratios1980.wt
cp energies.wt energies1980.wt
cp Ewave.wt Ewave1980.wt
mv DDzzduz1980.wt Out1980.wt Ratios1980.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe1980.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt - Começando 100 execução"
echo "$dt - Começando 100 execução" >> Tempo.txt
./a.out

cp Fe.wt Fe2000.wt
cp IL.wt IL2000.wt
cp IS.wt IS2000.wt
cp IT.wt IT2000.wt
cp DDzzduz.wt DDzzduz2000.wt
cp Out.wt Out2000.wt
cp Ratios.wt Ratios2000.wt
cp energies.wt energies2000.wt
cp Ewave.wt Ewave2000.wt
mv DDzzduz2000.wt Out2000.wt Ratios2000.wt Files 
mv Out.wt Ini.wt

if grep -q "NaN" Fe2000.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi

#---------------------------------------------

#paplay /usr/share/sounds/LinuxMint/stereo/desktop-login.ogg
shutdown now