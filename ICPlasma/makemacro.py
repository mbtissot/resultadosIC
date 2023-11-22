#!/usr/bin/env python3

import sys

def makestart(dtau: float, tauadd: int):
	names = ['DTau', 'TauAdd     ! Final time: Tau2= Tau1+TauAdd', 'TimeEvol', 'OneDfiles',
	'TwoDfiles', 'Onecolumnfiles', 'ADcoefficients', 'NewEffectsfiles']
	yesorno = ' ("Yes" or "No ")'
	values = ["{:.1e}".format(dtau).replace('e', 'E'), "{:.1e}".format(tauadd).replace('e', 'E'), '"Yes"', '"No "',
	'"Yes"', '"No "', '"No "', '"No "']
	with open('./Start.wt', 'w') as file:
		for i in range(len(names)):
			file.write(values[i].ljust(13) + names[i] + yesorno*(i>1) + '\n')

def makemacro(step: int, finish: int):
	import os
	digs = len(str(finish))
	sep = '\n#---------------------------------------------\n'
	howmany = 0
	counter = 0
	dt = "dt=$(date '+%d/%m/%Y %H:%M:%S');"
	os.system('cp ~/Desktop/Development/ICPlasma/macrosource.sh ./macro.sh')
	with open('macro.sh', 'a') as macro:
		while howmany<finish:
			counter=counter+1
			num = str(howmany+step).rjust(digs, '0')
			macro.write('\n'+dt+'\n')
			macro.write(f'echo "$dt - Começando {counter} execução"\n')
			macro.write(f'echo "$dt - Começando {counter} execução" >> Tempo.txt\n')
			macro.write("./a.out\n")
			macro.write(f'''\ncp Fe.wt Fe{num}.wt
cp IL.wt IL{num}.wt
cp IS.wt IS{num}.wt
cp IT.wt IT{num}.wt
cp DDzzduz.wt DDzzduz{num}.wt
cp Out.wt Out{num}.wt
cp Ratios.wt Ratios{num}.wt
mv IT{num}.wt DDzzduz{num}.wt Out{num}.wt Ratios{num}.wt Files 
mv Out.wt Ini.wt\n\n''')
			macro.write(check(num))
			macro.write(sep)
			howmany = howmany + step
		macro.write('''\npaplay /usr/share/sounds/LinuxMint/stereo/desktop-login.ogg\n#shutdown now''')

def check(num):
	stringIF = '''if grep -q "NaN" Fe'''+num+'''.wt;
then echo 'Achei um NaN na função distribuição e terminei a execução.' && return
fi\n'''
	return stringIF

if __name__ == '__main__':
	usage = '''USAGE: makemacro.py DTau (float or int) TauAdd (int) FinalTime (int).
		Make sure that FinalTime>TauAdd>DTau.

		It will create a Start.wt file with DTau and TauAdd, and a macro file that
		compiles the program, and runs it until FinalTime.
		Eg.: makemacro.py 0.1 10 100
		will create a "Start.wt" file with DTau=0.1 and TauAdd=10, and a macro file that compiles,
		runs and saves the results of 10 consecutive runs of the program.
		Saving them on files "Fe010.wt", "Fe020.wt" and so on until "Fe100.wt".'''
	if len(sys.argv) != 4:
		print(usage)
		sys.exit(1)
	else:
		try:
			dtau      = float(sys.argv[1])
			tauadd    = int(sys.argv[2])
			finaltime = int(sys.argv[3])
			if finaltime>=tauadd and tauadd>dtau:
				makestart(dtau, tauadd)
				makemacro(tauadd, finaltime)
			else:
				print('''You entered the right types, but either your DTau was bigger than TauAdd, or TauAdd was bigger than FinalTime.
					
					Please run "macro" on your terminal to see how to use it.
					''')
		except Exception:
			print('You inserted one wrong type of the values.\nPlease type "macro" on your terminal to see usage' )