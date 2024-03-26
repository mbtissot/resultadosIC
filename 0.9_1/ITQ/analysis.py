import matplotlib.pyplot as plt
import os

maxesPeak1 = []
maxesPeak2 = []
maxesPeak3 = []

path = os.getcwd()
pathList = path.split("/")
folder = pathList[-2]
#print(folder)

currentDir = os.listdir()

files = [file for file in currentDir if file[0:4]=="ITQ1"]

nFiles = len(files)

valores = [i*100 for i in range(1, nFiles+1)]

def openFile(path: str)->list[str]:
	with open(path, "r") as f:
		vals = f.readlines()
		return vals

for i in range(1,nFiles+1):
	file = "ITQ1" + "{:04d}".format(100*i)+".wt"

	rawValues = openFile(file)
	values = [float(rawValues[i].strip(" ")) for i in range(len(rawValues))]

	maxesPeak1.append(max(values[:10]))
	maxesPeak2.append(max(values[11:16]))
	maxesPeak3.append(max(values[16:23]))


with open("fs/"+folder+"_maxes.txt", "w") as f:
	f.write("{: ^24}{: ^24}{: ^24}{: ^24}\n".format("t", "Max1", "Max2", "Max3"))
	for i in range(len(valores)):
		f.write(f"{valores[i]}".ljust(24) + str(maxesPeak1[i]).ljust(24) + str(maxesPeak2[i]).ljust(24) + str(maxesPeak3[i]).ljust(24) + "\n")



ratios12 = [maxesPeak1[i]/maxesPeak2[i] for i in range(nFiles)]
ratios23 = [maxesPeak2[i]/maxesPeak3[i] for i in range(nFiles)]

with open("fs/"+folder+"_ratios.txt", "w") as f:
	f.write("{: ^24}{: ^24}{: ^24}\n".format("t","Max1/Max2","Max2/Max3"))
	for i in range(len(valores)):
		f.write(f"{valores[i]}".ljust(24) + str(ratios12[i]).ljust(24) + str(ratios23[i]).ljust(24)+"\n")
#print(maxesPeak2)

# print(valores)

# plt.vlines(valores[11], 1e-7, 1e-3, color="grey")

# for i in range(nFiles):
# 	plt.plot(valores[i], maxesPeak1[i],"bo")
# 	plt.plot(valores[i], maxesPeak2[i],"go")
# 	plt.plot(valores[i], maxesPeak3[i],"ro")
# 	plt.vlines(valores[i], maxesPeak2[i], maxesPeak1[i], "#007d7d80")
# 	plt.vlines(valores[i], maxesPeak3[i], maxesPeak2[i], "#9bc80080")
# plt.yscale("log")
# plt.ylabel("Valor Pico")
# plt.xlabel("Tempo")
# plt.xticks(valores[::4])
# plt.savefig(folder+"_maxes.png")
# plt.show()

# ratios12 = [maxesPeak1[i]/maxesPeak2[i] for i in range(nFiles)]
# ratios23 = [maxesPeak2[i]/maxesPeak3[i] for i in range(nFiles)]
# plt.plot(valores, ratios12, 'b', label="Pico1/Pico2")
# plt.plot(valores, ratios23, 'r', label="Pico2/Pico3")
# # plt.vlines(valores[11],0, 75, color="grey")
# plt.xticks(valores[::4])
# plt.xlabel("Tempo")
# plt.ylabel("Pico$_n$/Pico$_{n+1}$")
# plt.legend()
# plt.savefig(folder+"_ratios.png")
# plt.show()