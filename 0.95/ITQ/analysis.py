import matplotlib.pyplot as plt
import os

maxesPeak1 = []
maxesPeak2 = []

path = os.getcwd()
pathlist = path.split("/")
folder = pathlist[-2]

try:
	os.mkdir("maxes")
except Exception:
	print("diretorio maxes already exists")

filelist = list(os.listdir())

numfiles = int((len(filelist)-1)/2)
print(numfiles)

valores = [i*100 for i in range(1, numfiles)]

def openFile(path: str)->list[str]:
	with open(path, "r") as f:
		vals = f.readlines()
		return vals

for i in range(1,numfiles):
	file = "ITQ1" + "{:04d}".format(100*i)+".wt"#"ITQ11500.wt" 

	rawValues = openFile(file)
	values = [float(rawValues[i].strip(" ")) for i in range(len(rawValues))]

	maxesPeak1.append(max(values[:10]))
	maxesPeak2.append(max(values[11:16]))

# print(maxesPeak1[numfiles-1], maxesPeak2[numfiles-1])

for i in range(numfiles-1):
	plt.vlines(valores[i], maxesPeak2[i], maxesPeak1[i])
	plt.plot(valores[i], maxesPeak1[i],"bo")
	plt.plot(valores[i], maxesPeak2[i],"go")
plt.yscale("log")
plt.ylabel("Valor Pico")
plt.xlabel("Tempo")
plt.xticks(valores[::4])
plt.savefig("./"+folder+"_maxes.png")
# plt.show()

with open("./maxes/"+folder+"maxes.txt", "w") as f:
	for i in range(len(valores)):
		f.write(f"{valores[i]:5d}\t\t {maxesPeak1[i]:.2E}\t\t {maxesPeak2[i]:.2E}\n")

ratios = [maxesPeak2[i]/maxesPeak1[i] for i in range(numfiles-1)]
plt.plot(valores, ratios)
# plt.vlines(valores[11],0, 75, color="grey")
plt.xticks(valores[::4])
plt.xlabel("Tempo")
plt.ylabel("Pico2/Pico1")
plt.savefig("maxes/"+folder+"_ratios.png")
# plt.show()

with open("./maxes/"+folder+"_ratios.txt", "w") as f:
	for i in range(len(valores)):
		f.write(f"{valores[i]:5d}\t\t {ratios[i]:.2E}\n")