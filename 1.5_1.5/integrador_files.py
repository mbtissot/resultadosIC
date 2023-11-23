#!/bin/python3
import sys
import os
import re

def fixlines(lines: list[str])->list[str]:
	newLines = []
	for line in lines:
		line = " ".join(line.lstrip(" ").split(" ")).split()
		newLines.append(line)
	return newLines

def openAndRead(address: str):
	with open(address,"r") as f:
		lines = f.readlines()
		lines = fixlines(lines)

	return lines

def separateXandZ(lines: list[str]) -> list[list[float]]:
	values = []
	for line in lines:
		if len(line)>1 and float(line[0])>0:
			values.append([float(line[0]),float(line[2])])
	return values

def integrateOnX(XandZ):
	newDict = dict()
	for item in XandZ:
		if item[0] not in newDict:
			newDict[item[0]] = item[1]
		else:
			newDict[item[0]] += item[1]

	lista = list(newDict.items())
	return lista

def makeFileFromIntegral(listIntegral, address):
	with open('I'+address, 'w') as f:
		for item in listIntegral:
			# f.write(str(item[0]) + '  ' + str(item[1]) + '\n')
			f.write('{0:.18f}'.format(item[0]) + '   ' + '{0:.18f}'.format(item[1]) + '\n')



# def main():
# 	if len(sys.argv)!=2:
# 		print("""INVALID USE.    USAGE: integrate.py <filename>""")
# 		exit()
# 	filename = sys.argv[1]
# 	XandZ = separateXandZ(openAndRead(filename))
# 	integrado = integrateOnX(XandZ)
# 	makeFileFromIntegral(integrado, filename)

def main():
	folder = os.listdir()
	folder.sort()
	for file in folder:
		if re.search("^IT+([0-9]+)+\.+wt$", file):
			XandZ = separateXandZ(openAndRead(file))
			integrado = integrateOnX(XandZ)
			makeFileFromIntegral(integrado, file)


if __name__ == '__main__':
	main()