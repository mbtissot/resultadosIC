filename = 'output.txt'

# open the file for reading. stores its lines in a list of lines
with open(filename, 'r') as f:
	lines = f.readlines()

#grabs length of list of lines
length = len(lines)

# finds line where the start of 7th execution starts
for line in range(length):
	if lines[line].find("Começando 7")!=-1:
		linhacritica = line

# loops from linhacritica to the end of file writing to cleanoutput.txt
with open('cleanoutput.txt', 'w') as f:
	for i in range(linhacritica, len(lines)):
		f.write(lines[i])

# with clean output, read lines
with open('cleanoutput.txt', 'r') as f:
	cleanlines = f.readlines()

# writes values of 'E-XXX' to the file 'magnitudes.txt'
with open('magnitudes.txt', 'w') as f:
	for line in range(len(cleanlines)):
		index = cleanlines[line].find("E-")
		f.write(cleanlines[line][index:index+5])
		f.write('\n')