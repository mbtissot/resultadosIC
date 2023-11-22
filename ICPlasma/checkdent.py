#!/usr/bin/env python3

import sys

if __name__=="__main__":
    if len(sys.argv)!=2:
        print("Usage:")
        sys.exit(1)
    else:
        param1 = sys.argv[1]

        try:
            index = param1.rindex('/')
            filepath = param1[:index]
            filename = param1[index+1:]
        except:
            index = -1
            filepath = ''
            filename = param1

        print("Filepath: ", filepath)
        print("Filename: ", filename)

        output = filepath + '/'*(filepath!='') +'C' + filename

        print("Output: ", output)

        empty = " "

        f=open(filename,"r")
        lines=f.readlines()
        f.close()

        # Valores: Ux = [0], Uz = [1], Valor = [2]

        # # # Removing the Uz line of unwanted lines
        for lineidx in range(len(lines)):
            try:
                if float(lines[lineidx].split("     ")[1]) < 5.5 or float(lines[lineidx].split("     ")[1])>7.5:
                    lines[lineidx] = empty
            except:
                continue

        # # Removing the Ux line of unwanted lines
        for lineidx in range(len(lines)):
            try:
                if float(lines[lineidx].split("      ")[0]) < -1.5 or float(lines[lineidx].split("      ")[0])>1.5:
                    lines[lineidx] = empty
            except:
                continue



        with open(output, "w") as f:
            f.writelines(lines)