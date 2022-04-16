import math
import os

Nf = 2e4
Ncase = 20
Ndrop = 2
workdir = input("Please input working directory:")
os.chdir(workdir)

for i in range(1,Ncase+1):
    os.chdir(str(i))
    with open('JZIF3.dat','r+') as fid:
        fid2 = open('J.dat','w')
        for j in range(0,int(Ndrop+Nf+1)):
            line = fid.readline()
            fid2.write(line)
        fid2.close()
    os.chdir(workdir)