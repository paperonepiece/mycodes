import numpy as np
import os

Nt = int(2e3)

K    = np.zeros([5,Nt])
StdK = np.zeros([1,5])
workdir = input("Please input working directory:")
os.chdir(workdir)

with open('K_total_python.dat','r+') as fid:
        stringlist = fid.readlines()
        linenum = 0
        for line in stringlist:
            linelist = line.split()
            K[0,linenum] = linelist[0]
            K[1,linenum] = linelist[1]
            K[2,linenum] = linelist[2]
            K[3,linenum]= linelist[3]
            K[4,linenum]= linelist[4]
            linenum += 1
            if linenum >= Nt:
                break
StdT    = np.std(K[0,int(0.3*Nt):int(0.8*Nt)])
StdKV   = np.std(K[1,int(0.3*Nt):int(0.8*Nt)])
StdKK   = np.std(K[2,int(0.3*Nt):int(0.8*Nt)])
StdKC   = np.std(K[3,int(0.3*Nt):int(0.8*Nt)])
StdK    = np.std(K[4,int(0.3*Nt):int(0.8*Nt)])
with open('Std.txt','w') as fidw:
    fidw.write(format(StdT,'e')+' '+format(StdKV,'e')+' '+format(StdKK,'e')
               +' '+format(StdKC,'e')+' '+format(StdK,'e'))
