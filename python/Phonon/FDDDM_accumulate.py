import numpy as np
import os
import math


workdir = input("Please input working directory:")
os.chdir(workdir)
read1 = open("Heat_Flux_Spectrum_L2R.dat", "r+")
write1 = open("accumu.dat", "w")
ev = 1.60217662*1e-19  # ev to J
tconv = 1e-12          # ps to s
lconv = 1e-10          # A  to m
S =  1962.5*1e-20      # cross sectional area--m^2
length1 =  2e-9         # lenth of control volume--m
length2 =  2e-8            # lenth of heat transport region--m
deltaT = 1.0*1e9            # temperature gradient--K/m
allline = read1.readlines()
content = allline[1:]
accumu = np.zeros([len(content),5])

for i in range(0,len(content)):
    accumu[i,0] = float(content[i].split()[0])
    accumu[i,1] = float(content[i].split()[3])*1e6
    accumu[i,2] = accumu[i,1]*length1/deltaT/length1
for i in range(0,len(content)):
    if i < 1:
        accumu[i,3] = accumu[i,1]/np.sum(accumu,axis=0)[1]
        accumu[i,4] = accumu[i,3]*np.sum(accumu,axis=0)[2]
    else:
        accumu[i,3] = accumu[i,1]/np.sum(accumu,axis=0)[1]+accumu[i-1,3]
        accumu[i,4] = accumu[i,3]*np.sum(accumu,axis=0)[2]

write1.write('Frequency--THz  HeatCurrent--J/m^2/s  TC--W/m/k  AccumuTC--W/m/k  Normalized Accumu--1'+'\n')
for i in range(0,len(content)):
    write1.write(str(accumu[i,0])+'  '+str(accumu[i,1])+'  '
                 +str(accumu[i,2])+'  '+str(accumu[i,4])+'  '
                 +str(accumu[i,3])+'\n')
read1.close()
write1.close()
print('Calculation is done. Congratulations!')