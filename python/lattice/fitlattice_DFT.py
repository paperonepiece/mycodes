import math
import os
import numpy as np

olddir = os.getcwd()
os.chdir(olddir)
read = open("energy.dat", "r+")
alllines = read.readlines()
read.close()
count = len(alllines)
aE = np.zeros([count,2])
for i in range(0,count):
    aE[i,0] = float(alllines[i].split()[0])
    aE[i,1] = float(alllines[i].split()[1])
    
x = (aE[:,0]*5.468728)**(-2)
p = np.polyfit(x, aE[:,1], 3)
c0 = p[3]
c1 = p[2]
c2 = p[1]
c3 = p[0]

x1 = (math.sqrt(4*c2**2-12*c1*c3)-2*c2)/(6*c3)
final = 1/math.sqrt(x1)
print('The final lattice constant is '+str(final))