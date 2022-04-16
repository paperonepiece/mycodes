import math
import os
import numpy as np

olddir = os.getcwd()
os.chdir(olddir)
read1 = open("energy.dat", "r+")
read2 = open("POSCAR", "r+")
alllines = read1.readlines()
struinfo = read2.readlines()
read1.close()
read2.close()
count = len(alllines)
a0 =  float(struinfo[2].split()[0])
aE = np.zeros([count,2])
numlist = struinfo[6].split()
atomnum = sum([int(i) for i in numlist])
for i in range(0,count):
    aE[i,0] = float(alllines[i].split()[0])
    aE[i,1] = float(alllines[i].split()[1])
a1 = aE[:,0]*a0
Eper = aE[:,1]/atomnum

x = (a1)**(-2)
p = np.polyfit(x, aE[:,1], 3)
c0 = p[3]
c1 = p[2]
c2 = p[1]
c3 = p[0]
x1 = (math.sqrt(4*c2**2-12*c1*c3)-2*c2)/(6*c3)
final = 1/math.sqrt(x1)
finalE = (c0+c1*x1+c2*x1**2+c3*x1**3)/atomnum
aEout = np.vstack((a1,Eper))

output = open("out.dat", "w")
output.write('The final lattice constant is '+str(final)+'\n'+
             'The final cohesive energy is '+str(finalE)+'\n')
for i in range(0,count):
    #output.write(str(aEout[0,i])+' '+str(aEout[1,i])+'\n')
    output.write(format(aEout[0,i],'<.6f')+' '+format(aEout[1,i],'<.6f')+'\n')
output.close()
print('The final lattice constant is '+str(final)+'\n'+
      'The final cohesive energy is '+str(finalE))