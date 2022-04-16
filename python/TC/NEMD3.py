#NEMD method for TC caculation#
#------Import packages------#
import numpy as np
from scipy import integrate,fft,optimize
import os
from math import floor
#------Parameters setting------#
Nf = int(5e5)   #total steps
Nt = int(5e4)   #calculate temperature steps
N  = int(Nf//Nt)
Nt1 = int(1e3)  #thermo steps
Nt2 = int(50)  #layer number
Ndrop = 3       #tmp.profile
Ndrop1 = 1      #tmp.profile1
Ndrop2 = 2      #heat.dat
kb = 1.38064852e-23
eV = 1.60217662e-19
length = 200.5*1e-10
area = 353.44*1e-20
dt = 0.5e-15
scale_temp = 1/length
scale_heat = eV/dt/area
Ncase = 1
limt1 = 0.3
limt2 = 0.7    # linear fitting interval
Tc   = np.zeros([Ncase,N])
heat = np.zeros([2,int(Nf/Nt1)])
temp = np.zeros([Nt2,N,2])
temperature=np.zeros([2,Nt2])
#------Enter the path------#
olddir = os.getcwd()
print("Now we are in:"+str(olddir))
workdir = input("Please input working directory:")
os.chdir(workdir)
#------Reading and caculating------#
for i in range(1,Ncase+1):
    os.chdir(str(i))
    
    with open('GY_heatflux.dat','r+') as fid1:
        for j in range(1,Ndrop2+1):
            next(fid1)
        stringlist = fid1.readlines()
        linenum = 0
        for line in stringlist:
            linelist = line.split()
            heat[0,linenum] = linelist[0]
            heat[1,linenum] = linelist[1]
            linenum += 1
            if linenum >= int(Nf/Nt1):
                break
            
    with open('GYtmp.profile','r+') as fid2:
        for k in range(1,Ndrop+Ndrop1+1):
            next(fid2)
        for m in range(0,Nt2):
            string1 = fid2.readline()
            temp[m,0,0] = string1.split()[1]
            temp[m,0,1] = string1.split()[3]
        for n in range(1,N):
            next(fid2)
            for r in range(0,Nt2):
                string2 = fid2.readline()
                temp[r,n,0] = string2.split()[1]
                temp[r,n,1] = string2.split()[3]
    
    def f1(x, A, B):
        return A*x + B
    for j in range(1,N+1):
        Kh = optimize.curve_fit(f1,heat[0,(j-1)*int(Nt/Nt1):j*int(Nt/Nt1)],
                                heat[1,(j-1)*int(Nt/Nt1):j*int(Nt/Nt1)])[0]
        Kheat=Kh[0]
        temperature[0,:]=temp[:,j-1,0]
        temperature[1,:]=temp[:,j-1,1]
        Kt = optimize.curve_fit(f1,temperature[0,floor(limt1*Nt2):floor(limt2*Nt2)],
                   temperature[1,floor(limt1*Nt2):floor(limt2*Nt2)])[0]
        Ktemp= Kt[0]
        Tc[i-1,j-1]=Kheat*scale_heat/Ktemp/scale_temp
    os.chdir(workdir)
#------Output results------#     
AveTc = np.mean(Tc,axis=0)
StdTc = np.std(Tc,axis=0)
UltTc = np.mean(AveTc[int(0.5*N):])
with open('TC_python.dat','w') as fidw:
    fidw.write('Final TC= '+format(UltTc,'^.6f')+' W/mk'+'\n')
    for i in range(1,N+1):
        fidw.write(format(i*Nt*dt*1e12,'^.6f'))
        fidw.write(format(AveTc[i-1],'^12.6f')+format(StdTc[i-1],'^12.6f'))
        fidw.write('\n')
print("NEMD TC caculation is done!")
#---------------------Done-------------------------#

