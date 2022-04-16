#DOS caculation
#Applicable to python 3.6 and later version
#------Import libraries------#
import numpy as np
from scipy import integrate,fft,optimize
import os
from math import floor
#------Parameters setting------#
totalsteps = 1e5
thermsteps = 10
timesteps = 4e-15
times = int(totalsteps/thermsteps)
Ndrop0 = 3
Ndrop = 6
Ndrop1 = 10
xcorrvelocity = 0
ycorrvelocity = 0
zcorrvelocity = 0
vx = np.zeros([1,2*times-1])
vy = np.zeros([1,2*times-1])
vz = np.zeros([1,2*times-1])
vx_coe = np.zeros([1,2*times-1])
vy_coe = np.zeros([1,2*times-1])
vz_coe = np.zeros([1,2*times-1])
#------Enter the path------#
olddir = os.getcwd()
#print("Now we are in:"+str(olddir))
#workdir = input("Please input working directory, example like 'D:\\xx\\yy':")
os.chdir(olddir)
#------Reading and caculating------#
with open('velocityAlO.txt','r+') as fid1:
        print("Reading velocities......")
        for i in range(1,Ndrop0+1):
            next(fid1)
        num = int(fid1.readline().split()[0])
        xvelocity=np.zeros([num,times])
        yvelocity=np.zeros([num,times])
        zvelocity=np.zeros([num,times])
        for i in range(1,Ndrop):
            next(fid1)
        for i in range(0,num):
            string1 = fid1.readline()
            xvelocity[i,0] = string1.split()[0]
            yvelocity[i,0] = string1.split()[1]
            zvelocity[i,0] = string1.split()[2]
        for j in range(1,times):
            for i in range(1,Ndrop1):
                next(fid1)
            for i in range(0,num):
                string1 = fid1.readline()
                xvelocity[i,j] = string1.split()[0]
                yvelocity[i,j] = string1.split()[1]
                zvelocity[i,j] = string1.split()[2]
print("DOS caculation wolud take some time, please be patient~")
for i in range(0,num):
    vx = np.correlate(xvelocity[i,:],xvelocity[i,:],'full')
    vy = np.correlate(yvelocity[i,:],yvelocity[i,:],'full')
    vz = np.correlate(zvelocity[i,:],zvelocity[i,:],'full')
    vx_coe = vx/vx[times-1]
    vy_coe = vy/vy[times-1]
    vz_coe = vz/vz[times-1]
    xcorrvelocity=xcorrvelocity+vx_coe
    ycorrvelocity=ycorrvelocity+vy_coe
    zcorrvelocity=zcorrvelocity+vz_coe
 
xDOS = fft.fft(xcorrvelocity)
yDOS = fft.fft(ycorrvelocity)
zDOS = fft.fft(zcorrvelocity)
xPDOS=abs(xDOS[0:times])*2/num
yPDOS=abs(yDOS[0:times])*2/num
zPDOS=abs(zDOS[0:times])*2/num
PDOS = np.zeros([1,times])
FinalPDOS = np.zeros([times,5])

for i in range(0,times):
    PDOS[0,i]=(xPDOS[i]+yPDOS[i]+zPDOS[i])/3       
for i in range(0,times):    
    FinalPDOS[i,0] = 1/timesteps/1e12/2*(i)/(times-1);
    FinalPDOS[i,1] = xPDOS[i]
    FinalPDOS[i,2] = yPDOS[i]
    FinalPDOS[i,3] = zPDOS[i]
    FinalPDOS[i,4] = PDOS[0,i]
#------Output results------#
with open('PDOSAlO_python.txt','w') as fid2:
    fid2.write(' 1-Frecquency '+' 2-xPDOS '+' 3-yPDOS '
               +' 4-zPDOS '+' 5-TotalPDOS '+'\n')
    for i in range(0,times):
        fid2.write(format(FinalPDOS[i,0],'<14.10f')+
                   format(FinalPDOS[i,1],'^14.10f')+
                   format(FinalPDOS[i,2],'^14.10f')+
                   format(FinalPDOS[i,3],'^14.10f')+
                   format(FinalPDOS[i,4],'^14.10f')+'\n')
print("DOS caculation is done!")            
#---------------------Done-------------------------
            
            
