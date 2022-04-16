
########################################################################################
##this code bases on the old grammar of python2##
##be careful when using python3, some modifications are needed##
#import numpy as np
g     = 1.0e-3           # kg
A     = 1.0e-10          # m
nm    = 1.0e-9           # m
eV    = 1.6021766208e-19 # J
kJ    = 1.0e3            # J
kcal  = 4.1868e3         # J
mol   = 6.022140857e23   #
kB    = 1.38064852e-23   # Boltzmann constant
h     = 6.62606957e-34   # Planck    constant
fs    = 1.0e-15          # s
ps    = 1.0e-12          # s
ns    = 1.0e-9           # s
thz   = 1.0e12           # Hz
T     = 300.0            # K
dT    = 20               # K
bar   = 1.0e5
K = 1.0 # K

########################################################################################
num = 40       #the number of lines for average calculation of E1 and E2
natoms1 = 1000 #Pb
natoms2 = 3000 #I
natoms3 = 6000 #H
natoms4 = 1000 #N
natoms5 = 1000 #C
mass   = 207.2*g/mol*natoms1+126.9*g/mol*natoms2+1.0079*g/mol*natoms3+14*g/mol*natoms4+12.01*g/mol*natoms5 # kg
volume = 248553*A*A*A
density = mass/volume
start = []
final = []
#============Reading data===========#
fid = open("thermo_Cv.dat","r+")
linelist = fid.readlines()
for i in range(2,num+2):
    start.append(float(linelist[i].split()[2]))
for i in range(len(linelist)-num,len(linelist)):
    final.append(float(linelist[i].split()[2]))
fid.close()
#============Calculation============#
E1 = sum(start)/len(start)
E2 = sum(final)/len(final)
dE = E2 - E1   #kcal/mol
dE_dT  = dE/dT             # kcal/mol/K

Cv_Cal = dE_dT*kcal/mol/K/mass   # J/kgK
Cv_Cal_V = dE_dT*kcal/mol/K/mass*density # J/m3K
Cv_Ref = 966               # J/kgK Ref: https://www.periodic-table.org/copper-specific-heat/
#==============Output===============#
print()
print('Cv_Cal = ', Cv_Cal, ' J/kgK')
print('Cv_Cal_V = ', Cv_Cal_V, ' J/m3K')
print('Cv_Ref = ', Cv_Ref, ' J/kgK','#this ref stands for crystalline SiO2')
print()
with open('Cvout.dat','w') as fid:
    fid.write('Cv_Cal = '+format(Cv_Cal,'^.4f')+' J/kgK'+'\n')
    fid.write('Cv_Cal_V = '+format(Cv_Cal_V,'^.4f')+' J/m3K')
    
#==============Done===============#