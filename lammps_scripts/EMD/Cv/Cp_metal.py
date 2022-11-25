
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
num = 100       #the number of lines for average calculation of E1 and E2
natoms1 = 1728
natoms2 = 3456
mass   = 28.09*g/mol*natoms1+16*g/mol*natoms2 # kg
volume = 78935*A*A*A  # m3
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
dE = E2 - E1   #ev
dE_dT  = dE/dT             # eV/K

Cp_Cal = dE_dT*eV/K/mass   # J/kgK
Cp_Cal_V = dE_dT*eV/K/mass*density # J/m3K
Cp_Ref = 966               # J/kgK Ref: https://www.periodic-table.org/copper-specific-heat/
#==============Output===============#
print()
print('Cp_Cal = ', Cp_Cal, ' J/kgK')
print('Cp_Cal_V = ', Cp_Cal_V, ' J/kgK')
print('Cp_Ref = ', Cp_Ref, ' J/kgK','#this ref stands for crystalline SiO2')
print()
with open('Cp.dat','w') as fid:
    fid.write('Cv_Cal = '+format(Cv_Cal,'^.4f')+' J/kgK'+'\n')
    fid.write('Cv_Cal_V = '+format(Cv_Cal_V,'^.4f')+' J/m3K')
