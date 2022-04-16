
########################################################################################
##this code bases on the old grammar of python2##
##be careful when using python3, some modifications are needed##

g     = 1.0e-3           # kg
A     = 1.0e-10          # m
nm    = 1.0e-9           # m
um    = 1.0e-6           # m
m     = 1.0              # m
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
K     = 1.0 # K

########################################################################################
num = 20      #the number of lines for average calculation of L1 and L2
L0list = []
start  = []
final  = []
#============Reading data===========#
fid1   = open("thermo_L0.dat","r+")
fid2   = open("thermo_CLTE.dat","r+")
L0line = fid1.readlines()
for i in range(int(0.5*len(L0line)),len(L0line)):
    L0list.append(float(L0line[i].split()[1]))  #use the latter half of Lz data to calculate L0
fid1.close()
linelist = fid2.readlines()
for i in range(2,num+2):
    start.append(float(linelist[i].split()[2]))
for i in range(len(linelist)-num,len(linelist)):
    final.append(float(linelist[i].split()[2]))
fid2.close()
#============Calculation============#
L0    = sum(L0list)/len(L0list)     # A
L1    = sum(start)/len(start)       # A
L2    = sum(final)/len(final)       # A
dL    = L2 - L1  #A
dL_dT = dL/dT    # A/K

CTE_Cal = 1.0/L0*dL_dT/(um/m/K) # um/m/K
CTE_Ref = 0.5                   # um/m/K Ref: www.periodic-table.org/copper-thermal-expansion/
#==============Output===============#
print()
print('CTE_Cal = ', CTE_Cal, ' um/m/K')
print('CTE_Ref = ', CTE_Ref, ' um/m/K','#this ref stands for crystalline SiO2')
print()
with open('CTEout.dat','w') as fid:
    fid.write('CTE_Cal = '+format(CTE_Cal,'^.4f')+' 1e-6/K'+'\n')
#==============Done===============#