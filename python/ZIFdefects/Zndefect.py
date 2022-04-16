##Generate random Zn node vacancy defects in ZIF & Not complete
##Only applicable for creating defect with atoms have certain digit ordinal number(4 digit) 
import numpy as np
import random
import os

Ndrop = 16      # lines above atom position informations
totalnum = 2208 # total number of all types of atoms
abovenum = 2112 # the atoms which ordinal number smaller than defect atoms
Znnum    = 96   # total number of defect atoms before get deleted
belownum = totalnum-abovenum-Znnum  # atoms which ordinal number bigger than defect atoms
lendig = 4
dellist = [] # number list of deleted atoms

olddir = os.getcwd()
print("Now we are in:"+str(olddir))
workdir = input("Please input working directory, example like 'D:\\xx\\yy'：")
os.chdir(workdir)

#  defect concentration  0.05/0.1/0.15....
delnum = round(Znnum*eval(input("Please input the concentration of Zn node defect you need：")))
newnum = totalnum - delnum # total number of all types of atoms after defect generation
nodef = open("ZIF8.dat", "r+") 
defec = open("0.15ZIF8-2.dat", "w")

while(len(dellist)<delnum):
    x = random.randint(1,Znnum)
    if x not in dellist:
        dellist.append(x)
        
for k in range(0,Ndrop):
    line = nodef.readline()
    if k !=1:
        defec.write(line)
    else:
        line = line.replace(str(totalnum),str(newnum))
        defec.write(line)
for k in range(0,abovenum):
    line = nodef.readline()
    l = line.split()
    defec.write(line)

linenum = 0
linechange = 0
for k in range(abovenum+1,abovenum+Znnum+1):
    linenum += 1
    line = nodef.readline()
    if linenum not in dellist:
        line1 = line[0:lendig].replace(str(k),str(k-linechange))  ##four digit number
        line2 = line[lendig:]                                     ##change it as your wish
        line3 = line1 + line2
        defec.write(line3)
    else:
        linechange += 1
        continue
for k in range(totalnum-belownum+1,totalnum+1):
    leftline = nodef.readline()
    leftline = leftline.replace(str(k),str(k-linechange))
    defec.write(leftline)
defec.close()
nodef.close()
print('Zn node defects generation has done, '+
      'defect concentration is '+str(format(delnum/Znnum*100,'.0f'))+'%')



# with open('ZIF80.05.dat','r+') as fid:
#     next(fid)
#     linelist = fid.readline()
#     for k in range(0,Ndrop):    
#     stringlist = fid.readlines()
#     linenum = 0
#     for line in stringlist:
#         linelist = line.split()
#         if linelist[1] == 4:
