##Generate random Zn node vacancy defects in ZIF & Not complete
##Only applicable for creating defect with atoms have fixed digit ordinal number(4 digit) 
# Control the distance between delected atoms to prevent ligand missing
import numpy as np
import random
import os

Ndrop = 16      # lines above atom position informations
totalnum = 2208 # total number of all types of atoms
abovenum = 2112 # the atoms which ordinal number smaller than defect atoms
Znnum    = 96   # total number of defect atoms before get deleted
belownum = totalnum-abovenum-Znnum  # atoms which ordinal numbers are bigger than defect atoms
lendig = 4
Zndis = float(6.5) # Angstorm
dellist = [] # number list of deleted atoms
poslist = [] # position parameters list for the deleted atoms

olddir = os.getcwd()
print("Now we are in:"+str(olddir))
workdir = input("Please input working directory, example like 'D:\\xx\\yy'：")
os.chdir(workdir)

#  defect concentration  0.05/0.1/0.15....
delnum = round(Znnum*eval(input("Please input the concentration of Zn node defect you need：")))
newnum = totalnum - delnum # total number of all types of atoms after defect generation
nodef = open("ZIF8.dat", "r+") 
defec = open("0.03ZIF8-2.dat", "w")

#  generate defect nodes with enough distance to prevent ligand missing 
while(len(dellist)<delnum):
    break_flag = False
    nodef.seek(0,0)
    x = random.randint(1,Znnum)
    if x not in dellist:
        for i in range(0,Ndrop+abovenum+x-1):
            nodef.readline()
        line = nodef.readline().split()
        #pos = float(line[3])**2+float(line[4])**2+float(line[5])**2
        if len(dellist) < 1:
            poslist.append(line[3])
            poslist.append(line[4])
            poslist.append(line[5])
            dellist.append(x)         
        else:
            for k in range(int(len(poslist)/3)):
                if abs((float(line[3])-float(poslist[3*k]))**2+
                       (float(line[4])-float(poslist[3*k+1]))**2+
                       (float(line[5])-float(poslist[3*k+2]))**2)**0.5 > Zndis:
                    continue
                else:
                    print('Close atoms were detected. The numbers are '
                          +str(x)+','+str(dellist[k]))
                    break_flag = True
                    break
            if break_flag == True:
                continue
            else:
                poslist.append(line[3])
                poslist.append(line[4])
                poslist.append(line[5])
                dellist.append(x) 
    else:
        continue
nodef.seek(0,0)

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
