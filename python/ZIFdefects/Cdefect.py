## Generate node vacancy defects for one type of atom
## in ZIF/complex structural material & Not complete 
## Control the distance between delected atoms for some purpose
## Delete nearest single-bonded atoms at the same time if you wish
import numpy as np
import random
import os


#  locating work path
olddir = os.getcwd()
print("Now we are in:"+str(olddir))
workdir = input("Please input working directory, example like 'D:\\xx\\yy'：")
os.chdir(workdir)

#  seting parameters
nodef = open("ZIF4_6type.dat", "r+") 
defec = open("ZIF4-C1C2-0.40.dat", "w")
allline = nodef.readlines()
Ndrop = 19      # lines above atom position informations 16/17/19 for ZIF
totalnum = int(allline[1].split()[0]) # total number of atoms
typenum  = int(allline[2].split()[0]) # total number of atoms type
setdis   = float(1.6) # Angstorm atoms within this distance will not be deleted
neardis  = float(1.3) # Angstorm neighbor atoms within this distance will be deleted
Clist  = [] # number list of deleted C atoms
Hlist  = []
Cplist = [] # position parameters list for the deleted C atoms
Hplist = []
posinf = allline[Ndrop:Ndrop+totalnum]

#  determine the number of each type of atoms
typelist = [] # list of atoms number of each type
typename = []
for i in range(1,totalnum):
    if int(posinf[i].split()[1]) == int(posinf[i-1].split()[1]):
        continue
    else:
        if len(typelist) <= 0:
            typelist.append(i)
            typename.append(int(posinf[i-1].split()[1]))
            typename.append(int(posinf[i].split()[1]))
        elif len(typelist) < typenum-1:
            typelist.append(i-sum(typelist[0:len(typelist)]))
            typename.append(int(posinf[i].split()[1]))
        if len(typelist) == typenum-1:
            typelist.append(totalnum-sum(typelist[0:len(typelist)]))
            break

#  determine the type and number of the deleted atoms
deltype  = int(input("Please input the type of defect atoms："))
nebtype  = int(input("Please input the type of neareast atoms around defect atoms(zero for none)："))
Cnum = typelist[typename.index(deltype)] # total number of defect atoms before get deleted
abovenum = sum(typelist[0:typename.index(deltype)]) # atoms with number smaller than defect atoms
belownum = totalnum-abovenum-Cnum # atoms with number bigger than defect atoms
if nebtype != 0:
    Hnum = typelist[typename.index(nebtype)]
    aboveneb = sum(typelist[0:typename.index(nebtype)])
    belowneb = totalnum-aboveneb-Hnum
    if typename.index(deltype) > typename.index(nebtype):
        betwenHC = (sum(typelist[typename.index(nebtype):typename.index(deltype)])-
                typelist[typename.index(nebtype)])
        sequence = [aboveneb,belowneb,Hnum,abovenum,belownum,Cnum] 
    else:
        betwenHC = (sum(typelist[typename.index(deltype):typename.index(nebtype)])-
                typelist[typename.index(deltype)])
        sequence = [abovenum,belownum,Cnum,aboveneb,belowneb,Hnum] 

#  defect concentration  0.05/0.1/0.15....
delnum = round(Cnum*eval(input("Please input the concentration of defect you need：")))

#  generate defect nodes with long enough distance for special situation
closelist = []
while(len(Clist)<delnum):
    break_flag = False
    nodef.seek(0,0)
    list1 = list(range(1,Cnum+1))
    for i in Clist + closelist:
        list1.remove(i)
    x = random.choice(list1)
    if x not in Clist:
        for i in range(0,Ndrop+abovenum+x-1):
            nodef.readline()
        line = nodef.readline().split()
        if len(Clist) < 1:
            Cplist.append(line[3])
            Cplist.append(line[4])
            Cplist.append(line[5])
            Clist.append(x)         
        else:
            for k in range(int(len(Cplist)/3)):
                if abs((float(line[3])-float(Cplist[3*k]))**2+
                    (float(line[4])-float(Cplist[3*k+1]))**2+
                    (float(line[5])-float(Cplist[3*k+2]))**2)**0.5 > setdis:
                    continue
                else:
                    print('Close atoms were detected. The numbers are '
                        +str(x)+','+str(Clist[k]))
                    closelist.append(x) 
                    break_flag = True
                    break
            if break_flag == True:
                continue
            else:
                Cplist.append(line[3])
                Cplist.append(line[4])
                Cplist.append(line[5])
                Clist.append(x)
    else:
        continue
nodef.seek(0,0)

# seek for neighbor H atoms near the deleted C atoms
if nebtype != 0:
    print('Some extra neighbor atoms need to be deleted, '+
        'but this module is not perfect so far')
    cremove = []
    for i in range(0,delnum):
        Cp1 = float(Cplist[3*i])
        Cp2 = float(Cplist[3*i+1])
        Cp3 = float(Cplist[3*i+2])
        nearcount = 0       
        for x in range(0,Hnum):
            Hp1 = float(allline[Ndrop+aboveneb+x].split()[3])
            Hp2 = float(allline[Ndrop+aboveneb+x].split()[4])
            Hp3 = float(allline[Ndrop+aboveneb+x].split()[5])
            if abs((Cp1-Hp1)**2+(Cp2-Hp2)**2+
                    (Cp3-Hp3)**2)**0.5 < neardis:
                nearcount += 1
                Hlist.append(x+1)
                Hplist.append(Hp1)
                Hplist.append(Hp2)
                Hplist.append(Hp3)
            else:
                if x == Hnum-1 and len(Hlist) < i+1-len(cremove):
                    print('The neareast atom of defect atom '+str(Clist[i])
                            +' is not found, probably due to periodic boundary!')
                    cremove.append(Clist[i])
                else:
                    continue
    for j in cremove:
        Clist.remove(j)
    for k in range(0,Ndrop):
        if k !=1:
            defec.write(allline[k])
        else:
            line = allline[k].replace(str(totalnum),str(totalnum-len(Clist)-len(Hlist)))
            defec.write(line)
            print('The real delnum is '+str(len(Clist))+', and the real neardel is '
                +str(len(Hlist))+'. The real newnum is '
                +str(totalnum-len(Clist)-len(Hlist)))
    for k in range(0,sequence[0]):
        defec.write(allline[Ndrop+k])
    linechange= np.zeros([1,4])  
    for k in range(sequence[0],sequence[0]+sequence[2]):
        linechange[0,0] += 1
        sp = allline[Ndrop+k].split()
        if int(linechange[0,0]) not in Hlist:
            line1 = sp[0].replace(str(k+1),str(k+1-int(linechange[0,1])))
            line2 = allline[Ndrop+k][len(sp[0]):]
            line3 = line1 + line2
            defec.write(line3)
        else:
            linechange[0,1] += 1
            continue        
    for k in range(sequence[0]+sequence[2],sequence[0]+sequence[2]+betwenHC):
        linechange[0,0] += 1
        sp2 = allline[Ndrop+k].split()
        line4 = sp2[0].replace(str(k+1),str(k+1-int(linechange[0,1])))
        line5 = allline[Ndrop+k][len(sp2[0]):]
        line6 = line4 + line5
        defec.write(line6)
    for k in range(sequence[3],sequence[3]+sequence[5]):
        linechange[0,2] += 1
        sp3 = allline[Ndrop+k].split()
        if int(linechange[0,2]) not in Clist:
            line7 = sp3[0].replace(str(k+1),str(k+1-int(linechange[0,1]+
                                                linechange[0,3])))
            line8 = allline[Ndrop+k][len(sp3[0]):]
            line9 = line7 + line8
            defec.write(line9)
        else:
            linechange[0,3] += 1
            continue        
    for k in range(totalnum-sequence[4],totalnum):
        sp4 = allline[Ndrop+k].split()
        line10 = sp4[0].replace(str(k+1),str(k+1-int(linechange[0,1]+
                                                linechange[0,3])))
        line11 = allline[Ndrop+k][len(sp4[0]):]
        line12 = line10 + line11
        defec.write(line12)
        
# no neighbor atoms need to be deleted
else:
    for k in range(0,Ndrop):
        line = nodef.readline()
        if k !=1:
            defec.write(line)
        else:
            line = line.replace(str(totalnum),str(totalnum-len(Clist)))
            defec.write(line)
    for k in range(0,abovenum):
        line = nodef.readline()
        l = line.split()
        defec.write(line)
    linenum = 0
    linechange = 0
    for k in range(abovenum+1,abovenum+Cnum+1):
        linenum += 1
        line = nodef.readline()
        sp = line.split()
        if linenum not in Clist:
            line1 = sp[0].replace(str(k),str(k-linechange))
            line2 = line[len(sp[0]):]
            line3 = line1 + line2
            defec.write(line3)
        else:
            linechange += 1
            continue
    for k in range(totalnum-belownum+1,totalnum+1):
        leftline = nodef.readline()
        sp1 = line.split()
        line4 = sp1[0].replace(str(k),str(k-linechange))
        line5 = line[len(sp1[0]):]
        line6 = line4 + line5
        defec.write(line6)
defec.close()
nodef.close()
print('Node defects generation has done, '+
    'defect concentration is '+str(format(delnum/Cnum*100,'.0f'))+'%')

