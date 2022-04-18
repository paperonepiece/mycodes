## Written by Zheyi.
##E-Mail: Zhangzheyi@mail.dlut.edu.cn
## For converting different datafile to lammps datafile
import os
import random
import dpdata
import numpy as np
from ase import io
from ase.io import read,write


def hit_dat(hitfile,output):
#=======================convert .hit to lammps.dat=======================#
#==============================Basic Setting=============================#
    #hitfile= "amorph.hit"       ## inputfile name
    #output = "lammpsdata1.txt"  ## outputfile name
    typenum = []                ## type number and sign
    typenam = []                ## sign of all atoms with all types
    atominfo= []                ## all info of atoms
    cell = []                   ## cell size info
    diag = []
    upper= []
#===========================Cell and Atom infos==========================#
    hit = open(hitfile, "r+")
    while True:
        info = hit.readline()
        if 'ni-cell' in info:       
            break                   
    for i in range(3):
        info = hit.readline().split()
        cell.append(float(info[0]))
        cell.append(float(info[1]))
        cell.append(float(info[2]))
        diag.append(float(cell[4*i]))
    [upper.append(x) for x in cell if not x in diag]
    upper = np.array(upper)
    while True:
        info = hit.readline()
        if 'otal number' in info:
            natom = int(hit.readline())
            break
    while True:
        info = hit.readline()
        if 'tom site' in info:
            break
    for i in range(natom):
        line = hit.readline()
        if not line:
            print('Severe Error! Some Atoms info Not Found!')
            break
        info = line.split()
        if info[0] not in typenam:
            typenum.append(info[0])
        atominfo.append(info)    
        typenam.append(info[0])
    hit.close()
#==============================Output Files=============================#        
    if not upper.any():
        dat = open(output, "w")
        dat.write('# LAMMPS data file written by python\n')
        dat.write('%d atoms\n' %(natom))
        dat.write('%d atom types\n' %(len(typenum)))
        dat.write('0.0 %f xlo xhi\n0.0 %f ylo yhi\n0.0 %f zlo zhi\n'
                %(diag[0], diag[1], diag[2]))
        dat.write('%f %f %f xy xz yz\n' %(cell[1], cell[2], cell[5]))
        dat.write('\n'+'Atoms  # atomic\n\n')
        for i in range(natom):
            dat.write('%d %d %s %s %s  # %s\n' %(i+1, typenum.index(typenam[i])+1,
                    atominfo[i][1], atominfo[i][2], atominfo[i][3], typenam[i] ))
        dat.close()
    else:
        print('Need extra calculation for cell size conversion\n')
        print('Using ase to convert .vasp file\n')
        vasp = open("convert.vasp", "w")
        for i in range(len(typenum)):
            vasp.write(typenum[i]+' ')
        vasp.write('\n')
        vasp.write('  '+'1.0\n')  ##scale factor
        vasp.write('    %f  %f  %f\n' %(cell[0], cell[1], cell[2]))
        vasp.write('    %f  %f  %f\n' %(cell[3], cell[4], cell[5]))
        vasp.write('    %f  %f  %f\n' %(cell[6], cell[7], cell[8]))
        for i in range(len(typenum)):
            vasp.write('  '+typenum[i])
        vasp.write('\n')
        for i in range(len(typenum)):
            vasp.write('  '+str(typenam.count(typenum[i])))
        vasp.write('\n')
        vasp.write('Cartesian'+'\n')
        for i in range(len(typenum)):
            for j in range(natom):
                if typenum[i] == atominfo[j][0]:
                    vasp.write('   %s   %s   %s  # %s\n'
                    %(atominfo[j][1], atominfo[j][2], atominfo[j][3], typenum[i]))
        vasp.close()
        atoms = read('convert.vasp')
        atoms.write(output,format='lammps-data')
        os.remove('convert.vasp')
    print('Congrats! Data conversion has done ~')
#===============================Done=================================#

def car_npy(num):
#==================Convert outcar to npy for dpmd====================#
    num = 0
    f = open('OUTCAR','r+')
    dat = f.readlines()
    for line in dat:
        if 'ELECTRON-ION-THERMOSTAT' in line.split('\n')[0]:
            num += 1
    index = list(range(num))
    random.shuffle(index)
    train_set = index[:int(num*0.8)]    ##change the ratio as u want
    valid_set = index[int(num*0.8):int(num*0.9)]
    test_set  = index[int(num*0.9):]
    d_outcar=dpdata.LabeledSystem('OUTCAR')
    d_outcar.sub_system(train_set).to_deepmd_npy('data_train')
    d_outcar.sub_system(valid_set).to_deepmd_npy('data_valid')
    d_outcar.sub_system(test_set).to_deepmd_npy('data_test')
    print('Congrats! OUTCAR extraction has done ~')
#===============================Done=================================#
    
def gala_dat(galaname, dataname, full=True, manybody=True):
#==================Convert galamost to lammps.dat====================#
#==============================Basic info============================#
    fxml = open(galaname,'r+')
    if full:
        charge = False
    if manybody:
        bond = False
        angle= False
        dihe = False
        impo = False
    while manybody:
        info = fxml.readline()
        if 'natoms' in info:
            natoms = int(info.split('"')[-2])
            break
    while manybody:
        info = fxml.readline()
        if 'box' in info:
            lx = float(info.split('"')[1])*10 # nm to A
            ly = float(info.split('"')[3])*10 # nm
            lz = float(info.split('"')[5])*10 # nm
            break
#=========================Atom's type and positions==================#
    while manybody:
        info = fxml.readline()
        if 'position' in info:
            break
    for i in range(natoms):
        info = fxml.readline()
        if not info:
            print('Severe Error! Some Atoms info Not Found!')
            break
        info = info.split()
        line = np.array([[10*float(x) for x in info]])   # nm to A
        if i < 1:
            positions = line 
        else:
            positions = np.append(positions,line,axis=0)   

    index1  = []
    atoms = []
    atomtype = []
    atomtypenum = 0
    while manybody:
        info = fxml.readline()
        if 'type' in info:
            break
    for i in range(natoms):
        line = fxml.readline().split()
        info = line
        if not info:
            print('Severe Error! Some Atoms info Not Found!')
            break
        if info[0] not in atomtype:
            atomtypenum += 1
            index1.append(i)
            info.append(str(atomtypenum))
            info.append(str(i+1))
        else:
            info.append(atoms[atomtype.index(info[0])][1])
            info.append(str(i+1))
        atomtype.append(info[0])
        atoms.append(info)
#=====================Velocities masses and charges==================#
    while True:
        info = fxml.readline()
        if 'velocity' in info:
            vel=True
            break
        if not info:
            vel=False
            fxml.seek(0,0)
            break
    if vel:    
        for i in range(natoms):
            info = fxml.readline().split()
            line = np.array([[0.01*float(x) for x in info]]) # nm/ps to A/fs--real
            if i < 1:
                velocities = line 
            else:
                velocities = np.append(velocities,line,axis=0)
    else:print("Velocity Info Not Found! Velocities will not appear in datafile")

    masses = []
    while manybody:
        info = fxml.readline()
        if 'mass' in info:
            mass = True
            break
        if not info:
            fxml.seek(0,0)
            mass = False
            break
    if mass:
        for i in range(natoms):
            info = float(fxml.readline())
            if i in index1:
                masses.append(info)

    if full:
        charges = []
        while manybody:
            info = fxml.readline()
            if 'charge' in info:
                charge = True
                break
            if not info:
                fxml.seek(0,0)
                break
        if charge:
            for i in range(natoms):
                info = float(fxml.readline())
                charges.append(info)
        else:
            print("Charge Info Not Found! Charges will not appear in datafile")
            full = False
#=========================Bond and manybody params===================#
    if manybody:
        while manybody:
            info = fxml.readline()
            if 'bond' in info:
                bondnum = int(info.split('"')[-2])
                bond = True
                break
            if not info:
                fxml.seek(0,0)
                break
        if bond:
            bonds = []
            bondtype = []
            bondtypenum = 0
            for i in range(bondnum):
                line = fxml.readline().split()
                info = [line[0]]  
                if info[0] not in bondtype:
                    bondtypenum +=1
                    info.append(str(bondtypenum))
                else:
                    info.append(bonds[bondtype.index(info[0])][1])
                length = len(line)-1
                for l in range(length):
                    info.append(str(int(line[l-length])+1))
                bonds.append(info)
                bondtype.append(info[0])
        else:print("Bond Info Not Found! Bonds will not appear in datafile")
            
    if manybody:
        while manybody:
            info = fxml.readline()
            if 'angle' in info:
                anglenum = int(info.split('"')[-2])
                angle = True
                break
            if not info:
                fxml.seek(0,0)
                break
        if angle:
            angles = []
            angletype = []
            angletypenum = 0
            for i in range(anglenum):
                line = fxml.readline().split()
                info = [line[0]]
                if info[0] not in angletype:
                    angletypenum +=1
                    info.append(str(angletypenum))
                else:
                    info.append(angles[angletype.index(info[0])][1])
                length = len(line)-1
                for l in range(length):
                    info.append(str(int(line[l-length])+1))
                angles.append(info)
                angletype.append(info[0])
        else:print("Angle Info Not Found! Angles will not appear in datafile")

    if manybody:
        while manybody:
            info = fxml.readline()
            if 'dihedral' in info:
                dihenum = int(info.split('"')[-2])
                dihe = True
                break
            if not info:
                fxml.seek(0,0)
                break
        if dihe:
            dihes = []
            dihetype = []
            dihetypenum = 0
            for i in range(dihenum):
                line = fxml.readline().split()
                info = [line[0]]
                if info[0] not in dihetype:
                    dihetypenum +=1
                    info.append(str(dihetypenum))
                else:
                    info.append(dihes[dihetype.index(info[0])][1])
                length = len(line)-1
                for l in range(length):
                    info.append(str(int(line[l-length])+1))
                dihes.append(info)
                dihetype.append(info[0])
        else:print("Dihedral Info Not Found! Dihedrals will not appear in datafile")
    if manybody:
        while manybody:
            info = fxml.readline()
            if 'impro' in info:
                imponum = int(info.split('"')[-2])
                impo = True
                break
            if not info:
                fxml.seek(0,0)
                break
        if impo:
            pass
        else:print("Improper Info Not Found! Impropers will not appear in datafile")
    fxml.close()
    print("Reading info from galamost file has done.")
#=============================Write datafile=======================#
    fdat = open(dataname,'w')
    fdat.write('Lammps data file generated by python code from ZhangZheyi\n')
    fdat.write('%d atoms\n' %(natoms))
    if bond:
        fdat.write('%d bonds\n' %(bondnum))
    if angle:
        fdat.write('%d angles\n' %(anglenum))
    if dihe:
        fdat.write('%d dihedrals\n' %(dihenum))
    fdat.write('\n')
    fdat.write('%d atom types\n' %(atomtypenum))
    if bond:
        fdat.write('%d bond types\n' %(bondtypenum))
    if angle:
        fdat.write('%d angle types\n' %(angletypenum))
    if dihe:
        fdat.write('%d dihedral types\n' %(dihetypenum))
    fdat.write('\n')
    fdat.write('%f  %f xlo xhi\n%f  %f ylo yhi\n%f  %f zlo zhi\n' %
            (-lx/2, lx/2, -ly/2, ly/2, -lz/2, lz/2))
    fdat.write('0.000 0.000 0.000 xy xz yz\n\n')
    if masses:
        fdat.write('Masses\n\n')
        for i in range(len(masses)):
            fdat.write('%d %f # %s\n' %(i+1,masses[i],atoms[index1[i]][0]))
        fdat.write('\n')
    if full:
        molecule = 1
        fdat.write('Atoms # full\n\n')
        for i in range(natoms):
            fdat.write('%d %d %s %f %f %f %f # %s\n' % (i+1,molecule,
                        atoms[i][1],charges[i],positions[i][0],positions[i][1],
                        positions[i][2],atomtype[i]))
    else:
        fdat.write('Atoms # atomic\n\n')
        for i in range(natoms):
            fdat.write('%d %s %f %f %f # %s\n' % (i+1,atoms[i][1],
                    positions[i][0],positions[i][1],positions[i][2],atomtype[i]))
    if vel:
        fdat.write('\nVelocities\n\n')
        for i in range(natoms):
            fdat.write('%d %f %f %f # %s\n' % (i+1,
                velocities[i][0],velocities[i][1],velocities[i][2],atomtype[i]))    
    if manybody:
        if bond:
            fdat.write('\nBonds\n\n')
            for i in range(bondnum):
                fdat.write('%d %s %s %s # %s\n' % (i+1,bonds[i][1],
                            bonds[i][2],bonds[i][3],bonds[i][0]))
        if angle:
            fdat.write('\nAngles\n\n')
            for i in range(anglenum):
                fdat.write('%d %s %s %s %s # %s\n' % (i+1,angles[i][1],
                            angles[i][2],angles[i][3],angles[i][4],angles[i][0]))
        if dihe:
            fdat.write('\nDihedrals\n\n')
            for i in range(dihenum):
                fdat.write('%d %s %s %s %s %s # %s\n' % (i+1,dihes[i][1],
                            dihes[i][2],dihes[i][3],dihes[i][4],dihes[i][5],dihes[i][0]))
        if impo:pass
    fdat.close()
    print("Writing data to lammps-datafile has done.")
#==============================Done==============================#  