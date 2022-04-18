## Written by Zheyi
## For converting hit datafile to lammps datafile
import numpy as np
from ase.io import read,write
import os

def main():
#============================Basic Setting===========================#
    hitfile= "amorph_50_1.hit"       ## inputfile name
    output = "lammpsdata1.txt"  ## outputfile name
    typenum = []                ## type number and sign
    typenam = []                ## sign of all atoms with all types
    atominfo= []                ## all info of atoms
    cell = []                   ## cell size info
    diag = []
    upper= []
#=========================Cell and Atom infos========================#
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
#============================Output Files===========================#        
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
#===============================Done================================#
if __name__ == '__main__':
    main()