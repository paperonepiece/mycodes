# By Peng Hu

from ase.db import connect
import numpy as np
import os

db_name = 'sset_NiAl_10_energy.db'
db = connect(db_name)
element = {'Ni': 0, 'Al': 1}
selection = 'natoms<=7'
dir_name = '../data/'

for row in db.select(selection=selection):
    pos = np.array([row.positions.flatten()])
    box = np.array([row.cell.flatten()])
    ener = np.array([row.total_energy])
    t_map = row.symbols
    t = np.array([element[i] for i in t_map])

    # print("pos:", pos.shape)
    # print("box:", box.shape)
    # print("ener:", ener.shape)

    _dir_name = dir_name + str(row.id) + '/set.000'
    os.makedirs(_dirname)
    np.save(_dir_name + '/coord.npy', pos)
    np.save(_dir_name + '/box.npy', box)
    np.save(_dir_name + '/energy.npy', ener)
    with open(dir_name + str(row.id) + '/type_map.raw', 'w') as f:
        for k in element.keys():
            f.writelines(k + '\n')
    np.savetxt(dir_name + str(row.id) + '/type.raw', t, fmt="%d")
