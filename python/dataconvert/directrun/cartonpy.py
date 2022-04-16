import random
import os
import math
import dpdata

def main():
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
    
if __name__ == '__main__':
    main()