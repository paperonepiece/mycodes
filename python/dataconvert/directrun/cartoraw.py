import random
import os
import math
import dpdata

def main():
    num = int(input("How many steps were run in AIMD: "))
    index = list(range(num))
    random.shuffle(index)
    train_set = index[:int(num*0.8)]
    valid_set = index[int(num*0.8):int(num*0.9)]
    test_set  = index[int(num*0.9):]
    #os.mkdir("Train")
    #os.mkdir("Valid")
    #os.mkdir("Test")
    d_outcar=dpdata.LabeledSystem('OUTCAR')
    d_outcar.sub_system(train_set).to_deepmd_raw('train_raw')
    d_outcar.sub_system(valid_set).to_deepmd_raw('valid_raw')
    d_outcar.sub_system(test_set).to_deepmd_raw('test_raw')
    
if __name__ == '__main__':
    main()