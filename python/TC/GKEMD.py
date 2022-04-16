#GK—EMD method for TC caculation
#引入计算库
import numpy as np
from scipy import integrate,fft,optimize
import os

#参数设置
Nf = int(1e5)
Nt = int(1e4)
N = int(Nf//Nt)
T = 300      
dt = 1e-15*10
kb = 1.38064852e-23
Vo =78935.105*1e-30
Ndrop = 2
Ncase = 20
interval = 1
dump_t = Nt/interval

#定义空矩阵，KK、KV同理
K  = np.zeros([Ncase,Nt])
Kx = np.zeros([Ncase,Nt])
Ky = np.zeros([Ncase,Nt])
Kz = np.zeros([Ncase,Nt])
# KK = np.zeros([Ncase,Nt]) 对流项
# KKx= np.zeros([Ncase,Nt])
# KKy= np.zeros([Ncase,Nt])
# KKz= np.zeros([Ncase,Nt])
J  = np.zeros([3,Nf])
JK = np.zeros([3,Nf])
J1 = np.zeros([1,Nt])                   
Dt = (np.arange(1,Nt+1)*dt) 
CJ0= np.zeros([1,2*Nt-1])
CJ = np.zeros([3,Nt]) 
CJK= np.zeros([3,Nt])

#显示当前路径并修改工作路径
olddir = os.getcwd()
print("Now we are in:"+str(olddir))
workdir = input("Please input working directory, example like 'D:\\xx\\yy'：")
#直接复制工作文件夹路径地址--粘贴--回车
os.chdir(workdir)

for i in range(1,Ncase+1):
    os.chdir(str(i))
    with open('J.dat','r+') as fid:
        for k in range(1,Ndrop+1):
            next(fid)
        stringlist = fid.readlines()
        linenum = 0
        for line in stringlist:
            linelist = line.split()
            J[0,linenum] = linelist[2]
            J[1,linenum] = linelist[3]
            J[2,linenum] = linelist[4]
#             JK[0,linenum]= linelist[5]
#             JK[1,linenum]= linelist[6]
#             JK[2,linenum]= linelist[7]
            linenum += 1
            if linenum >= Nf:
                break
    os.chdir(workdir)
    for j in range(1,N+1):
        J1  = (J[0,((j-1)*Nt):j*Nt])*1.6e-17/Vo#*0.04345*1000 可直接用于metal
        J2  = (J[1,((j-1)*Nt):j*Nt])*1.6e-17/Vo#*0.04345*1000 real需要取消注释
        J3  = (J[2,((j-1)*Nt):j*Nt])*1.6e-17/Vo#*0.04345*1000
#         JK1 = (JK[0,((j-1)*Nt):j*Nt])*1.6e-17/Vo*0.04345*1000
#         JK2 = (JK[1,((j-1)*Nt):j*Nt])*1.6e-17/Vo*0.04345*1000
#         JK3 = (JK[2,((j-1)*Nt):j*Nt])*1.6e-17/Vo*0.04345*1000
        
#互相关+无偏化计算
        CJ0 = np.correlate(J1,J1, "full")
        CJ0_un = np.concatenate([CJ0[0:Nt]/(len(J1)+np.arange((-Nt+1),1,1)),
                                 CJ0[Nt:(2*Nt)]/(len(J1)-np.arange(1,Nt,1))])
        CJ1 = np.correlate(J2,J2, "full")
        CJ1_un = np.concatenate([CJ1[0:Nt]/(len(J2)+np.arange((-Nt+1),1,1)),
                                 CJ1[Nt:(2*Nt)]/(len(J2)-np.arange(1,Nt,1))])
        CJ2 = np.correlate(J3,J3, "full")
        CJ2_un = np.concatenate([CJ2[0:Nt]/(len(J3)+np.arange((-Nt+1),1,1)),
                                 CJ2[Nt:(2*Nt)]/(len(J3)-np.arange(1,Nt,1))])
#         CJK0 = np.correlate(JK1,JK1, "full")
#         CJK0_un = np.concatenate([CJK0[0:Nt]/(len(JK1)+np.arange((-Nt+1),1,1)),
#                                   CJK0[Nt:(2*Nt)]/(len(JK1)-np.arange(1,Nt,1))])
#         CJK1 = np.correlate(JK2,JK2, "full")
#         CJK1_un = np.concatenate([CJK1[0:Nt]/(len(JK2)+np.arange((-Nt+1),1,1)),
#                                   CJK1[Nt:(2*Nt)]/(len(JK2)-np.arange(1,Nt,1))])
#         CJK2 = np.correlate(JK3,JK3, "full")
#         CJK2_un = np.concatenate([CJK2[0:Nt]/(len(JK3)+np.arange((-Nt+1),1,1)),
#                                   CJK2[Nt:(2*Nt)]/(len(JK3)-np.arange(1,Nt,1))])
        CJ[0,:] = CJ0_un[Nt-1:2*Nt]
        CJ[1,:] = CJ1_un[Nt-1:2*Nt]
        CJ[2,:] = CJ2_un[Nt-1:2*Nt]
#         CJK[0,:] = CJK0_un[Nt-1:2*Nt]
#         CJK[1,:] = CJK1_un[Nt-1:2*Nt]
#         CJK[2,:] = CJK2_un[Nt-1:2*Nt]
        
#梯形法求积分
        K[i-1,:] = K[i-1,:] +(integrate.cumtrapz(CJ[0,:], Dt, initial=0)
                             +(integrate.cumtrapz(CJ[1,:], Dt, initial=0))
                             +(integrate.cumtrapz(CJ[2,:], Dt, initial=0)))*Vo/kb/T/T/N/3
        Kx[i-1,:]= Kx[i-1,:] +(integrate.cumtrapz(CJ[0,:], Dt, initial=0))*Vo/kb/T/T/N
        Ky[i-1,:]= Ky[i-1,:] +(integrate.cumtrapz(CJ[1,:], Dt, initial=0))*Vo/kb/T/T/N
        Kz[i-1,:]= Kz[i-1,:] +(integrate.cumtrapz(CJ[2,:], Dt, initial=0))*Vo/kb/T/T/N
#         KK[i-1,:]= KK[i-1,:] +(integrate.cumtrapz(CJK[0,:], Dt, initial=0)
#                               +(integrate.cumtrapz(CJK[1,:], Dt, initial=0))
#                               +(integrate.cumtrapz(CJK[2,:], Dt, initial=0)))*Vo/kb/T/T/N/3 
#         KKx[i-1,:]= KKx[i-1,:] +(integrate.cumtrapz(CJK[0,:], Dt, initial=0))*Vo/kb/T/T/N
#         KKy[i-1,:]= KKy[i-1,:] +(integrate.cumtrapz(CJK[1,:], Dt, initial=0))*Vo/kb/T/T/N
#         KKz[i-1,:]= KKz[i-1,:] +(integrate.cumtrapz(CJK[2,:], Dt, initial=0))*Vo/kb/T/T/N
#平均值
AveK   = np.mean(K, axis=0)
AveKx  = np.mean(Kx, axis=0)
AveKy  = np.mean(Ky, axis=0)
AveKz  = np.mean(Kz, axis=0)
# AveKK  = np.mean(KK, axis=0)
# AveKKx = np.mean(KKx, axis=0)
# AveKKy = np.mean(KKy, axis=0)
# AveKKz = np.mean(KKz, axis=0)
StdK   = np.std(K, axis=0)
# StdKK  = np.std(KK, axis=0)
Tc     = np.mean(AveK[int(0.2*Nt):int(0.8*Nt)])
Tcx    = np.mean(AveKx[int(0.2*Nt):int(0.8*Nt)])
Tcy    = np.mean(AveKy[int(0.2*Nt):int(0.8*Nt)])
Tcz    = np.mean(AveKz[int(0.2*Nt):int(0.8*Nt)])
# TcK    = np.mean(AveKK[int(0.2*Nt):int(0.8*Nt)])
# TcKx   = np.mean(AveKKx[int(0.2*Nt):int(0.8*Nt)])
# TcKy   = np.mean(AveKKy[int(0.2*Nt):int(0.8*Nt)])
# TcKz   = np.mean(AveKKz[int(0.2*Nt):int(0.8*Nt)])
with open('K_correlation.dat','w') as fidw,\
     open('pyTC.dat','w') as fidtc:
    for m in range(1,int(dump_t+1)):
        fidw.write(format(m*interval*dt*1e12,'<10.5f'))
        fidw.write(format(AveK[(m-1)*interval],'^12.6f')+
#                    format(AveKK[(m-1)*interval],'^12.6f')+
                   format(StdK[(m-1)*interval],'^12.6f'))
        fidw.write('\n')
    fidtc.write('Tc='+format(Tc,'^.4f')+'\n')
    fidtc.write('Tcx='+format(Tcx,'^.4f')+'\n')
    fidtc.write('Tcy='+format(Tcy,'^.4f')+'\n')
    fidtc.write('Tcz='+format(Tcz,'^.4f')+'\n')
#     fidtc.write('TcK='+format(TcK,'^.4f')+'\n')
#     fidtc.write('TcKx='+format(TcKx,'^.4f')+'\n')
#     fidtc.write('TcKy='+format(TcKy,'^.4f')+'\n')
#     fidtc.write('TcKz='+format(TcKz,'^.4f')+'\n')
print("GKEMD TC caculation is done!")
#---------------------Done-------------------------
