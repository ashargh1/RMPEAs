import numpy as np

def TH(T,c): #delS,delH,Omega,Eth,K1
    A=np.loadtxt('MixingEnthalpy.txt')
    B=np.loadtxt('InterMetallicEnthalpy.txt')
    R=8.3144598
    delS=0
    delH=0
    delHf=[]
    delHIM=0
    Omega=0
    
    for i in range(c.shape[0]):
        if c[i]!=0:
            delS=delS+c[i]*np.log(c[i])
    delS=-R*delS 
    
        
    for i in range(c.shape[0]-1): #
        for j in range(i+1,c.shape[0]):
            delH=delH+4*A[i,j]*c[i]*c[j]  
    Omega=(T*delS)/(1000*np.abs(delH)) #kJ into J convertion with 1000 
       
    for i in range(c.shape[0]-1):
        for j in range(i+1,c.shape[0]):
            if c[i]!=0 and c[j]!=0:
                delHf=np.append(delHf,B[i,j])
    Eth=(-0.8*T*delS)/(96.4916*np.abs(np.min(delHf))) #mev/atom into J/mol conversion with 96.4916
    
        
    for i in range(c.shape[0]-1):
        for j in range(i+1,c.shape[0]):
            delHIM=delHIM+4*B[i,j]*c[i]*c[j] 
            
          
    K1=(1+(0.4*T*delS/(1000*np.abs(delH))))/(96.4916*delHIM/(1000*delH))
    Output=[]
    Output=np.append(Output,[delS,delH,Omega,Eth,K1])
    Output=np.reshape(Output,(1,5))  
    return Output

def HR(c): #delta,E2/E0,VEC,delX
    A=np.loadtxt('Element.txt') #r VEC X Tm
    delta=0
    E=0
    VEC=0
    delX=0
    rbar=0
    Xbar=0
    for i in range(c.shape[0]):
        rbar=rbar+c[i]*A[i,0]
        Xbar=Xbar+c[i]*A[i,2]
        
    for i in range(c.shape[0]):
        delta=delta+c[i]*(np.power((1-(A[i,0]/rbar)),2))
        VEC=VEC+c[i]*A[i,1]
        delX=delX+c[i]*np.power((A[i,2]-Xbar),2)
    delta=(np.sqrt(delta))*100
    delX=np.sqrt(delX)
    
    for i in range(c.shape[0]):
        for j in range(i,c.shape[0]):
            E=E+((c[i]*c[j])*(np.power((A[i,0]+A[j,0]-2*rbar),2)))/(np.power(2*rbar,2))
        
    Output=[]
    Output=np.append(Output,[delta,E,VEC,delX])
    Output=np.reshape(Output,(1,4))  
    return Output

def PD(T,c): #PFP-A1,PFP-A2,PFP-A3,PFP-B2,PFP-Laves,PFP-Sigma,PSP with c=[0,0,0,0,0,0] for Ti Fe Al V Ni Nb Zr as binary
    Tpf=0.8*T # i.e. Phase formation temperature
    #print('Tpf=',Tpf)
    A1x=[]
    PFP_A1=0
    A2x=[]
    PFP_A2=0
    A3x=[]
    PFP_A3=0
    B2x=[]
    PFP_B2=0
    Lvx=[]
    PFP_Lv=0
    Sigmax=[]
    PFP_Sigma=0
    PSP=0
    
    if c[0] !=0 and c[1] !=0: #Ti-Fe only has A2, B2 and Laves phases
        In=np.loadtxt('BinPD/TiFe-B2.txt')
        B2x=np.append(B2x,IntSec(In,Tpf)[0])
        if np.shape(B2x)[0]==0:
            PFP_B2=0
        elif np.shape(B2x)[0]==2: 
            PFP_B2=np.abs(B2x[1]-B2x[0])
            
        In=np.loadtxt('BinPD/TiFe-A2.txt')
        C1=In[0:864][:] #First region
        A2x=np.append(A2x,IntSec(C1,Tpf)[0])
        C1=In[864:1479][:] #Second region
        A2x=np.append(A2x,IntSec(C1,Tpf)[0])   
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])
        elif np.shape(A2x)[0]==4:
            PFP_A2=np.abs(A2x[1]-A2x[0])+np.abs(A2x[3]-A2x[2])
            
        In=np.loadtxt('BinPD/TiFe-Laves.txt')
        Lvx=np.append(Lvx,IntSec(In,Tpf)[0])
        if np.shape(Lvx)[0]==0:
            PFP_Lv=0
        elif np.shape(Lvx)[0]==2: 
            PFP_Lv=np.abs(Lvx[1]-Lvx[0])
        
        
    if c[0] !=0 and c[2] !=0: #Ti-Al only has A1, A3, A2, Sigma and Laves phases
        In=np.loadtxt('BinPD/TiAl-Al.txt')
        A1x=np.append(A1x,IntSec(In,Tpf)[0])
        if np.shape(A1x)[0]==0:
            PFP_A1=0
        elif np.shape(A1x)[0]==2: 
            PFP_A1=np.abs(A1x[1]-A1x[0])
        
        In=np.loadtxt('BinPD/TiAl-Laves.txt')
        Lvx=np.append(Lvx,IntSec(In,Tpf)[0])
        if np.shape(Lvx)[0]==0:
            PFP_Lv=0
        elif np.shape(Lvx)[0]==2: 
            PFP_Lv=np.abs(Lvx[1]-Lvx[0])
        
        In=np.loadtxt('BinPD/TiAl-A2.txt')
        A2x=np.append(A2x,IntSec(In,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2: 
            PFP_A2=np.abs(A2x[1]-A2x[0])
           
        In=np.loadtxt('BinPD/TiAl-A3.txt') 
        if Tpf<=450:
            C1=In[0:656][:] #x<0.66
            A3x=np.append(A3x,IntSec(C1,Tpf)[0])
            C1=In[656:1133][:]
            A3x=np.append(A3x,IntSec(C1,Tpf)[0])
            if np.shape(A3x)[0]==4:
                PFP_A3=np.abs(A3x[1]-A3x[0])+np.abs(A3x[3]-A3x[2])
        elif Tpf>=845 and Tpf<=1140:
            C1=In[0:209][:] #x<0.28
            A3x=np.append(A3x,IntSec(C1,Tpf)[0])
            C1=In[210:1133][:]
            A3x=np.append(A3x,IntSec(C1,Tpf)[0])
            if np.shape(A3x)[0]==4:
                PFP_A3=np.abs(A3x[1]-A3x[0])+np.abs(A3x[3]-A3x[2])    
        else:
            A3x=np.append(A3x,IntSec(In,Tpf)[0])
            if np.shape(A3x)[0]==0:
                PFP_A3=0
            elif np.shape(A3x)[0]==2: 
                PFP_A3=np.abs(A3x[1]-A3x[0])
      
        In=np.loadtxt('BinPD/TiAl-Sigma.txt')
        Sigmax=np.append(Sigmax,IntSec(In,Tpf)[0])
        if np.shape(Sigmax)[0]==0:
            PFP_Sigma=0
        elif np.shape(Sigmax)[0]==2: 
            PFP_Sigma=np.abs(Sigmax[1]-Sigmax[0])    
        
        
    if c[0] !=0 and c[3] !=0: #Ti-V only has A2 and A3 phases
        In=np.loadtxt('BinPD/TiV-A2.txt')
        if Tpf>=1864 and Tpf<=1940:
            C1=In[0:878][:] #x<0.63
            A2x=np.append(A2x,IntSec(C1,Tpf)[0])
            C1=In[879:1283][:]
            A2x=np.append(A2x,IntSec(C1,Tpf)[0])
            if np.shape(A2x)[0]==4:
                PFP_A2=np.abs(A2x[1]-A2x[0])+np.abs(A2x[3]-A2x[2])        
        else:
            A2x=np.append(A2x,IntSec(In,Tpf)[0])
            if np.shape(A2x)[0]==0:
                PFP_A2=0
            elif np.shape(A2x)[0]==2: 
                PFP_A2=np.abs(A2x[1]-A2x[0])
                
        In=np.loadtxt('BinPD/TiV-A3.txt')
        A3x=np.append(A3x,IntSec(In,Tpf)[0])
        if np.shape(A3x)[0]==0:
            PFP_A3=0
        elif np.shape(A3x)[0]==2: 
            PFP_A3=np.abs(A3x[1]-A3x[0]) 
        
        
    if c[0] !=0 and c[4] !=0: #Ti-Ni only has A1, A2 and B2 phases
        In=np.loadtxt('BinPD/TiNi-B2.txt')
        B2x=np.append(B2x,IntSec(In,Tpf)[0])
        if np.shape(B2x)[0]==0:
            PFP_B2=0
        elif np.shape(B2x)[0]==2: 
            PFP_B2=np.abs(B2x[1]-B2x[0]) 
            
        In=np.loadtxt('BinPD/TiNi-A2.txt')
        A2x=np.append(A2x,IntSec(In,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2: 
            PFP_A2=np.abs(A2x[1]-A2x[0])     
  
        In=np.loadtxt('BinPD/TiNi-A1.txt') 
        if Tpf<=580:
            C1=In[0:638][:] #x<0.011
            A1x=np.append(A1x,IntSec(C1,Tpf)[0])
            C1=In[638:1671][:]
            A1x=np.append(A1x,IntSec(C1,Tpf)[0])
            if np.shape(A1x)[0]==4:
                PFP_A1=np.abs(A1x[1]-A1x[0])+np.abs(A1x[3]-A1x[2])        
        else:
            A1x=np.append(A1x,IntSec(In,Tpf)[0])
            if np.shape(A1x)[0]==0:
                PFP_A3=0
            elif np.shape(A1x)[0]==2: 
                PFP_A1=np.abs(A1x[1]-A1x[0])

        
    if c[0] !=0 and c[5] !=0: #Ti-Nb only has A2 and A3 phases
        In=np.loadtxt('BinPD/TiNb-A2.txt')
        A2x=np.append(A2x,IntSec(In,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])
            
        In=np.loadtxt('BinPD/TiNb-A3.txt')
        A3x=np.append(A3x,IntSec(In,Tpf)[0])
        if np.shape(A3x)[0]==0:
            PFP_A3=0
        elif np.shape(A3x)[0]==2:
            PFP_A3=np.abs(A3x[1]-A3x[0])
     
    
    if c[0] !=0 and c[6] !=0: #Ti-Zr only has A2 and A3 phases
        In=np.loadtxt('BinPD/TiZr-A2.txt')
        if Tpf>=1827 and Tpf<=1940:
            C1=In[0:545][:] #x<0.66
            A2x=np.append(A2x,IntSec(C1,Tpf)[0])
            C1=In[545:927][:]
            A2x=np.append(A2x,IntSec(C1,Tpf)[0])
        else:
            A2x=np.append(A2x,IntSec(In,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])
        elif np.shape(A2x)[0]==4:
            PFP_A2=np.abs(A2x[1]-A2x[0])+np.abs(A2x[3]-A2x[2])
            
        In=np.loadtxt('BinPD/TiZr-A3.txt')
        if Tpf>=874 and Tpf<=1142:
            C1=In[0:407][:] #x<0.5
            A3x=np.append(A3x,IntSec(C1,Tpf)[0])
            C1=In[407:819][:]
            A3x=np.append(A3x,IntSec(C1,Tpf)[0])
        else:
            A3x=np.append(A3x,IntSec(In,Tpf)[0])
        if np.shape(A3x)[0]==0:
            PFP_A3=0
        elif np.shape(A3x)[0]==2:
            PFP_A3=np.abs(A3x[1]-A3x[0])
        elif np.shape(A3x)[0]==4:
            PFP_A3=np.abs(A3x[1]-A3x[0])+np.abs(A3x[3]-A3x[2])
            
            
    if c[1] !=0 and c[2] !=0: #Fe-Al only has A1, A2 and B2 phases
        In=np.loadtxt('BinPD/FeAl-B2.txt')
        B2x=np.append(B2x,IntSec(In,Tpf)[0])
        if np.shape(B2x)[0]==0:
            PFP_B2=0
        elif np.shape(B2x)[0]==2:
            PFP_B2=np.abs(B2x[1]-B2x[0])
            
        In=np.loadtxt('BinPD/FeAl-A2.txt')
        A2x=np.append(A2x,IntSec(In,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])    
              
        In=np.loadtxt('BinPD/FeAl-A1.txt')
        A1x=np.append(A1x,IntSec(In,Tpf)[0])
        if np.shape(A1x)[0]==0:
            PFP_A1=0
        elif np.shape(A1x)[0]==2:
            PFP_A1=np.abs(A1x[1]-A1x[0]) 

        
    if c[1] !=0 and c[3] !=0: #Fe-V only has A1, Sigma and A2 phases
        In=np.loadtxt('BinPD/FeV-A1.txt')
        A1x=np.append(A1x,IntSec(In,Tpf)[0])
        if np.shape(A1x)[0]==0:
            PFP_A1=0
        elif np.shape(A1x)[0]==2:
            PFP_A1=np.abs(A1x[1]-A1x[0])
        
        In=np.loadtxt('BinPD/FeV-Sigma.txt')
        Sigmax=np.append(Sigmax,IntSec(In,Tpf)[0])
        if np.shape(Sigmax)[0]==0:
            PFP_Sigma=0
        elif np.shape(Sigmax)[0]==2:
            PFP_Sigma=np.abs(Sigmax[1]-Sigmax[0]) 
            
        In=np.loadtxt('BinPD/FeV-A2.txt') 
        if Tpf<=1535:
            C1=In[0:1101][:] #x<0.52
            A2x=np.append(A2x,IntSec(C1,Tpf)[0])
            C1=In[1101:1968][:]
            A2x=np.append(A2x,IntSec(C1,Tpf)[0])
        elif Tpf>=1740 and Tpf<=1811:
            C1=In[0:1124][:] #x<0.6
            A2x=np.append(A2x,IntSec(C1,Tpf)[0])
            C1=In[1124:1968][:]
            A2x=np.append(A2x,IntSec(C1,Tpf)[0])
        else:
            A2x=np.append(A2x,IntSec(In,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])
        elif np.shape(A2x)[0]==4:
            PFP_A2=np.abs(A2x[1]-A2x[0])+np.abs(A2x[3]-A2x[2])

                    
    if c[1] !=0 and c[4] !=0: #Fe-Ni only has A1 and A2 phases
        In=np.loadtxt('BinPD/FeNi-A1.txt')
        A1x=np.append(A1x,IntSec(In,Tpf)[0])
        if np.shape(A1x)[0]==0:
            PFP_A1=0
        elif np.shape(A1x)[0]==2:
            PFP_A1=np.abs(A1x[1]-A1x[0])
            
        In=np.loadtxt('BinPD/FeNi-A2.txt')
        A2x=np.append(A2x,IntSec(In,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])
    
    
    if c[1] !=0 and c[5] !=0: #Fe-Nb only has A1, Laves and A2 phases
        In=np.loadtxt('BinPD/FeNb-A1.txt')
        A1x=np.append(A1x,IntSec(In,Tpf)[0])
        if np.shape(A1x)[0]==0:
            PFP_A1=0
        elif np.shape(A1x)[0]==2:
            PFP_A1=np.abs(A1x[1]-A1x[0])
        
        In=np.loadtxt('BinPD/FeNb-Laves.txt')
        Lvx=np.append(Lvx,IntSec(In,Tpf)[0])
        if np.shape(Lvx)[0]==0:
            PFP_Lv=0
        elif np.shape(Lvx)[0]==2:
            PFP_Lv=np.abs(Lvx[1]-Lvx[0])
              
        In=np.loadtxt('BinPD/FeNb-A2.txt')
        C1=In[0:1169][:] #First region
        A2x=np.append(A2x,IntSec(C1,Tpf)[0])
        C1=In[1169:1982][:] #Second region
        A2x=np.append(A2x,IntSec(C1,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])
        elif np.shape(A2x)[0]==4:
            PFP_A2=np.abs(A2x[1]-A2x[0])+np.abs(A2x[3]-A2x[2]) 
        
    if c[1] !=0 and c[6] !=0: #Fe-Zr only has A1, Laves and A2 phases
        In=np.loadtxt('BinPD/FeZr-A1.txt')
        A1x=np.append(A1x,IntSec(In,Tpf)[0])
        if np.shape(A1x)[0]==0:
            PFP_A1=0
        elif np.shape(A1x)[0]==2:
            PFP_A1=np.abs(A1x[1]-A1x[0])
        
        In=np.loadtxt('BinPD/FeZr-Laves.txt')
        Lvx=np.append(Lvx,IntSec(In,Tpf)[0])
        if np.shape(Lvx)[0]==0:
            PFP_Lv=0
        elif np.shape(Lvx)[0]==2:
            PFP_Lv=np.abs(Lvx[1]-Lvx[0])
        
        In=np.loadtxt('BinPD/FeZr-A2.txt')
        C1=In[0:649][:] #First region
        A2x=np.append(A2x,IntSec(C1,Tpf)[0])
        C1=In[649:777][:] #Second region
        A2x=np.append(A2x,IntSec(C1,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])
        elif np.shape(A2x)[0]==4:
            PFP_A2=np.abs(A2x[1]-A2x[0])+np.abs(A2x[3]-A2x[2]) 
    
    
    if c[2] !=0 and c[3] !=0: #Al-V only has A1 and A2 phases
        In=np.loadtxt('BinPD/AlV-A1.txt')
        A1x=np.append(A1x,IntSec(In,Tpf)[0])
        if np.shape(A1x)[0]==0:
            PFP_A1=0
        elif np.shape(A1x)[0]==2:
            PFP_A1=np.abs(A1x[1]-A1x[0])
        
        In=np.loadtxt('BinPD/AlV-A2.txt')
        A2x=np.append(A2x,IntSec(In,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])

    
    if c[2] !=0 and c[4] !=0: #Al-Ni only has A1 and B2 phases
        In=np.loadtxt('BinPD/AlNi-B2.txt')
        B2x=np.append(B2x,IntSec(In,Tpf)[0])
        if np.shape(B2x)[0]==0:
            PFP_B2=0
        elif np.shape(B2x)[0]==2:
            PFP_B2=np.abs(B2x[1]-B2x[0])
        
        In=np.loadtxt('BinPD/AlNi-A1.txt')
        C1=In[0:890][:] #First region
        A1x=np.append(A1x,IntSec(C1,Tpf)[0])
        C1=In[890:1062][:] #Second region
        A1x=np.append(A1x,IntSec(C1,Tpf)[0])
        if np.shape(A1x)[0]==0:
            PFP_A1=0
        elif np.shape(A1x)[0]==2:
            PFP_A1=np.abs(A1x[1]-A1x[0])
        elif np.shape(A1x)[0]==4:
            PFP_A1=np.abs(A1x[1]-A1x[0])+np.abs(A1x[3]-A1x[2])            
            
        
    if c[2] !=0 and c[5] !=0: #Al-Nb only has A1, Sigma, Laves, A2 and B2 phases
        In=np.loadtxt('BinPD/AlNb-A1.txt')
        A1x=np.append(A1x,IntSec(In,Tpf)[0])
        if np.shape(A1x)[0]==0:
            PFP_A1=0
        elif np.shape(A1x)[0]==2:
            PFP_A1=np.abs(A1x[1]-A1x[0])  
            
        In=np.loadtxt('BinPD/AlNb-Sigma.txt')
        Sigmax=np.append(Sigmax,IntSec(In,Tpf)[0])
        if np.shape(Sigmax)[0]==0:
            PFP_Sigma=0
        elif np.shape(Sigmax)[0]==2:
            PFP_Sigma=np.abs(Sigmax[1]-Sigmax[0]) 
            
        In=np.loadtxt('BinPD/AlNb-Laves.txt')
        Lvx=np.append(Lvx,IntSec(In,Tpf)[0])
        if np.shape(Lvx)[0]==0:
            PFP_Lv=0
        elif np.shape(Lvx)[0]==2:
            PFP_Lv=np.abs(Lvx[1]-Lvx[0]) 
        
        In=np.loadtxt('BinPD/AlNb-B2.txt')
        B2x=np.append(B2x,IntSec(In,Tpf)[0])
        if np.shape(B2x)[0]==0:
            PFP_B2=0
        elif np.shape(B2x)[0]==2:
            PFP_B2=np.abs(B2x[1]-B2x[0]) 
        
        In=np.loadtxt('BinPD/AlNb-A2.txt')
        A2x=np.append(A2x,IntSec(In,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0]) 
                        
            
    if c[2] !=0 and c[6] !=0: #Al-Zr only has A1, A3, laves and A2 phases
        In=np.loadtxt('BinPD/AlZr-A2.txt')
        A2x=np.append(A2x,IntSec(In,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])
        
        In=np.loadtxt('BinPD/AlZr-A1.txt')
        A1x=np.append(A1x,IntSec(In,Tpf)[0])
        if np.shape(A1x)[0]==0:
            PFP_A1=0
        elif np.shape(A1x)[0]==2:
            PFP_A1=np.abs(A1x[1]-A1x[0])
            
        In=np.loadtxt('BinPD/AlZr-A3.txt')
        A3x=np.append(A3x,IntSec(In,Tpf)[0])
        if np.shape(A3x)[0]==0:
            PFP_A3=0
        elif np.shape(A3x)[0]==2:
            PFP_A3=np.abs(A3x[1]-A3x[0])
        
        In=np.loadtxt('BinPD/AlZr-Laves.txt')
        Lvx=np.append(Lvx,IntSec(In,Tpf)[0])
        if np.shape(Lvx)[0]==0:
            PFP_Lv=0
        elif np.shape(Lvx)[0]==2:
            PFP_Lv=np.abs(Lvx[1]-Lvx[0])
        
        
    if c[3] !=0 and c[4] !=0: #V-Ni only has A1, Sigma and A2 phases
        In=np.loadtxt('BinPD/VNi-A2.txt')
        A2x=np.append(A2x,IntSec(In,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])
        
        In=np.loadtxt('BinPD/VNi-A1.txt')
        A1x=np.append(A1x,IntSec(In,Tpf)[0])
        if np.shape(A1x)[0]==0:
            PFP_A1=0
        elif np.shape(A1x)[0]==2:
            PFP_A1=np.abs(A1x[1]-A1x[0])
        
        In=np.loadtxt('BinPD/VNi-Sigma.txt')
        Sigmax=np.append(Sigmax,IntSec(In,Tpf)[0])
        if np.shape(Sigmax)[0]==0:
            PFP_Sigma=0
        elif np.shape(Sigmax)[0]==2:
            PFP_Sigma=np.abs(Sigmax[1]-Sigmax[0])  
       
    
    if c[3] !=0 and c[5] !=0: #V-Nb only has A2 phase            
        In=np.loadtxt('BinPD/VNb-A2.txt') 
        if Tpf>=2130 and Tpf<=2188:
            C1=In[0:716][:] #x<0.75
            A2x=np.append(A2x,IntSec(C1,Tpf)[0])
            C1=In[716:1412][:]
            A2x=np.append(A2x,IntSec(C1,Tpf)[0])
        else:
            A2x=np.append(A2x,IntSec(In,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])
        elif np.shape(A2x)[0]==4:
            PFP_A2=np.abs(A2x[1]-A2x[0])+np.abs(A2x[3]-A2x[2])
        
        
    if c[3] !=0 and c[6] !=0: #V-Zr only has A3, Laves and A2 phases
        In=np.loadtxt('BinPD/VZr-A3.txt')
        A3x=np.append(A3x,IntSec(In,Tpf)[0])
        if np.shape(A3x)[0]==0:
            PFP_A3=0
        elif np.shape(A3x)[0]==2:
            PFP_A3=np.abs(A3x[1]-A3x[0])
        
        In=np.loadtxt('BinPD/VZr-Laves.txt')
        Lvx=np.append(Lvx,IntSec(In,Tpf)[0])
        if np.shape(Lvx)[0]==0:
            PFP_Lv=0
        elif np.shape(Lvx)[0]==2:
            PFP_Lv=np.abs(Lvx[1]-Lvx[0])
        
        In=np.loadtxt('BinPD/VZr-A2.txt')
        C1=In[0:621][:] #First region
        A2x=np.append(A2x,IntSec(C1,Tpf)[0])
        C1=In[621:1761][:] #Second region
        A2x=np.append(A2x,IntSec(C1,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])
        elif np.shape(A2x)[0]==4:
            PFP_A2=np.abs(A2x[1]-A2x[0])+np.abs(A2x[3]-A2x[2])
        
        
    if c[4] !=0 and c[5] !=0: #Ni-Nb only has A1, A3, A2 and Laves phase
        In=np.loadtxt('BinPD/NiNb-A2.txt')
        A2x=np.append(A2x,IntSec(In,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])
        
        In=np.loadtxt('BinPD/NiNb-A1.txt')
        C1=In[0:106][:] #First region
        A1x=np.append(A1x,IntSec(C1,Tpf)[0])
        C1=In[106:1005][:] #Second region
        A1x=np.append(A1x,IntSec(C1,Tpf)[0])
        if np.shape(A1x)[0]==0:
            PFP_A1=0
        elif np.shape(A1x)[0]==2:
            PFP_A1=np.abs(A1x[1]-A1x[0])
        elif np.shape(A1x)[0]==4:
            PFP_A1=np.abs(A1x[1]-A1x[0])+np.abs(A1x[3]-A1x[2])
            
        In=np.loadtxt('BinPD/NiNb-Laves.txt')
        Lvx=np.append(Lvx,IntSec(In,Tpf)[0])
        if np.shape(Lvx)[0]==0:
            PFP_Lv=0
        elif np.shape(Lvx)[0]==2:
            PFP_Lv=np.abs(Lvx[1]-Lvx[0])

        In=np.loadtxt('BinPD/NiNb-A3.txt') 
        if Tpf>=375 and Tpf<=395:
            C1=In[0:371][:] #x<0.65
            A3x=np.append(A3x,IntSec(C1,Tpf)[0])
            C1=In[371:450][:]
            A3x=np.append(A3x,IntSec(C1,Tpf)[0])
        else:
            A3x=np.append(A3x,IntSec(In,Tpf)[0])
        if np.shape(A3x)[0]==0:
            PFP_A3=0
        elif np.shape(A3x)[0]==2:
            PFP_A3=np.abs(A3x[1]-A3x[0])
        elif np.shape(A3x)[0]==4:
            PFP_A3=np.abs(A3x[1]-A3x[0])+np.abs(A3x[3]-A3x[2])
    
        
    if c[4] !=0 and c[6] !=0: #Ni-Zr only has A1, A2 and B2 phase
        In=np.loadtxt('BinPD/NiZr-A1.txt')
        A1x=np.append(A1x,IntSec(In,Tpf)[0])
        if np.shape(A1x)[0]==0:
            PFP_A1=0
        elif np.shape(A1x)[0]==2:
            PFP_A1=np.abs(A1x[1]-A1x[0])
        
        In=np.loadtxt('BinPD/NiZr-B2.txt')
        B2x=np.append(B2x,IntSec(In,Tpf)[0])
        if np.shape(B2x)[0]==0:
            PFP_B2=0
        elif np.shape(B2x)[0]==2:
            PFP_B2=np.abs(B2x[1]-B2x[0])
            
        In=np.loadtxt('BinPD/NiZr-A2.txt')
        A2x=np.append(A2x,IntSec(In,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])
    
        
    if c[5] !=0 and c[6] !=0: #Nb-Zr only has A3 and A2 phase
        In=np.loadtxt('BinPD/NbZr-A3.txt')
        A3x=np.append(A3x,IntSec(In,Tpf)[0])
        if np.shape(A3x)[0]==0:
            PFP_A3=0
        elif np.shape(A3x)[0]==2:
            PFP_A3=np.abs(A3x[1]-A3x[0])
        
        In=np.loadtxt('BinPD/NbZr-A2.txt') 
        if Tpf>=2013 and Tpf<=2132:
            C1=In[0:464][:] #x<0.22
            A2x=np.append(A2x,IntSec(C1,Tpf)[0])
            C1=In[464:1378][:]
            A2x=np.append(A2x,IntSec(C1,Tpf)[0])
        else:
            A2x=np.append(A2x,IntSec(In,Tpf)[0])
        if np.shape(A2x)[0]==0:
            PFP_A2=0
        elif np.shape(A2x)[0]==2:
            PFP_A2=np.abs(A2x[1]-A2x[0])
        elif np.shape(A2x)[0]==4:
            PFP_A2=np.abs(A2x[1]-A2x[0])+np.abs(A2x[3]-A2x[2])
            
        
    PSP=1-PFP_A1-PFP_A2-PFP_A3-PFP_B2-PFP_Lv-PFP_Sigma
    return [PFP_A1,PFP_A2,PFP_A3,PFP_B2,PFP_Lv,PFP_Sigma,PSP]
        

def IntSec(C1,Tpf): #It only gives 0 or 2 points 
    Xint=[] #Initial output
    Yint=[] 
    Xint0=[] #Final output
    Yint0=[]
    y1=Tpf*np.ones(np.shape(C1)[0])
    x=C1[:,1]
    y=C1[:,2]
    A=y1 - y
    for i in range(np.shape(A)[0]):
        if np.abs(A[i])<=10: #Adjustable parameter
            Xint=np.append(Xint,x[i])
            Yint=np.append(Yint,y[i])
    if np.shape(Xint)[0] !=0: 
        A0=np.abs(Yint-Tpf)
        Xint0=np.append(Xint0,Xint[np.argmin(A0)]) #First intersection point
        Yint0=np.append(Yint0,Yint[np.argmin(A0)])
        Yint=np.delete(Yint,np.argmin(A0))
        Xint=np.delete(Xint,np.argmin(A0))
        
        while np.shape(Xint0)[0] <2 and np.shape(Xint)[0] !=0 :
            A0=np.abs(Yint-Tpf)
            if np.abs(Xint[np.argmin(A0)]-Xint0) <0.000001: 
                Xint=np.delete(Xint,np.argmin(A0))
                Yint=np.delete(Yint,np.argmin(A0))
            else:
                Xint0=np.append(Xint0,Xint[np.argmin(A0)])
                Yint0=np.append(Yint0,Yint[np.argmin(A0)])
        if np.shape(Xint0)[0] !=2:
            Xint0=[]
            Yint0=[]
    return [Xint0,Yint0]

def PFP_PSP(T,c):
    PFP_T_A1=0
    PFP_D_A1=0
    PFP_T_A2=0
    PFP_D_A2=0
    PFP_T_A3=0
    PFP_D_A3=0
    PFP_T_B2=0
    PFP_D_B2=0
    PFP_T_Lv=0
    PFP_D_Lv=0
    PFP_T_Sigma=0
    PFP_D_Sigma=0
    PSP_T=0
    PSP_D=0

    for i in range(np.shape(c)[0]):
        for j in range(i+1,np.shape(c)[0]):
            if c[i] !=0 and c[j] !=0:  
                a=[]
                for m in range(np.shape(c)[0]):
                    a=np.append(a,0)
                a[i]=1
                a[j]=1

                PFP_T_A1=PFP_T_A1+(PD(T,a)[0])*(c[i]*c[j])
                PFP_D_A1=PFP_D_A1+c[i]*c[j]
                PFP_T_A2=PFP_T_A2+(PD(T,a)[1])*(c[i]*c[j])
                PFP_D_A2=PFP_D_A2+c[i]*c[j]
                PFP_T_A3=PFP_T_A3+(PD(T,a)[2])*(c[i]*c[j])
                PFP_D_A3=PFP_D_A3+c[i]*c[j]
                PFP_T_B2=PFP_T_B2+(PD(T,a)[3])*(c[i]*c[j])
                PFP_D_B2=PFP_D_B2+c[i]*c[j]
                PFP_T_Lv=PFP_T_Lv+(PD(T,a)[4])*(c[i]*c[j])
                PFP_D_Lv=PFP_D_Lv+c[i]*c[j]
                PFP_T_Sigma=PFP_T_Sigma+(PD(T,a)[5])*(c[i]*c[j])
                PFP_D_Sigma=PFP_D_Sigma+c[i]*c[j]
                PSP_T=PSP_T+(PD(T,a)[6])*(c[i]*c[j])
                PSP_D=PSP_D+(1-PD(T,a)[6])*(c[i]*c[j])


    if PFP_D_A1 !=0:
        PFP_A1=PFP_T_A1/PFP_D_A1
    else:
        PFP_A1=0

    if PFP_D_A2 !=0:
        PFP_A2=PFP_T_A2/PFP_D_A2
    else:
        PFP_A2=0

    if PFP_D_A3 !=0:
        PFP_A3=PFP_T_A3/PFP_D_A3
    else:
        PFP_A3=0

    if PFP_D_B2 !=0:
        PFP_B2=PFP_T_B2/PFP_D_B2
    else:
        PFP_B2=0

    if PFP_D_Lv !=0:
        PFP_Lv=PFP_T_Lv/PFP_D_Lv
    else:
        PFP_Lv=0

    if PFP_D_Sigma !=0:
        PFP_Sigma=PFP_T_Sigma/PFP_D_Sigma
    else:
        PFP_Sigma=0

    if PSP_D !=0:
        PSP=PSP_T/PSP_D
    else:
        PSP=0

    Output=[]
    Output=np.append(Output,[PFP_A1,PFP_A2,PFP_A3,PFP_B2,PFP_Lv,PFP_Sigma,PSP])
    Output=np.reshape(Output,(1,7)) 
    
    return Output 