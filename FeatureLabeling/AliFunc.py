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