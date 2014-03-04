#!/usr/bin/python

"""Image Creator EM reverberation chamber
Visit http://http://projets.ietr.fr/imagetheorymodel/
"""
from __future__ import division

__author__ = "Emmanuel Amador (emmanuel.amador@insa-rennes.fr)"
__version__ = "$Revision: 0.2 $"
__date__ = "$Date: 2014/03/02$"
__copyright__ = "Copyright (c) 2014 E. Amador"
__license__ = "NDA"


from numpy import *

c = 299792458.0

def IC(Tm,l,p,h,X,Y,Z,tilt,azimut): #Computes images' postions and image's orientation
	dims = array([l,p,h])
	Tp = Tm+3./c*dims.max()
	dmax=c*Tp
	order = round(dmax/dims.min())
	POSp = array([X,Y,Z, 0, tilt, azimut]);
 
	#1D
	for i in range(1,int(order)):
		POSp = vstack([POSp,[2*i*l-X, Y, Z, abs(2*i)-1, pi-tilt, 2*pi-azimut],[2*i*l+X, Y, Z, abs(2*i), tilt, azimut]])
	POSi=POSp.copy();
	POSi[:,1] = array([p*ones((1,len(POSp)))-POSi[:,1]])
	POSi[:,4] = array([pi*ones((1,len(POSp)))-POSi[:,4]])
	POSi[:,5] = array([pi*ones((1,len(POSp)))-POSi[:,5]])
	
	#2D
	POSpp = POSp.copy();
	POSpp[:,1] = array([p*ones((1,len(POSp)))]+POSp[:,1])
	POSpp[:,3] = array([ones((1,len(POSp)))]+POSp[:,3])
	POSii=POSi.copy()
	POSii[:,1] = array([p*ones((1,len(POSi)))]+POSi[:,1])
	POSii[:,3] = array([ones((1,len(POSi)))]+POSi[:,3])
	POS1 = POSpp.copy()
	POS1 = concatenate((POS1,POSii),axis=0)
	dist = (POS1[:,0]**2+POS1[:,1]**2+POS1[:,2]**2)**.5
	POS1 = concatenate((POS1,dist.reshape(-1,1)),axis=1)
	POS1 = POS1[POS1[:,6].argsort(),]
	if POS1[:,6].max()>c*Tp:
		U = where(POS1[:,6]>c*Tp)
		POS1=delete(POS1, s_[min(min(U)):], axis=0)
	
	POS1 = delete(POS1, s_[6], axis=1)
	POSP=POS1.copy()
	
	for jj in range (2,int(order),2):
		POSpp = POSp.copy()
		POSpp[:,1] = array([jj*p*ones((1,len(POSp)))]+POSp[:,1])
		POSpp[:,3] = array([abs(jj)*ones((1,len(POSp)))]+POSp[:,3])
		POSii=POSi.copy()
		POSii[:,1] = array([(jj+1)*p*ones((1,len(POSi)))]+POSi[:,1])
		POSii[:,3] = array([abs(jj+1)*ones((1,len(POSi)))]+POSi[:,3])
		POS2 = POSpp.copy()
		POS2 = concatenate((POS2,POSii),axis=0)
		dist = (POS2[:,0]**2+POS2[:,1]**2+POS2[:,2]**2)**.5
		POS2 = concatenate((POS2,dist.reshape(-1,1)),axis=1)
		POS2 = POS2[POS2[:,6].argsort(),]
		if POS2[:,6].max()>c*Tp:
			V = where(POS2[:,6]>c*Tp)
			POS2=delete(POS2, s_[min(min(V)):], axis=0)
		POS2 = delete(POS2, s_[6], axis=1)
		POSP = concatenate((POSP,POS2),axis=0)
		
	
	POSI = POSP.copy();
	POSI[:,2] = array([h*ones((1,len(POSI)))-POSI[:,2]])
	POSI[:,5] = mod(array([pi*ones((1,len(POSI)))+POSI[:,5]]),2*pi)
	
	#3D
	POSPP = POSP.copy()
	POSII=POSI.copy()
	POSII[:,2] = array([h*ones((1,len(POSI)))]+POSI[:,2])
	POSII[:,3] = array([ones((1,len(POSI)))]+POSI[:,3])
	POS3 = POSPP.copy()
	POS3 = concatenate((POS3,POSII),axis=0)
	dist = (POS3[:,0]**2+POS3[:,1]**2+POS3[:,2]**2)**.5
	POS3 = concatenate((POS3,dist.reshape(-1,1)),axis=1)
	POS3 = POS3[POS3[:,6].argsort(),]
	if POS3[:,6].max()>c*Tp:
		W = where(POS3[:,6]>c*Tp)
		POS3=delete(POS3, s_[min(min(W)):], axis=0)
	
	POS3 = delete(POS3, s_[6], axis=1)
	POS=POS3.copy()		
	
	for k in range (2,int(order),2):
		POSPP = POSP.copy()
		POSPP[:,2] = array([k*h*ones((1,len(POSP)))]+POSP[:,2])
		POSPP[:,3] = array([(k)*ones((1,len(POSP)))]+POSP[:,3])
		POSII=POSI.copy()
		POSII[:,2] = array([(k+1)*h*ones((1,len(POSI)))]+POSI[:,2])
		POSII[:,3] = array([(k+1)*ones((1,len(POSI)))]+POSI[:,3])
		POS4 = POSPP.copy()
		POS4 = concatenate((POS4,POSII),axis=0)
		dist = (POS4[:,0]**2+POS4[:,1]**2+POS4[:,2]**2)**.5
		POS4 = concatenate((POS4,dist.reshape(-1,1)),axis=1)
		POS4 = POS4[POS4[:,6].argsort(),]
		if POS4[:,6].max()>c*Tp:
			Q = where(POS4[:,6]>c*Tp)
			POS4=delete(POS4, s_[min(min(Q)):], axis=0)
		POS4 = delete(POS4, s_[6], axis=1)
		POS = concatenate((POS,POS4),axis=0)
	
	POS=vstack((POSp[0,:],POS))
	return POS

if __name__ == '__main__':
    #Image Creation
    Lt=1e-6 #length of the time window in s
    
    #dimensions of the reverb chamber in m
    l=8.7   
    p=3.7
    h=2.9
    
    #position of the emitter and angular orientation	
    X = 1.
    Y = 2.
    Z = 1.
    tilt = pi/2-math.acos(sqrt(2./3));
    azimut = pi/4;
    print('Lt= %2.3f mus simulation in a  %2.3f x %2.3f x %2.3f m3 chamber\nEmitter at (%2.3f,%2.3f,%2.3f)\nCreating images...' %(Lt/1e-6,l,p,h,X,Y,Z))    
    POS=IC(Lt,l,p,h,X,Y,Z,tilt,azimut)
    
    savez('POS.npz',Lt=Lt,POS=POS,l=l,p=p,h=h,X=X,Y=Y,Z=Z)    


