#!/usr/bin/python

"""Channel impulse response in an EM reverberation chamber
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

def CIR(a,b,c,POS,R): #change the POS matrix t get the whole amount of images
	#1/8
	Sx1,Sy1,Sz1,t1 = CIR8th(a,b,c,POS,R)
	2/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Sx2,Sy2,Sz2,t2 = CIR8th(a,b,c,POS,R)
	POS[:,2]=-POS[:,2]
 	#3/8
	POS[:,1]=-POS[:,1]
	POS[:,4]=array([pi*ones((1,len(POS)))-POS[:,4]])
	Sx3,Sy3,Sz3,t3 = CIR8th(a,b,c,POS,R)
	#4/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Sx4,Sy4,Sz4,t4 = CIR8th(a,b,c,POS,R)
	POS[:,1]=-POS[:,1]
	POS[:,2]=-POS[:,2]
	#5/8
	POS[:,3]=array([POS[:,3]-ones((1,len(POS)))])
	POS[:,0]=-POS[:,0]
	Sx5,Sy5,Sz5,t5 = CIR8th(a,b,c,POS,R)
	#6/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Sx6,Sy6,Sz6,t6 = CIR8th(a,b,c,POS,R)
	POS[:,2]=-POS[:,2]
	POS[:,4]=array([pi*ones((1,len(POS)))-POS[:,4]])
	#7/8
	POS[:,1]=-POS[:,1]
	POS[:,5]=-POS[:,5]
	Sx7,Sy7,Sz7,t7 = CIR8th(a,b,c,POS,R)
	#8/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Sx8,Sy8,Sz8,t8 = CIR8th(a,b,c,POS,R)
	POS[:,0]=-POS[:,0]
	POS[:,1]=-POS[:,1]
	POS[:,2]=-POS[:,2]
	POS[:,3]=array([-3.*ones((1,len(POS)))+POS[:,3]])
	POS[:,5]=mod(POS[:,5],2*pi)
	Sx=concatenate((Sx1,Sx2,Sx3,Sx4,Sx5,Sx6,Sx7,Sx8),axis=1)
	Sy=concatenate((Sy1,Sy2,Sy3,Sy4,Sy5,Sy6,Sy7,Sy8),axis=1)
	Sz=concatenate((Sz1,Sz2,Sz3,Sz4,Sz5,Sz6,Sz7,Sz8),axis=1)
	t=concatenate((t1,t2,t3,t4,t5,t6,t7,t8),axis=1)
	return Sx,Sy,Sz,t

def CIR8th(q,s,u,POS,R):  #computes the E-field
	Sx8th=zeros((len(q),len(POS)))
	Sy8th=zeros((len(q),len(POS)))
	Sz8th=zeros((len(q),len(POS)))
	delay=zeros((len(q),len(POS)))
	for i in range(0,len(q)):
		DX = q[i]-POS[:,0]
		DY = s[i]-POS[:,1]
		DZ = u[i]-POS[:,2]
		r = sqrt(DX**2+DY**2+DZ**2)
		delay[i,:]=r/c
		ca    = cos(POS[:,4])
		sa    = sin(POS[:,4])
		cb    = cos(POS[:,5])
		sb    = sin(POS[:,5])
		rx = ((-sb)**2+(1-(-sb)**2)*ca)*DX+(-sb*cb*(1-ca))*DY+(cb*sa)*DZ
		ry = (-sb*cb*(1-ca))*DX+((cb)**2+(1-cb**2)*ca)*DY+(sb*sa)*DZ
		rz = (-cb*sa)*DX+(-sb*sa)*DY+ca*DZ
		rxy = sqrt(rx**2+ry**2)
		costheta = rz/r
		sintheta = rxy/r
		cosphi   = rx/rxy
		sinphi   = ry/rxy
		L =R**(POS[:,3])*1/r      
		Sx8th[i,:] = L*(((-sb)**2+(1-(-sb)**2)*ca)*(-sintheta*costheta*cosphi)+(-sb*cb*(1-ca))*(-sintheta*costheta*sinphi)+(-cb*sa)*(-sintheta*(-sintheta)));
		Sy8th[i,:] = L*((-sb*cb*(1-ca))*(-sintheta*costheta*cosphi)+((cb)**2+(1-(cb)**2)*ca)*(-sintheta*costheta*sinphi)+(-sb*sa)*(-sintheta*(-sintheta)));
		Sz8th[i,:] = L*((cb*sa)*(-sintheta*costheta*cosphi)+(sb*sa)*(-sintheta*costheta*sinphi)+ca*(-sintheta*(-sintheta)));
	return Sx8th,Sy8th,Sz8th,delay
