#!/usr/bin/python

"""Frequency response in an EM reverberation chamber
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

def FR(a,b,c,POS,R,f):
	#1/8
	Fx1,Fy1,Fz1 = FR8th(a,b,c,POS,R,f)
	2/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Fx2,Fy2,Fz2 = FR8th(a,b,c,POS,R,f)
	POS[:,2]=-POS[:,2]
 	3/8
	POS[:,1]=-POS[:,1]
	POS[:,4]=array([pi*ones((1,len(POS)))-POS[:,4]])
	Fx3,Fy3,Fz3 = FR8th(a,b,c,POS,R,f)
	4/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Fx4,Fy4,Fz4 = FR8th(a,b,c,POS,R,f)
	POS[:,1]=-POS[:,1]
	POS[:,2]=-POS[:,2]
	5/8
	POS[:,3]=array([POS[:,3]-ones((1,len(POS)))])
	POS[:,0]=-POS[:,0]
	Fx5,Fy5,Fz5 = FR8th(a,b,c,POS,R,f)
	6/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Fx6,Fy6,Fz6 = FR8th(a,b,c,POS,R,f)
	POS[:,2]=-POS[:,2]
	POS[:,4]=array([pi*ones((1,len(POS)))-POS[:,4]])
	7/8
	POS[:,1]=-POS[:,1]
	POS[:,5]=-POS[:,5]
	Fx7,Fy7,Fz7 = FR8th(a,b,c,POS,R,f)
	8/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Fx8,Fy8,Fz8 = FR8th(a,b,c,POS,R,f)
	POS[:,0]=-POS[:,0]
	POS[:,1]=-POS[:,1]
	POS[:,2]=-POS[:,2]
	POS[:,3]=array([-3.*ones((1,len(POS)))+POS[:,3]])
	POS[:,5]=mod(POS[:,5],2*pi)
	Fx=array([Fx1]+[Fx2]+[Fx3]+[Fx4]+[Fx5]+[Fx6]+[Fx7]+[Fx8]).sum(axis=0)
	Fy=array([Fy1]+[Fy2]+[Fy3]+[Fy4]+[Fy5]+[Fy6]+[Fy7]+[Fy8]).sum(axis=0)
	Fz=array([Fz1]+[Fz2]+[Fz3]+[Fz4]+[Fz5]+[Fz6]+[Fz7]+[Fz8]).sum(axis=0)
	return Fx,Fy,Fz


def FR8th(r,s,t,POS,R,f):
	Ex8th=zeros((len(r),len(f)),'complex')
	Ey8th=zeros((len(r),len(f)),'complex')
	Ez8th=zeros((len(r),len(f)),'complex')
	for i in range(0,len(r)):
		DX = r[i]-POS[:,0]
		DY = s[i]-POS[:,1]
		DZ = t[i]-POS[:,2]
		dist = sqrt(DX**2+DY**2+DZ**2)
		dp=tile(dist, (len(f),1))
		fp=tile(f,(len(dist),1))
		phase=2*pi*dp*fp.T/c
		ca    = cos(POS[:,4])
		sa    = sin(POS[:,4])
		cb    = cos(POS[:,5])
		sb    = sin(POS[:,5])
		distx = ((-sb)**2+(1-(-sb)**2)*ca)*DX+(-sb*cb*(1-ca))*DY+(cb*sa)*DZ
		disty = (-sb*cb*(1-ca))*DX+((cb)**2+(1-cb**2)*ca)*DY+(sb*sa)*DZ
		distz = (-cb*sa)*DX+(-sb*sa)*DY+ca*DZ
		DXY=sqrt(DX**2+DY**2)
		distxy = sqrt(distx**2+disty**2)
		costheta = distz/dist
		sintheta = distxy/dist
		cosphi   = distx/distxy
		sinphi   = disty/distxy
		L =tile(R**(POS[:,3]),(len(f),1))*1/dp #Amplitude & free space attenuation
		Ex8th[i,:] = sum(exp(1j*phase)*L*tile(((((-sb)**2+(1-(-sb)**2)*ca)*(-sintheta*costheta*cosphi)+(-sb*cb*(1-ca))*(-sintheta*costheta*sinphi)+(-cb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
		Ey8th[i,:] = sum(exp(1j*phase)*L*tile((((-sb*cb*(1-ca))*(-sintheta*costheta*cosphi)+((cb)**2+(1-(cb)**2)*ca)*(-sintheta*costheta*sinphi)+(-sb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
		Ez8th[i,:] = sum(exp(1j*phase)*L*tile((((cb*sa)*(-sintheta*costheta*cosphi)+(sb*sa)*(-sintheta*costheta*sinphi)+ca*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)           
	return Ex8th,Ey8th,Ez8th
