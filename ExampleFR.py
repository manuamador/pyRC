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
from pylab import *
import os

from FR import *

c = 299792458.0

if os.path.isfile('POS.npz')==0:
    from Image_Creator import *
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

    POS=IC(Lt,l,p,h,X,Y,Z,tilt,azimut)
    savez('POS.npz',Lt=Lt,POS=POS,l=l,p=p,h=h,X=X,Y=Y,Z=Z)
else:
    MAT=load('POS.npz')
    POS=MAT['POS']
    l=MAT['l']
    p=MAT['p']
    h=MAT['h']
    X=MAT['X']
    Y=MAT['Y']
    Z=MAT['Z']
    Lt=MAT['Lt']

print('Lt= %2.3f mus, %2.3f x %2.3f x %2.3f m3' %(Lt/1e-6,l,p,h))

#Frequency response calculation
f = array(arange(10e6,500e6,1e6))

R=0.998 #loss coefficient

#measurement point
X_1=array([4.5])
Y_1=array([3])
Z_1=array([1])

Fx,Fy,Fz=FR(X_1,Y_1,Z_1,POS,R,f)

for i in range(0,len(Fz[:,0])):
    figure(i)
    subplot(311)
    l = plot(f/1e6,20*log10(abs(Fx[i,:])))
    grid(True)
    title('Frequency response')
    ylabel('$E_x$')
    subplot(312)
    l = plot(f/1e6,20*log10(abs(Fy[i,:])))
    grid(True)
    ylabel('$E_y$')
    subplot(313)
    l = plot(f/1e6,20*log10(abs(Fz[i,:])))
    ylabel('$E_z$')
    grid(True)
    xlabel('frequency (MHz)')

show()
