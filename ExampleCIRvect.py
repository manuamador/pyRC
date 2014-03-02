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
from pylab import *
import os

from CIR import *


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

#Channel impulse response calculation
R=0.998 #loss coefficient

#measurement point
X_1=array([4.5])
Y_1=array([3])
Z_1=array([1])

Sx,Sy,Sz,t=CIR(X_1,Y_1,Z_1,POS,R)


#POST TREATMENT
#sampling
f0=1e9
N=int(5*Lt*f0)
ts=Lt/N*arange(1,N+1)
Sxs=zeros((len(Sz[:,0]),N))
Sys=zeros((len(Sz[:,0]),N))
Szs=zeros((len(Sz[:,0]),N))

for i in range(0,len(Sz[:,0])):
	for j in range (0,len(Sz[0,:])):
		u=int(t[i,j]/Lt*N)
		if u<N:
			Sxs[i,u]=Sxs[i,u]+Sx[i,j]
			Sys[i,u]=Sys[i,u]+Sy[i,j]
			Szs[i,u]=Szs[i,u]+Sz[i,j]

for i in range(0,len(Sz[:,0])):
	figure(i)
	subplot(311)
	l = plot(ts/1e-6,(Sxs[i,:]))
	grid(True)
	title('Channel impulse response')
	ylabel('$E_x$')	
	subplot(312)
	l = plot(ts/1e-6,(Sys[i,:]))
	grid(True)
	ylabel('$E_y$')
	subplot(313)
	l = plot(ts/1e-6,(Szs[i,:]))
	ylabel('$E_z$')
	grid(True)
	xlabel('time ($\mu$s)')

def nextpow2(v):
    v -= 1
    v |= v >> 1
    v |= v >> 2
    v |= v >> 4
    v |= v >> 8
    v |= v >> 16
    return v + 1

Fs = N/Lt #sampling frequency
periode = 1/Fs                     # sample length
L = N                  #Number of points
tt = array(arange(0,L-1,1))*periode                # time...
NFFT = 2^nextpow2(N) # Next power of 2 from length of y

FFTx=zeros((len(Szs[:,0]),NFFT/2-1),'complex')
FFTy=zeros((len(Szs[:,0]),NFFT/2-1),'complex')
FFTz=zeros((len(Szs[:,0]),NFFT/2-1),'complex')

for i in range(0,len(Sz[:,0])):
	Yx = fft.fft(Sxs[i,:],NFFT)
	Yy = fft.fft(Sys[i,:],NFFT)
	Yz = fft.fft(Szs[i,:],NFFT)
	f = Fs/2*linspace(0,1,NFFT/2)
	FFTx[i,:] = Yx[1:NFFT/2]
	FFTy[i,:] = Yy[1:NFFT/2]
	FFTz[i,:] = Yz[1:NFFT/2]
	freq = f

for i in range(0,len(Sz[:,0])):
	figure(i+2)
	subplot(311)
	plot(freq[0:len(FFTx[i,:])]/1e6,20*log(abs(FFTx[i,:])))
	grid(True)
	title('Frequency response')
	ylabel('$FFT_x$')	
	xlim( 0, 500 ) 
	subplot(312)
	plot(freq[0:len(FFTy[i,:])]/1e6,20*log(abs(FFTy[i,:])))
	grid(True)
	ylabel('$FFT_y$')
	xlim( 0, 500 ) 
	subplot(313)
	plot(freq[0:len(FFTz[i,:])]/1e6,20*log(abs(FFTz[i,:])))
	ylabel('$FFT_z$')
	grid(True)
	xlabel('frequency (MHz)')
	xlim( 0, 500 )  

show()
