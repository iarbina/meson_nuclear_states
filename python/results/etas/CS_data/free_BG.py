#!usr/bin/python
# coding: latin-1

import os, glob
import math
import numpy as np
import matplotlib.pyplot as plt

def main():
	
	
	def reading(datafile,n):
		
		x1 = []
		x2 = []
		x3 = []
		x4 = []
		x5 = []
		p  = []
		
		myfile = open(datafile,'r')
		lines = myfile.readlines()[1:]
		myfile.close()
		
		for line in lines:
			p = line.split()
			x1.append(float(p[0]))
			x2.append(float(p[1]))
			x3.append(float(p[2]))
			x4.append(float(p[3]))
			if n == 5:
				x5.append(float(p[4]))
					
		a = np.array(x1)
		b = np.array(x2)
		c = np.array(x3)
		d = np.array(x4)
		if n == 5:
			e = np.array(x5)
			return a, b, c, d, e
		else:
			return a, b, c, d
		
#----------------------------------------------------------------------

	datafile1 = "../CS_old/fCS_mares.dat"
	Z, A, B, G = reading(datafile1, 4)
	
	datafile2 = "free_CS.dat"
	Z2, A2, B2, G2, SQRTS = reading(datafile2, 5)

#----------------------------------------------------------------------
	
	pi = 3.14159265358979323846260
	hbarc = 197.3269710 # MeV fm
	rhoc = 0.16
	alpha = 0.27
	
	etamass = 547.853
	hmass = 493.70   # Mass of the hadron (in this case the kaon)
	nmass = 939.0    # Average nucleon mass
	mu = hmass*nmass*A/(hmass+nmass*A) # Reduced mass

	Bn = 8.5
	Eth = etamass + nmass

	AA = []
	BB = []
	GG = []
	AA2 = []
	BB2 = []
	GG2 = []
	for i in range(len(A)):
		AA.append(A[i]**(2./3))
		BB.append(B[i])
		GG.append(G[i])
		AA2.append(A2[i]**(2./3))
		BB2.append(B2[i])
		GG2.append(G2[i])


# Plot----------------------------------------------------------------	
	
	plt.rc('text',usetex=True)
	font = {'family':'serif','size':15}
	plt.rc('font',**font)
	
	xx = np.linspace(0,12,50)

	ystep = 2.0
	xstep = -1.0

	fig = plt.figure(facecolor='w',figsize=(5.5,4.5))
	ax = fig.add_subplot(111)
	ax.plot(AA, BB, '--', markersize=15, linewidth=0.5, label=r'Free CS ($\sqrt{s}$ Mares et. al.)')
	ax.plot(AA2, BB2, '--', markersize=15, linewidth=0.5, label=r'Free CS ($\sqrt{s}$ ours)')
	for i in range(len(A)):
		ax.scatter(AA[i], BB[i])
		ax.scatter(AA2[i], BB2[i])
	ax.set_xlabel(r"$A^{2/3}$")
	ax.set_ylabel(r"$B_{\eta}$ (MeV)")
	ax.set_xlim(0,40)
	ax.set_ylim(0,20)
	ax.tick_params(which='both',direction='in',top=True,right=True)
	ax.legend()
	plt.text(AA[0]+xstep, BB2[0]+ystep, r'C')
	plt.text(AA[1]+xstep-0.5, BB[1]+ystep, r'Mg')
	plt.text(AA[2]+xstep, BB[2]+ystep, r'S')
	plt.text(AA[3]+xstep, BB[3]+ystep, r'Ca')
	plt.text(AA[4]+xstep, BB2[4]+ystep, r'Pb')
	plt.locator_params(axis='x', nbins=10)
	plt.locator_params(axis='y', nbins=10)
	plt.tight_layout()
	
	
	fig2 = plt.figure(facecolor='w',figsize=(5.5,4.5))
	ax2 = fig2.add_subplot(111)
	ax2.plot(AA, GG, '--', markersize=15, linewidth=0.5, label=r'Free CS ($\sqrt{s}$ Mares et. al.)')
	ax2.plot(AA2, GG2, '--', markersize=15, linewidth=0.5, label=r'Free CS ($\sqrt{s}$ ours)')
	for i in range(len(A)):
		ax2.scatter(AA[i], GG[i])
		ax2.scatter(AA2[i], GG2[i])
	ax2.set_xlabel(r"$A^{2/3}$")
	ax2.set_ylabel(r"$\Gamma_{\eta}$ (MeV)")
	ax2.set_xlim(0,40)
	ax2.set_ylim(0,20)
	ax2.tick_params(which='both',direction='in',top=True,right=True)
	ax2.legend()
	plt.text(AA[0]+xstep, GG2[0]+ystep, r'C')
	plt.text(AA[1]+xstep-0.85, GG[1]+ystep, r'Mg')
	plt.text(AA[2]+xstep, GG[2]+ystep, r'S')
	plt.text(AA[3]+xstep, GG[3]+ystep, r'Ca')
	plt.text(AA[4]+xstep, GG[4]+ystep, r'Pb')
	plt.locator_params(axis='x', nbins=10)
	plt.locator_params(axis='y', nbins=10)
	plt.tight_layout()
	
	
	plt.show()
	
main()
