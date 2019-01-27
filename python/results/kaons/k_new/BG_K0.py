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

	datafile1 = "antiK.dat"
	Z1, A1, B1, G1, SQRTS1 = reading(datafile1, 5)
	
#----------------------------------------------------------------------
	
	pi = 3.14159265358979323846260
	hbarc = 197.3269710 # MeV fm
	rhoc = 0.16
	alpha = 0.27
	
	etamass = 547.853
	hmass = 493.70   # Mass of the hadron (in this case the kaon)
	nmass = 939.0    # Average nucleon mass

	Bn = 8.5
	Eth = etamass + nmass

	AA1 = []
	BB1 = []
	GG1 = []
	for i in range(len(A1)):
		AA1.append(A1[i]**(2./3))
		BB1.append(B1[i])
		GG1.append(G1[i])


# Plot----------------------------------------------------------------	
	
	plt.rc('text',usetex=True)
	font = {'family':'serif','size':15}
	plt.rc('font',**font)
	
	xx = np.linspace(0,12,50)

	ystep = 1.0
	xstep = -1.0

	indices = [1, 2, 4, 6, 7]
	fig = plt.figure(facecolor='w',figsize=(5.5,4.5))
	ax = fig.add_subplot(111)
	AA = [AA1[1], AA1[2], AA1[4], AA1[6], AA1[7]]
	BB = [BB1[1], BB1[2], BB1[4], BB1[6], BB1[7]]
	ax.plot(AA, BB, '--', markersize=15, linewidth=0.5, label=r'Free CS')
	#for i in (1, 2, 4, 6, 7):
	for i in range(len(A1)):
		ax.scatter(AA1[i], BB1[i])
	ax.set_xlabel(r"$A^{2/3}$")
	ax.set_ylabel(r"$B_{\eta}$ (MeV)")
	#ax.set_xlim(0,40)
	#ax.set_ylim(0,20)
	ax.tick_params(which='both',direction='in',top=True,right=True)
	ax.legend()
	#plt.text(AA1[0]+xstep, BB1[0]+ystep, r'B')
	plt.text(AA1[1]+xstep, BB1[1]+ystep, r'C')
	plt.text(AA1[2]+xstep, BB1[2]+ystep, r'O')
	plt.text(AA1[3]+xstep-1, BB1[3]+ystep, r'Mg')
	plt.text(AA1[4]+xstep, BB1[4]+ystep, r'Si')
	plt.text(AA1[5]+xstep, BB1[5]+ystep, r'S')
	plt.text(AA1[6]+xstep, BB1[6]+ystep, r'Ca')
	plt.text(AA1[7]+xstep, BB1[7]-ystep, r'Pb')
	plt.locator_params(axis='x', nbins=10)
	plt.locator_params(axis='y', nbins=10)
	plt.tight_layout()
	
	
	fig2 = plt.figure(facecolor='w',figsize=(5.5,4.5))
	ax2 = fig2.add_subplot(111)
	ax2.plot(AA1, GG1, '--', markersize=15, linewidth=0.5, label=r'Free CS')
	for i in range(len(A1)):
		ax2.scatter(AA1[i], GG1[i])
	ax2.set_xlabel(r"$A^{2/3}$")
	ax2.set_ylabel(r"$\Gamma_{\eta}$ (MeV)")
	#ax2.set_xlim(0,40)
	#ax2.set_ylim(0,20)
	ax2.tick_params(which='both',direction='in',top=True,right=True)
	ax2.legend()
	plt.text(AA1[0]+xstep, GG1[0]+ystep, r'B')
	plt.text(AA1[1]+xstep, GG1[1]+ystep, r'C')
	plt.text(AA1[2]+xstep, GG1[2]+ystep, r'O')
	plt.text(AA1[3]+xstep, GG1[3]+ystep, r'Mg')
	plt.text(AA1[4]+xstep, GG1[4]+ystep, r'Si')
	plt.text(AA1[5]+xstep, GG1[5]+ystep, r'S')
	plt.text(AA1[6]+xstep, GG1[6]+ystep, r'Ca')
	plt.text(AA1[7]+xstep, GG1[7]+ystep, r'Pb')
	plt.locator_params(axis='x', nbins=10)
	plt.locator_params(axis='y', nbins=10)
	plt.tight_layout()
	
	
	plt.show()
	
main()
