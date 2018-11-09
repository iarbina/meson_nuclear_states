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
					
		a = np.array(x1)
		b = np.array(x2)
		c = np.array(x3)
		d = np.array(x4)
		return a, b, c, d
		
#----------------------------------------------------------------------

	datafile1 = "antiK_mares.dat"
	Z1, A1, B1, G1 = reading(datafile1, 5)
	
	datafile2 = "antiK_ours.dat"
	Z2, A2, B2, G2 = reading(datafile2, 5)

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
	AA2 = []
	BB2 = []
	GG2 = []
	for i in range(len(A1)):
		AA1.append(A1[i]**(2./3))
		BB1.append(B1[i])
		GG1.append(G1[i])
		AA2.append(A2[i]**(2./3))
		BB2.append(B2[i])
		GG2.append(G2[i])


# Plot----------------------------------------------------------------	
	
	plt.rc('text',usetex=True)
	font = {'family':'serif','size':17}
	plt.rc('font',**font)
	
	xx = np.linspace(0,12,50)

	ystep = 2.0
	xstep = -1.0

	fig = plt.figure(facecolor='w',figsize=(5.5,4.5))
	ax = fig.add_subplot(111)
	ax.plot(AA1, BB1, '--', markersize=15, linewidth=0.5, label=r'$\sqrt{s}$ Mares et. al.')
	ax.plot(AA2, BB2, '--', markersize=15, linewidth=0.5, label=r'$\sqrt{s}$ ours')
	for i in range(len(A1)):
		ax.scatter(AA1[i], BB1[i], c='C' + str(i))
		ax.scatter(AA2[i], BB2[i], c='C' + str(i))
	ax.set_xlabel(r"$A^{2/3}$")
	ax.set_ylabel(r"$B_{K^-}$ (MeV)")
	ax.set_xlim(0,40)
	ax.set_ylim(0,40)
	ax.tick_params(which='both',direction='in',top=True,right=True)
	ax.legend()
	plt.text(AA1[0]+xstep, BB1[0]+ystep, r'C')
	plt.text(AA1[1]+xstep-1.5, BB1[1]+ystep, r'Mg')
	plt.text(AA1[2]+xstep, BB1[2]+ystep, r'S')
	plt.text(AA1[3]+xstep, BB1[3]+ystep, r'Ca')
	plt.text(AA1[4]+xstep, BB1[4]-ystep-1.9, r'Pb')
	plt.locator_params(axis='x', nbins=10)
	plt.locator_params(axis='y', nbins=10)
	plt.tight_layout()
	
	
	fig2 = plt.figure(facecolor='w',figsize=(5.5,4.5))
	ax2 = fig2.add_subplot(111)
	ax2.plot(AA1, GG1, '--', markersize=15, linewidth=0.5, label=r'$\sqrt{s}$ Mares et. al.')
	ax2.plot(AA2, GG2, '--', markersize=15, linewidth=0.5, label=r'$\sqrt{s}$ ours')
	for i in range(len(A1)):
		ax2.scatter(AA1[i], GG1[i], c='C' + str(i))
		ax2.scatter(AA2[i], GG2[i], c='C' + str(i))
	ax2.set_xlabel(r"$A^{2/3}$")
	ax2.set_ylabel(r"$\Gamma_{K^-}$ (MeV)")
	ax2.set_xlim(0,40)
	ax2.set_ylim(60,90)
	ax2.tick_params(which='both',direction='in',top=True,right=True)
	ax2.legend()
	plt.text(AA1[0]+xstep, GG2[0]+ystep, r'C')
	plt.text(AA1[1]+xstep-1.0, GG1[1]+ystep+2.0, r'Mg')
	plt.text(AA1[2]+xstep, GG1[2]+ystep, r'S')
	plt.text(AA1[3]+xstep, GG1[3]+ystep, r'Ca')
	plt.text(AA1[4]+xstep, GG1[4]+ystep, r'Pb')
	plt.locator_params(axis='x', nbins=10)
	plt.locator_params(axis='y', nbins=10)
	plt.tight_layout()
	
	
	plt.show()
	
main()
