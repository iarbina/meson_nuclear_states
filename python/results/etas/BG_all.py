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

	datafile1 = "CS_data/free_CS.dat"
	Z1, A1, B1, G1, SQRTS1 = reading(datafile1, 5)
	
	datafile2 = "GW_data/free_GW.dat"
	Z2, A2, B2, G2, SQRTS2 = reading(datafile2, 5)

	datafile3 = "IOV_data/free_IOV.dat"
	Z3, A3, B3, G3, SQRTS3 = reading(datafile3, 5)
	
	datafile4 = "KSW_data/free_KSW.dat"
	Z4, A4, B4, G4, SQRTS4 = reading(datafile4, 5)

	datafile5 = "M2_data/free_M2.dat"
	Z5, A5, B5, G5, SQRTS5 = reading(datafile5, 5)
	
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
	AA3 = []
	BB3 = []
	GG3 = []
	AA4 = []
	BB4 = []
	GG4 = []
	AA5 = []
	BB5 = []
	GG5 = []
	for i in range(len(A1)):
		AA1.append(A1[i]**(2./3))
		BB1.append(B1[i])
		GG1.append(G1[i])
		AA2.append(A2[i]**(2./3))
		BB2.append(B2[i])
		GG2.append(G2[i])
		AA3.append(A3[i]**(2./3))
		BB3.append(B3[i])
		GG3.append(G3[i])
		AA4.append(A4[i]**(2./3))
		BB4.append(B4[i])
		GG4.append(G4[i])
		AA5.append(A5[i]**(2./3))
		BB5.append(B5[i])
		GG5.append(G5[i])


# Plot----------------------------------------------------------------	
	
	plt.rc('text',usetex=True)
	font = {'family':'serif','size':15}
	plt.rc('font',**font)
	
	xx = np.linspace(0,12,50)

	ystep = 2.0
	xstep = -1.0

	xs = -15.0

	fig = plt.figure(facecolor='w',figsize=(5.5,4.5))
	ax = fig.add_subplot(111)
	ax.plot(AA1, BB1, '--', markersize=15, linewidth=0.5, label=r'Free CS')
	ax.plot(AA2, BB2, '--', markersize=15, linewidth=0.5, label=r'Free GW')
	ax.plot(AA3, BB3, '--', markersize=15, linewidth=0.5, label=r'Free IOV')
	ax.plot(AA4, BB4, '--', markersize=15, linewidth=0.5, label=r'Free KWS')
	ax.plot(AA5, BB5, '--', markersize=15, linewidth=0.5, label=r'Free M2')
	plt.text(AA1[-1]+xs, BB1[-1]-1, r'CS')
	plt.text(AA2[-1]+xs, BB2[-1]-2, r'GW')
	plt.text(AA3[-1]+xs+5, BB3[-1]-4, r'IOV')
	plt.text(AA4[-1]+xs, BB4[-1]-4.5, r'KWS')
	plt.text(AA5[-1]+xs+5, BB5[-1]-0.5, r' M2')
	for i in range(len(A1)):
		ax.scatter(AA1[i], BB1[i], marker='s', color='C'+str(i))
		ax.scatter(AA2[i], BB2[i], marker='v', color='C'+str(i))
		ax.scatter(AA3[i], BB3[i], marker='o', color='C'+str(i))
		ax.scatter(AA4[i], BB4[i], marker='*', color='C'+str(i))
		ax.scatter(AA5[i], BB5[i], marker='d', color='C'+str(i))
	ax.set_xlabel(r"$A^{2/3}$")
	ax.set_ylabel(r"$B_{\eta}$ (MeV)")
	#ax.set_xlim(0,40)
	ax.set_ylim(0,35)
	ax.tick_params(which='both',direction='in',top=True,right=True)
	#ax.legend(frameon=False)
	plt.text(AA1[0]+xstep, BB2[0]+ystep, r'C')
	plt.text(AA1[1]+xstep-0.5, BB2[1]+ystep, r'Mg')
	plt.text(AA1[2]+xstep, BB2[2]+ystep, r'S')
	plt.text(AA1[3]+xstep, BB2[3]+ystep, r'Ca')
	plt.text(AA1[4]+xstep, BB2[4]+ystep, r'Pb')
	plt.locator_params(axis='x', nbins=10)
	plt.locator_params(axis='y', nbins=10)
	plt.tight_layout()
	
	
	fig2 = plt.figure(facecolor='w',figsize=(5.5,4.5))
	ax2 = fig2.add_subplot(111)
	ax2.plot(AA1, GG1, '--', markersize=15, linewidth=0.5, label=r'Free CS')
	ax2.plot(AA2, GG2, '--', markersize=15, linewidth=0.5, label=r'Free GW')
	ax2.plot(AA3, GG3, '--', markersize=15, linewidth=0.5, label=r'Free IOV')
	ax2.plot(AA4, GG4, '--', markersize=15, linewidth=0.5, label=r'Free KWS')
	ax2.plot(AA5, GG5, '--', markersize=15, linewidth=0.5, label=r'Free M2')
	plt.text(AA1[-1]+xs+5,  GG1[-1]+0.5,  r'CS')
	plt.text(AA2[-1]+xs+5,  GG2[-1]+0.5,  r'GW')
	plt.text(AA3[-1]+xs,    GG3[-1]-1,    r'IOV')
	plt.text(AA4[-1]+xs,    GG4[-1]+0.5,r'KWS')
	plt.text(AA5[-1]+xs,    GG5[-1]-2,  r' M2')
	for i in range(len(A1)):
		ax2.scatter(AA1[i], GG1[i], marker='s', color='C'+str(i))
		ax2.scatter(AA2[i], GG2[i], marker='v', color='C'+str(i))
		ax2.scatter(AA3[i], GG3[i], marker='o', color='C'+str(i))
		ax2.scatter(AA4[i], GG4[i], marker='*', color='C'+str(i))
		ax2.scatter(AA5[i], GG5[i], marker='d', color='C'+str(i))
	ax2.set_xlabel(r"$A^{2/3}$")
	ax2.set_ylabel(r"$\Gamma_{\eta}$ (MeV)")
	#ax2.set_xlim(0,40)
	ax2.set_ylim(0,15)
	ax2.tick_params(which='both',direction='in',top=True,right=True)
	#ax2.legend(frameon=False)
	plt.text(AA1[0]+xstep, GG4[0]+ystep-1, r'C')
	plt.text(AA1[1]+xstep-0.2, GG4[1]+ystep, r'Mg')
	plt.text(AA1[2]+xstep, GG4[2]+ystep, r'Si')
	plt.text(AA1[3]+xstep, GG4[3]+ystep, r'Ca')
	plt.text(AA1[4]+xstep, GG4[4]+ystep, r'Pb')
	plt.locator_params(axis='x', nbins=10)
	plt.locator_params(axis='y', nbins=10)
	plt.tight_layout()
	
	
	plt.show()
	
main()
