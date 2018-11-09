#!usr/bin/python
# coding: latin-1

import os, glob
import numpy as np
import matplotlib.pyplot as plt

def main():
	
	
	def reading(datafile,n):
		
		x1 = []
		x2 = []
		x3 = []
		x4 = []
		x5 = []
		x6 = []
		p  = []
		
		myfile = open(datafile,'r')
		if n == 8:
			lines = myfile.readlines()[0]
		else:
			lines = myfile.readlines()
			
		myfile.close()
			
		for line in lines:
			p = line.split()
			x1.append(float(p[0]))
			x2.append(float(p[1]))
			x3.append(float(p[2]))
			x4.append(float(p[3]))
			x5.append(float(p[4]))
			x6.append(float(p[5]))
					
		a = np.array(x1)
		b = np.array(x2)
		c = np.array(x3)
		d = np.array(x4)
		e = np.array(x5)
		f = np.array(x6)
		
		return a, b, c, d, e, f
		
#----------------------------------------------------------------------
	# read data from files
	
	file_path = "../../inp/nuclei/"

	datafile1 = file_path + "input.inp"
	x1, x2, x3, x4, x5, x6 = reading(datafile1, 6)
	Z = x1[0]
	A = x2[0]
	Rn = x3[0]
	a0 = x4[0]
	pqn = x5[0]
	aqn = x6[0]
	print Z, A, Rn, a0, pqn, aqn
	
	datafile2 = "../../out/scatt_rho.dat"
	sqrts, rho, tpr, tpi, tnr, tni = reading(datafile2, 7)
	
#----------------------------------------------------------------------
	# Construct Re and Im parts of the scattering amplitudes

	def tKN(tp, tn):
		return 0.5 * (tp + tn)
	
	def TKN2(tpr, tpi, tnr, tni):
		return 0.25 * ( (tpr + tnr)**2 + (tpi + tni)**2 )

#----------------------------------------------------------------------
	# some nuclear parameters and functions

	pi = 3.14159265358979323846260
	hbarc = 197.3269710 # MeV fm
	rhoc = 0.16
	alpha = 0.27
	
	A = 2
	hmass = 493.70   # Mass of the hadron (in this case the kaon)
	nmass = 939.0    # Average nucleon mass
	mu = hmass*nmass*A/(hmass+nmass*A) # Reduced mass
	Eth = hmass + nmass
	dx = 1e-2
   	
   	def DenSuma():
   		suma = 0.0
		N = 10 # times the nuclear radius
		x = 0.0
		if (A <= 16):
			suma = (Rn*np.sqrt(pi))**3*(1.0+1.50*a0)
		elif (A > 16):
			while (x < N*Rn):
				suma = suma + x*x*dx / (1.0 + np.exp((x-Rn)/a0))
				x = x + dx
			suma = 4*pi*suma
		DenSuma = suma
		return DenSuma
   	
   	rho0=A/DenSuma() 
   	
   	def Density(x):
   		if (A <= 16):
			Density = rho0*(1.0+a0*(x/Rn)**2)*np.exp(-x*x/(Rn*Rn))
		elif (A > 16):
			Density = rho0/(1.0+np.exp((x-Rn)/a0))
		return Density
   	

#----------------------------------------------------------------------
	# Plot 1	
	
	rows = 201
	cols = 6

	plt.rc('text',usetex=True)
	font = {'family':'serif','size':15}
	plt.rc('font',**font)
	
	sval = 150

	rhos = np.unique(rho)
	ww = np.where(sqrts == sqrts[sval])
	tprs = tpr[ww]	
	tpis = tpi[ww]	
	tnrs = tnr[ww]	
	tnis = tni[ww]

	ww2 = np.where((sqrts == sqrts[sval]) & ((rho == 0.00) | (rho == 0.04) | (rho == 0.08) | (rho == 0.12) | (rho == 0.16)))
	rhos2 = rho[ww2]
	tpr0 = tpr[ww2]	
	tpi0 = tpi[ww2]	
	tnr0 = tnr[ww2]	
	tni0 = tni[ww2]
	print 'rhos2 =', rhos2


	pt1 = '.'
	ps1 = 4


	fig = plt.figure(facecolor='w',figsize=(7,5))
	ax = fig.add_subplot(111)
	ax.plot(rhos, tKN(tprs, tnrs), pt1, markersize = ps1, label = r"Cubic spline")
	ax.plot(rhos2, tKN(tpr0, tnr0), 'o-', label="Original data")
	ax.set_xlabel(r"$\rho$ ($\mathrm{fm^{-3}}$)")
	ax.set_ylabel(r"$\mathrm{Re}[T_{KN}]$ ($\mathrm{MeV}^{-1}$)")
	#ax.set_xlim(0,6)
	#ax.set_ylim(-210,10)
	ax.tick_params(which='both',direction='in',top=True,right=True)
	plt.title(r"$\sqrt{s} =$" + str(sqrts[sval]))
	ax.legend()
	#plt.axvline(x=0, linestyle='--', linewidth=0.4, color='grey')
	plt.tight_layout()


	# Plot 2

	pt2 = '.'
	ps2 = 4
	
	fig2 = plt.figure(facecolor='w',figsize=(7,5))
	ax2 = fig2.add_subplot(111)
	ax2.plot(rhos, tKN(tpis, tnis), pt2, markersize = ps2, label=r"Cubic spline")
	ax2.plot(rhos2, tKN(tpi0, tni0), 'o-', label="Original data")
	ax2.set_xlabel(r"$\rho$ ($\mathrm{fm^{-3}}$)")
	ax2.set_ylabel(r"$\mathrm{Im}[T_{KN}]$ ($\mathrm{MeV}^{-1}$)")
	#ax2.set_xlim(0,6)
	#ax2.set_ylim(-210,10)
	ax2.tick_params(which='both',direction='in',top=True,right=True)
	plt.title(r"$\sqrt{s} =$" + str(sqrts[sval]))
	ax2.legend()
	#plt.axvline(x=0, linestyle='--', linewidth=0.4, color='grey')
	plt.tight_layout()
	
    # Plot 3

	pt3 = '.'
	ps3 = 4
	
	fig3 = plt.figure(facecolor='w',figsize=(7,5))
	ax3 = fig3.add_subplot(111)
	ax3.plot(rhos, TKN2(tprs, tpis, tnrs, tnis), pt2, markersize = ps2, label = r"Cubic spline")
	ax3.plot(rhos2, TKN2(tpr0, tpi0, tnr0, tni0), 'o-', label="Original data")
	ax3.set_xlabel(r"$\rho$ ($\mathrm{fm^{-3}}$)")
	ax3.set_ylabel(r"$|T_{KN}|^2$ ($\mathrm{MeV}^{-2}$)")
	#ax3.set_xlim(0,6)
	#ax3.set_ylim(-210,10)
	ax3.tick_params(which='both',direction='in',top=True,right=True)
	plt.title(r"$\sqrt{s} =$" + str(sqrts[sval]))
	ax3.legend()
	#plt.axvline(x=0, linestyle='--', linewidth=0.4, color='grey')
	plt.tight_layout()
	
	plt.show()
	
main()
