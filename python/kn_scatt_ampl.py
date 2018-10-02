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
		if n == 6:
			lines = myfile.readlines()
		else:
			lines = myfile.readlines()[1:]
		myfile.close()
			
		for line in lines:
			p = line.split()
			x1.append(float(p[0]))
			x2.append(float(p[1]))
			x3.append(float(p[2]))
			x4.append(float(p[3]))
			x5.append(float(p[4]))
			if n == 6:
				x6.append(float(p[5]))
					
		a = np.array(x1)
		b = np.array(x2)
		c = np.array(x3)
		d = np.array(x4)
		e = np.array(x5)
		if n == 6:
			f = np.array(x6)
			return a, b, c, d, e, f
		else:
			return a, b, c, d, e	
		
#----------------------------------------------------------------------
	# read data from files
	
	datafile1 = "../inp/input.inp"
	Z, A, Rn, a0, pqn, aqn = reading(datafile1, 6)
	print Z, A, Rn, a0, pqn, aqn
	
	datafile2 = "../inp/amplKN00.inp"
	sqrts00, tpr00, tpi00, tnr00, tni00 = reading(datafile2, 5)
	
	datafile3 = "../inp/amplKN025.inp"
	sqrts025, tpr025, tpi025, tnr025, tni025 = reading(datafile3, 5)
	
	datafile4 = "../inp/amplKN05.inp"
	sqrts05, tpr05, tpi05, tnr05, tni05 = reading(datafile4, 5)
	
	datafile5 = "../inp/amplKN075.inp"
	sqrts075, tpr075, tpi075, tnr075, tni075 = reading(datafile5, 5)
	
	datafile6 = "../inp/amplKN10.inp"
	sqrts10, tpr10, tpi10, tnr10, tni10 = reading(datafile6, 5)
	

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
	
	plt.rc('text',usetex=True)
	font = {'family':'serif','size':15}
	plt.rc('font',**font)
	
	xx = np.linspace(0,12,50)

	pt1 = '.'
	ps1 = 4

	fig = plt.figure(facecolor='w',figsize=(7,5))
	ax = fig.add_subplot(111)
	ax.plot(sqrts00-0*Eth, tKN(tpr00, tnr00), pt1, markersize = ps1, label = r"$0.00 \rho_0$")
	ax.plot(sqrts025-0*Eth, tKN(tpr025, tnr025), pt1, markersize = ps1, label = r"$0.25 \rho_0$")
	ax.plot(sqrts05-0*Eth, tKN(tpr05, tnr05), pt1, markersize = ps1, label = r"$0.50 \rho_0$")
	ax.plot(sqrts075-0*Eth, tKN(tpr075, tnr075), pt1, markersize = ps1, label=r"$0.75 \rho_0$")
	ax.plot(sqrts10-0*Eth, tKN(tpr10, tnr10), pt1, markersize = ps1, label=r"$1.00 \rho_0$")
	ax.set_xlabel(r"$\sqrt{s}$ (MeV)")
	ax.set_ylabel(r"$\mathrm{Re}[T_{KN}]$ ($\mathrm{MeV}^{-1}$)")
	#ax.set_xlim(0,6)
	#ax.set_ylim(-210,10)
	ax.tick_params(which='both',direction='in',top=True,right=True)
	ax.legend()
	#plt.axvline(x=0, linestyle='--', linewidth=0.4, color='grey')
	plt.tight_layout()


	# Plot 2

	pt2 = '.'
	ps2 = 4
	
	fig2 = plt.figure(facecolor='w',figsize=(7,5))
	ax2 = fig2.add_subplot(111)
	ax2.plot(sqrts00-0*Eth, tKN(tpi00, tni00), pt2, markersize = ps2, label = r"$0.00 \rho_0$")
	ax2.plot(sqrts025-0*Eth, tKN(tpi025, tni025), pt2, markersize = ps2, label = r"$0.25 \rho_0$")
	ax2.plot(sqrts05-0*Eth, tKN(tpi05, tni05), pt2, markersize = ps2, label = r"$0.50 \rho_0$")
	ax2.plot(sqrts075-0*Eth, tKN(tpi075, tni075), pt2, markersize = ps2, label = r"$0.75 \rho_0$")
	ax2.plot(sqrts10-0*Eth, tKN(tpi10, tni10), pt2, markersize = ps2, label = r"$1.00 \rho_0$")
	ax2.set_xlabel(r"$\sqrt{s}$ (MeV)")
	ax2.set_ylabel(r"$\mathrm{Im}[T_{KN}]$ ($\mathrm{MeV}^{-1}$)")
	#ax2.set_xlim(0,6)
	#ax2.set_ylim(-210,10)
	ax2.tick_params(which='both',direction='in',top=True,right=True)
	ax2.legend()
	#plt.axvline(x=0, linestyle='--', linewidth=0.4, color='grey')
	plt.tight_layout()
	
    # Plot 3

	pt3 = '.'
	ps3 = 4
	
	fig3 = plt.figure(facecolor='w',figsize=(7,5))
	ax3 = fig3.add_subplot(111)
	ax3.plot(sqrts00-0*Eth, TKN2(tpr00, tpi00, tnr00, tni00), pt2, markersize = ps2, label = r"$0.00 \rho_0$")
	ax3.plot(sqrts025-0*Eth, TKN2(tpr025, tpi025, tnr025, tni025), pt2, markersize = ps2, label = r"$0.25 \rho_0$")
	ax3.plot(sqrts05-0*Eth, TKN2(tpr05, tpi05, tnr05, tni05), pt2, markersize = ps2, label = r"$0.50 \rho_0$")
	ax3.plot(sqrts075-0*Eth, TKN2(tpr075, tpi075, tnr075, tni075), pt2, markersize = ps2, label = r"$0.75 \rho_0$")
	ax3.plot(sqrts10-0*Eth, TKN2(tpr10, tpi10, tnr10, tni10), pt2, markersize = ps2, label = r"$1.00 \rho_0$")
	ax3.set_xlabel(r"$\sqrt{s}$ (MeV)")
	ax3.set_ylabel(r"$|T_{KN}|^2$ ($\mathrm{MeV}^{-2}$)")
	#ax3.set_xlim(0,6)
	#ax3.set_ylim(-210,10)
	ax3.tick_params(which='both',direction='in',top=True,right=True)
	ax3.legend()
	#plt.axvline(x=0, linestyle='--', linewidth=0.4, color='grey')
	plt.tight_layout()
	
	plt.show()
	
main()
