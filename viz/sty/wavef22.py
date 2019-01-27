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
		lines = myfile.readlines()
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
	
	datafile1 = "../inp/nuclei/input.inp"
	x1, x2, x3, x4, x5, x6 = reading(datafile1, 6)
	Z = x1[0]
	A = x2[0]
	Rn = x3[0]
	a0 = x4[0]
	pqn = x5[0]
	aqn = x6[0]
	print Z, A, Rn, a0, pqn, aqn
	
	datafile2 = "lwf.dat"
	li, lx, lrwf, liwf, lxturn = reading(datafile2,5)
	
	datafile3 = "rwf.dat"
	ri, rx, rrwf, riwf, rxturn = reading(datafile3,5)
	

#----------------------------------------------------------------------
	
	pi = 3.14159265358979323846260
	hbarc = 197.3269710 # MeV fm
	rhoc = 0.16
	alpha = 0.27
	
	A = 2
	hmass = 493.70   # Mass of the hadron (in this case the kaon)
	nmass = 939.0    # Average nucleon mass
	mu = hmass*nmass*A/(hmass+nmass*A) # Reduced mass
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
   	

# Plot	
	
	plt.rc('text',usetex=True)
	font = {'family':'serif','size':17}
	plt.rc('font',**font)
	
	xx = np.linspace(0,12,50)

	fig = plt.figure(facecolor='w',figsize=(5.5,4.5))
	ax = fig.add_subplot(111)
	ax.plot(lx,lrwf,'*',label=r"Forward int.")
	ax.plot(rx,rrwf,'.',label=r"Backward int.")
	ax.set_xlabel(r"$r$ (fm)")
	ax.set_ylabel(r"$\mathrm{Re}(u)$ ($\mathrm{fm}^{-1/2}$)")
	#ax.set_xlim(0,6)
	#ax.set_ylim(-3,5)
	ax.tick_params(which='both',direction='in',top=True,right=True)
	ax.legend(loc=1, prop={'size': 15}, frameon=False, handletextpad=0.1, markerfirst=False)
	plt.axvline(x=lx[-1], linestyle='--', linewidth=0.5, color='grey')
	plt.locator_params(axis='x', nbins=12)
	plt.locator_params(axis='y', nbins=10)
	plt.tight_layout()
	
	fig2 = plt.figure(facecolor='w',figsize=(5.5,4.5))
	ax2 = fig2.add_subplot(111)
	ax2.plot(lx,liwf,'*',label=r"Forward int.")
	ax2.plot(rx,riwf,'.',label=r"Backward int.")
	ax2.set_xlabel(r"$r$ (fm)")
	ax2.set_ylabel(r"$\mathrm{Im}(u)$ ($\mathrm{fm}^{-1/2}$)")
	#ax2.set_xlim(0,6)
	#ax2.set_ylim(-1,13)
	ax2.tick_params(which='both',direction='in',top=True,right=True)
	ax2.legend(loc=1, prop={'size': 15}, frameon=False, handletextpad=0.1, markerfirst=False)
	plt.axvline(x=lx[-1], linestyle='--', linewidth=0.5, color='grey')
	ax2.annotate('continuity point',
			xy=(rx[0]+0.1,riwf[0]),  # theta, radius
			xytext=(0.6, 0.55),    # fraction, fraction
			textcoords='figure fraction',
			arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3"),
			horizontalalignment='left',
			verticalalignment='bottom')
	plt.locator_params(axis='x', nbins=12)
	plt.locator_params(axis='y', nbins=10)
	plt.tight_layout()
	
	plt.show()
	
main()
