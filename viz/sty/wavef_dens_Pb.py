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
			if n == 6:
				x5.append(float(p[4]))
				x6.append(float(p[5]))
					
		a = np.array(x1)
		b = np.array(x2)
		c = np.array(x3)
		d = np.array(x4)
		if n == 6:
			e = np.array(x5)
			f = np.array(x6)
			return a, b, c, d, e, f
		else:
			return a, b, c, d	
		
#----------------------------------------------------------------------
	
	datafile1 = "../inp/nuclei/input.inp"
	x1, x2, x3, x4, x5, x6 = reading(datafile1, 6)
	nuc = 7
	Z = x1[nuc]
	A = x2[nuc]
	Rn = x3[nuc]
	a0 = x4[nuc]
	pqn = x5[nuc]
	aqn = x6[nuc]
	print Z, A, Rn, a0, pqn, aqn
	
	datafile2 = "wf2_Pb_CS.dat"
	i, x, swf, xturn = reading(datafile2,4)
	
	datafile3 = "wf2_Pb_inCS.dat"
	i3, x3, swf3, xturn3 = reading(datafile3,4)

#----------------------------------------------------------------------
	
	pi = 3.14159265358979323846260
	hbarc = 197.3269710 # MeV fm
	rhoc = 0.16
	alpha = 0.27
	
	A = 82
	hmass = 493.70   # Mass of the hadron (in this case the kaon)
	nmass = 939.0    # Average nucleon mass
	mu = hmass*nmass*A/(hmass+nmass*A) # Reduced mass
	dx = 1e-3
   	
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
   
	def Norma(r,wf):
		rr = np.array(r)
		ww = np.array(wf)
		return np.sqrt(sum(rr*ww))

	print Norma(x, swf)

# Plot	
	
	plt.rc('text',usetex=True)
	font = {'family':'serif','size':15}
	plt.rc('font',**font)
	
	xx = np.linspace(0,12,50)

	fig = plt.figure(facecolor='w',figsize=(5.5,4.5))
	ax = fig.add_subplot(111)
	ax.plot(x,swf/Norma(x, swf),'-',label=r"$1s_{K^-}$")
	ax.plot(xx,1.5*Density(xx),'-',label=r"$\rho(r)$")
	ax.set_xlabel(r"$r$ (fm)")
	ax.set_ylabel(r"$|u|^2$ ($\mathrm{fm}^{-1}$)")
	ax.set_xlim(0,12)
	ax.set_ylim(0,0.2)
	ax.tick_params(which='both',direction='in',top=True,right=True)
	ax.legend(frameon=False)
	plt.text(10,0.05,r'${}^{208}{\textrm{Pb}}$')
	plt.locator_params(axis='x', nbins=10)
	plt.locator_params(axis='y', nbins=10)
	plt.tight_layout()
	plt.show()
	
main()
