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

	datafile1 = "../../inp/nuclei/input.inp"
	x1, x2, x3, x4, x5, x6 = reading(datafile1, 6)
	Z = x1[0]
	A = x2[0]
	Rn = x3[0]
	a0 = x4[0]
	pqn = x5[0]
	aqn = x6[0]
	print Z, A, Rn, a0, pqn, aqn
	
	datafile2 = "../../inp/potentials/vpot.pika"
	r, rrho0, pf, VReal, VImag = reading(datafile2,5)
	#print r, rrho0, pf, VReal, VImag

#----------------------------------------------------------------------
	
	pi = 3.14159265358979323846260
	hbarc = 197.3269710 # MeV fm
	rhoc = 0.16
	alpha = 0.27
	
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
   	
   	def Vc(x):
   		if (x.any() <= Rn):
			Vc = -0.50*Z*hbarc*(3.0-(x/Rn)**2)/(Rn*137.0359991)
		else:
			Vc = -Z*hbarc/(x*137.0359991)
		return Vc
		
	def RVopt(x):
		b0 = 0.520
		Vopt = -2.0*pi*(1.0+mu/nmass)*b0*Density(x)*hbarc**2/mu
		return Vopt
	
	def IVopt(x):
		b0 = 0.8
		Vopt = -2.0*pi*(1.0+mu/nmass)*b0*Density(x)*hbarc**2/mu
		return Vopt
		
	def RDDVopt(x):
		b0exp = -0.15
		B0 = 1.62
		DDVopt = -2.*pi*(1.+mu/nmass)*Density(x)*hbarc**2/mu \
		         *(b0exp+B0*(Density(x)/rhoc)**(alpha))
		return DDVopt
	
	def IDDVopt(x):
		b0exp = 0.62
		B0 = -0.028
		DDVopt = -2.*pi*(1.+mu/nmass)*Density(x)*hbarc**2/mu \
		         *(b0exp+B0*(Density(x)/rhoc)**(alpha))
		return DDVopt

# Plot	
	
	plt.rc('text',usetex=True)
	font = {'family':'serif','size':15}
	plt.rc('font',**font)
	
	xx = np.linspace(0,12,50)

	fig = plt.figure(facecolor='w',figsize=(7,5))
	ax = fig.add_subplot(111)
	ax.plot(r,VReal*mu/hmass,label=r"$V_{\mathrm{opt}}^{(1)}$")
	ax.plot(xx,RVopt(xx),label=r"$V_{\mathrm{opt}}^{(2)}$")
	ax.plot(xx,RDDVopt(xx),label=r"$V_{\mathrm{opt}}^{(2\mathrm{DD})}$")
	ax.plot(xx,Vc(xx),label=r"$V_{\mathrm{C}}$")
	ax.set_xlabel(r"$r$ (fm)")
	ax.set_ylabel(r"Re[$V_{\mathrm{opt}}$] (MeV)")
	ax.set_xlim(0,6)
	ax.set_ylim(-210,10)
	ax.tick_params(which='both',direction='in',top=True,right=True)
	ax.legend()
	plt.tight_layout()
	plt.show()
	
main()
