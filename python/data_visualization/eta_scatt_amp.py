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
		lines = myfile.readlines()[1:]
		myfile.close()
			
		for line in lines:
			p = line.split()
			x1.append(float(p[0]))
			x2.append(float(p[1]))
					
		a = np.array(x1)
		b = np.array(x2)
		return a, b	
		
#----------------------------------------------------------------------
	# read data from files
	
	file_path = "../../inp/eta_amplitudes/"

	datafile2 = file_path + "CS_model/free_Re_CS.dat"
	ReS_CS, ReF_CS = reading(datafile2, 5)
	
	datafile3 = file_path + "CS_model/free_Im_CS.dat"
	ImS_CS, ImF_CS = reading(datafile3, 5)
	
	datafile4 = file_path + "GW_model/free_Re_GW.dat"
	ReS_GW, ReF_GW = reading(datafile4, 5)
	
	datafile5 = file_path + "GW_model/free_Im_GW.dat"
	ImS_GW, ImF_GW = reading(datafile5, 5)
	
	datafile6 = file_path + "IOV_model/free_Re_IOV.dat"
	ReS_IOV, ReF_IOV = reading(datafile6, 5)
	
	datafile7 = file_path + "IOV_model/free_Im_IOV.dat"
	ImS_IOV, ImF_IOV = reading(datafile7, 5)
	
	datafile8 = file_path + "KSW_model/free_Re_KSW.dat"
	ReS_KSW, ReF_KSW = reading(datafile8, 5)
	
	datafile9 = file_path + "KSW_model/free_Im_KSW.dat"
	ImS_KSW, ImF_KSW = reading(datafile9, 5)
	
	datafile10 = file_path + "M2_model/free_Re_M2.dat"
	ReS_M2, ReF_M2 = reading(datafile10, 5)
	
	datafile11 = file_path + "M2_model/free_Im_M2.dat"
	ImS_M2, ImF_M2 = reading(datafile11, 5)
	
	

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
	
	hmass = 548.0   # Mass of the hadron (in this case the kaon)
	nmass = 939.0    # Average nucleon mass
	Eth = hmass + nmass
   	

#----------------------------------------------------------------------
	# Plot 1	
	
	plt.rc('text',usetex=True)
	font = {'family':'serif','size':17}
	plt.rc('font',**font)
	
	xx = np.linspace(0,12,50)

	pt1 = '.'
	ps1 = 4

	fig = plt.figure(facecolor='w',figsize=(5.5,4.5))
	ax = fig.add_subplot(111)
	ax.plot(ReS_CS, ReF_CS, ':', markersize=ps1, label=r'Free CS')
	ax.plot(ReS_GW, ReF_GW, '--', markersize=ps1, label=r'Free GW')
	ax.plot(ReS_IOV, ReF_IOV, '-.', markersize=ps1, label=r'Free IOV')
	ax.plot(ReS_KSW, ReF_KSW, '--', dashes=(5,1), markersize=ps1, label=r'Free KSW')
	ax.plot(ReS_M2, ReF_M2, '-', markersize=ps1, label=r'Free M2')
	ax.set_xlabel(r"$\sqrt{s}$ (MeV)")
	ax.set_ylabel(r"$\mathrm{Re}[F_{\eta N}]$ (fm)")
	ax.set_xlim(1400,1600)
	#ax.set_ylim(-210,10)
	ax.tick_params(which='both',direction='in',top=True,right=True)
	#ax.legend(loc=1, prop={'size': 15}, frameon=False, handletextpad=0.1, markerfirst=False)
	plt.axvline(x=Eth, linestyle='--', linewidth=0.4, color='grey')
	plt.tight_layout()


	# Plot 2

	pt2 = '.'
	ps2 = 4
	
	fig2 = plt.figure(facecolor='w',figsize=(5.5,4.5))
	ax2 = fig2.add_subplot(111)
	ax2.plot(ImS_CS, ImF_CS, ':', markersize=ps1, label=r'Free CS')
	ax2.plot(ImS_GW, ImF_GW, '--', markersize=ps1, label=r'Free GW')
	ax2.plot(ImS_IOV, ImF_IOV, '-.', markersize=ps1, label=r'Free IOV')
	ax2.plot(ImS_KSW, ImF_KSW, '--', dashes=(5,1), markersize=ps1, label=r'Free KSW')
	ax2.plot(ImS_M2, ImF_M2, '-', markersize=ps1, label=r'Free M2')
	ax2.set_xlabel(r"$\sqrt{s}$ (MeV)")
	ax2.set_ylabel(r"$\mathrm{Im}[F_{\eta N}]$ (fm)")
	ax2.set_xlim(1400,1600)
	#ax2.set_ylim(-210,10)
	ax2.tick_params(which='both',direction='in',top=True,right=True)
	ax2.legend(loc=2, prop={'size': 15}, frameon=False, handletextpad=0.1, markerfirst=True)
	plt.axvline(x=Eth, linestyle='--', linewidth=0.4, color='grey')
	plt.tight_layout()
	
	
	plt.show()
	
main()
