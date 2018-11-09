#!usr/bin/python
# coding: latin-1

import os, glob
import numpy as np
import matplotlib.pyplot as plt
from math import factorial

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

	file_path = "../../inp/amplitudes/"

	datafile2 = file_path + "amplKN00.inp"
	sqrts00, tpr00, tpi00, tnr00, tni00 = reading(datafile2, 5)
	
	datafile3 = file_path + "amplKN025.inp"
	sqrts025, tpr025, tpi025, tnr025, tni025 = reading(datafile3, 5)
	
	datafile4 = file_path + "amplKN05.inp"
	sqrts05, tpr05, tpi05, tnr05, tni05 = reading(datafile4, 5)
	
	datafile5 = file_path + "amplKN075.inp"
	sqrts075, tpr075, tpi075, tnr075, tni075 = reading(datafile5, 5)
	
	datafile6 = file_path + "amplKN10.inp"
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
	
	hmass = 493.70   # Mass of the hadron (in this case the kaon)
	nmass = 939.0    # Average nucleon mass
	Eth = hmass + nmass
	dx = 1e-2
   	
#----------------------------------------------------------------------

	

	def savitzky_golay(y, window_size, order, deriv=0, rate=1):
		r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
		The Savitzky-Golay filter removes high frequency noise from data.
		It has the advantage of preserving the original shape and
		features of the signal better than other types of filtering
		approaches, such as moving averages techniques.
		Parameters
		----------
		y : array_like, shape (N,)
			the values of the time history of the signal.
		window_size : int
			the length of the window. Must be an odd integer number.
		order : int
			the order of the polynomial used in the filtering.
			Must be less then `window_size` - 1.
		deriv: int
		    the order of the derivative to compute (default = 0 means only smoothing)
		Returns
		-------
		ys : ndarray, shape (N)
			the smoothed signal (or it's n-th derivative).
		Notes
		-----
		The Savitzky-Golay is a type of low-pass filter, particularly
		suited for smoothing noisy data. The main idea behind this
		approach is to make for each point a least-square fit with a
		polynomial of high order over a odd-sized window centered at
		the point.
		Examples
		--------
		t = np.linspace(-4, 4, 500)
		y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
		ysg = savitzky_golay(y, window_size=31, order=4)
		import matplotlib.pyplot as plt
		plt.plot(t, y, label='Noisy signal')
		plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
		plt.plot(t, ysg, 'r', label='Filtered signal')
		plt.legend()
		plt.show()
		References
		----------
		.. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
			Data by Simplified Least Squares Procedures. Analytical
			Chemistry, 1964, 36 (8), pp 1627-1639.
		.. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
			W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
			Cambridge University Press ISBN-13: 9780521880688
		"""

		try:
			window_size = np.abs(np.int(window_size))
			order = np.abs(np.int(order))
		except ValueError, msg:
			raise ValueError("window_size and order have to be of type int")
		if window_size % 2 != 1 or window_size < 1:
			raise TypeError("window_size size must be a positive odd number")
		if window_size < order + 2:
			raise TypeError("window_size is too small for the polynomials order")
		order_range = range(order+1)
		half_window = (window_size -1) // 2
		# precompute coefficients
		b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
		m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
		# pad the signal at the extremes with
		# values taken from the signal itself
		firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
		lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
		y = np.concatenate((firstvals, y, lastvals))
		return np.convolve( m[::-1], y, mode='valid')


# Smoothing of the whole set of data points

	w025 = 101
	o025 = 5
	tpr025_smooth = savitzky_golay(tpr025, w025, o025)
	tpi025_smooth = savitzky_golay(tpi025, w025, o025)
	tnr025_smooth = savitzky_golay(tnr025, w025, o025)
	tni025_smooth = savitzky_golay(tni025, w025, o025)

	w05 = 101
	o05 = 5
	tpr05_smooth = savitzky_golay(tpr05, w05, o05)
	tpi05_smooth = savitzky_golay(tpi05, w05, o05)
	tnr05_smooth = savitzky_golay(tnr05, w05, o05)
	tni05_smooth = savitzky_golay(tni05, w05, o05)

	w075 = 101
	o075 = 5
	tpr075_smooth = savitzky_golay(tpr075, w075, o075)
	tpi075_smooth = savitzky_golay(tpi075, w075, o075)
	tnr075_smooth = savitzky_golay(tnr075, w075, o075)
	tni075_smooth = savitzky_golay(tni075, w075, o075)


# Smoothing of just the needed part of the data set
	
	# First we select the data

	ww025 = np.where(sqrts025 < 1430)
	sqrts025_s = sqrts025[ww025]
	tpr025_s = tpr025[ww025]
	tpi025_s = tpi025[ww025]
	tnr025_s = tnr025[ww025]
	tni025_s = tni025[ww025]

	ww05 = np.where(sqrts05 < 1420)
	sqrts05_s = sqrts05[ww05]
	tpr05_s = tpr05[ww05]
	tpi05_s = tpi05[ww05]
	tnr05_s = tnr05[ww05]
	tni05_s = tni05[ww05]

	ww075 = np.where(sqrts075 < 1420)
	sqrts075_s = sqrts075[ww075]
	tpr075_s = tpr075[ww075]
	tpi075_s = tpi075[ww075]
	tnr075_s = tnr075[ww075]
	tni075_s = tni075[ww075]


	# Second: smooth the selected data

	w025 = 71
	o025 = 4
	tpr025_sm = savitzky_golay(tpr025_s, w025, o025)
	tpi025_sm = savitzky_golay(tpi025_s, w025, o025)
	tnr025_sm = savitzky_golay(tnr025_s, w025, o025)
	tni025_sm = savitzky_golay(tni025_s, w025, o025)

	w05 = 71
	o05 = 5
	tpr05_sm = savitzky_golay(tpr05_s, w05, o05)
	tpi05_sm = savitzky_golay(tpi05_s, w05, o05)
	tnr05_sm = savitzky_golay(tnr05_s, w05, o05)
	tni05_sm = savitzky_golay(tni05_s, w05, o05)

	w075 = 71
	o075 = 5
	tpr075_sm = savitzky_golay(tpr075_s, w075, o075)
	tpi075_sm = savitzky_golay(tpi075_s, w075, o075)
	tnr075_sm = savitzky_golay(tnr075_s, w075, o075)
	tni075_sm = savitzky_golay(tni075_s, w075, o075)


	np.savetxt('amplKN025_smooth.inp', np.array([sqrts025, tpr025_smooth, tpi025_smooth, tnr025_smooth, tni025_smooth]).T, fmt='%.6f',
				header='s12               (k-p)              (k-n)')
	np.savetxt('amplKN05_smooth.inp', np.array([sqrts05, tpr05_smooth, tpi05_smooth, tnr05_smooth, tni05_smooth]).T, fmt='%.6f',
				header='s12               (k-p)              (k-n)')
	np.savetxt('amplKN075_smooth.inp', np.array([sqrts075, tpr075_smooth, tpi075_smooth, tnr075_smooth, tni075_smooth]).T, fmt='%.6f',
				header='s12               (k-p)              (k-n)')


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
	ax.plot(sqrts025_s-0*Eth, tKN(tpr025_sm, tnr025_sm), '-', markersize = ps1, label = r"$0.25 \rho_0$ (smooth)")
	ax.plot(sqrts05-0*Eth, tKN(tpr05, tnr05), pt1, markersize = ps1, label = r"$0.50 \rho_0$")
	ax.plot(sqrts05_s-0*Eth, tKN(tpr05_sm, tnr05_sm), '-', markersize = ps1, label = r"$0.50 \rho_0$ (smooth)")
	ax.plot(sqrts075-0*Eth, tKN(tpr075, tnr075), pt1, markersize = ps1, label=r"$0.75 \rho_0$")
	ax.plot(sqrts075_s-0*Eth, tKN(tpr075_sm, tnr075_sm), '-', markersize = ps1, label=r"$0.75 \rho_0$ (smooth)")
	ax.plot(sqrts10-0*Eth, tKN(tpr10, tnr10), pt1, markersize = ps1, label=r"$1.00 \rho_0$")
	ax.set_xlabel(r"$\sqrt{s}$ (MeV)")
	ax.set_ylabel(r"$\mathrm{Re}[T_{KN}]$ ($\mathrm{MeV}^{-1}$)")
	#ax.set_xlim(0,6)
	#ax.set_ylim(-210,10)
	ax.tick_params(which='both',direction='in',top=True,right=True)
	ax.legend()
	plt.title('Smooth selected')
	#plt.axvline(x=0, linestyle='--', linewidth=0.4, color='grey')
	plt.tight_layout()


	# Plot 2

	pt2 = '.'
	ps2 = 4
	
	fig2 = plt.figure(facecolor='w',figsize=(7,5))
	ax2 = fig2.add_subplot(111)
	ax2.plot(sqrts00-0*Eth, tKN(tpi00, tni00), pt2, markersize = ps2, label = r"$0.00 \rho_0$")
	ax2.plot(sqrts025-0*Eth, tKN(tpi025, tni025), pt2, markersize = ps2, label = r"$0.25 \rho_0$")
	ax2.plot(sqrts025_s-0*Eth, tKN(tpi025_sm, tni025_sm), '-', markersize = ps2, label = r"$0.25 \rho_0$ (smooth)")
	ax2.plot(sqrts05-0*Eth, tKN(tpi05, tni05), pt2, markersize = ps2, label = r"$0.50 \rho_0$")
	ax2.plot(sqrts05_s-0*Eth, tKN(tpi05_sm, tni05_sm), '-', markersize = ps2, label = r"$0.50 \rho_0$ (smooth)")
	ax2.plot(sqrts075-0*Eth, tKN(tpi075, tni075), pt2, markersize = ps2, label = r"$0.75 \rho_0$")
	ax2.plot(sqrts075_s-0*Eth, tKN(tpi075_sm, tni075_sm), '-', markersize = ps2, label = r"$0.75 \rho_0$ (smooth)")
	ax2.plot(sqrts10-0*Eth, tKN(tpi10, tni10), pt2, markersize = ps2, label = r"$1.00 \rho_0$")
	ax2.set_xlabel(r"$\sqrt{s}$ (MeV)")
	ax2.set_ylabel(r"$\mathrm{Im}[T_{KN}]$ ($\mathrm{MeV}^{-1}$)")
	#ax2.set_xlim(0,6)
	#ax2.set_ylim(-210,10)
	ax2.tick_params(which='both',direction='in',top=True,right=True)
	ax2.legend()
	plt.title('Smooth selected')
	#plt.axvline(x=0, linestyle='--', linewidth=0.4, color='grey')
	plt.tight_layout()
	
    # Plot 3

	pt3 = '.'
	ps3 = 4
	
	fig3 = plt.figure(facecolor='w',figsize=(7,5))
	ax3 = fig3.add_subplot(111)
	ax3.plot(sqrts00-0*Eth, TKN2(tpr00, tpi00, tnr00, tni00), pt2, markersize = ps2, label = r"$0.00 \rho_0$")
	ax3.plot(sqrts025-0*Eth, TKN2(tpr025, tpi025, tnr025, tni025), pt2, markersize = ps2, label = r"$0.25 \rho_0$")
	ax3.plot(sqrts025_s-0*Eth, TKN2(tpr025_sm, tpi025_sm, tnr025_sm, tni025_sm), '-', markersize = ps3, label = r"$0.25 \rho_0$ (smooth)")
	ax3.plot(sqrts05-0*Eth, TKN2(tpr05, tpi05, tnr05, tni05), pt2, markersize = ps2, label = r"$0.50 \rho_0$")
	ax3.plot(sqrts05_s-0*Eth, TKN2(tpr05_sm, tpi05_sm, tnr05_sm, tni05_sm), '-', markersize = ps3, label = r"$0.25 \rho_0$ (smooth)")
	ax3.plot(sqrts075-0*Eth, TKN2(tpr075, tpi075, tnr075, tni075), pt2, markersize = ps2, label = r"$0.75 \rho_0$")
	ax3.plot(sqrts075_s-0*Eth, TKN2(tpr075_sm, tpi075_sm, tnr075_sm, tni075_sm), '-', markersize = ps3, label = r"$0.25 \rho_0$ (smooth)")
	ax3.plot(sqrts10-0*Eth, TKN2(tpr10, tpi10, tnr10, tni10), pt2, markersize = ps2, label = r"$1.00 \rho_0$")
	ax3.set_xlabel(r"$\sqrt{s}$ (MeV)")
	ax3.set_ylabel(r"$|T_{KN}|^2$ ($\mathrm{MeV}^{-2}$)")
	#ax3.set_xlim(0,6)
	#ax3.set_ylim(-210,10)
	ax3.tick_params(which='both',direction='in',top=True,right=True)
	ax3.legend()
	plt.title('Smooth selected')
	#plt.axvline(x=0, linestyle='--', linewidth=0.4, color='grey')
	plt.tight_layout()
	
	
	# Plot 4
	
	pt4 = '.'
	ps4 = 4
	pt4s = '-'
	ps4s = 4

	fig4 = plt.figure(facecolor='w',figsize=(7,5))
	ax4 = fig4.add_subplot(111)
	ax4.plot(sqrts00-0*Eth, tKN(tpr00, tnr00), pt4, markersize = ps4, label = r"$0.00 \rho_0$")
	ax4.plot(sqrts025-0*Eth, tKN(tpr025, tnr025), pt4, markersize = ps4, label = r"$0.25 \rho_0$")
	ax4.plot(sqrts025-0*Eth, tKN(tpr025_smooth, tnr025_smooth), pt4s, markersize = ps4, label = r"$0.25 \rho_0$ (smooth)")
	ax4.plot(sqrts05-0*Eth, tKN(tpr05, tnr05), pt1, markersize = ps4, label = r"$0.25 \rho_0$")
	ax4.plot(sqrts05-0*Eth, tKN(tpr05_smooth, tnr05_smooth), pt4s, markersize = ps4, label = r"$0.50 \rho_0$ (smooth)")
	ax4.plot(sqrts075-0*Eth, tKN(tpr075, tnr075), pt4, markersize = ps4, label = r"$0.25 \rho_0$")
	ax4.plot(sqrts075-0*Eth, tKN(tpr075_smooth, tnr075_smooth), pt4s, markersize = ps4, label=r"$0.75 \rho_0$ (smooth)")
	ax4.plot(sqrts10-0*Eth, tKN(tpr10, tnr10), pt4, markersize = ps4, label=r"$1.00 \rho_0$")
	ax4.set_xlabel(r"$\sqrt{s}$ (MeV)")
	ax4.set_ylabel(r"$\mathrm{Re}[T_{KN}]$ ($\mathrm{MeV}^{-1}$)")
	#ax4.set_xlim(0,6)
	#ax4.set_ylim(-210,10)
	ax4.tick_params(which='both',direction='in',top=True,right=True)
	ax4.legend()
	plt.title('Smooth all')
	#plt.axvline(x=0, linestyle='--', linewidth=0.4, color='grey')
	plt.tight_layout()


	# Plot 5

	pt5 = '.'
	ps5 = 4
	pt5s = '-'
	ps5s = 5
	
	fig5 = plt.figure(facecolor='w',figsize=(7,5))
	ax5 = fig5.add_subplot(111)
	ax5.plot(sqrts00-0*Eth, tKN(tpi00, tni00), pt5, markersize = ps5, label = r"$0.00 \rho_0$")
	ax5.plot(sqrts025-0*Eth, tKN(tpi025, tni025), pt5, markersize = ps5, label = r"$0.25 \rho_0$")
	ax5.plot(sqrts025-0*Eth, tKN(tpi025_smooth, tni025_smooth), pt5s, markersize = ps5, label = r"$0.25 \rho_0$ (smooth)")
	ax5.plot(sqrts025-0*Eth, tKN(tpi05, tni05), pt4, markersize = ps5, label = r"$0.25 \rho_0$")
	ax5.plot(sqrts05-0*Eth, tKN(tpi05_smooth, tni05_smooth), pt5s, markersize = ps5, label = r"$0.50 \rho_0$ (smooth)")
	ax5.plot(sqrts025-0*Eth, tKN(tpi075, tni075), pt5, markersize = ps5, label = r"$0.25 \rho_0$")
	ax5.plot(sqrts075-0*Eth, tKN(tpi075_smooth, tni075_smooth), pt5s, markersize = ps5, label = r"$0.75 \rho_0$ (smooth)")
	ax5.plot(sqrts10-0*Eth, tKN(tpi10, tni10), pt5, markersize = ps5, label = r"$1.00 \rho_0$")
	ax5.set_xlabel(r"$\sqrt{s}$ (MeV)")
	ax5.set_ylabel(r"$\mathrm{Im}[T_{KN}]$ ($\mathrm{MeV}^{-1}$)")
	#ax5.set_xlim(0,6)
	#ax5.set_ylim(-210,10)
	ax5.tick_params(which='both',direction='in',top=True,right=True)
	ax5.legend()
	plt.title('Smooth all')
	#plt.axvline(x=0, linestyle='--', linewidth=0.4, color='grey')
	plt.tight_layout()
	
    
	# Plot 6

	pt6 = '.'
	ps6 = 4
	
	fig6 = plt.figure(facecolor='w',figsize=(7,5))
	ax6 = fig6.add_subplot(111)
	ax6.plot(sqrts00-0*Eth, TKN2(tpr00, tpi00, tnr00, tni00), pt6, markersize = ps6, label = r"$0.00 \rho_0$")
	ax6.plot(sqrts025-0*Eth, TKN2(tpr025, tpi025, tnr025, tni025), pt6, markersize = ps6, label = r"$0.25 \rho_0$")
	ax6.plot(sqrts025-0*Eth, TKN2(tpr025_smooth, tpi025_smooth, tnr025_smooth, tni025_smooth), '-', markersize = ps6, label = r"$0.25 \rho_0$ (smooth)")
	ax6.plot(sqrts05-0*Eth, TKN2(tpr05, tpi05, tnr05, tni05), pt2, markersize = ps6, label = r"$0.50 \rho_0$")
	ax6.plot(sqrts05-0*Eth, TKN2(tpr05_smooth, tpi05_smooth, tnr05_smooth, tni05_smooth), '-', markersize = ps6, label = r"$0.25 \rho_0$ (smooth)")
	ax6.plot(sqrts075-0*Eth, TKN2(tpr075, tpi075, tnr075, tni075), pt6, markersize = ps6, label = r"$0.75 \rho_0$")
	ax6.plot(sqrts075-0*Eth, TKN2(tpr075_smooth, tpi075_smooth, tnr075_smooth, tni075_smooth), '-', markersize = ps6, label = r"$0.25 \rho_0$ (smooth)")
	ax6.plot(sqrts10-0*Eth, TKN2(tpr10, tpi10, tnr10, tni10), pt6, markersize = ps6, label = r"$1.00 \rho_0$")
	ax6.set_xlabel(r"$\sqrt{s}$ (MeV)")
	ax6.set_ylabel(r"$|T_{KN}|^2$ ($\mathrm{MeV}^{-2}$)")
	#ax6.set_xlim(0,6)
	#ax6.set_ylim(-210,10)
	ax6.tick_params(which='both',direction='in',top=True,right=True)
	ax6.legend()
	plt.title('Smooth all')
	#plt.axvline(x=0, linestyle='--', linewidth=0.4, color='grey')
	plt.tight_layout()
	
	
	plt.show()
	
main()
