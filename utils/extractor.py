#!/usr/bin/env python
# -*- coding: utf-8 -*-

#  Copyright 2017 Bastien Rouzé
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  




def SAextractor_slice(data_dictionnary):
	print 'Importing matplotlib...'
	import matplotlib as mpl	
	print("current backend is %s" % mpl.get_backend())
	mpl.use('TkAgg')
	print("switched backend is %s" % mpl.get_backend())
	
	
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	from matplotlib import cm
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	from matplotlib.figure import Figure
		
	
		
	print 'Importing numpy...'
	import numpy as np
	print 'Importing astropy...'
	from astropy.io import fits
	from astropy import units as u
	from astropy import coordinates as coord
	from astropy.wcs import WCS
	from astropy.io import fits
	from astropy.utils.data import get_pkg_data_filename
	from astropy.modeling import models, fitting
	from astropy.stats import sigma_clip
	
	print 'Importing scipy...'
	from scipy.interpolate import interp1d
	from scipy.optimize import curve_fit
	from scipy.interpolate import UnivariateSpline

	print 'Importing sys and tk...'
	import sys
	import commands
	import math
	import csv
	from functools import partial
		
	import tkMessageBox
	import tkFileDialog
	
	# Unpack the dictionnary--------------------------------------------
		
	fits_path=data_dictionnary['fits_path']
	slice_wavelength=data_dictionnary['slice_wavelength']
	slice_mode=data_dictionnary['slice_mode']
	print slice_mode
	Hori_or_Verti=data_dictionnary['Hori_or_Verti']
	SA_is_running=data_dictionnary['SAextrator_running']
	
	print 'Starting SA extractor_slice...'
	print 'Opening fits...'
	fits_main = fits.open(fits_path)
	fits_main.info()
	fits_name = fits_path.split('/')
	fits_name= fits_name[-1]
	print 'Getting data...'
	image_data = fits.getdata(fits_path, ext=0)
	if Hori_or_Verti == 'Hori':
		print 'Horizontal spectra'
		print 'Getting data from header...'
		crval1 = fits_main[0].header['CRVAL1']
		cdelt1 = fits_main[0].header['CDELT1']
		naxis1 = int(fits_main[0].header['NAXIS1'])
		cmax1 = crval1 + naxis1*cdelt1
		crval2 = int(fits_main[0].header['CRVAL2'])
		cdelt2 = int(fits_main[0].header['CDELT2'])
		naxis2 = int(fits_main[0].header['NAXIS2'])
		cmax2 = crval2 + naxis2*cdelt2 -1
	
	Instrument = fits_main[0].header['INSTRUME']
	InstruType = fits_main[0].header['HIERARCH ESO TPL NAME']
	
	print '-----------------'
	print 'CRVAL1 = ', crval1
	print 'CDELT1 = ',cdelt1
	print 'Maximum CMAX1 =', cmax1
	print 'CRVAL2 = ', crval2
	print 'CDELT2 = ', cdelt2
	print 'Maximum CMAX2 =', cmax2
	print '-----------------'
	print 'Instrument = ',Instrument
	print 'Type = ',InstruType
	print '-----------------'
	slice_wavelength = int((float(slice_wavelength)-crval1)/cdelt1)
	x_table = range(int(cmax2))
	g_init = models.Gaussian1D(amplitude=image_data[int(cmax2/2),slice_wavelength], mean=cmax2/2, stddev=1.)
	fit_g = fitting.LevMarLSQFitter()
	g = fit_g(g_init, x_table, image_data[:,slice_wavelength])
	FWHM = (g.stddev.value)*2*np.sqrt(2*np.log(2))
	plt.ion()
	if SA_is_running:
		plt.figure()
	if slice_mode == "Sing":
		print 'Mode one by one'
		plt.figure()
	plt.title('Gaussian fitting of spatial profile')
	plt.plot(x_table, image_data[:,slice_wavelength], '+',label='Data')
	plt.plot(x_table, g(x_table), '-',label='Gaussian Fit')
	plt.plot([g.mean.value,g.mean.value], [0,g(g.mean.value)], '--',label='Mean = '+str(g.mean.value),marker='o')
	plt.plot([x_table[0],x_table[-1]],[g(g.mean.value)/2,g(g.mean.value)/2], 'k--',label='FWHM = '+str(FWHM))
	plt.xlabel('Spatial profile (Pixel)')
	plt.ylabel('Flux (/)')
	#plt.legend(loc=2)
	

def SAextractor(data_dictionnary):


	print 'Importing matplotlib...'
	import matplotlib as mpl	
	print("current backend is %s" % mpl.get_backend())
	mpl.use('TkAgg')
	print("switched backend is %s" % mpl.get_backend())
	
	
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	from matplotlib import cm
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	from matplotlib.figure import Figure
		
	
		
	print 'Importing numpy...'
	import numpy as np
	print 'Importing astropy...'
	from astropy.io import fits
	from astropy import units as u
	from astropy import coordinates as coord
	from astropy.wcs import WCS
	from astropy.io import fits
	from astropy.utils.data import get_pkg_data_filename
	from astropy.modeling import models, fitting
	from astropy.stats import sigma_clip
	
	print 'Importing scipy...'
	from scipy.interpolate import interp1d
	from scipy.optimize import curve_fit
	from scipy.interpolate import UnivariateSpline

	print 'Importing sys and tk...'
	import sys
	import commands
	import math
	import csv
	from functools import partial
		
	import tkMessageBox
	import tkFileDialog
	
	# Unpack the dictionnary--------------------------------------------
		
	fits_path = data_dictionnary['fits_path']
	up_roiP = data_dictionnary['up_roiP']
	low_roiP = data_dictionnary['low_roiP']
	up_roiW = data_dictionnary['up_roiW']
	low_roiW = data_dictionnary['low_roiW']
	Wave_or_Velo=  data_dictionnary['Wave_or_Velo']
	Hori_or_Verti = data_dictionnary['Hori_or_Verti']
	issetnoise = data_dictionnary['issetnoise']
	spin_selector_noise = data_dictionnary['spin_selector_noise']
	spatial_offset_choice = data_dictionnary['spatial_offset_choice']
	set3Dplot = data_dictionnary['set3Dplot']
	
		
	# Kernel Layer------------------------------------------------------
	
	def NoiseSetter():
		if (spin_selector_noise== '0' or spin_selector_noise==''):
				print 'No noise regions. Computing the sigma of the full signal...'
				issetnoise = False
				noisearrays = None
				
		else:
			tkMessageBox.showinfo('Noise','You will be asked to fill the noise samples along the spectra. Check the terminal.')
			issetnoise = True
			samplenumber = int(spin_selector_noise)
			noisearrays = range(samplenumber)
			
			for i in range(samplenumber):
				print 'Setting ',i+1,'/',samplenumber, 'noise region.'
				print '/!\ you must restart if you make a mistake. This will be improved soon.'
				lownoiserange = raw_input('Lower bound : ')
				lownoiserange = int((float(lownoiserange)-crval1)/cdelt1)
				
				upnoiserange = raw_input('Upper bound : ')
				upnoiserange = int((float(upnoiserange)-crval1)/cdelt1)
				
				noisearrays[i] = range(lownoiserange,upnoiserange)
		return issetnoise, noisearrays			
	
	print 'Starting SA extractor...'
	print 'Opening fits...'
	fits_main = fits.open(fits_path)
	fits_main.info()
	fits_name = fits_path.split('/')
	fits_name= fits_name[-1]
	
	
	
	print 'Getting data...'
	image_data = fits.getdata(fits_path, ext=0)
	if (int(fits_main[0].header['NAXIS1'])) > (int(fits_main[0].header['NAXIS2'])):
		print 'Horizontal spectra'
		print 'Getting data from header...'
		crval1 = fits_main[0].header['CRVAL1']
		cdelt1 = fits_main[0].header['CDELT1']
		naxis1 = int(fits_main[0].header['NAXIS1'])
		cmax1 = crval1 + naxis1*cdelt1
		crval2 = int(fits_main[0].header['CRVAL2'])
		cdelt2 = int(fits_main[0].header['CDELT2'])
		naxis2 = int(fits_main[0].header['NAXIS2'])
		cmax2 = crval2 + naxis2*cdelt2 -1
	else :
		print 'Vertical spectra'
		print 'Getting data from header...'
		crval1 = fits_main[0].header['CRVAL2']
		cdelt1 = fits_main[0].header['CDELT2']
		naxis1 = int(fits_main[0].header['NAXIS2'])
		cmax1 = crval1 + naxis1*cdelt1
		crval2 = int(fits_main[0].header['CRVAL1'])
		cdelt2 = int(fits_main[0].header['CDELT1'])
		naxis2 = int(fits_main[0].header['NAXIS1'])
		cmax2 = crval2 + naxis2*cdelt2 -1
		image_data=image_data.T
		
	Instrument = fits_main[0].header['INSTRUME']
	InstruType = fits_main[0].header['HIERARCH ESO TPL NAME']
	
	print '-----------------'
	print 'CRVAL1 = ', crval1
	print 'CDELT1 = ',cdelt1
	print 'Maximum CMAX1 =', cmax1
	print 'CRVAL2 = ', crval2
	print 'CDELT2 = ', cdelt2
	print 'Maximum CMAX2 =', cmax2
	print '-----------------'
	print 'Instrument = ',Instrument
	print 'Type = ',InstruType
	print '-----------------'
	
	# Google ESO + instrument name + pixel scale to find it
	# As all data can be processed if I don't know the keyword
	# I will put 1.0. Feel free to add keywords !
	# Further use for displaying with matplotlib
	yc_text = 'Centroid (Arcsec)'
	yf_text = 'FWHM (Arcsec)'
	
	if Instrument == 'UVES':
		if 'Red' in InstruType:
			pixel2angle = 0.16
		if 'Blue' in InstruType:
			pixel2angle = 0.16
		else:
			# Unknow issue I met with red data without keyword
			pixel2angle = 0.16
						
	elif Instrument == 'XSHOOTER':
		pixel2angle = 0.17
			
	else:
		pixel2angle = 1.0
		# If we don't know, it will display in pixel scale
		yc_text = 'Centroid (Pixel)'
		yf_text = 'FWHM (Pixel)'
	
	
	# Setting the wavelength/velocity arrays
	# Both are always computed because other tasks can call
	# A velocity displaying while you wanted here a wavelength
	# displaying for extraction the SA signal
			
	wavelength = []
	velocity = []
	lightspeed=299792.458 # rare class attribute
	refline = 6562.8 # This is the H alpha line in angstroms
	l = crval1 
	
	while l < cmax1 :
		wavelength.append(l)
		velocity.append(((l - refline)/refline)*lightspeed)
		l = l + cdelt1
	
	
	# misc------------
	pixels = []
	p = crval2
	
	while p < cmax2 :
		pixels.append(p)
		p = p + cdelt2
	#-----------------
	
	# Get the spatial pixel region of interest	
	if up_roiP !='':
		P_roi_max = int(up_roiP)
	else:
		P_roi_max = naxis2
	if low_roiP !='':
		P_roi_min = int(low_roiP)
	else:
		P_roi_min = 0
	
	# Compute into width and the datatable		
	P_roi_width = P_roi_max - P_roi_min
	x_table = range(P_roi_min,P_roi_max)
	
	
	# Get the wavelength of region of interest
	if up_roiW !='':
		W_roi_max = float(up_roiW)			
	else:
		W_roi_max = cmax1
	
	if low_roiW !='':
		W_roi_min = float(low_roiW)
	else:
		W_roi_min = crval1
	
	# Convert it into pixel (at lower int value) to compute the image
	W_roi_max = int((float(W_roi_max)-crval1)/cdelt1)
	W_roi_min = int((float(W_roi_min)-crval1)/cdelt1)
		
	print 'Wavelength pixel start = ',W_roi_min
	print 'Wavelength pixel stop = ',W_roi_max
		
	W_roi_width = W_roi_max - W_roi_min 
	sline = range(W_roi_width)  # Spectral LINE array
		
	# Initialisation of the mean (centroid) x position
	ROI = range(W_roi_min,W_roi_max)
	X_centroid = range(W_roi_width)
	FWHM = range(W_roi_width)
	G = range(W_roi_width)
				
	#g_init = models.Gaussian1D(amplitude=1., mean=cmax2/2, stddev=1.)
	fit_g = fitting.LevMarLSQFitter()
	
	#g_init2 = models.Gaussian1D(amplitude=1., mean=cmax2/2, stddev=1.)
	#fit = fitting.LevMarLSQFitter()
	#or_fit = fitting.LevMarLSQFitter()
	#x_table2 = np.array(x_table)
	#sline2 = []
	
	#X2 = range(W_roi_width)
	
	
	# ------------------- Plot 2D Spectrum map -------------------------
	plt.ion()
	
	imdisp = plt.figure()
	plt.subplot(211)
	plt.imshow(image_data[P_roi_min:P_roi_max,(W_roi_min-100):(W_roi_max+100)], extent=[W_roi_min-100,W_roi_max+100,P_roi_min,P_roi_max],aspect='auto',origin='lower',interpolation='none')
	plt.colorbar()
	plt.plot([W_roi_min, W_roi_min],[P_roi_min, P_roi_max], 'k-',linewidth=2.5)
	plt.plot([W_roi_max, W_roi_max],[P_roi_min, P_roi_max], 'k-',linewidth=2.5)
	plt.draw()
	

	#SAextractor core
	# -------------- Astropy method + 'American' method ----------------
	# Astropy / Scipy inspired :: Jets from Young Stars II
	#                             The Spectro-Astrometry
	#
	# American :: https://arxiv.org/pdf/0805.3314.pdf
	Xestimate = []
	Festimate = []
	x_o = 32    # maybe unused	
	j=0 # Count for sline (spectral - line)
	
	
	if set3Dplot == 'Y':
		fig = plt.figure()
		ax = fig.gca(projection='3d')
	
	for l in ROI:
		sys.stdout.write ('Processing pixel :    ' + str(l) + chr(13))
		sys.stdout.flush()
		
		sline[j] = image_data[P_roi_min:P_roi_max,l]
		
		Xe_to_append = 0
		Fe_to_append = 0
		Fe_square_weight_sum = 0
		for i in range(len(x_table)):
			Xe_to_append = Xe_to_append + (x_table[i])*sline[j][i]
			
		Xestimate.append(Xe_to_append/sum(sline[j]))
		
		for i in range(len(x_table)):
			Fe_to_append = Fe_to_append + sline[j][i]*(abs(x_table[i]-Xestimate[j])**2)
			Fe_square_weight_sum = Fe_square_weight_sum + sline[j][i]**2
		
		Fe_correction_coefficient = sum(sline[j])/((sum(sline[j])**2)-Fe_square_weight_sum)
		
		Festimate.append((64/63)*Fe_to_append/(sum(sline[j])))
		
		g_init = models.Gaussian1D(amplitude=sline[j][int(P_roi_width/2)], mean=int(naxis2/2), stddev=1.)
		#fit_g = fitting.LevMarLSQFitter()
		g = fit_g(g_init, x_table, sline[j])
		X_centroid[j] = g.mean.value
						
		if math.isnan(X_centroid[j]):
			# NaN values are replaced by j-1 value
			# Feel free to change this !
			print 'NaN found at ', j,' removed by former value.'
			X_centroid[j] = X_centroid[j-1]
			
		if (X_centroid[j]>P_roi_max or X_centroid[j]<P_roi_min):
			print 'Corrupted data found at ',j,' corresponding to ',l,' with value ',X_centroid[j],' removed by former value ',X_centroid[j-1]
			X_centroid[j] = X_centroid[j-1]
			#plt.figure()
			#plt.plot(x_table,g(x_table), 'b-',label='Corrupted'+str(l))
			#plt.plot(x_table, sline[j], 'gx', label='Data'+str(l))
			#plt.legend(loc=2)
		
		FWHM[j] = (g.stddev.value)*2*np.sqrt(2*np.log(2))
		Festimate[j] = np.sqrt(Festimate[j])*2*np.sqrt(2*np.log(2))
		G[j]=g(x_table)
		
		if set3Dplot == 'Y':
			i=0
			tabwave = []
		#print g(X_centroid[j])
			for x in x_table:
				tabwave.append(l)
		#ax.plot(np.array(x_table),np.array(tabwave),G[j],'b-')
		#ax.plot([X_centroid[j]],[np.array(l)],[g(X_centroid[j])],'r+')
			plt.draw()
		#for x in x_table:
			#print [x], np.array(self.l), self.G[self.j][i]
			#ax.plot([x],np.array(l),G[j][i],'r+')
			#i = i + 1
		j=j + 1
		#loop
	#end of for
	
	
	if set3Dplot == 'Y':
		X,Y = np.meshgrid(x_table,ROI)
		X = np.reshape(X,len(x_table)*len(ROI))
		Y = np.reshape(Y,len(x_table)*len(ROI))
		print "\n Length of X = ",len(X)
		print "\n Length of Y = ",len(Y)
			
		G = np.array(G)
		print "\n Length of G = ",len(G)
		print "\n Length of G[0] = ",len(G[0])
		G = np.reshape(G,len(x_table)*len(ROI))
		ax.plot_trisurf(X,Y,G,cmap=cm.jet)
		ax.plot(X_centroid,tabwave,g(X_centroid))

	
	print '\nSA method Astropy done\nSA method bis done'
	#
	if set3Dplot == 'Y':
		blabla = raw_input('Enter')
	#------------- Scipy method ----------------------------------------
	def func(x, a, b, c):
		return a*np.exp(-(x-b)**2/(2*c**2))
	
	yaxis = naxis2
	central_ap=int(yaxis/2)
	image_data2 = np.transpose(image_data)
	i=0
	centroid = []
	fwhm = []		
	y=np.asarray(range(yaxis))
	
	while i < len(ROI):
		try:
			popt, pcov = curve_fit(func, y-central_ap, image_data2[W_roi_min+i])
			centroid.append(popt[1]+central_ap)
			calc=2*(np.sqrt(2*np.log(2))*popt[2])
			fwhm.append(abs(calc))
		
		except RuntimeError:
			print("Error - curve_fit failed - data replaced by i-1 value")
			centroid.append(centroid[-1])
			fwhm.append(abs(calc))
		
		i += 1
	
	print 'SA method scipy gaussian fitting done\n'
	
	plt.plot(ROI,X_centroid, 'c-',label='Gaussian Kernel ')
	plt.plot(ROI,Xestimate, 'k-',label='Statistical Mean')
	#plt.plot(ROI,centroid, 'g-',label='Scipy Gaussian')
	plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, borderaxespad=0.)
	

	mean_array = []
	sigma_array = []
	

	
	i = 0
	image_centroid = []
	
		
	plt.subplot(212)
	plt.plot(wavelength[W_roi_min:W_roi_max], X_centroid, '-', label='Gaussian Kernel')	

	issetnoise, noisearrays = NoiseSetter()
	if issetnoise:
		Noise_signal = []
		for tab in noisearrays :
			Noise_signal.extend(X_centroid[(tab[0]-W_roi_min):(tab[-1]-W_roi_min)])		
	
	
		
	tiltorder = 1	
	if spatial_offset_choice=='Noise':	
		Tilt = np.poly1d(np.polyfit(ROI[0:len(Noise_signal)], Noise_signal, tiltorder))(ROI) 
		Tilte = np.poly1d(np.polyfit(ROI[0:len(Noise_signal)], Noise_signal, tiltorder))(ROI)      #ROI
	elif spatial_offset_choice=='Middle':
		Tilt = []
		Tilte = []
		for x in X_centroid:
			Tilt.append(int((P_roi_min+P_roi_max)/2))
			Tilte.append(int((P_roi_min+P_roi_max)/2))
	else:
		Tilt = np.poly1d(np.polyfit(ROI, X_centroid, tiltorder))(ROI)
		Tilte = np.poly1d(np.polyfit(ROI, Xestimate, tiltorder))(ROI)
			
	
	
	SA_signal = []
	SA_estimated = []
	
	for i in range(len(X_centroid)):
		SA_signal.append((X_centroid[i]-Tilt[i]))
		SA_estimated.append((Xestimate[i]-Tilte[i]))
		SA_estimated[i] = SA_estimated[i]*pixel2angle
		SA_signal[i] = SA_signal[i]*pixel2angle
		FWHM[i] = FWHM[i]*pixel2angle
		Festimate[i] = Festimate[i]*pixel2angle
		
		if (issetnoise and i<len(Noise_signal)):
			Noise_signal[i] = (Noise_signal[i])*pixel2angle

	SA_mean = sum(SA_signal)/float(len(SA_signal))
	SA_sigma = np.sqrt(sum([(x-SA_mean)**2 for x in SA_signal])/float(len([(x-SA_mean)**2 for x in SA_signal])))
	print "SA mean value = ",SA_mean,
	print "SA sigma = ",SA_sigma
	Median_SA = np.median(SA_signal)
	Median_FWHM = np.median(FWHM)
		
			
	if issetnoise:
		Noisemean = sum(Noise_signal)/float(len(Noise_signal))
		Noisesigma = np.sqrt(sum([(x-Noisemean)**2 for x in Noise_signal])/float(len([(x-Noisemean)**2 for x in Noise_signal])))
		print "Noise mean value = ",Noisemean,
		print "Noise sigma = ",Noisesigma		
		for l in X_centroid:
			mean_array.append(Noisemean)
			sigma_array.append(Noisesigma)
	else:
		for l in SA_signal:
			mean_array.append(SA_mean)
			sigma_array.append(SA_sigma)
	
	
	
	print 'Plotting...'
	

	if issetnoise:
		print "SAsigma/Noisesigma = ", SA_sigma/Noisesigma
	
	print 'Plotting SA Signal from ', W_roi_min, ' to ', W_roi_min+len(SA_signal)
	
	plt.plot(wavelength[W_roi_min:W_roi_max], Tilt, '--', label='Tilt')
	plt.plot(wavelength[W_roi_min:W_roi_max], Xestimate, 'k-', label='Statistical Mean')
	#plt.plot(wavelength[W_roi_min:W_roi_max], centroid, 'g-', label='X method Scipy')
	plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, borderaxespad=0.)	
	plt.legend(loc=2)
	
	i=0
	for x in X_centroid:
		#print 'Centroid ', x
		#print 'Flux ', image_data[int(x),W_roi_min+i]
		#print 'At coordinate ', W_roi_min+i
		image_centroid.append(image_data[int(x),W_roi_min+i])
		i = i+1
		#raw_input('Enter')
	
	#print image_data[47,43669]
	#print image_data[47,43670]	
	
	
	main_fig, (fig1, fig2, fig3) = plt.subplots(nrows = 3, ncols=1)
	plt.suptitle('Spectro-astrometric extracted signal :\n'+fits_name)	
	if Wave_or_Velo=='V':
		#fig1 = plt.subplot(311)
		fig1.plot(velocity[W_roi_min:W_roi_max], image_centroid, '-', label='Flux')
		#fig1.title('Spectro-astrometric Signal')		
		fig1.set_ylabel('Flux(/)')
			
		#fig2 = plt.subplot(312)
		fig2.plot(velocity[W_roi_min:W_roi_max], SA_signal, '-')
		fig2.plot(velocity[W_roi_min:W_roi_max], sigma_array,linestyle='--', c='k',label='$\sigma \simeq$'+str(round(sigma_array[0],5)))
		fig2.plot(velocity[W_roi_min:W_roi_max], [i*-1 for i in sigma_array],linestyle='--', c='k')
		fig2.plot(velocity[W_roi_min:W_roi_max], [i*3 for i in sigma_array],linestyle='--', c='y',label='3$\sigma$ level')
		fig2.plot(velocity[W_roi_min:W_roi_max], [i*-3 for i in sigma_array],linestyle='--', c='y')
		fig2.set_xlabel('Velocity (km/s)')
		fig2.set_ylabel(yc_text)
		fig2.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, borderaxespad=0.)
		#fig3 = plt.subplot(313)
		fig3.plot(velocity[W_roi_min:W_roi_max], FWHM, '-', color='green')
		fig3.set_xlabel('Velocity (km/s)')
		fig3.set_ylabel(yf_text)
		
	else:
		#fig1 = plt.subplot(311)
		fig1.plot(wavelength[W_roi_min:W_roi_max], image_centroid, '-', label='Flux')
		#fig1.title('Spectro-astrometric Signal')
		fig1.set_ylabel('Flux(/)')
			
		#fig2 = plt.subplot(312)
		fig2.plot(wavelength[W_roi_min:W_roi_max], SA_signal, '-')
		fig2.plot(wavelength[W_roi_min:W_roi_max], sigma_array,linestyle='--', c='k',label='$\sigma \simeq$'+str(round(sigma_array[0],5)))
		fig2.plot(wavelength[W_roi_min:W_roi_max], [i*-1 for i in sigma_array],linestyle='--', c='k')
		fig2.plot(wavelength[W_roi_min:W_roi_max], [i*3 for i in sigma_array],linestyle='--', c='y',label='3$\sigma$ level')
		fig2.plot(wavelength[W_roi_min:W_roi_max], [i*-3 for i in sigma_array],linestyle='--', c='y')
		fig2.set_xlabel('Wavelength $\lambda$ (Angstroms)')
		fig2.set_ylabel(yc_text)
		fig2.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, borderaxespad=0.)
	
		
		#fig3 = plt.subplot(313)
		fig3.plot(wavelength[W_roi_min:W_roi_max], FWHM, '-', color='green')
		fig3.set_xlabel('Wavelength $\lambda$ (Angstroms)')
		fig3.set_ylabel(yf_text)
	
	
	
	plt.figure()
	plt.subplot(211)
	plt.plot(wavelength[W_roi_min:W_roi_max], SA_signal, '-',label='Gaussian Kernel')
	plt.xlabel('Wavelength $\lambda$ (Angstroms)')
	plt.ylabel(yc_text)
	plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, borderaxespad=0.)
	plt.subplot(212)
	plt.plot(wavelength[W_roi_min:W_roi_max], SA_estimated, '-',label='Statistical Mean')
	plt.xlabel('Wavelength $\lambda$ (Angstroms)')
	plt.ylabel(yc_text)
	plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, borderaxespad=0.)
	
	
	
	plt.figure()
	plt.subplot(211)
	plt.plot(wavelength[W_roi_min:W_roi_max], SA_signal, '-',label='Gaussian Kernel')
	plt.plot(wavelength[W_roi_min:W_roi_max], SA_estimated, 'r-',label='Statistical Mean')
	plt.xlabel('Wavelength $\lambda$ (Angstroms)')
	plt.ylabel(yc_text)
	plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, borderaxespad=0.)
           
	#plt.subplot(312)
	#plt.plot(wavelength[W_roi_min:W_roi_max], FWHM, '-', color='green')
	#plt.plot(wavelength[W_roi_min:W_roi_max], Festimate, 'k-')
	#plt.xlabel('Wavelength $\lambda$ (Angstroms)')
	#plt.ylabel(yf_text)
	
	plt.subplot(212)
	plt.plot(wavelength[W_roi_min:W_roi_max], np.array(SA_signal)-np.array(SA_estimated), '-', color='green',label='Numerical substraction difference')
	plt.xlabel('Wavelength $\lambda$ (Angstroms)')
	plt.ylabel('Substraction (Arcsec)')
	plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, borderaxespad=0.)
	plt.draw()
	
	
		
	if tkMessageBox.askyesno('Information','Save your work ?'):
		print 'Creating the file'
		#u=commands.getoutput('touch'+' '+fits_path+'.txt')	
		#filesave = open(fits_path+'.txt','w')
		filesave = open(tkFileDialog.asksaveasfilename(title='spectroastrometry_data')+'.txt','w')
		filesave.write('Velocity'+'	'+'Wavelength'+'	'+'Flux'+'	'+'SA_signal'+'	'+'FHWM\n')
		i=0
		print 'Saving...'
		for l in wavelength[W_roi_min:W_roi_max]: 
			filesave.write(str(velocity[W_roi_min+i])+'	'+str(l)+'	'+str(image_centroid[i])+'	'+str(SA_signal[i])+'	'+str(FWHM[i])+'\n')
			i = i+1
		filesave.close()
		tkMessageBox.showinfo('Information','Data saved !')
	else:
		print 'end'
		
	#if tkMessageBox.askyesno('Information','Move to data editor ?'):
	#	root = Tk()
	#	GUISAC = GUI_SAdataeditor(root)
	

	
	
	#Instrument = fits_main[0].header['INSTRUME']
	#InstruType = fits_main[0].header['HIERARCH ESO TPL NAME']
	
	
	
