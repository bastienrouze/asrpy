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

def SADatacube(data_dictionnary):
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
	from astropy.coordinates import SkyCoord
	
	from astropy.wcs import WCS
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
	Wave_or_Velo = data_dictionnary['Wave_or_Velo']
	central_wavelength = data_dictionnary['central_wavelength']
	range_wavelength = data_dictionnary['range_wavelength']
	range_decra = data_dictionnary['range_decra']
	refline = data_dictionnary['refline']
	line_name = data_dictionnary['line_name']
	print 'Starting SA extractor...'
	print 'Opening fits...'
	fits_main = fits.open(fits_path, mode='update')
	fits_main.info()
	fits_name = fits_path.split('/')
	fits_name = fits_name[-1]
	
	print 'Getting data...'
	image_data = fits.getdata(fits_path, ext=0)
	
	# ------------------------------------------------------------------
	
	if range_decra%2 != 0:
		range_decra = range_decra - 1
	
	print 'Horizontal spectra'
	print 'Getting data from header...'
	coord_radec_fk5 = SkyCoord(ra=fits_main[0].header['RA']*u.degree, dec=fits_main[0].header['DEC']*u.degree, frame='fk5')
	coord_galactic = coord_radec_fk5.galactic
	print coord_galactic
	axis_data = {}
	axis_data[1] = {}
	axis_data[1]['ctype'] = fits_main[0].header['CTYPE1']
	if 'WAVE'  in axis_data[1]['ctype']:
		wave_axis = 1
	if 'RA'  in axis_data[1]['ctype']:
		ra_axis = 1
	if 'DEC'  in axis_data[1]['ctype']:
		dec_axis = 1
	axis_data[1]['crpix'] = int(fits_main[0].header['CRPIX1'])
	axis_data[1]['crval'] = float(fits_main[0].header['CRVAL1'])
	axis_data[1]['cdelt'] = float(fits_main[0].header['CDELT1'])
	axis_data[1]['cunit'] = fits_main[0].header['CUNIT1']
	axis_data[1]['naxis'] = int(fits_main[0].header['NAXIS1'])
		
	axis_data[2] = {}
	axis_data[2]['ctype'] = fits_main[0].header['CTYPE2']
	if 'WAVE'  in axis_data[2]['ctype']:
		wave_axis = 2
	if 'RA'  in axis_data[2]['ctype']:
		ra_axis = 2
	if 'DEC'  in axis_data[2]['ctype']:
		dec_axis = 2
	axis_data[2]['crpix'] = int(fits_main[0].header['CRPIX2'])
	axis_data[2]['crval'] = float(fits_main[0].header['CRVAL2'])
	axis_data[2]['cdelt'] = float(fits_main[0].header['CDELT2'])
	axis_data[2]['cunit'] = fits_main[0].header['CUNIT2']
	axis_data[2]['naxis'] = int(fits_main[0].header['NAXIS2'])
	
	axis_data[3] = {}
	axis_data[3]['ctype'] = fits_main[0].header['CTYPE3']
	if 'WAVE'  in axis_data[3]['ctype']:
		wave_axis = 3
	if 'RA'  in axis_data[3]['ctype']:
		ra_axis = 3
	if 'DEC'  in axis_data[3]['ctype']:
		dec_axis = 3
	axis_data[3]['crpix'] = int(fits_main[0].header['CRPIX3'])
	axis_data[3]['crval'] = float(fits_main[0].header['CRVAL3'])
	axis_data[3]['cdelt'] = float(fits_main[0].header['CDELT3'])
	axis_data[3]['cunit'] = fits_main[0].header['CUNIT3']
	axis_data[3]['naxis'] = int(fits_main[0].header['NAXIS3'])
	
	print axis_data
	
	
	print 'Reference pixel X (Projected Right Ascension) = ', axis_data[2]['crpix']
	print 'Reference pixel Y (Projected Declinaison)     = ', axis_data[1]['crpix']
	print 'X pix to deg scale = ', axis_data[2]['cdelt'] 
	print 'Y pix to deg scale = ', axis_data[1]['cdelt'] 
	
	print 'Reference pixels flux = ', image_data[central_wavelength-1][axis_data[2]['crpix']-1][axis_data[1]['crpix']-1]
	
	
	
	g_init = models.Gaussian2D(amplitude = image_data[central_wavelength,axis_data[2]['crpix'],axis_data[1]['crpix']], x_mean = axis_data[ra_axis]['crpix'], y_mean = axis_data[dec_axis]['crpix'],x_stddev = 1.0,y_stddev = 1.0)
	fit_g = fitting.LevMarLSQFitter()
	
	G = []
	X_centroid = []
	Y_centroid = []
	x_array = []
	y_array = []
	z_array = []
	
	for i in range(-int(range_decra/2),int(range_decra/2)+1) :
		x_array.append(axis_data[1]['crpix']-1+i)
	
	for i in range(-int(range_decra/2),int(range_decra/2)+1)  :
		y_array.append(axis_data[2]['crpix']-1+i)
	
	for i in range(-int(range_wavelength/2),int(range_wavelength/2)+1) :
		z_array.append(central_wavelength-1+i)
		
	print len(z_array)
	xlow = axis_data[1]['crpix']-1-int(range_decra/2)
	xup = axis_data[1]['crpix']-1+int(range_decra/2)
	
	ylow = axis_data[2]['crpix']-1-int(range_decra/2)
	yup = axis_data[2]['crpix']-1+int(range_decra/2)
	
	zlow = central_wavelength-1-int(range_wavelength/2)
	zup = central_wavelength-1+int(range_wavelength/2)+1 #will be explained
	
	print '       | px | px |'
	print 'x-axis |', xlow,'|', xup,'|'
	print 'y-axis |', ylow,'|', yup,'|'
	print 'z-axis |', zlow,'|', zup,'|'
	
	x,y = np.mgrid[xlow:xup,ylow:yup]
	
	
	
	
	wavelength = range(axis_data[3]['naxis'])
	wavelength[axis_data[3]['crpix']-1]=axis_data[3]['crval'] 
	i = (axis_data[3]['crpix']-1)-1
	while i >= 0 :
		wavelength[i]=wavelength[i+1]-axis_data[3]['cdelt'] 
		i = i-1
	i = (axis_data[3]['crpix']-1)+1
	while i < axis_data[3]['naxis']:
		wavelength[i]=wavelength[i-1]+axis_data[3]['cdelt'] 
		i = i+1
	
	
	
	
	#y, x = np.mgrid[axis_data[2]['crpix']-1-15:axis_data[2]['crpix']-1+16, axis_data[1]['crpix']-1-15:axis_data[1]['crpix']-1+16]
	#print x
	
		#if set3Dplot == 'Y':
		#X,Y = np.meshgrid(x_table,ROI)
		#X = np.reshape(X,len(x_table)*len(ROI))
		#Y = np.reshape(Y,len(x_table)*len(ROI))
		#print "\n Length of X = ",len(X)
		#print "\n Length of Y = ",len(Y)
			
		#G = np.array(G)
		#print "\n Length of G = ",len(G)
		#print "\n Length of G[0] = ",len(G[0])
		#G = np.reshape(G,len(x_table)*len(ROI))
		#ax.plot_trisurf(X,Y,G,cmap=cm.jet)
		#ax.plot(X_centroid,tabwave,g(X_centroid))
	
	
	
	
	Flux = []
	
	#--------------------------------------SPECTROASTROMETRY------------
	
	lightspeed=299792.458 # c
	velocity = []
	w0 = wavelength[0]
	RVC_B = fits_main[0].header['HIERARCH ESO QC VRAD BARYCOR']  #Radial Velocity Correction Barycentric
	RVC_H = fits_main[0].header['HIERARCH ESO QC VRAD HELICOR']  #Radial Velocity Correction Heliocentric
	l = coord_galactic.l.radian#*(math.pi/180)   # Galactic Longitude of Y*O Oph1
	b = coord_galactic.b.radian#*(math.pi/180)   # Galactic Lattitude of Y*O Oph1
	RVC_lb = 9*math.cos(l)*math.cos(b)+12*math.sin(l)*math.cos(b)+7*math.sin(b)   #Radial Velocity Correction LSR
	RVC = RVC_H + RVC_lb # Total correction
	while w0 < wavelength[-1] :
		velocity.append(((w0 - refline)/refline)*lightspeed+RVC)
		w0 = w0 + axis_data[3]['cdelt'] 
	
	print 'Central wavelength = ', wavelength[central_wavelength-1], axis_data[3]['cunit']
	print 'Reference wavelength = ', refline, axis_data[3]['cunit']
	
	print'Radial velocity correction (RVC) = ', RVC, ' km/s'
	
	declinaison = range(axis_data[1]['naxis'])
	declinaison[axis_data[1]['crpix']-1]=axis_data[1]['crval'] 
	i = (axis_data[1]['crpix']-1)-1
	while i >= 0 :
		declinaison[i]=declinaison[i+1]-axis_data[1]['cdelt'] 
		i = i-1
	i = (axis_data[1]['crpix']-1)+1
	while i < axis_data[1]['naxis']:
		declinaison[i]=declinaison[i-1]+axis_data[1]['cdelt'] 
		i = i+1
	
	rightascension = range(axis_data[2]['naxis'])
	rightascension[axis_data[2]['crpix']-1]=axis_data[2]['crval'] 
	i = (axis_data[2]['crpix']-1)-1
	while i >= 0 :
		rightascension[i]=rightascension[i+1]-axis_data[2]['cdelt'] 
		i = i-1
	i = (axis_data[2]['crpix']-1)+1
	while i < axis_data[2]['naxis']:
		rightascension[i]=rightascension[i-1]+axis_data[2]['cdelt'] 
		i = i+1
	
	
	
	
	for z in z_array :
		Zdata = []
		for i in range(range_decra):
			Zdata.append(image_data[z][ylow+i][xlow:xup])
		g = fit_g(g_init,x,y,Zdata)
		Flux.append(image_data[z][axis_data[2]['crpix']-1][axis_data[1]['crpix']-1])
		X_centroid.append(g.x_mean.value)
		Y_centroid.append(g.y_mean.value)
	
	#print len(X_centroid)
	#print len(wavelength[zlow:zup])
	
	XTilt = np.poly1d(np.polyfit(z_array, X_centroid, 1))(z_array)
	YTilt = np.poly1d(np.polyfit(z_array, Y_centroid, 1))(z_array)
	
	X_centroid2 , Y_centroid2 = [], []
	
	for i in range(len(X_centroid)):
		X_centroid2.append(-1*(X_centroid[i] - XTilt[i])*axis_data[1]['cdelt']*3600*1000) # pixel to deg to arcsec  #(-1) is here just because of the strange way of reading data...
		Y_centroid2.append(-1*(Y_centroid[i] - YTilt[i])*axis_data[2]['cdelt']*3600*1000)
	
	
	Xhalflength = int(len(X_centroid2)/2)
	Yhalflength = int(len(Y_centroid2)/2)
	Noise_arrayX1 = X_centroid2[0:Xhalflength-15]
	Noise_arrayY1 = Y_centroid2[0:Yhalflength-15]
	Noise_arrayX1.extend(X_centroid2[Xhalflength+14:])
	Noise_arrayY1.extend(Y_centroid2[Yhalflength+14:])
	
	Noise_meanX = sum(Noise_arrayX1)/float(len(Noise_arrayX1))
	Noise_sigmaX = np.sqrt(sum([(x-Noise_meanX)**2 for x in Noise_arrayX1])/float(len([(x-Noise_meanX)**2 for x in Noise_arrayX1])))
	
	Noise_meanY = sum(Noise_arrayY1)/float(len(Noise_arrayY1))
	Noise_sigmaY = np.sqrt(sum([(x-Noise_meanY)**2 for x in Noise_arrayY1])/float(len([(x-Noise_meanY)**2 for x in Noise_arrayY1])))
	#print len(velocity[zlow:zup])
	#print len(X_centroid2)

	
	
	
	
	SA_meanX = sum(X_centroid2)/float(len(X_centroid2))
	SA_sigmaX = np.sqrt(sum([(x-SA_meanX)**2 for x in X_centroid2])/float(len([(x-SA_meanX)**2 for x in X_centroid2])))
	
	SA_meanY = sum(Y_centroid2)/float(len(Y_centroid2))
	SA_sigmaY = np.sqrt(sum([(x-SA_meanY)**2 for x in Y_centroid2])/float(len([(x-SA_meanY)**2 for x in Y_centroid2])))
	
	
	SA_sigmaR2 = np.sqrt((SA_sigmaX**2)+(SA_sigmaY**2))
	SA_sigmaR = np.sqrt((Noise_sigmaX**2)+(Noise_sigmaY**2))
	
	countX = 0
	countY = 0
	countR = 0
	for i in range(len(X_centroid2)):
		#print X_centroid2[i]
		if abs(X_centroid2[i]) <= Noise_sigmaX:
			countX = countX+1
		if abs(Y_centroid2[i]) <= Noise_sigmaY:
			countY = countY+1
		if (X_centroid2[i]**2)+(Y_centroid2[i]**2) <= SA_sigmaR**2:
			countR  = countR +1
			
	resX = float(countX)/len(X_centroid2)
	resY = float(countY)/len(Y_centroid2)
	resR = float(countR)/len(Y_centroid2)
	print 'Ratio (X<SigmaX , Y<SigmaY , R^2<SigmaR^2) ', resX, resY, resR
		
	#print 'Mean (X,Y)    ',SA_meanX, SA_meanY
	#print 'Sigma (X,Y)    ',SA_sigmaX, SA_sigmaY
	print 'Sigma Noise (X,Y) ',Noise_sigmaX,Noise_sigmaY
	print 'Sigma (Noise) (R)    ',SA_sigmaR
	
	
	
	
	cm = plt.cm.get_cmap('RdYlBu')
	plt.ion()
	#-------------------------------------------------------------------
	
	plt.figure()
	plt.imshow(image_data[central_wavelength-1],aspect='auto',origin='lower',interpolation='none')
	#plt.contourf(image_data[central_wavelength-1],origin='lower',extent=[-32, 48, -33, 41])
	#plt.contourf(image_data[central_wavelength-1],origin='lower',extent=[declinaison[0]-axis_data[1]['crval'], declinaison[-1]-axis_data[1]['crval'], rightascension[0]-axis_data[2]['crval'], rightascension[-1]-axis_data[2]['crval']])
	plt.colorbar()
	#-------------------------------------------------------------------
	
	x_axis = wavelength
	x_axis_label = 'Wavelength ('+axis_data[wave_axis]['cunit']+')'
	if Wave_or_Velo == 'V':
		x_axis = velocity
		x_axis_label = 'LSR Velocity (km/s)'
	
	main_fig, (fig1, fig2, fig3) = plt.subplots(nrows = 3, ncols=1)
	main_fig.subplots_adjust(hspace=1.0)
	plt.suptitle('Spectro-astrometry of '+line_name)
	
	fig3.plot(x_axis[zlow:zup],Flux,'k',label='Flux')
	fig3.set_xlabel(x_axis_label)
	fig3.set_ylabel('Flux (Arbitrary Units)')
	
	fig1.plot(x_axis[zlow:zup],X_centroid,'b-',label='DEC')
	fig1.plot(wavelength[zlow:zup],XTilt,'k--',label='DEC-Continuum')
	fig1.set_xlabel(x_axis_label)
	fig1.set_ylabel('DEC centroids (pixel)')
	fig1.legend(loc=2)
	
	fig2.plot(x_axis[zlow:zup],Y_centroid,'g-',label='RA')
	fig2.plot(x_axis[zlow:zup],YTilt,'k--',label='RA-Continuum')
	fig2.set_xlabel(x_axis_label)
	fig2.set_ylabel('RA centroids (pixel)')
	fig2.legend(loc=2)
	
	#-------------------------------------------------------------------
	
	main_fig2, (fig21, fig23) = plt.subplots(nrows = 2, ncols=1)
	main_fig2.subplots_adjust(hspace=0.5)
	plt.suptitle('Spectro-astrometry of '+line_name)
	
	fig21.plot(x_axis[zlow:zup],Flux,'k',label='Flux')
	fig21.set_xlabel(x_axis_label)
	fig21.set_ylabel('Flux (Arbitrary Units)')
	
	fig23.plot(x_axis[zlow:zup],X_centroid2,'b-',label='DEC ')
	fig23.plot(x_axis[zlow:zup],Y_centroid2,'g-',label='RA')
	fig23.set_xlabel(x_axis_label)
	fig23.set_ylabel('Offset (mas)')
	fig23.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, borderaxespad=0.)
	
	#fig21.plot(velocity[zlow:zup],X_centroid2,'b-',label='DEC')
	#fig21.set_xlabel('LSR Velocity (km/s)')
	#fig21.set_ylabel('DEC centroids (mas)')
	#fig21.legend(loc=2)
	
	#fig22.plot(velocity[zlow:zup],Y_centroid2,'g-',label='RA')
	
	#fig22.set_xlabel('LSR Velocity (km/s)')
	#fig22.set_ylabel('RA centroids (mas)')
	#fig22.legend(loc=2)
	
	#-------------------------------------------------------------------
		
	plt.figure()	
	
	plt.title('Spectro-astrometric map of '+line_name)
	plt.scatter(Y_centroid2,X_centroid2, c=x_axis[zlow:zup])
		
	plt.gca().invert_xaxis()
	
	xarrow = input('Please enter x arrow East-North')
	yarrow = input('Please enter y arrow East-North')
	
	darrow = input('Please enter dx (dy) for arrow East-North')
	
	plt.axis('equal')
	plt.arrow(xarrow,yarrow,darrow,0,width=0.0001,head_width=0.1)
	plt.text(xarrow+darrow+0.6,yarrow-0.1,'E',fontdict=None, withdash=False)
	plt.arrow(xarrow,yarrow,0,darrow,width=0.0001,head_width=0.1)
	plt.text(xarrow+0.1,yarrow+darrow+0.5,'N',fontdict=None, withdash=False)
	plt.colorbar(label=x_axis_label)
	plt.ylabel('DEC (mas)')
	plt.xlabel('RA (mas)')
	fig = plt.gcf()
	ax = fig.gca()
	circ = plt.Circle((0, 0), radius=SA_sigmaR, color='g',fill= False)
	ax.add_artist(circ)
		
	
	
			
		
		
		

		
		
		
	
