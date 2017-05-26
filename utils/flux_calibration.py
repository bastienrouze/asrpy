#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  flux_calibration.py
#  
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


def flux_calibration(fits_path_standard,fits_path_object):
	print 'Flux calibration routine'
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
	from astropy.modeling import models, fitting
	from astropy.stats import sigma_clip
	from astropy import units as u
	from astropy.analytic_functions import blackbody_lambda, blackbody_nu
		
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
	
	
	fits_standard = fits.open(fits_path_standard)
	fits_standard.info
	standard_data = fits_standard[0].data
	
	fits_object = fits.open(fits_path_object)
	fits_object.info
	object_data = fits_object[0].data
	
	
	crval_standard = fits_standard[0].header['CRVAL1']
	cdelt_standard = fits_standard[0].header['CDELT1']
	maxpix_standard = len(standard_data)
	
	count = 0
	wavelengths = [crval_standard]
	while count<maxpix_standard-1:
		wavelengths.append(wavelengths[count]+cdelt_standard)
		count = count+1
	print maxpix_standard
	print len(wavelengths)
	plt.ion()
	plt.figure()
	plt.plot(wavelengths,standard_data)
	plt.plot(wavelengths,object_data)
	plt.xlabel('Wavelength (microns)')
	plt.ylabel('Flux (Counts)')
	
	
	# useless----------
	def blackbody_lam(lam, T):
		""" Blackbody as a function of wavelength (m) and temperature (K).

		returns units of W/(cm2 . um)
		"""
		from scipy.constants import h,k,c
		return 10000*2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))
	
	def func(wa,T):
		return blackbody_lam(wa,T)
	# end of useless------------
	
	# magnitude in JHK
	
	mJ = 7.203
	mH = 7.249
	mK = 7.200
		
	
	# flux0(J,H,K) are in W / (cm2 . um)
		
	flux0J = 3.129E-13
	flux0H = 1.133E-13
	flux0K = 4.283E-14
	
	# F(STD) = F0 * 10 ^ (-m / 2.5)
	
	fluxJ = flux0J*10**(-mJ/2.5)
	fluxH = flux0H*10**(-mH/2.5)
	fluxK = flux0K*10**(-mK/2.5)
	
		
	fluxlist = [np.log10(fluxJ),np.log10(fluxH),np.log10(fluxK)]
	
	
	# wavelengths are in um
	
	lambdaJ = 1.235
	lambdaH = 1.662
	lambdaK = 2.159
	
	lambdalist = [np.log10(lambdaJ),np.log10(lambdaH),np.log10(lambdaK)]
		
	p = np.polyfit(lambdalist,fluxlist,1)  
	
	print p  # polynomial y = p[0]*x + p[1] = ax+b
			
	logF_interpolation = np.poly1d(p)(lambdalist)
	
	plt.figure()
	plt.subplot(211)
	plt.plot(lambdalist,fluxlist,'+')
	plt.plot(lambdalist,logF_interpolation)
	plt.xlabel('log$_{10}(\lambda)$')
	plt.ylabel('log$_{10}(Flux)$')
	
	# log10 F = a log10 w + b
	# F = 10**b . w**a
	
	powconstant = 10**p[1]
	blackbody_fit =  []
	f_test = []
	w_range = []
	i = 1.2
	while i<2.4:
		w_range.append(i)
		i=i+0.01
	for w in w_range:
		f_test.append((1e6)*powconstant*(w**p[0]))
	
	
	for w in wavelengths:
		blackbody_fit.append((1e6)*powconstant*(w**p[0]))
	
		
	
	plt.subplot(212)
	plt.plot([lambdaJ,lambdaH,lambdaK],[fluxJ,fluxH,fluxK],'+')
	plt.plot(w_range,f_test,':')
	plt.plot(wavelengths,blackbody_fit)
	plt.xlabel('Wavelength (microns)')
	plt.ylabel('Flux (erg/s/cm2/$\mu$m)')
		
	calib_factor= []
	
	for i in range(maxpix_standard):
		# This loop gives the calibration factor for each wavelength of the 1D spec
		# calib_factor is in Flux/Count => W/(cm2 um count)
		# 
		calib_factor.append(blackbody_fit[i]/standard_data[i])
		
	object_calibrated = []
	
	# 1 W/(cm2 um) = 1 (W/cm2)*(1/um) = 10 000 000 ((erg/s) / cm2)*(1/um) 
	
	for i in range(maxpix_standard):
		object_calibrated.append(calib_factor[i]*object_data[i]*(fits_standard[0].header['EXPTIME']/fits_object[0].header['EXPTIME']))
	
	#
	plt.figure()
	plt.plot(wavelengths,object_calibrated)
	plt.plot(2.159,2.899E-11,'k+')
	plt.plot(1.662,1.774E-11,'k+')
	plt.plot(1.235,0.57E-11,'k+')
	plt.title('Flux calibration')
	plt.xlabel('Wavelength (microns)')
	plt.ylabel('Intensity (erg/s/cm$^2$/$\mu$m)')
	
	hdu = fits.PrimaryHDU()
	hdu.data = calib_factor
		
	hdu.header['BITPIX']=-32
	hdu.header['NAXIS']=1
	hdu.header['NAXIS1']=len(calib_factor)
	
	hdu.header['CRPIX1']=1
	hdu.header['CRVAL1']=wavelengths[0]
	hdu.header['CDELT1']=cdelt_standard
	hdu.header['CD1_1']=cdelt_standard
	hdu.header['CUNIT1']=fits_object[0].header['CUNIT1']
	hdu.header['DISPAXIS']=1

	hdu.header['STANDARD']=fits_standard[0].header['HIERARCH OBSERVED OBJECT NAME']
	hdu.writeto('./products/flux_calibration/fluxcalib_response.fits',output_verify='exception', overwrite=True, checksum=True)
	
	hdu = None 
	
	hdu = fits.PrimaryHDU()
	hdu.data = object_calibrated
		
	hdu.header['BITPIX']=-32
	hdu.header['NAXIS']=1
	hdu.header['NAXIS1']=len(object_calibrated)
	
	hdu.header['CRPIX1']=1
	hdu.header['CRVAL1']=wavelengths[0]
	hdu.header['CDELT1']=cdelt_standard
	hdu.header['CD1_1']=cdelt_standard
	hdu.header['CUNIT1']=fits_object[0].header['CUNIT1']
	hdu.header['DISPAXIS']=1
	
	
	hdu.header['AIRMASS OBJECT']=fits_object[0].header['AIRMASS']
	hdu.header['AIRMASS_STANDARD']=fits_standard[0].header['AIRMASS']
	hdu.header['EXPTIME OBJECT']=fits_object[0].header['EXPTIME']
	hdu.header['EXPTIME STANDARD']=fits_standard[0].header['EXPTIME']
	
	hdu.header['OBJECT']=fits_object[0].header['HIERARCH OBSERVED OBJECT NAME']
	hdu.header['STANDARD']=fits_standard[0].header['HIERARCH OBSERVED OBJECT NAME']
	hdu.writeto('./products/flux_calibration/fluxcalib_output.fits',output_verify='exception', overwrite=True, checksum=True)
	
	print 'end'
	
	
	
	
	
	

    


