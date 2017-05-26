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

def FLUXextraction_3D(data_dictionnary):
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
	
	
	print 'Starting FLUXextraction_datacube...'
	# Unpack the dictionnary--------------------------------------------
		
	fits_path = data_dictionnary['fits_path']
	central_wavelength = data_dictionnary['central_wavelength']
	range_wavelength = data_dictionnary['range_wavelength']
	range_decra = data_dictionnary['range_decra']
	line_name = data_dictionnary['name']
		
	print 'Opening fits...'
	fits_main = fits.open(fits_path, mode='update')
	fits_main.info()
	fits_name = fits_path.split('/')
	fits_name = fits_name[-1]
	
	print 'Getting data...'
	image_data = np.nan_to_num(fits.getdata(fits_path, ext=0))
	
	# ------------------------------------------------------------------
	
	if range_decra%2 != 0:
		range_decra = range_decra - 1
	
	print 'Horizontal spectra'
	print 'Getting data from header...'
	
	axis_data = {}
	axis_data[1] = {}
	axis_data[1]['ctype'] = fits_main[0].header['CTYPE1']
	axis_data[1]['crpix'] = int(fits_main[0].header['CRPIX1'])
	axis_data[1]['crval'] = float(fits_main[0].header['CRVAL1'])
	axis_data[1]['cdelt'] = float(fits_main[0].header['CDELT1'])
	axis_data[1]['cunit'] = fits_main[0].header['CUNIT1']
	axis_data[1]['naxis'] = int(fits_main[0].header['NAXIS1'])
		
	axis_data[2] = {}
	axis_data[2]['ctype'] = fits_main[0].header['CTYPE2']
	axis_data[2]['crpix'] = int(fits_main[0].header['CRPIX2'])
	axis_data[2]['crval'] = float(fits_main[0].header['CRVAL2'])
	axis_data[2]['cdelt'] = float(fits_main[0].header['CDELT2'])
	axis_data[2]['cunit'] = fits_main[0].header['CUNIT2']
	axis_data[2]['naxis'] = int(fits_main[0].header['NAXIS2'])

	axis_data[3] = {}
	axis_data[3]['ctype'] = fits_main[0].header['CTYPE3']
	axis_data[3]['crpix'] = int(fits_main[0].header['CRPIX3'])
	axis_data[3]['crval'] = float(fits_main[0].header['CRVAL3'])
	axis_data[3]['cdelt'] = float(fits_main[0].header['CDELT3'])
	axis_data[3]['cunit'] = fits_main[0].header['CUNIT3']
	axis_data[3]['naxis'] = int(fits_main[0].header['NAXIS3'])
	
	print axis_data
		
	print 'Reference pixel X (Projected Right Ascension) = ', axis_data[2]['crpix']
	print 'Reference pixel Y (Projected Declination)     = ', axis_data[1]['crpix']
	print 'Reference pixels flux = ', image_data[central_wavelength-1][axis_data[2]['crpix']-1][axis_data[1]['crpix']-1]
		

	z_array = []
	
		
	for i in range(-int(range_wavelength/2),int(range_wavelength/2)+1) :
		z_array.append(central_wavelength-1+i)
		
	print len(z_array)
	
	zlow = central_wavelength-1-int(range_wavelength/2)
	zup = central_wavelength-1+int(range_wavelength/2)+1 
	
	print 'z-axis |', zlow,'|', zup,'|'
	
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
	
	#--------------------------------------FLUX EXTRACTION--------------	
	Flux = []
	#image_data = fits_main[0].data
	Flux_ref = []
	
	for z in z_array :
		arrayflux = []
		totalflux = []
		for i in range(range_decra):
			arrayflux.extend(image_data[z][(axis_data[2]['crpix']-1-int(range_decra/2)):(axis_data[2]['crpix']-1+int(range_decra/2))][i][(axis_data[1]['crpix']-1-int(range_decra/2)):(axis_data[1]['crpix']-1+int(range_decra/2))])
			totalflux.extend(image_data[z][5:50][i][5:50])
		#print arrayflux
		Flux_ref.append(np.sum(totalflux))
		Flux.append(np.sum(arrayflux))
		#input('test')
		#print Flux
		
	for f in Flux_ref:
		if math.isnan(f):
			# NaN values are replaced by 0 value
			# Feel free to change this !
			f=0
	
	plt.ion()
	
	main_fig, (fig1, fig2) = plt.subplots(nrows = 2, ncols=1)
	main_fig.subplots_adjust(hspace=0.5)
	plt.suptitle('Flux calibration')
	
	fig1.plot(wavelength[zlow:zup],Flux,'k',label='Flux')
	fig1.set_xlabel('Wavelength ($\mu m$)')
	fig1.set_ylabel('Flux (Arbitrary Units)')
	
	fig2.imshow(np.sum(image_data,axis=0),cmap='gist_heat',origin='lower',aspect='auto',interpolation='none')
		
	hdu = fits.PrimaryHDU()
	hdu.data = Flux
	
	hdu.header['BITPIX']=-32
	hdu.header['NAXIS']=1
	hdu.header['NAXIS1']=len(Flux)
	
	hdu.header['CRPIX1']=1
	hdu.header['CRVAL1']=wavelength[zlow]
	hdu.header['CDELT1']=axis_data[3]['cdelt']
	hdu.header['CD1_1']=axis_data[3]['cdelt']
	hdu.header['CUNIT1']=axis_data[3]['cunit']
	hdu.header['DISPAXIS']=1
	hdu.header['AIRMASS']=fits_main[0].header['HIERARCH ESO TEL AIRM START']
	hdu.header['EXPTIME']=float(fits_main[0].header['HIERARCH ESO DET NDIT'])*float(fits_main[0].header['HIERARCH ESO DET DIT'])
	#
	hdu.header['OBSERVED OBJECT NAME']=fits_main[0].header['HIERARCH ESO OBS TARG NAME']
	hdu.writeto('./products/flux_extraction/output.fits',output_verify='exception', overwrite=True, checksum=True)
		
		
		
		

		
		
		
	
