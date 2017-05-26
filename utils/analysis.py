def SAanalysis(data_dictionnary):
	
	
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
	
	path1 = data_dictionnary['path1'] 
	Dbtn = data_dictionnary['Sing_or_Comp']
	Wbtn = data_dictionnary['Wave_or_velo'] 
	
	
	print 'Opening file...'
	SAfile1 = open(path1,'r')
	SAfile1.readline()
	reader = csv.reader(SAfile1, delimiter = '\t')	
	SA_signal1 = []
	Wavelength1 = []
	Velocity1 = []
	FWHM1 = []
	image_data1 = []
	print 'Reading file...'
	for line in reader:
		#print line
		Velocity1.append(float(line[0]))
		Wavelength1.append(float(line[1]))
		image_data1.append(float(line[2]))
		SA_signal1.append(float(line[3]))
		FWHM1.append(float(line[4]))
	print 'Done'
	print len(Wavelength1)
	SAfile1.close()
	if Dbtn=='Comp':
		path2 = data_dictionnary['path2'] 
		print 'Opening file 2...'
		SAfile2 = open(path2,'r')
		SAfile2.readline()
		reader = csv.reader(SAfile2, delimiter = '\t')	
		SA_signal2 = []
		Wavelength2 = []
		Velocity2 = []
		FWHM2 = []
		image_data2 = []
		print 'Reading file...'
		for line in reader:
			#print line
			Velocity2.append(float(line[0]))
			Wavelength2.append(float(line[1]))
			image_data2.append(float(line[2]))
			SA_signal2.append(float(line[3]))
			FWHM2.append(float(line[4]))
		SAfile2.close()
		print 'Done'
		print len(Wavelength2)
	plt.ion()
	plt.figure()
	if Dbtn == 'Comp':
		plt.suptitle('Spectro-astrometric extracted signal :\n'+path1.split('/')[-1]+'\n'+path2.split('/')[-1])
	else:
		plt.suptitle('Spectro-astrometric extracted signal :\n'+path1.split('/')[-1])	
	if Wbtn == 'V':
		fig = plt.subplot(311)
		#plt.title('Spectro-astrometric Signal')
		plt.plot(Velocity1,image_data1,label='Flux 1')
		plt.ylabel('Flux')
		if Dbtn=='Comp':
			plt.plot(Velocity2,image_data2,label='Flux 2')
		plt.legend(loc=2)
		fig = plt.subplot(312)
		plt.plot(Velocity1,SA_signal1,label='Signal 1')
		if Dbtn=='Comp':
			plt.plot(Velocity2,SA_signal2,label='Signal 2')
		plt.xlabel('Velocity (km/s)')
		plt.ylabel('Centroid (Arcsec)')
		plt.legend(loc=2)
		fig = plt.subplot(313)
		plt.plot(Velocity1,FWHM1,label='FWHM 1')
		plt.xlabel('Velocity (km/s)')
		plt.ylabel('FWHM (Arcsec)')
		if Dbtn=='Comp':
			plt.plot(Velocity2,FWHM2,label='FWHM 2')
		plt.legend(loc=2)
	else:
		fig = plt.subplot(311)
		#plt.title('Spectro-astrometric Signal')
		plt.plot(Wavelength1,image_data1,label='Flux 1')
		if Dbtn=='Comp':
			plt.plot(Wavelength2,image_data2,label='Flux 2')
		plt.ylabel('Flux')
		plt.legend(loc=2)
		fig = plt.subplot(312)
		plt.plot(Wavelength1,SA_signal1,label='Signal 1')
		if Dbtn=='Comp':
			plt.plot(Wavelength2,SA_signal2,label='Signal 2')
		plt.xlabel('Wavelength (Angstroms)')
		plt.ylabel('Centroid (Arcsec)')
		plt.legend(loc=2)
		fig = plt.subplot(313)
		plt.plot(Wavelength1,FWHM1,label='FWHM 1')
		plt.xlabel('Wavelength (Angstroms)')
		plt.ylabel('FWHM (Arcsec)')
		if Dbtn=='Comp':
			plt.plot(Wavelength2,FWHM2,label='FWHM 2')
		plt.legend(loc=2)	
			
		

