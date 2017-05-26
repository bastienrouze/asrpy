#!/usr/bin/env python
# -*- coding: utf-8 -*-

# All import------------------------------------------------------------
from Tkinter import *
from tkMessageBox import *
import subprocess
import tkFileDialog
import Tix
import ttk
import tkColorChooser
import tkCommonDialog
import tkSimpleDialog
import tkFont
import Tkdnd
import ScrolledText
import sys
import commands
import math
import csv
from functools import partial
import re
import matplotlib
# ----------------------------------------------------------------------
from utils.extractor import SAextractor
from utils.analysis import SAanalysis
from utils.extractor3D import SADatacube
from utils.flux_calibration import flux_calibration
# ----------------------------------------------------------------------


#-----------------------------------------------------------------------		
# CLASS DECLARATION

class GUI_SAextractor(Frame):
	''' Class for 2D spectro-astrometry (applied to a 2D spectra):
	An instance of this class opens a GUI to help the user to 
	perform a spectro-astrometry study with differents parameters.
	Spectro-astrometry : for each wavelength in the array the centroid 
	of the spatial profile is extracted by a chosen method.
	Methods : Gaussian Kernel+FHWM / Statistical Mean+Std
	Customization : user can change
	-	the spectral window (if he does not 
		want to extract the full spectra, which is often useless) 
		without	trimming ahead the file he wanted to work with.
	-	the pixel (aperture in IRAF) window
	-	the displaying output (velocity or wavelength by redshift)
	'''
	
	
	def __init__(self, parent):
		
		# __init__
		Frame.__init__(self, parent)
		self.parent = parent  
		# Initialzsation of the FITS viewer (SaO ds9 by default)
		self.fitsviewer = "/home/postgrad/bin/ds9"	
		# Initialization of the FITS path (None by default)
		self.fits_path = ''
		abouttext = 'Released in Python 2.7.11\nAuthor : Bastien Rouzé\n Contact : rouzebastien(at)gmail.com'
		self.firstuse = True
		self.isrunning = False
		self.noisefirstsetup = True
		self.issetnoise = False
		self.set3Dplot = StringVar()
		self.q = 0 #Noise sample counter
		self.spatial_offset_choice = StringVar()
		# Initialization of the user interface   
		self.initUI()
		
	def initUI(self):
		# Window title
		self.parent.title("SA Extractor")
		self.pack()
		# Definition of the top menubar
		self.menubar = Menu(self.parent)
		self.parent.config(menu=self.menubar)
		# Creation of commands under the 'File' menu
		self.fileMenu = Menu(self.menubar)
		self.fileMenu.add_command(label="Open", command=self.onOpen)
		# Creation of commands under the 'View' menu
		self.displayMenu = Menu(self.menubar)
		self.displayMenu.add_command(label="Display", command=self.displayFits)
		self.displayMenu.add_command(label="Change FITS viewer", command=self.changeViewer)
		# Creation of commands under the 'Addon' menu
		self.addonMenu = Menu(self.menubar)
		self.addonMenu.add_command(label="Offset / Tilt", command=self.chooseTilt)
		self.addonMenu.add_command(label="Spatial profile", command=self.slice_selector)
		self.addonMenu.add_command(label="3D plot", command=self.make3D)
		# Creation of commands under the 'Addon Datacube' menu
		# Creation of commands under the 'Help' menu
		self.helpMenu = Menu(self.menubar)
		self.helpMenu.add_command(label="Help", command=self.Help)
		self.helpMenu.add_command(label="About...", command=self.About)
		# Definition of menulabels
		self.menubar.add_cascade(label="File", menu=self.fileMenu)
		self.menubar.add_cascade(label="View", menu=self.displayMenu)
		self.menubar.add_cascade(label="2D Add-on", menu=self.addonMenu)
		self.menubar.add_cascade(label="Help", menu=self.helpMenu)
		# Presentation Frame - Optional
		self.l1 = LabelFrame(self, text="Spectro-Astrometric", padx=20, pady=10)
		self.l1.pack(side=TOP) 
		Label(self.l1, text="Spectro-astrometry Signal Extractor Task").pack()
		
		# Definition of the Displaying choice --------------------------	
		self.radioframe = LabelFrame(self, text="Displaying", padx=20, pady=10)
		self.radioframe.pack(fill=X)
		
		self.Wave_or_Velo_button = StringVar()
		Radiobutton(self.radioframe,text="Wavelength",value = "W",variable = self.Wave_or_Velo_button).pack(side=LEFT)
		Radiobutton(self.radioframe,text="Velocity",value = "V",variable = self.Wave_or_Velo_button).pack(side=LEFT)
		
		self.radioframe2 = LabelFrame(self, text="Dispersion Axis", padx=20, pady=10)
		self.radioframe2.pack(fill=X)
		self.Hori_or_Verti_button = StringVar()				
		Radiobutton(self.radioframe2,text="Horizontal",value = "Hori",variable = self.Hori_or_Verti_button).pack(side=LEFT)
		Radiobutton(self.radioframe2,text="Vertical",value = "Vert",variable = self.Hori_or_Verti_button).pack(side=LEFT)
		# --------------------------------------------------------------
		
		# Set the pixel and lambda ROIs --------------------------------
		# first is for spectral direction 
		# ...2 is for spatial one.
		
		self.roiframe = LabelFrame(self, text="Set Regions of Interest", padx=20, pady=10)
		self.roiframe.pack(fill=X)
		self.lowroi = Frame(self.roiframe)
		self.lowroi.pack(fill=X)
		self.lowlabel = Label(self.lowroi, text="Wavelength Start :  ")
		self.lowlabel.pack(side=LEFT)
		self.entry_low_roiW = Entry(self.lowroi)
		self.entry_low_roiW.pack(side=LEFT)
				
		self.uproi = Frame(self.roiframe)
		self.uproi.pack(fill=X)
		self.uplabel = Label(self.uproi, text="Wavelength Stop :  ")
		self.uplabel.pack(side=LEFT)
		self.entry_up_roiW = Entry(self.uproi)
		self.entry_up_roiW.pack(side=RIGHT)
		
		Label(self.roiframe, text='\n').pack()
		
		self.lowroi2 = Frame(self.roiframe)
		self.lowroi2.pack(fill=X)
		self.lowlabel2 = Label(self.lowroi2, text='Spatial Start :')
		self.lowlabel2.pack(side=LEFT)
		self.entry_low_roiP = Entry(self.lowroi2)
		self.entry_low_roiP.pack(side=RIGHT)
		
		self.uproi2 = Frame(self.roiframe)
		self.uproi2.pack(fill=X)
		self.uplabel2 = Label(self.uproi2, text="Spatial Stop  : ")
		self.uplabel2.pack(side=LEFT)
		self.entry_up_roiP = Entry(self.uproi2)
		self.entry_up_roiP.pack(side=RIGHT)
		
		Label(self.roiframe, text='\n').pack()
		# --------------------------------------------------------------
		self.noiseframe = LabelFrame(self,text='Noise samples',padx=20, pady=10)
		self.noiseframe.pack(fill=X)
		self.spin = Spinbox(self.noiseframe, from_=0, to=10)
		self.spin.pack()
		# --------------------------------------------------------------
		Button(self, text='Display', command=self.displayFits).pack()
		self.frame1 = LabelFrame(self,text="FITS name", height=100, width=100)
		self.frame1.pack()
		Button(self, text='Run', command=self.SAcheck).pack(side=BOTTOM)
		
		self.pack()
	
	# Out Layer---------------------------------------------------------
	def onOpen(self):
		'''# File -> Open is the path to this function
		# It will open a file selection window
		# The user must choose a file or nothing happens'''
		self.choice = tkFileDialog.askopenfilename()
		#print choice
		#print type(choice)
		if str(type(self.choice)) != '<type \'tuple\'>':
			self.fits_path = self.choice
			self.choice = self.choice.split('/')
			if self.firstuse == False:
				self.pathlabel.pack_forget()
			self.pathlabel = Label(self.frame1, text=self.choice[-1])
			self.pathlabel.pack(padx=10, pady=10)
			self.firstuse = False
			
	# GUI Layer---------------------------------------------------------	
	def changeViewer(self):
		'''# Allow the user to change to another fits displayer system
		# Default is mine from my own laptop
		# Future version will save the fits viewer into a log or txt'''
		def retrieve_input():
			# 'Get' button action from changeviewer to save the new fits viewer path into 
			# the class instance.
			self.fitsviewer = self.e.get()
			print self.fitsviewer
			self.change.destroy()
		self.change = Tk()
		self.change.geometry("220x50")
		self.change.title("Change the fits viewer")
		self.e = Entry(self.change)
		self.e.pack()
		b = Button(self.change, text="Get", width=10, command=retrieve_input).pack()
		
	def displayFits(self):	
		'''# Start the displayer (for example it can be SaO ds9)
		# With a fits file if loaded'''
		print 'Opening ', self.fits_path, ' with ',self.fitsviewer	
		print self.fitsviewer,self.fits_path
		subprocess.Popen([self.fitsviewer+' '+self.fits_path], shell=True)

	def About(self):
		'''# Display an about window. Nothing special here'''
		self.about = Tk()
		self.about.title("About")
		abouttext = 'Spectro Astrometric GUI\nSignal Extractor\nReleased in Python 2.7.11\nAuthor : Bastien Rouzé\n Contact : rouzebastien(at)gmail.com'
		Label(self.about, text=abouttext).pack()

	def Help(self):
		'''# Display a little help for the task
		# Future version will include an help.txt + this function.'''
		self.thehelp = Tk()
		self.thehelp.title("Help")
		helptext = 'Help for SA Extractor\nSA Extractor is a code that extract for you the spectro-astrometric signal of a 2D spectrum\n\n(1) Choose your final data recording Wavelength/Velocity. No choice will put Wavelength by default.\n(2) Click on File->Load or on Run. Choose your fits file.\n(3)'
		Label(self.thehelp, text=helptext).pack()
	
	def chooseTilt(self):
		''' #The spectrum may be tilted on spatial axis (especially if you do not use the UVES pipelin..
		This function just allow you to choose to compute offset from noise sample, whole signal or just an integer'''
		self.addoffsetTilt = Tk()
		self.addoffsetTilt.title("Addon - Offset / Tilt")
		self.addoffsetTilt.geometry("300x100")
		Label(self.addoffsetTilt, text='/!\ Default is Signal').pack(side=TOP)
		self.spatial_offset_choice = StringVar()
		Radiobutton(self.addoffsetTilt,text="Compute from noises samples",value = "Noise",variable = self.spatial_offset_choice).pack()
		Radiobutton(self.addoffsetTilt,text="Compute from whole signal",value = "Signal",variable = self.spatial_offset_choice).pack()
		Radiobutton(self.addoffsetTilt,text="Compute from middle of spatial axis",value = "Middle",variable = self.spatial_offset_choice).pack()
	
	def slice_selector(self):
		self.slice_tk = Tk()
		self.slice_tk.title("Addon - Spatial profile")
		#self.slice_tk.geometry("300x100")
		Label(self.slice_tk, text='Select a spectral line (in wavelength)').pack(side=TOP)
		self.entry_slice = Entry(self.slice_tk)
		self.entry_slice.pack()
		self.slice_mode = StringVar()
		Label(self.slice_tk, text='Tick if you want to extract others slices').pack()
		Radiobutton(self.slice_tk,text="One-by-one",value = "Sing",variable = self.slice_mode).pack()
		Radiobutton(self.slice_tk,text="All in one",value = "Comp",variable = self.slice_mode).pack()
		Label(self.slice_tk, text='\n').pack(fill=X)
		Button(self.slice_tk, text='Display the full FITS spectrum', command=self.displayFits).pack()
		Button(self.slice_tk, text='Run', command=self.SAslice_check).pack(side=BOTTOM)
	
	def make3D(self):
		self.make3d_tk = Tk()
		self.make3d_tk.title("Addon - 3D plot options")
		Label(self.make3d_tk, text='/!\ Default is no').pack(side=TOP)
		self.set3Dplot = StringVar()
		Radiobutton(self.make3d_tk,text="Yes",value = "Y",variable = self.set3Dplot).pack()
		Radiobutton(self.make3d_tk,text="No",value = "N",variable = self.set3Dplot).pack()
		
		
		
	# Kernel Layer------------------------------------------------------
	
	def SAcheck(self):
		# Will check if all is ok to run the SA extractor
		if self.fits_path=='':
			self.onOpen()
		if self.Wave_or_Velo_button.get()=='':
			showwarning('Warning','No choice between Wavelength and Velocity.\nWavelength is set by default.')
			self.Wave_or_Velo = 'W'
			self.Wave_or_Velo_button = StringVar()
		else:
			self.Wave_or_Velo = self.Wave_or_Velo_button.get()
		if self.Hori_or_Verti_button.get()=='':
			showwarning('Warning','No choice between Horizontal or Vertical spectrum.\nHorizontal is set by default.')
			self.Hori_or_Verti = 'Hori'
			self.Hori_or_Verti_button = StringVar()
		else:
			self.Hori_or_Verti=self.Hori_or_Verti_button.get()
		# Future version will include other parameters.	
		# Last parameter to check is the ROI
		if self.set3Dplot.get() != '':
			self.set3Dplotv = self.set3Dplot.get()
		else :
			self.set3Dplotv = 'N'
		
		if (self.entry_low_roiW.get()=='' or self.entry_up_roiW.get()==''):
			if askyesno('Information','ROI is empty or half filled. Continue?'):
				self.isrunning = True 
		else:
			self.isrunning = True 
		
		if self.isrunning	:
			print 'Loading dictionnary...'
			self.data_dictionnary ={}
			self.data_dictionnary['fits_path']=self.fits_path
			self.data_dictionnary['up_roiP']=self.entry_up_roiP.get()
			self.data_dictionnary['low_roiP']=self.entry_low_roiP.get()
			self.data_dictionnary['up_roiW']=self.entry_up_roiW.get()
			self.data_dictionnary['low_roiW']=self.entry_low_roiW.get()
			self.data_dictionnary['Wave_or_Velo']=self.Wave_or_Velo
			self.data_dictionnary['Hori_or_Verti']=self.Hori_or_Verti
			self.data_dictionnary['issetnoise']=self.issetnoise
			self.data_dictionnary['spin_selector_noise']=self.spin.get()
			self.data_dictionnary['spatial_offset_choice']=self.spatial_offset_choice.get()
			self.data_dictionnary['set3Dplot']=self.set3Dplotv
			print 'done'
			SAextractor(self.data_dictionnary)
			
			
			
	def SAslice_check(self):
		if self.fits_path=='':
			self.onOpen()
		if self.entry_slice.get()!='':
			print 'Loading dictionnary...'
			self.data_dictionnary ={}
			self.data_dictionnary['fits_path']=self.fits_path
			self.data_dictionnary['slice_wavelength']=self.entry_slice.get()
			self.data_dictionnary['slice_mode']=self.slice_mode
			self.data_dictionnary['Hori_or_Verti']='Hori'
			self.data_dictionnary['SAextrator_running']=self.isrunning
			print 'done'
			from utils.extractor import SAextractor_slice
			SAextractor_slice(self.data_dictionnary)
		else :
			showerror('Error','You must choose a wavelength !')
	
		
