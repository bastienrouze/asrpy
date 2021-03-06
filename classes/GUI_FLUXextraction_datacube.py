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

class GUI_FLUXextraction_datacube(Frame):
	
	
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
		self.parent.title("SA Datacube")
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

		# Creation of commands under the 'Addon Datacube' menu
		# Creation of commands under the 'Help' menu
		self.helpMenu = Menu(self.menubar)
		self.helpMenu.add_command(label="Help", command=self.Help)
		self.helpMenu.add_command(label="About...", command=self.About)
		# Definition of menulabels
		self.menubar.add_cascade(label="File", menu=self.fileMenu)
		self.menubar.add_cascade(label="View", menu=self.displayMenu)
		self.menubar.add_cascade(label="Options", menu=self.addonMenu)
		self.menubar.add_cascade(label="Help", menu=self.helpMenu)
		# Presentation Frame - Optional
		self.l1 = LabelFrame(self, text="Spectro-Astrometric", padx=20, pady=10)
		self.l1.pack(side=TOP) 
		Label(self.l1, text="Spectro-astrometry Datacube Tasks").pack()
		
		# Definition of the Displaying choice --------------------------	
		self.displayframe = LabelFrame(self, text="Displaying [Optionnal]", padx=20, pady=10)
		self.displayframe.pack(fill=X)
		
		self.radioframe = Frame(self.displayframe)
		self.radioframe.pack(fill=X)
								
		self.namelabel = Label(self.displayframe, text="Caption : ")
		self.namelabel.pack(side=LEFT)
		self.entry_name = Entry(self.displayframe)
		self.entry_name.pack(side=RIGHT)
		
		# Definition of the Job choice --------------------------	
		
				
		# --------------------------------------------------------------
		
		# Set the pixel and lambda ROIs --------------------------------
		# first is for spectral direction 
		# ...2 is for spatial one.
		
		self.roiframe = LabelFrame(self, text="Set Ranges of Interest [in pixels]", padx=20, pady=10)
		self.roiframe.pack(fill=X)
		
		self.winterest = Frame(self.roiframe)
		self.winterest.pack(fill=X)
		self.winlabel = Label(self.winterest, text="Wavelength of interest slice:  ")
		self.winlabel.pack(side=LEFT)
		self.entry_W = Entry(self.winterest)
		self.entry_W.pack(side=RIGHT)
		
		self.wroi = Frame(self.roiframe)
		self.wroi.pack(fill=X)
		self.wlabel = Label(self.wroi, text="Wavelength range length :  ")
		self.wlabel.pack(side=LEFT)
		self.entry_roiW = Entry(self.wroi)
		self.entry_roiW.pack(side=RIGHT)
				
		self.decraroi = Frame(self.roiframe)
		self.decraroi.pack(fill=X)
		self.decralabel = Label(self.decraroi, text="DEC / RA range length : ")
		self.decralabel.pack(side=LEFT)
		self.entry_decra_roi = Entry(self.decraroi)
		self.entry_decra_roi.pack(side=RIGHT)

		# --------------------------------------------------------------
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
	

		
		
	# Kernel Layer------------------------------------------------------
	
	def SAcheck(self):
		# Will check if all is ok to run the SA extractor
		if self.fits_path=='':
			self.onOpen()
		if self.entry_W.get()=='':
			showerror('Error','Wavelength of interest is empty')
			self.isrunning = False
		else :
			self.isrunning = True
		
		if self.isrunning	:
			print 'Loading dictionnary...'
			self.data_dictionnary ={}
			self.data_dictionnary['fits_path']=self.fits_path
			self.data_dictionnary['central_wavelength']=int(self.entry_W.get())
			self.data_dictionnary['range_wavelength']=int(self.entry_roiW.get())
			self.data_dictionnary['range_decra']=int(self.entry_decra_roi.get())
			self.data_dictionnary['name']=self.entry_name.get()
			print 'done'
			from utils.flux3D import FLUXextraction_3D
			FLUXextraction_3D(self.data_dictionnary)
