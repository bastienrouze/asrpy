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

class GUI_SAanalysis(Frame):
	
	
	def __init__(self, parent):
		Frame.__init__(self, parent)
		self.parent = parent        
		self.initUI()
		self.path1 = ''
		self.path2 = ''
		self.firstuse = True
		self.fitsviewer = '/home/postgrad/bin/ds9'
		abouttext = 'Released in Python 2.7.11\nAuthor : Bastien Rouzé\n Contact : rouzebastien(at)gmail.com'
		
	def initUI(self):
		self.parent.title("SA Analysis")
		self.pack()
		self.menubar = Menu(self.parent)
		self.parent.config(menu=self.menubar)
		self.fileMenu = Menu(self.menubar)
		self.fileMenu.add_command(label="Open", command=self.onOpen)
		self.displayMenu = Menu(self.menubar)
		self.displayMenu.add_command(label="Run FITS viewer", command=self.displayFits)
		self.displayMenu.add_command(label="Change FITS viewer", command=self.changeviewer)
		self.helpMenu = Menu(self.menubar)
		self.helpMenu.add_command(label="Help", command=self.Help)
		self.helpMenu.add_command(label="About...", command=self.About)
		
		self.menubar.add_cascade(label="File", menu=self.fileMenu)
		self.menubar.add_cascade(label="View", menu=self.displayMenu)
		self.menubar.add_cascade(label="Help", menu=self.helpMenu)
		
		self.l1 = LabelFrame(self, text="Spectro-Astrometric", padx=20, pady=10)
		self.l1.pack(side=TOP) 
		Label(self.l1, text="Spectro-astrometry Signal Analysis Task").pack()
		
		
		self.radioframe = LabelFrame(self, text="Displaying", padx=20, pady=10)
		self.radioframe.pack(fill=X)
		self.Wbtn = StringVar()
				
		Radiobutton(self.radioframe,text="Wavelength",value = "W",variable = self.Wbtn).pack(side=LEFT)
		Radiobutton(self.radioframe,text="Velocity",value = "V",variable = self.Wbtn).pack(side=LEFT)
		
		self.radioframe2 = LabelFrame(self, text="Mode", padx=20, pady=10)
		self.radioframe2.pack(fill=X)
		self.Dbtn = StringVar()
				
		Radiobutton(self.radioframe2,text="Single",value = "Sing",variable = self.Dbtn).pack(side=LEFT)
		Radiobutton(self.radioframe2,text="Comparison",value = "Comp",variable = self.Dbtn).pack(side=LEFT)
		
		Button(self, text='Run', command=self.SAcheck).pack(side=BOTTOM)
		self.frame1 = LabelFrame(self,text="TXT path", height=100, width=100)
		self.frame1.pack(side=BOTTOM)
		self.pack()


	def SAcheck(self):
		print 'Making dictionnary...'
		self.dictionnary = {}
		if self.path1 == '':
			self.path1 = tkFileDialog.askopenfilename()
		if self.Dbtn.get()=='Comp':
			self.path2 = tkFileDialog.askopenfilename()
			self.forComparison = self.Dbtn.get()
			self.dictionnary['path2'] = self.path2
		else :
			self.forComparison = 'Sing'
		
		if self.Wbtn.get()=='':
			print '/!\ No choice between Wavelength and Velocity. Wavelength is set by default.'
			showwarning('Warning','No choice between Wavelength and Velocity.\n Wavelength is set by default.')
			self.WaveorVelo = 'W'
		elif self.Wbtn.get()=='V':
			self.WaveorVelo = 'V'
		else :
			self.WaveorVelo = 'W'
		
		self.dictionnary['path1'] = self.path1
		self.dictionnary['Sing_or_Comp']= self.forComparison 
		self.dictionnary['Wave_or_velo'] = self.WaveorVelo
		print 'done'
		SAanalysis(self.dictionnary)
		# Reset the path
		self.path1 = ''
		self.path2 = ''		
					

	def retrieve_input(self):
		self.fitsviewer = self.e.get()
		print self.fitsviewer
		self.change.destroy()
		
	def changeviewer(self):
		self.change = Tk()
		self.change.geometry("220x50")
		self.change.title("Change the fits viewer")
		self.e = Entry(self.change)
		self.e.pack()
		b = Button(self.change, text="Get", width=10, command=self.retrieve_input).pack()
		
		
	def displayFits(self):	
		#commands.getoutput('ds9'+' '+self.fits_path)
		self.fits_path = tkFileDialog.askopenfilename()
		if str(type(self.fits_path)) != '<type \'tuple\'>':
			print 'Opening ', self.fits_path, ' with ',self.fitsviewer	
			print self.fitsviewer,self.fits_path
			subprocess.Popen([self.fitsviewer+' '+self.fits_path], shell=True)

	def About(self):
		self.about = Tk()
		self.about.title("About")
		Label(self.about, text=abouttext).pack()

	def Help(self):
		self.thehelp = Tk()
		self.thehelp.title("Help")
		helptext = 'Help for SA Analysis\n'
		Label(self.thehelp, text=helptext).pack()
		
	def onOpen(self):
		self.choice = tkFileDialog.askopenfilename()
		if choice != '':			
			self.self.fits_path1 = choice
			Label(self.frame1, text=self.fits_path1).pack(padx=10, pady=10)
		if self.Dbtn.get()=='Comp':
			self.choice = tkFileDialog.askopenfilename()
			if self.choice != '':			
				self.fits_path2 = self.choice
				Label(self.frame1, text=self.fits_path2).pack(padx=10, pady=10)	
