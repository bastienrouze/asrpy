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
from classes.GUI_SAextractor import GUI_SAextractor
from classes.GUI_SAanalysis import GUI_SAanalysis
from classes.GUI_SAdatacube import GUI_SAdatacube
from classes.GUI_FLUXextraction_datacube import GUI_FLUXextraction_datacube
from classes.GUI_FLUXcalibration import GUI_FLUXcalibration

# ----------------------------------------------------------------------


# 2 methods for main----------------------------------------------------
# Have to be improved / removed
def fetch():
	'''Run the selected task from main actor task selector
	by destroying et rebuilding the window'''

	choice = listbox.get(ACTIVE)
	print choice
	destroyer('main')
	#The if decides which class instance invoke.
	if choice == 'Spectroastrometry - extract/2dspec':
		root = Tk()
		GUISAE = GUI_SAextractor(root)
	elif choice == 'Spectroastrometry - analysis/2dspec':
		root = Tk()		
		GUISAA = GUI_SAanalysis(root)
	elif choice == 'Spectroastrometry - extract/3dspec':
		root = Tk()
		GUISAC = GUI_SAdatacube(root)
	elif choice == 'Flux - fluxextract/3dspec':
		root = Tk()
		GUIMISCFLUXE = GUI_FLUXextraction_datacube(root)
	elif choice == 'Flux - fluxcalib':
		root = Tk()
		GUIMISCFLUXC = GUI_FLUXcalibration(root)		
	else:
		print 'This task does not exist ! Please Retry...'

def destroyer(who):
	'''Destroy a tkinter window'''
	if who=='main':
		#root.quit()
		root.destroy()	
	if who=='extractor':
		root.quit()
		
#-----------------------------------------------------------------------

def main(args):
	# The main contains the first window which give the task list to
	# the current user. The user can select one task and click on run
	# but he can also decide to exit the tool
	# In future version main will be replaced by a class instance
	program_name = 'ASPpy - Astrospec\nAstronomical Spectroscopy Package python\nGraphical User Interface for data processing\n2017'
	#----------------
	global root
	root = Tk()
	root.wm_title("Astrospec")
	global listbox
	frame1 = Frame(root, bg='gray')
	frame1.pack(fill=X)
	Label(frame1, text=program_name).pack(padx=100,pady=10)
	#photo = PhotoImage(file="MU_logo.png")

	#canvas = Canvas(frame1,width=350, height=70)
	#canvas.create_image(240, 70/2, anchor=CENTER, image=photo)
	#canvas.pack()
	frame2 = Frame(root)
	frame2.pack(fill=X)
	Label(frame2, text='--Task selector--').pack(padx=100,pady=10)
	
	listbox = Listbox(root)
	listbox.pack(fill=X)
	Button(root, text='Run', command=fetch).pack()
	Button(root, text='Exit', command=root.quit).pack()
	listbox.insert(1, 'Spectroastrometry - extract/2dspec')
	listbox.insert(2, 'Spectroastrometry - analysis/2dspec')
	listbox.insert(3, 'Spectroastrometry - extract/3dspec')
	listbox.insert(4, 'Flux - fluxextract/3dspec')
	listbox.insert(5, 'Flux - fluxcalib')
	root.mainloop()
	#----------------
	print 'end'
	return 0


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
