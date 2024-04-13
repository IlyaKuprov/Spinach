# python script to directly modify definitions defined
# in the current experiment - designed primarily to
# be used with py_run.m
# 
# The input should be pairwise strings, 
# (var1_name var1_value var2_name var2_value...),
# being the variable name (as named in the .def file)
# and the variable value - the name and value should be
# input as strings, separated by spaces.
#
# NOTE:	The correct location of the XeprAPI module should
# 	be inserted into the code below
#
# david.goodwin@kit.edu

import os,sys,numpy; 
from numpy import genfromtxt

# xepr import
try:
	sys.path.insert(0, os.popen("Xepr --apipath").read());  # this locates the XeprAPI module
	sys.path.append('/opt/Bruker/xepr/sharedProDeL/Examples/XeprAPI-examples');
	sys.path.append('/opt/Bruker/xepr/sharedProDeL/Standard/XeprAPI');

	import XeprAPI      # load the Xepr API module
	Xepr=XeprAPI.Xepr()
except:
	print "xepr_import_error"
	sys.exit(1)

try:
	currentExp = Xepr.XeprExperiment()
	hiddenExp = Xepr.XeprExperiment("AcqHidden")
except:
	print("No Experiment in Primary: "+"No Experiment has been selected in the Primary Viewport of Xepr.")
	sys.exit(2)

# get the start value for the variable
# search for our desired variable
# first get the text of the full PulseSpel def
fullDefs = currentExp.getParam("PlsSPELGlbTxt").value
# need to check if fullDefs is empty and exit cause pulsespel not being loaded
fullDefs = fullDefs.split("\n")

no_defs=(len(sys.argv)-1)/2

try:
	for value in fullDefs:
		for index in range(1,no_defs+1):
			if str(sys.argv[(2*index)-1]) in value:
				cmdStr = str(sys.argv[(2*index)-1])+" = "+str(sys.argv[2*index])
				currentExp["ftEPR.PlsSPELSetVar"].value = cmdStr
except:
	print "error changing pulseSPEL defs"
	sys.exit(4)

