# python script to get data from spectrometer - designed
# primarily to be used with py_run.m
# 
# First input should be "TM", "Signal", or "RM"
# 
# Second input should be a string containing the name of the 
# experiment to run, usually "AWGTransient" for the Signal or
# "AWG_1pulseSHP" for the TM. 
# 
# Third input should be a string containing
# the name of the save files (without extensions).
#
# NOTE:	The correct location of the XeprAPI module should
# 	be inserted into the code below
#
# david.goodwin@kit.edu

import os,sys,numpy,time; 
from numpy import genfromtxt
from pprint import pprint as p

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
	currentExp=Xepr.XeprExperiment()
	hiddenExp = Xepr.XeprExperiment("AcqHidden")
except:
	print("No Experiment in Primary: "+"No Experiment has been selected in the Primary Viewport of Xepr.")
	sys.exit(2)

# change detection mode and experiement select
try:
	if hiddenExp["ftBridge.Detection"].value != str(sys.argv[1]):
		hiddenExp["ftBridge.Detection"].value = str(sys.argv[1])
	if currentExp["ftEpr.PlsSPELEXPSlct"].value != str(sys.argv[2]):
		currentExp["ftEpr.PlsSPELEXPSlct"].value = str(sys.argv[2])
except:
	print "error changing detection mode to TM, or experiment selection"
	sys.exit(3)

# run new experiment
try:
    currentExp.aqExpRunAndWait()
except:
    print "error running current experiment"
    sys.exit(4)

#retrieve data
dset = Xepr.XeprDataset()

# no data? -- try to run current experiment (if possible)
if not dset.datasetAvailable():
    try:
        print "Trying to run current experiment to create some data..."
        Xepr.XeprExperiment().aqExpRunAndWait()
    except XeprAPI.ExperimentError:
        print "No dataset available and no (working) experiment to run...giving up..."
        sys.exit(5)
if not dset.datasetAvailable():
    print "(Still) no dataset available...giving up..."
    sys.exit(6)

# save data
try:
    ordinate=dset.O
    numpy.savetxt((str(sys.argv[3]) + "_X.txt"),dset.X)
    numpy.savetxt((str(sys.argv[3]) + "_rY.txt"),ordinate.real)
    numpy.savetxt((str(sys.argv[3]) + "_iY.txt"),ordinate.imag)
except:
    print "Problem saving dataset to file"
    sys.exit(7)

