# python script to copy the current experiment and use 
# it to replace the current experiment - designed 
# primarily to be used with py_run.m
# 
# Should be used to reset the AWG after 114 sequential
# shape "load" and "run" events.
#
# NOTE:	The correct location of the XeprAPI module should
# 	be inserted into the code below
#
# david.goodwin@kit.edu

import os,subprocess,sys,numpy; 
from pprint import pprint as p

# Xepr import
try:
    sys.path.insert(0, os.popen("Xepr --apipath").read());  # this locates the XeprAPI module
    sys.path.append('/opt/Bruker/xepr/sharedProDeL/Examples/XeprAPI-examples');
    sys.path.append('/opt/Bruker/xepr/sharedProDeL/Standard/XeprAPI');

    import XeprAPI      # load the Xepr API module
    Xepr=XeprAPI.Xepr()
except:
    print "xepr_import_error"
    sys.exit(1)

# get current experiment name
curr_exp=Xepr.XeprExperiment()
expt_name=curr_exp.aqGetExpName()

# duplicate current experiment
print("replacing experiment <" + expt_name +"> with new instance of its copy...")
Xepr.XeprCmds.aqExpCut(expt_name)

print("new instance of experiment <" + expt_name +"> created...")
Xepr.XeprCmds.aqExpPaste()

# select new experiment
Xepr.XeprCmds.aqExpSelect(expt_name)
print("new instance of experiment <" + expt_name +"> selected...")

# activate new experiment
Xepr.XeprCmds.aqExpActivate(expt_name)
print("new instance of experiment <" + expt_name +"> activated.")

# open parameter panel
Xepr.XeprCmds.aqParOpen()

