# python script to load and compile a program defined
# by the single input to this function - designed
# primarily to be used with py_run.m
# 
# The single input should be a string containing the name
# of the program file (full path with .exp extension)
#
# NOTE:	The correct location of the XeprAPI module should
# 	be inserted into the code below
#
# david.goodwin@kit.edu

import os,subprocess,sys;

# xepr import
try:
    sys.path.insert(0, os.popen("Xepr --apipath").read());  # this locates the XeprAPI module
    sys.path.append('/opt/Bruker/xepr/sharedProDeL/Examples/XeprAPI-examples');
    sys.path.append('/opt/Bruker/xepr/sharedProDeL/Standard/XeprAPI');

    import XeprAPI      # load the Xepr API module
    Xepr=XeprAPI.Xepr()
except:
    print "XeprAPI_import_error"
    sys.exit(1)

# python inputs
try:
    program_file=(str(sys.argv[1]))
except:
    print "python_input_error"
    sys.exit(2)

# Attempt to change the shape file
try:
	Xepr.XeprCmds.aqPgLoad(program_file)
	Xepr.XeprCmds.aqPgShowPrg()
	Xepr.XeprCmds.aqPgCompValid()
	Xepr.XeprCmds.aqPgCompile()
except:
	print "program_load_compile_error"
	sys.exit(3)
