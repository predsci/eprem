#!/usr/bin/env python3

import argparse
import netCDF4 as nc
import numpy as np
import bottleneck as bn

#############################################################################
## argument parsing##########################################################
#############################################################################
def argParsing():

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('file_ref',
                        help='Reference netCDF EPREM output file')

    parser.add_argument('file_new',
                        help='New netCDF EPREM output file')

    parser.add_argument('-v',
                        help='Output the location and values of max errors',
                        dest='verbose',
                        action='store_true',
                        default=False,
                        required=False)

    parser.add_argument('--version', action='version', version='%(prog)s 1.2')

    return parser.parse_args()


#############################################################################
## read netCDF file in ######################################################
#############################################################################
def openNetCDF(ncFile):

    rootGrp = nc.Dataset(ncFile, 'r')

    class dataStruct:
        def __init__(self):
            ## variables
            self.Dist   = np.copy(rootGrp.variables['Dist'])

    dStruct = dataStruct()
    rootGrp.close()

    return dStruct

#############################################################################
## output builder ###########################################################
#############################################################################
def outputBuilder(args, data_0, data_1):

    ## error checking on the file dimension
    if (data_0.Dist.shape != data_1.Dist.shape):
        print('The dimensions of the files are not the same.\n')
        quit()

    np.seterr(divide='ignore', invalid='ignore')

    d1md0     = np.subtract(data_1.Dist,data_0.Dist)
    absrelerr = np.absolute(np.divide(d1md0,data_0.Dist))

    maxabs_index  = bn.nanargmax(np.absolute(d1md0))
    maxabs_index  = np.unravel_index(maxabs_index, d1md0.shape)
    maxrabs_index = bn.nanargmax(absrelerr)
    maxrabs_index = np.unravel_index(maxrabs_index, d1md0.shape)
    sumsqsdiff    = bn.ss(d1md0)       #np.sum(np.square(d1md0))
    sumsqs0       = bn.ss(data_0.Dist) #np.sum(np.square(data_0.Dist))
    N             = np.prod(data_0.Dist.shape)
    #Nnon0         = np.count_nonzero(data_0.Dist)

    maxabs     = abs(d1md0[maxabs_index])
    maxrabs    = absrelerr[maxrabs_index]
    meanrabs   = bn.nanmean(absrelerr)
    rmsd       = np.sqrt(sumsqsdiff/N)
    cvrmsd     = np.sqrt(sumsqsdiff/sumsqs0)


    if (args.verbose):
        print(' ')
        print('--------------------------------------------------------------------------------------------------------------------')
        print(' a: '+str(args.file_ref))
        print(' b: '+str(args.file_new))
        print('--------------------------------------------------------------------------------------------------------------------')
        print('%14s    %14s    %14s    %20s    %14s    %20s' % (' ','sqrt(ssd/ssa)','sqrt(ssd/N)','max|(b-a)/a|','mean|(b-a)/a|','max|b-a|'))
        print('--------------------------------------------------------------------------------------------------------------------')
        print('%14s    %14.8e    %14.8e    %20.8e    %14.8e    %20.8e' % ('Error value',cvrmsd, rmsd, maxrabs, meanrabs, maxabs))
        print('%14s    %14s    %14s    %20s    %14s    %20s' % ('Array Index',' ',' ',str(maxrabs_index),' ',str(maxabs_index)))
        print('%14s    %14s    %14s    %20.8e    %14s    %20.8e' % ('Data a value',' ',' ',data_0.Dist[maxrabs_index],' ', data_0.Dist[maxabs_index]))
        print('%14s    %14s    %14s    %20.8e    %14s    %20.8e' % ('Data b value',' ',' ',data_1.Dist[maxrabs_index],' ', data_1.Dist[maxabs_index]))
    else:
        print(' ')
        print('--------------------------------------------------------------------------------------')
        print(' a: '+str(args.file_ref))
        print(' b: '+str(args.file_new))
        print('--------------------------------------------------------------------------------------')
        print('%14s    %14s    %14s    %14s    %14s' % ('sqrt(ssd/ssa)','sqrt(ssd/N)','max|(b-a)/a|','mean|(b-a)/a|','max|b-a|'))
        print('--------------------------------------------------------------------------------------')
        print('%14.8e    %14.8e    %14.8e    %14.8e    %14.8e' % (cvrmsd, rmsd, maxrabs, meanrabs, maxabs))


    # Figure out is the comparison has "passed" or not.
    passfail = 0
    fptol = 1.0e-8

# Binary
    if (maxabs == 0):
        passfail=-1
#Floating point:
    if (cvrmsd > fptol):
        passfail=max(passfail,1)
    if (rmsd > fptol):
        passfail = max(passfail,1)
    if (maxabs > fptol):
        passfail = max(passfail,1)
#For relative diff, check for inf
    if (np.isfinite(maxrabs)):
        if (maxrabs > fptol):
            passfail = max(passfail,1)
    if (np.isfinite(meanrabs)):
        if (meanrabs > fptol):
            passfail = max(passfail,1)

    if (passfail == 0):
        print('PASS? (FP TOL 1e-11)')
    if (passfail == -1):
        print('PASS (BINARY)')
    if (passfail == 1):
        print('FAIL')

    return 0

#############################################################################
## main #####################################################################
#############################################################################

#Get the command line arguments:
args = argParsing()

#Read in the netCDF data:
data_0 = openNetCDF(args.file_ref)
data_1 = openNetCDF(args.file_new)

#Perform the comparions and output results:
outputBuilder(args, data_0, data_1)



