import h5py, seaborn, matplotlib.pyplot as plt, numpy as np
from matplotlib import colors
import numpy.ma as ma, pandas as pds
import gdal, gdalconst, csv, sys, h5py

def read_ascii_to_array(ascii_path, xmin, ymax, cols_wide, rows_high):
    """ Reads an ascii grid file into a numpy array and returns the array. This
        function first checks that the number of rows and columns and the upper
        left corner co-ordinates are as expected.
    """
    # Register drivers
    gdal.AllRegister()

    # Open the file with GDAL
    ds = gdal.Open(ascii_path, gdalconst.GA_ReadOnly)
    if ds is None:
        print('Could not open ' + ascii_path)
        sys.exit(1)

    # Get dataset properties
    geotransform = ds.GetGeoTransform()
    origin_x = geotransform[0]
    origin_y = geotransform[3]

    # Read the data to an array and return it
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray()

    # Get the no data value
    ndv = band.GetNoDataValue()

    # Check properties are correct
    assert origin_x == xmin
    assert origin_y == ymax
    assert data.shape == (rows_high, cols_wide)

    # Mask the array and return it
    masked = ma.masked_equal(data, ndv)
    return masked

def read_from_h5(hdf5_file, dset_path):
    """ Reads a masked array from the specified location in an HDF5 file.
    """
    # Open the HDF5 file
    h5file = h5py.File(hdf5_file, 'r')

    # Get the array
    dset = h5file[dset_path]
    array = dset[:].astype(float)

    # Read the array and mask no data
    array = ma.masked_equal(array, float(0))
    #array = ma.masked_equal(array, float(dset.attrs['No_Data_Value']))

    # Close the HDF5 file
    h5file.close()

    return array

###############################################################################
# User input
#h5_path = r'NIRAMS_Output_2020_ParamCombo/NIRAMS_2020_MultiPar.h5'
ascii_gwbodygroups = r'../code/data/NewData2017/gwbodygroups.txt'

st_yr = 2001  # Min = 2001; max = 2018 
end_yr = 2018 # Min = 2001; max = 2018
xmn = 0
nxs = 485
ymn = 0
nys = 715
runmin = 1      # Min = 1; max = 81
runmax = 81
#runmax = 1
grpmin = 1
#grpmin = 6
#grpmax = 6
grpmax = 60
#grpmax = 1
runavlen=6
###############################################################################

# Read the GWB groups into an array
print('Reading GW body group data...')
gwb_groups = read_ascii_to_array(ascii_gwbodygroups, 0, 1235000, 485, 715)

for targ_gwbgrp in range(grpmin,grpmax+1):
	filename = 'IndividualGWBGroupData_2020/GWBGroupData_%i' % (targ_gwbgrp)
	avfilename = 'IndividualGWBGroupData_2020/GWBGroupData_YearAvs_%i' % (targ_gwbgrp)
	celllistfile = 'IndividualGWBGroupData_2020/GWBGroupCellList_%i' % (targ_gwbgrp)
	target = open(filename, "w")
	avtarget = open(avfilename, "w")
	gwgrp = ma.masked_not_equal(gwb_groups, targ_gwbgrp)
	target.write(r'%i %i' % (gwgrp.count(),end_yr-st_yr+1))
	target.write("\n")
	avtarget.write(r'%i %i' % (gwgrp.count(),end_yr-st_yr+1))
	avtarget.write("\n")
	nonzlist=gwgrp.nonzero()
	celllisttarg=open(celllistfile, "w")
	celllisttarg.write(r'%i' % (gwgrp.count()))
	celllisttarg.write("\n")
	for i in range(gwgrp.count()):
		celllisttarg.write(r'%i %i' % (nonzlist[0][i],nonzlist[1][i]))
		celllisttarg.write("\n")
	celllisttarg.close()

	yeardata=np.zeros((gwgrp.count(),end_yr-st_yr+1))
	nyrdata=np.zeros((gwgrp.count(),end_yr-st_yr+1))
	yearav=np.zeros(gwgrp.count())
	nyrav=np.zeros(gwgrp.count())
	for run in range(runmin, runmax+1):
	    h5_path = r'NIRAMS_Output_2020_ParamCombo/temp/run_%03d.h5' % (run)
	    #for year in range(hms_df.index[0], hms_df.index[-1]+1):
	    #out_row = np.array(['Run{rval}'.format(rval=run)])
	    yeardata.fill(0)
	    nyrdata.fill(0)
	    yearav.fill(0)
	    nyrav.fill(0)
	    for year in range(2001,2019):
	        #np.append(out_row,['Year{yr}'.format(yr=year)])
	        print('    Currently processing: GWB group %i, run %03d, %s' % (targ_gwbgrp, run, year))
	
	        # Read the arrays for this year
	        drain = read_from_h5(h5_path,
	                             '/run_%03d/yearly_total_drainage_%s' % (run,year))
	        no3 = read_from_h5(h5_path,
	                           '/run_%03d/yearly_n_leached_%s' % (run,year))

	        for i in range(gwgrp.count()):
	                if ((drain[nonzlist[0][i]][nonzlist[1][i]] != -9999) and (no3[nonzlist[0][i]][nonzlist[1][i]] != -9999)):
	                        mod_n_conc = 100.*no3[nonzlist[0][i]][nonzlist[1][i]]/drain[nonzlist[0][i]][nonzlist[1][i]]
	                        target.write(r'%10.8e ' % (mod_n_conc))
	                        yeardata[i,year-st_yr]=mod_n_conc
	                        nyrdata[i,year-st_yr]=1 
	                else:
	                        target.write("-9999 ")
	                        nyrdata[i,year-st_yr]=0

	        # Mask drainage values not in current catchment
	        #drain[gwb_groups!=targ_gwbgrp] = ma.masked
	        #drain[no3==-9999] = ma.masked
	        # Mask N values not in current catchment
	        #no3[gwb_groups!=targ_gwbgrp] = ma.masked
	        #no3[no3==-9999] = ma.masked
		#mod_n_conc = 100.*drain/no3
		#for val in np.array(mod_n_conc[mod_n_conc.nonzero()]):
		#	target.write(r'%10.8e ' % (val))
		#np.append(out_row,np.array(mod_n_conc[mod_n_conc.nonzero()]))
		#target.write(out_row.tolist())
	        target.write("\n")
	        #drain.mask = ma.nomask
	        #no3.mask = ma.nomask
	    for j in range(runavlen):
	    	for i in range(gwgrp.count()):
	    		if (nyrdata[i,j]):
	    			yearav[i]+=yeardata[i,j]
	    			nyrav[i]+=1
	    for i in range(gwgrp.count()):
	    	avtarget.write(r'%10.8e ' % (yearav[i]/nyrav[i]))
	    avtarget.write("\n")
	    for j in range(runavlen,end_yr-st_yr+1):
	    	for i in range(gwgrp.count()):
	    		if (nyrdata[i,j-runavlen]):
	    			yearav[i]-=yeardata[i,j-runavlen]
	    			nyrav[i]-=1
	    		if (nyrdata[i,j]):
	    			yearav[i]+=yeardata[i,j]
	    			nyrav[i]+=1
	    	for i in range(gwgrp.count()):
	    		avtarget.write(r'%10.8e ' % (yearav[i]/nyrav[i]))
	    	avtarget.write("\n")


	target.close()
	avtarget.close()
	gwgrp.mask = ma.nomask
