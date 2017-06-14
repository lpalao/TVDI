# TVDI
# Author: Lorena Villano, Leo Kris Palao, and Jorrel Aunario
# Function that will calculate and produce the TVDI values and image, plus graph of the NDVI-LST triangle and scatter plot

"""
Author: Lorie Villano
Revised by: Leo Kris Palao
parameters in the function:
    ndviFile - filename of the ndvi image
    lstFile -filename of the LST image
    sPlot - filename for the output plot
    tvdiFile - filename of the output TVDI image """

#  import necessary libraries and packages
import gdal, gdalnumeric
from matplotlib import pyplot as plt
from numpy import polyfit
import numpy as np
import glob
import itertools, os

def array_to_raster(array,filename):
    dst_filename = filename
    proj=ds.GetProjection()
    x_pixels = ds.RasterXSize  # number of pixels in x
    y_pixels = ds.RasterYSize  # number of pixels in y
    PIXEL_SIZE = dsTrans[1]  # size of the pixel...        
    x_min = dsTrans[0] 
    y_max = dsTrans[3]  # x_min & y_max are like the "top left" corner.
    wkt_projection = proj
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(dst_filename,x_pixels,y_pixels,1,gdal.GDT_Float32, ) #GDT_Float32,gdal.GDT_Int16
    dataset.SetGeoTransform((x_min,PIXEL_SIZE,0,y_max,0,-PIXEL_SIZE))  
    dataset.SetProjection(wkt_projection)
    dataset.GetRasterBand(1).WriteArray(array)
    dataset.FlushCache()
    print os.path.basename(filename)
    return dataset

def calcTVDI(ndviFile, lstFile, sPlot, tvdiname):

    # image file to array
    arrNdvi = gdalnumeric.LoadFile(ndviFile) #.astype(np.dtype(float))

    # image file to array
    arrLst = gdalnumeric.LoadFile(lstFile)

    # mask out the temperature values <0 -erroneous values from converting temp from Kelvin to Celsius in ENVI
    maskLST = np.where((arrLst>=0), 1, 0)
    arrLSTm = np.choose(maskLST, (np.nan, arrLst))

    vNdvi = np.ravel(arrNdvi) # image file to 1D array a.k.a flatten array
    vLst = np.ravel(arrLst) # image file to 1D array a.k.a flatten array

    # Check if NDVI and LST have the same dimension
    lenNdvi = np.size(arrNdvi)
    lenLst = np.size(arrLst)
    dims = lenNdvi == lenLst

    if not dims == True:
        raise IOError("The dimension of NDVI and LST are not the same.\
        Please check number of cols and rows")

    """ Calculate the dry and wet edges """

    #Compute for slope and intercept from relationship of NDVI and LST
    # Compute for LSTmax
    mask = np.where((vNdvi >= 0.3) & (vNdvi <= 0.8), 1, 0) # set threshold for NDVI values for LSTmax calculation
    rNdvi = np.choose(mask, (np.nan, vNdvi)) # mask the selected NDVI
    Nmax = np.nanmax(rNdvi) # maximum NDVI value
    Nmin = np.nanmin(rNdvi) # minimum NDVI value

    print "NDVI threshold for LSTmax calculation : 0.3 to 0.8"

    # set the step interval for the edges
    step =0.025

    UpPercentile = 99 # the upper percentile of LST for calculating LSTmax per interval
    LoPercentile = 1 # the lower percentile of LST for calculating LSTmin per interval
    nObsMin = 20 # minimum no. of observation for every step interval

    print "Upper percentile of LST for calculating LSTmax: 99th percentile"
    print "Lower percentile of LST for calculating LSTmin: 1st percentile"
    print "Step interval used: 0.025"

    # get the no. of interval for the computation
    ni = int(round((Nmax - Nmin)/0.025))+1

    idxValue = []  # list for the mid NDVI value within an interval
    LSTmin=[]   # list of the 1st percentile LST for every NDVI interval
    LSTmax = [] # list of the 95th percentile LST for every NDVI interval

    # for every step within NDVI interval:
    for r in range(0,ni-1):
        rMin = r*step + Nmin  # minimum interval value
        rMax = r*step + Nmin + step # maximum interval value

        SelNdvi = np.where((vNdvi >= rMin) & (vNdvi <= rMax), 1, 0) #filter pixels with NDVI value within the interval

        temp_ = np.extract(SelNdvi, vLst) # get equivalent LST pixels of the filtered NDVIs
        temp = temp_[~np.isnan(temp_)]

        if np.size(temp)> nObsMin:  # check the number of pixels observations greater than 20
            LSTmaxpt = np.percentile(temp, UpPercentile)   #get the LST 99th percentile
            LSTminpt = np.percentile(temp, LoPercentile)    #get the LST 1st percentile
            LSTmax.append(LSTmaxpt)             # put the 99th percentile value on a list - LSTmax for the regression line calculation
            LSTmin.append(LSTminpt)             # put the 1st percentile value on a list - LST min for the wet edge calculation
            midNdvi = (rMax - rMin)*2 + rMax  # get the mid value of the NDVI interval
            idxValue.append(midNdvi)          # put mid value NDVI on a list - NDVI for the regression line calculation
        else: print "There is not enough LST observations in this interval."  # message if number of observations is less than 20

    # compute for the dry edge regression line
    # polyfit, with 1 degree polynomial computes for a linear fit, same with an ordinary least square regression OLS)
    b,a = polyfit(idxValue, LSTmax, 1)     # returns the slope (b) and the intercept (a), respectively
    print "Dry edge = {}".format(a), " + ({}".format(b), " * NDVI)"

    # Compute the dry edge from OLS regression
    Dry_edge = a + (b * np.array(vNdvi))

    # Compute for cold edge (LSTmin) - the minimum among the 1st percentile of LST
    Cold_edge =  np.min(LSTmin)    # Minimum of the 10th percentile LST for every NDVI interval
    print "Cold edge = %3f" % Cold_edge

    # mask the NDVI and the LST values for the scatter plot - showing the acceptable range
    Nmask = np.where((vNdvi > 0.1) & (vNdvi <= 1), 1, 0)
    Lmask = np.where((vLst > 0) & (vLst <= 70), 1, 0)
    Cmask = Nmask * Lmask
    NMasked = np.choose(Cmask, (np.nan, vNdvi))
    LMasked = np.choose(Cmask, (np.nan, vLst))

    # plot NDVI vs LST, dry edge and cold edge
    plt.scatter(NMasked, LMasked, c='green', alpha=0.3, s=1, zorder=1)
    plt.plot(vNdvi,Dry_edge, '-', color = 'red', zorder=2)
    #plt.hlines(Cold_edge, 0.3, midNdvi, color='blue')
    plt.hlines(Cold_edge, 0.0, 1.0, color='blue', zorder=3)
    plt.scatter(idxValue, LSTmax, c='red', marker='.', alpha=1, s=45, zorder=4) #- this will show the scatter plot for the dry edge points
    plt.scatter(idxValue, LSTmin, c='blue', marker='v', alpha=1, s=45, zorder=5) # - this will show the scatter plot for the wet edge points
    plt.minorticks_on()
    plt.grid(alpha=1)
    plt.gcf()
    plt.xlabel('Normalized Difference Vegetation Index')
    plt.xlim([0,1])
    plt.ylabel('Land Surface Temperature ($^\circ$C)')
    plt.ylim([10,55])
    plt.savefig(sPlot)
    #plt.show()
    plt.close()

    """ TVDI calculation  """
    #outData = np.empty((rows_Ndvi, cols_Ndvi), np.dtype(float))  # create an empty array for the output data
    outData = (arrLSTm - Cold_edge) / (a + b*arrNdvi - Cold_edge) # calculate TVDI   
    array_to_raster(outData,tvdiname)

##################################################################################################

Dir = "D:\z_To_delete\hh" # directory of NDVI and LST files
srcdir = os.path.abspath(Dir)  # setting the default directory as working directory
os.chdir(srcdir)
    
# search for ndvi and lst files
ndvifiles = glob.glob("*NDVI.tif")
lstfiles = glob.glob("*LST.tif")

ds = gdal.Open(ndvifiles[0])
dsTrans = ds.GetGeoTransform()

outDir = srcdir+"\\outTVDI"
if not os.path.exists(outDir):
    os.mkdir(outDir)

# calculate TVDI for all matching NDVI-LST input files
for n,l in itertools.izip(ndvifiles,lstfiles):
    MOD, yr, doy, tile, name = n.split("_") #MOD_yr2015_doy001_h30v08_LST.ti
    print "Processing TVDI for %s & %s" %(n,l)
    plotname = outDir + "%s_%s_%s_%s_plot.png" %(MOD, yr, doy, tile)
    tvdiname = outDir + "%s_%s_%s_%s_TVDI.tif" %(MOD, yr, doy, tile)
    calcTVDI(n, l, plotname,tvdiname)
