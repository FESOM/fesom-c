#!/usr/bin/env python
# coding: utf-8

# In[1]:

"""
short: Script regrids FESOM-C data to transects with custom settings.

full:
Script for regridding FESOM-C ocean model output. Key features:
- Regrids unstructured data to transects, defined by lines from shapefiles
  (lfile) drawable in QGIS(for example, python example below). Each line in lfile represents a separate transect.
- Interpolates data from sigma to z coordinates, preserving original data.
- User-defined parameters:
  - distance_delta: Desired horizontal resolution of lines.
  - influence: Radius for searching data for interpolation.
  - kk: Number of points used in interpolation.
  - var: List of variables to process.
  - tsteps: Time steps to process (e.g., tsteps=np.arange(1,6899,200)).
            An empty list means all time steps.
  - zlevs: Z levels (e.g., zlevs=np.arange(0,20,1)).
  - infile: Input file name.
  - outfile: Base name for output files (.nc will be added).
  - lfile: Shapefile with transect lines.

Ideal for oceanographers analyzing FESOM-C outputs, this script supports 
data conversion for detailed transect analysis with customizable resolutions 
and interpolation settings.

python script for line definition:
    point1 = (longitude1, latitude1)
    point2 = (longitude2, latitude2)
    line = LineString([point1, point2])
    gdf = gpd.GeoDataFrame(crs="EPSG:4326", geometry=[line])
    gdf.to_file("my_line_shapefile.shp")

@author: Ivan Kuznetsov (kuivi)
""" 

import numpy as np
from netCDF4 import Dataset
import pyproj
import os
from scipy.spatial import cKDTree

import geopandas as gpd
from shapely.geometry import LineString, Point
from shapely.ops import unary_union


path2='./'
infile = path2+'LE.nc'
outfile = path2+'lines_regular_'# .nc will be added
lfile = 'alongx.shp'
path2output = path2
path2png    = path2+'png/'
if not os.path.exists(path2png): os.makedirs(path2png)

#define resolution x,y and variable name
#resolution = [1000,1000]
distance_delta = 1000 #resolution of the lines
influence=10000 # sear Radi
kk = 10 # number of points to use for interpolation
var=['eta','u','v','temperature', 'salinity']    # variable name
tsteps=[]    #time steps indxes, []=all
#tsteps=np.array([-3,-2,-1])    #time steps indxes, []=all
tsteps=np.array([0])    #time steps indxes, []=all
tsteps=np.arange(0,2100,21)
#tindw=np.array([2])
#tindw=np.array([len(intime)-1])
zlevs = np.arange(0,20,0.5) #define z levels
nzlevs = len(zlevs)


# Read the shapefile
gdf = gpd.read_file(lfile)
crs = gdf.crs
# Get the CRS of the GeoDataFrame
gdfc = gdf.to_crs(crs=3995)
# Iterate through each line and generate points
lengths=[]
for idx, row in gdfc.iterrows():
    line = row['geometry']
    lengths.append(line.length)

multipoint = []
ldist = []
for idx, row in gdfc.iterrows():
    line = row['geometry']
    distances = np.arange(0, lengths[idx], distance_delta)
    points = [line.interpolate(distance) for distance in distances] #+ [line.boundary[1]]
    multipoint.append(unary_union(points))  # or new_line = LineString(points)
    ldist.append(distances)
# Create a GeoDataFrame
gdfmc = gpd.GeoDataFrame({'geometry': multipoint}, geometry='geometry', crs=3995)
gdfm = gdfmc.to_crs(crs=crs)


# Read the shapefile
#gdf = gpd.read_file(lfile)
linex=[]
liney=[]
nlines = len(gdfm['geometry'])
for l in gdfm['geometry']:
    #points = list(l)
    lx=[]
    ly=[]
    for p in l.geoms:
        lx.append(p.x)
        ly.append(p.y)
    linex.append(np.array(lx))    
    liney.append(np.array(ly))



#open infile
ncfin = Dataset(infile)
#read lon,lat
mesh_x2=ncfin.variables['lon'][:].data
mesh_y2=ncfin.variables['lat'][:].data
mesh_x2_e=ncfin.variables['lon_elem'][:].data
mesh_y2_e=ncfin.variables['lat_elem'][:].data




#check variables
dimvar=[]
node=[]
nz=[]
print('Variables, dimension, on nodes:')  
for i,v in enumerate(var):
    dimvar.append(len(ncfin.variables[v].shape))
    node.append(len(mesh_x2) in ncfin.variables[v].shape)
    print(v,dimvar[-1],node[-1])    
    
if (np.max(dimvar)>2):
    nsig=ncfin.dimensions['sigma'].size
    nsigm1=ncfin.dimensions['sigmam1'].size
    

sigma_lev = ncfin.variables['sigma_lev'][:].data.squeeze()
sigmam1_lev = (sigma_lev[:-1]+sigma_lev[1:])/2

#check variables
nz=[]
varsigm = []
for i,v in enumerate(var):
    if ('sigma' in ncfin.variables[v].dimensions):
        nz.append(nsig)
        varsigm.append(2)
    else:
        if ('sigmam1' in ncfin.variables[v].dimensions):
            nz.append(nsigm1)
            varsigm.append(1)
        else:
            nz.append(0)
            varsigm.append(0)


radius_of_influence=influence


print("number of lines: ", nlines)
print("points number in lines: ")
s=0
for i in range(nlines):
    s=s+len(linex[i])
    print(i,': ',len(linex[i]))
print("Total points: ",s)



#read time from infile
intime = ncfin.variables['time'][:].data

if (len(tsteps)==0):
    tindw=np.arange(0,len(intime),1)
else:
    tindw=tsteps
tot_time = tindw.size

print("total number of time steps: ",tot_time)
print("time steps: ", tindw)





#convert 2d lon,lat of regular grid to a stereographic coordinates 
pste = pyproj.Proj(proj="stere", errcheck="True",ellps='WGS84', lat_0=mesh_y2.mean(), lon_0=mesh_x2.mean())

    



#prepare output file
ifiles = []
for iline in range(nlines):
    filename=outfile+str(iline)+'.nc'
    if os.path.exists(filename):
        os.remove(filename)
    ifile = Dataset(filename,'w')

    #did_lon   =  ifile.createDimension('lon',linex[iline].size);
    #did_lat   =  ifile.createDimension('lat',liney[iline].size);
    did_x   =  ifile.createDimension('x',liney[iline].size);
    did_time  =  ifile.createDimension('time',tot_time);
    did_z  =  ifile.createDimension('z',nzlevs);

    did_sigma =  ifile.createDimension('sigma', nsig)#z.size )
    did_sigmam1 =  ifile.createDimension('sigmam1', nsigm1)#z.size )

    vid_lat  = ifile.createVariable('lat','f8',('x',));
    vid_lon  = ifile.createVariable('lon','f8',('x',));
    vid_time = ifile.createVariable('time','f8',('time',));
    vid_dist = ifile.createVariable('x','f8',('x',));
    vid_depth  = ifile.createVariable('z','f8',('z',));
    
    vid_sigma  = ifile.createVariable('sigma','f8',('sigma',));
    vid_sigmam1  = ifile.createVariable('sigmam1','f8',('sigmam1',));    
    
    vid_datavar=[]    
    vid_datavar_sigm=[]        
    for i,v in enumerate(var):
        if (dimvar[i]==2):
            vid_datavar.append(ifile.createVariable(v,'f8',('time','x',)))
            vid_datavar_sigm.append(vid_datavar[-1])
        elif (dimvar[i]==1):
            vid_datavar.append(ifile.createVariable(v,'f8',('x',)))
            vid_datavar_sigm.append(vid_datavar[-1])
        else:
            vid_datavar.append(ifile.createVariable(v,'f8',('time','z','x',)))            
            if (nz[i]==nsig):
                vid_datavar_sigm.append(ifile.createVariable(v+'_s','f8',('time','sigma','x',)))
            else:
                vid_datavar_sigm.append(ifile.createVariable(v+'_s','f8',('time','sigmam1','x',)))            
                
            

    vid_lon.axis  =  ''
    vid_lon.units =  'degrees_east'
    vid_lon.standard_name =  'longitude'
    vid_lon.long_name = 'longitude'
    vid_lon[:] = linex[iline][:]
    vid_lat.axis  =  ''
    vid_lat.units =  'degrees_north'
    vid_lat.standard_name =  'latitude'
    vid_lat.long_name = 'latitude'
    vid_lat[:] = liney[iline][:]
    vid_dist.axis  =  'X'
    vid_dist.units =  'm'
    vid_dist.standard_name =  'distance'
    vid_dist.long_name = 'distance'
    vid_dist[:] = ldist[iline][:]                #################3
    vid_depth.axis  =  'Z'
    vid_depth.units =  'm'
    vid_depth.standard_name =  'depth'
    vid_depth.long_name = 'depth'
    vid_depth[:] = zlevs

    vid_sigma.axis  =  'Z'
    vid_sigma.units =  'sigma_lev'
    vid_sigma.standard_name =  'sigma_lev'
    vid_sigma.long_name = 'sigma levels'
    vid_sigma[:] = sigma_lev[:]
    
    vid_sigmam1.axis  =  'Z'
    vid_sigmam1.units =  'sigmam1_lev'
    vid_sigmam1.standard_name =  'sigmam1_lev'
    vid_sigmam1.long_name = 'sigma levels'
    vid_sigmam1[:] = sigmam1_lev[:]    
    
    #use time from infile
    vid_time.units =  ncfin.variables['time'].units
    vid_time.standard_name =  'time'
    vid_time.long_name = 'time'
    vid_time.calendar = 'noleap'
    vid_time[:] = intime[tindw]

    for i,v in enumerate(var):
        vid_datavar[i].units =  ncfin.variables[v].units
        vid_datavar[i].standard_name =  ncfin.variables[v].standard_name
        if (dimvar[i]>2):
            vid_datavar_sigm[i].units =  ncfin.variables[v].units
            vid_datavar_sigm[i].standard_name =  ncfin.variables[v].standard_name
            
            


    #make inds and distances
    (xsn, ysn) = pste(mesh_x2,mesh_y2)
    (xt, yt) = pste(linex[iline], liney[iline])
    treen = cKDTree(list(zip(xsn, ysn)))
    distancesn, indsn = treen.query(list(zip(xt, yt)), k=kk, workers=-1)

    if any([not n for n in node]):  #if any var on elements
        (xse, yse) = pste(mesh_x2_e,mesh_y2_e)
        treee = cKDTree(list(zip(xse, yse)))
        distancese, indse = treee.query(list(zip(xt, yt)), k=kk, workers=-1)

    #if ('depth' in list(ncfin.variables)):
    inds=indsn.copy()
    distances=distancesn.copy()    
    depth = ncfin.variables['depth'][:].data.squeeze()
    distances_ma = np.ma.masked_greater(distances, radius_of_influence)
    w = 1.0 / distances_ma ** 2
    depth = np.ma.sum(w * depth[inds], axis=1) / np.ma.sum(w, axis=1)                
    depth = np.ma.masked_invalid(depth)
    depth=depth.data

    if  not ('eta' in list(ncfin.variables)):        
        print("No eta in input file, 0 for eta will be used.")
    
    
    for itime,time in enumerate(tindw):
        #if (itime%(tot_time%50)==0): print(round(itime/tot_time*100),"%")
        inds=indsn.copy()
        distances=distancesn.copy()
        if ('eta' in list(ncfin.variables)):
            eta = ncfin.variables['eta'][time,:].data.squeeze()
            distances_ma = np.ma.masked_greater(distances, radius_of_influence)
            w = 1.0 / distances_ma ** 2
            eta = np.ma.sum(w * eta[inds], axis=1) / np.ma.sum(w, axis=1)                
            eta = np.ma.masked_invalid(eta)
            eta=eta.data
        else:
            eta = depth*0
        fulldepth = depth+eta    
        s,d = np.meshgrid((1-sigma_lev),fulldepth)
        sm1,dm1 = np.meshgrid((1-sigmam1_lev),fulldepth)
        z_u = d*s #2d array with z levels at each node #sigma
        z_t = (z_u[:,0:-1]+z_u[:,1:])/2.0    #sigmam1
        
        print(itime,"/",tot_time)
        for i,v in enumerate(var):
            if (node[i]):
                inds=indsn.copy()
                distances=distancesn.copy()
            else:
                inds=indse.copy()
                distancesn=distancese.copy()
            #read data
            if (dimvar[i]==2):
                datain = ncfin.variables[v][time,:].data.squeeze()
                distances_ma = np.ma.masked_greater(distances, radius_of_influence)
                w = 1.0 / distances_ma ** 2
                datai = np.ma.sum(w * datain[inds], axis=1) / np.ma.sum(w, axis=1)                
                datai = np.ma.masked_invalid(datai)
                vid_datavar[i][itime,:]=datai#datai2    
            elif (dimvar[i]==1):
                datain = ncfin.variables[v][:].data.squeeze()
                distances_ma = np.ma.masked_greater(distances, radius_of_influence)
                w = 1.0 / distances_ma ** 2
                datai = np.ma.sum(w * datain[inds], axis=1) / np.ma.sum(w, axis=1)                
                datai = np.ma.masked_invalid(datai)
                vid_datavar[i][:]=datai#datai2    
            else:
                datain = ncfin.variables[v][time,:,:].data.squeeze()
                datasigma = np.zeros((nz[i],liney[iline].size))
                
                for iz in range(nz[i]):
                    distances_ma = np.ma.masked_greater(distances, radius_of_influence)
                    w = 1.0 / distances_ma ** 2
                    datai = np.ma.sum(w * datain[inds,iz], axis=1) / np.ma.sum(w, axis=1)                
                    datai = np.ma.masked_invalid(datai)
                    datasigma[iz,:]=datai
                
                if (nz[i]==nsig):
                    z_sig = z_u
                else:
                    z_sig = z_t

                data_intp = np.zeros((nzlevs,liney[iline].size))
                data_intp[:,:] = np.nan                    
                for j in range(liney[iline].size):
                    zind = np.where(zlevs < fulldepth[j])[0][-1]+2
                    data_intp[:zind,j] = np.interp(zlevs[:zind],z_sig[j,:],datasigma[:,j],left=datasigma[0,j],right=datasigma[-1,j])
                    
                vid_datavar[i][itime,:,:]=data_intp[:,:]
                
                #add on sigma, org data
                vid_datavar_sigm[i][itime,:,:]=datasigma[:,:]
                    
    ifile.close()
           


# In[17]:


ncfin.close()


# In[18]:





# In[ ]:




