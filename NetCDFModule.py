import numpy as np
import scipy.interpolate
import netCDF4

def process():

    '''
    Process function generates a full grid NetCDF file containing surface concentration
    data from PyGnome NetCDF output file named 'pygnome.nc'
    
    Parameters
    ----------
    None

    Returns
    -------
    None
    '''

    """
    Initializing Variables
    """
    surface_concentrations = []
    particle_count = []
    longitudes = []
    latitudes = []
    times = []

    """
    Reading Source File Data
    """
    pygnome_dataset = netCDF4.Dataset('pygnome.nc', 'r')
    surface_concentrations[:] = pygnome_dataset.variables['surface_concentration'][:]
    particle_count[:] = pygnome_dataset.variables['particle_count'][:]
    longitudes[:] = pygnome_dataset.variables['longitude'][:]
    latitudes[:] = pygnome_dataset.variables['latitude'][:]
    times[:] = pygnome_dataset.variables['time'][:]
    times_units = pygnome_dataset.variables['time'].units
    # times_start = dt.datetime.strptime(times_units, 'seconds since %Y-%m-%dT%H:%M:%S')
    # ts_size = times[1] - times[0]

    """
    Creating Destination Grid
    """
    lon_TimeStep_range = (min(longitudes), max(longitudes))
    lat_TimeStep_range = (min(latitudes), max(latitudes))
    lon_min = round(round(lon_TimeStep_range[0] / 0.001) * 0.001, 3) - 0.01 # rounding a variable to 0.001 decision 
    lon_max = round(round(lon_TimeStep_range[1] / 0.001) * 0.001, 3) + 0.01
    lons = np.round(np.arange(lon_min, lon_max+.001, 0.001), 3)
    lat_min = round(round(lat_TimeStep_range[0] / 0.001) * 0.001, 3) - 0.01
    lat_max = round(round(lat_TimeStep_range[1] / 0.001) * 0.001, 3) + 0.01
    lats = np.round(np.arange(lat_min, lat_max+.001, 0.001), 3)
    lon_grid,lat_grid = np.meshgrid(lons, lats)

    """
    Creating NetCDF File
    """
    # Creating NetCDF file
    output_file = netCDF4.Dataset('output.nc','w', format='NETCDF4')

    # Creating dimensions
    output_file.createDimension('lon', len(lons))
    output_file.createDimension('lat', len(lats))
    output_file.createDimension('time', None)

    # Creating global attributes
    output_file.title='PyGnome Post Analysis Output File'
    output_file.source='PyGnome 0.6.5'
    output_file.author='Farrokh A. Ghavanini'
    output_file.contact='ghavanini[at]gmail[dot]com'

    # Creating variables
    lon = output_file.createVariable('lon', np.float64, ('lon',))
    lon.units = 'degrees_east'
    lon.long_name = 'longitude'
    lon.point_spacing = 'even'
    lon.datatype
    lat = output_file.createVariable('lat', np.float64, ('lat',))
    lat.units = 'degrees_north'
    lat.long_name = 'latitude'
    lat.point_spacing = 'even'
    time = output_file.createVariable('time', np.float64, ('time',))
    time.units = times_units
    time.calendar = 'gregorian'
    time.long_name = 'time'
    surface_concentration = output_file.createVariable('surfcon', np.float64,('time','lat','lon'))
    surface_concentration.units = 'g m-2'
    surface_concentration.long_name = 'Surface Concentration of Pollutant'

    # Filling variables
    lon[:] = lons
    lat[:] = lats
    time[:] = times
    surface_concentration[:,:,:] = np.empty((0, len(lats), len(lons)))
    item = 0
    for timestep, particles in enumerate(particle_count):
        if particles == 0:
            surface_concentration[timestep,:,:] = np.NaN
            next
        else:
            lon_array = np.array(longitudes[item:item+particles])
            lat_array = np.array(latitudes[item:item+particles])
            sc_array = np.array(surface_concentrations[item:item+particles])
            zi = scipy.interpolate.griddata((lon_array,lat_array), sc_array, (lon_grid,lat_grid), method='linear')
            surface_concentration[timestep,:,:] = zi
            item = item + particles

    # Closing NetCDF file
    output_file.close()