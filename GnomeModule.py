# -*- coding: utf-8 -*-

import gnome.scripting as gs  # many modules are in this folder -init file.py
import datetime
import os
from gnome.utilities.projections import GeoProjection
from gnome.scripting import subsurface_plume_spill
from gnome.utilities import distributions
from gnome.movers import GridWindMover

# class GnomeProcess:    

def process (model_strt,                 
             model_timeRange,
             wind_file_path,
             wind_cnst_values,
             current_file_path,
             current_components,
             output_directory_path,
             release_strt,
             release_end,
             release_location,
             pollutant_amount,
             pollutant_unit):
    '''
    Process function calculates and generates pollutant trace plots and coressponding NetCDF file

    :param model_strt: model start time
    :type model_strt: datetime.datetime

    :param model_timeRange: model time range
    :type model_timeRange: datetime.timedelta

    :param wind_file_path: path to the user-provided wind NetCDF file
    :type wind_file_path: string

    :param wind_cnst_values: wind velocity components in term of u and v
    :type wind_cnst_values: 2-tuple of floats (u, v)

    :param current_file_path: path to the user-provided current NetCDF file (ROMS)
    :type current_file_path: string

    :param output_directory_path: user-provided directory for storing output NetCDF file
    :type output_directory_path: string

    :param release_strt: start time of spill
    :type release_strt: datetime.datetime

    :param release_end: end time of spill
    :type release_end: datetime.datetime   

    :param release_location: the location of spill in lon, lat and depth
    :type release_location: 3-tuple of floats (long, lat, z)

    :param pollutant_amount: amount of pollutant
    :type pollutant_amount: float

    :param pollutant_unit: amount unit of pollutant ('BBL', 'kg')
    :type pollutant_unit: string    

    :param pollutant_material: type of pollutant ('MEG', 'Natural Gas', 'Condensate')
    :type pollutant_material: string
    '''

    # model_strt            = datetime.datetime(1992, 5, 10, 12, 0, 0, 0)
    # model_timeRange       = datetime.timedelta(days = 2)
    # wind_file_path        = u'/media/iwf/Archive/Farrokh/Personal/Programming/Python_Projects/GUI/New_Project/withQFileDialog/wind_1992.nc'
    # wind_cnst_values      = None # (10, 90)  # (wind speed, wind direction from North)
    # current_file_path     = u'/media/iwf/Archive/Farrokh/Personal/Programming/Python_Projects/GUI/New_Project/withQFileDialog/his_00005.nc'
    # current_components    = None # (-1, -1.2, 0) # (u,v,w)
    # output_directory_path = u'/media/iwf/Archive/Farrokh/Personal/Programming/Python_Projects/GUI/New_Project/withQFileDialog'
    # release_strt          = datetime.datetime(1992, 5, 10, 12, 0, 0, 0)
    # release_end           = datetime.datetime(1992, 5, 10, 18, 0, 0, 0)
    # release_location      = (52.0, 27.0, 5.0)
    # pollutant_amount      = 1000
    # pollutant_unit        = 'bbl'
    # pollutant_material    = u'SHARJAH CONDENSATE'

    # 1- Model Definition    
    mymap = gs.MapFromBNA('coast.bna', refloat_halflife = 1)
    model = gs.Model(start_time = model_strt,
                     duration = model_timeRange,
                     time_step = gs.minutes(30))
    model.map = mymap

    # 2- Define spill
    ud = distributions.UniformDistribution(10e-6,300e-6) #droplets in the range 10-300 microns
    spill = gs.subsurface_plume_spill(release_time = release_strt,
                                      end_release_time = release_end,
                                      distribution_type = 'droplet_size',
                                      distribution = ud,
                                      start_position = (release_location[0], release_location[1], 0),  # fixing 0 as depth, for now
                                      num_elements = 1000,
                                      amount = pollutant_amount,
                                      units = pollutant_unit,
                                      substance = 'SHARJAH CONDENSATE')
    model.spills += spill

    # 3- Environment Objects
    water = gs.Water(temperature = 300.0,
                     salinity = 35.0,
                     sediment = 0.005,
                     wave_height = None,
                     fetch = None,
                     units = None,
                     name = 'Water')

    # 4- Movers
    #  4-1- Wind
    if wind_file_path != None:
        wind_mover = GridWindMover(filename = wind_file_path)
        model.movers += wind_mover
    else:
        uniform_wind = gs.constant_wind(wind_cnst_values[0],
                                        wind_cnst_values[1],
                                        units = 'm/s')
        uniform_wind_mover = gs.WindMover(uniform_wind)                                             
        model.movers += uniform_wind_mover

    #  4-2- Current
    grid_topo = {'node_lon': 'lon_rho',
                'node_lat': 'lat_rho'}
    if current_file_path != None:
        current_mover = gs.PyCurrentMover(filename = current_file_path,
                                            grid_topology = grid_topo)
        model.movers += current_mover
    else:
        current_velocity = (current_components[0],
                            current_components[1],
                            0) #(u, v, w) in m/s
        uniform_velocity_mover = gs.SimpleMover(current_velocity)
        model.movers += uniform_velocity_mover

    #  4-3- Diffusion
    model.movers += gs.RandomMover(diffusion_coef = 100000)

    # Weatherers
    # wind = gs.constant_wind(0, 0, 'knots')
    # model.weatherers += Evaporation(water, wind)

    #  4-4- Outputters
    renderer = gs.Renderer('coast.bna',
                           output_dir = 'image20',
                           image_size = (800, 600),
                           projection_class = GeoProjection,
                           output_timestep = gs.minutes(30))
    model.outputters += renderer

    netcdf_file = os.path.join(output_directory_path, 'pygnome.nc')
    gs.remove_netcdf(netcdf_file)
    nc_outputter = gs.NetCDFOutput(netcdf_file,
                                    which_data = 'most',
                                    output_timestep = gs.minutes(30))
    model.outputters += nc_outputter
    # model.outputters += gs.ShapeOutput('gnome_results',
    #                                    zip_output = True,
    #                                    output_timestep = gs.minutes(30))
    # model.outputters += gs.KMZOutput('gnome_results.kmz',
    #                                  output_timestep = gs.minutes(30))

    # 5- Run model
    x = []
    y = []
    for step in model:
        positions = model.get_spill_property('positions')
    x.append(positions[:, 0])
    y.append(positions[:, 1])
    model.list_spill_properties()
    model.full_run()