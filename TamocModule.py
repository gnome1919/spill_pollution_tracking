from tamoc import dispersed_phases
from tamoc import bent_plume_model
import numpy as np
import matplotlib.pyplot as plt
from tamoc import seawater, ambient, dbm ,dbm_utilities, MEG_bent_plume_model, bent_plume_model, MEG_lmp,lmp
from tamoc import dispersed_phases
import netCDF4 as nc4
from netCDF4 import date2num
from tamoc import single_bubble_model
from pylab import figure
from os import remove
import datetime as datetime

def process(method = None,
            current_file_path = None,
            release_location = None,
            release_timestep_index = None,
            release_timestep_datetime = None,
            release_time = None,
            release_pollutant = None,
            pollutant_amount = None,
            pollutant_unit = None,
            current_components = None,
            current_sal = None,
            current_temp = None):
    '''
    Process function calculates and generates time-series  and cross-sectinal plots.
    Returns x_final and y_final in case of 'Natural Gas' to use in AERMOD input file.

    Parameters
    ----------
    :param method: determines whether the calculations are for constant 
                   current conditions or user has provided a pre-generated
                   current file (ROMS) - ('file' or 'constant')
    :type method: string

    :param current_file_path: path to the user-provided current NetCDF file (ROMS)
                              (when method is 'file')
    :type current_file_path: string

    :param release_location: the location of spill in lon, lat and depth
    :type release_location: 3-tuple of floats (long, lat, z)

    :param release_timestep_index: time-step number of initial spill
    :type release_timestep_index: integer

    :param release_timestep_datetime: time of initial spill
    :type release_timestep_datetime: datetime.datetime

    :param release_time: total time of release in seconds
    :type release_time: float

    :param release_pollutant: type of pollutant ('MEG', 'Natural Gas', 'Condensate')
    :type release_pollutant: string

    :param pollutant_amount: amount of pollutant released
    :type pollutant_amount: float

    :param pollutant_unit: amount unit of pollutant released ('kg', 'bbl')
    :type pollutant_unit: string

    :param current_components: current velocity components in form u, v and w when
                               method is 'constant'
    :type current_components: 3-tuple of floats (u, v, w)

    :param current_sal: water salinity (PSU) when method is 'constant'
    :type current_sal: float

    :param current_temp: water temprature (C) when method is 'constant'
    :type current_temp: float

    Returns
    -------
    nc : `netCDF4.Dataset` object
        An object that contains the new netCDF dataset object/
    '''

    # method = 'file'
    # current_file_path = u'/media/iwf/Archive/Farrokh/Personal/Programming/Python_Projects/GUI/New_Project/withQFileDialog/his_00005.nc'
    # release_location = (52.0, 27.0, 25.0)
    # release_timestep_index = 0
    # release_timestep_datetime = datetime(1992, 5, 10, 12, 0, 0, 0)
    # release_pollutant = 'MEG'

    # method = 'constant'
    # release_location = (52.0, 27.0, 25.0)
    # release_strt = datetime(1992, 5, 10, 12, 0, 0, 0)
    # release_pollutant = 'MEG'
    # current_components = (-1, -1.2)
    # current_sal = 40
    # current_temp = 20

    x_final = None
    y_final = None

    # 1- Ambient profile
    #  1-1- Method1 - from ROMS output
    #   Read lon lat depth from netcdf
    if method == 'file':
        fn = nc4.Dataset(current_file_path)
        lon_rho = fn['lon_rho'][:]
        lat_rho = fn['lat_rho'][:]
        lon = release_location[0]        # longitude of release
        lat = release_location[1]        # latitude of release
        t_idx = release_timestep_index   # time of release
        j_idx,i_idx = np.where((np.abs(lon_rho[:] - lon) < 0.017)&(np.abs(lat_rho[:] - lat) < 0.017))
        j_idx = np.ndarray.tolist(j_idx)
        j_idx = j_idx[0]
        i_idx = np.ndarray.tolist(i_idx)
        i_idx = i_idx[0]

        (nc, current_file_path) = ambient.get_nc_db_from_roms(current_file_path, 'TAMOC.tmp', t_idx,j_idx, i_idx,['u_eastward','v_northward'])
        current_file_path.close()
        roms = ambient.Profile(nc)
        roms.nc.close()
        roms = ambient.Profile('TAMOC.tmp', chem_names = 'all')
        z = roms.nc.variables['z'][:]
        T = roms.nc.variables['temperature'][:]
        S = roms.nc.variables['salinity'][:]
        P = roms.nc.variables['pressure'][:]
        u = roms.nc.variables['u_eastward'][:]
        v = roms.nc.variables['v_northward'][:]
        ua = u
        va = v
        data1 = np.vstack((z, ua)).transpose()
        symbols1 = ['z', 'ua']
        units = ['m', 'm/s']
        comments = ['measured', 'arbitrary crossflow velocity']
        roms.append(data1, symbols1, units, comments, 0)
        data2 = np.vstack((z, va)).transpose()
        symbols2 = ['z', 'va']
        roms.append(data2, symbols2, units, comments, 0)
        roms.close_nc()
        z = roms.z
        z0 = release_location[2]

    #  1-2- Method2 - Ambient constant
    else:
        z = np.array([0.0, 5, 10, 20, 30])

        # # method1: if T and S are  default values
        #   T = np.array([21.0, 20.0, 19, 18, 17]) + 273.15
        #   S = np.array([40, 41, 42, 45, 48])
        # # method2: if T and S are constant from user:
        T = current_temp
        S = current_sal
        T = np.array([T, T, T, T, T]) + 273.15
        S = np.array([S, S, S, S, S])

        profile = np.vstack((z, T, S)).transpose()
        summary = 'Dataset'
        source = 'Synthetic data for idealized laboratory conditions'
        sea_name = 'Persian Gulf'
        p_lat = release_location[1]    # lon of release point
        p_lon = release_location[0]    # lat of release point
        p_time = date2num(release_timestep_datetime, units='seconds since 1970-01-01 00:00:00 0:00', calendar='julian')   # start time of release
        nc = ambient.create_nc_db('TAMOC.tmp', summary, source, sea_name, p_lat, p_lon, p_time)

        # Insert the CTD data into the netCDF dataset
        comments = ['synthetic'] * 3
        nc = ambient.fill_nc_db(nc, profile, ['z', 'temperature', 'salinity'], ['m', 'K', 'psu'], comments, 0)

        # Calculate and insert the pressure data
        z = nc.variables['z'][:]
        T = nc.variables['temperature'][:]
        S = nc.variables['salinity'][:]
        P = ambient.compute_pressure(z, T, S, 0)
        P_data = np.vstack((z, P)).transpose()
        nc = ambient.fill_nc_db(nc, P_data, ['z', 'pressure'], ['m', 'Pa'], ['measured', 'computed'], 0)

        u = current_components[0]   # east component of currnt profile - magnitude and direction values should have already been converted to u and v components
        v = current_components[1]   # north component of currnt profile - magnitude and direction values should have already been converted to u and v components
        nc = ambient.Profile(nc)

        ua = np.array([u, u, u, u, u])
        va = np.array([v, v, v, v, v])
        data1 = np.vstack((z, ua)).transpose()
        symbols1 = ['z', 'ua']
        units = ['m', 'm/s']
        comments = ['measured', 'arbitrary crossflow velocity']
        nc.append(data1, symbols1, units, comments, 0)
        data2 = np.vstack((z, va)).transpose()
        symbols2 = ['z', 'va']
        nc.append(data2, symbols2, units, comments, 0)
        nc.close_nc()
        roms = nc
        z = roms.z
        z0 = release_location[2]

    # 2- Simulation (By Materials)
    if pollutant_unit == 'bbl':
        pollutant_amount = pollutant_amount * 139.456  # 1 bbl = 139.456 kg
    # discharge = pollutant_amount / release_time
    mb0 = (pollutant_amount / release_time) / 2  # half of total mass flux in kg/s (mass / release time in s)    
    #  2-1- Natural Gas
    if release_pollutant != 'MEG':
        if release_pollutant == 'Natural Gas':
            Vj = 0  # Q=32 m^3/day or 4.7 m^3/s
            phi_0 = -np.pi / 2.
            theta_0 = 0.
            D = 0.8128  # 32 inch= 32*2.54= 81.28 cm
            Tj = 273.15 + 35.
            Sj = 35
            cj = 0.1
            chem_name = 'MEG'
            Ta, Sa, Pj = roms.get_values(z0, ['temperature', 'salinity','pressure'])
            # Create the stratified plume model object
            bpm = bent_plume_model.Model(roms)
            gas_compounds = ['methane', 'ethane', 'propane', 'isobutane', 'n-butane']
            gas_fractions = np.array([0.939, 0.042, 0.0184, 0.0003, 0.0003])
            bubble = dbm.FluidParticle(gas_compounds)
            disp_phases = []
            # Larger free gas bubbles
            # mb0 = 222.7  # total mass flux in kg/s (mass / release time in s)  1 bbl = 139.456 kg
            de = 0.005  # bubble diameter in m
            lambda_1 = 0.85
            (m0, T0, nb0, P, Sa, Ta) = dispersed_phases.initial_conditions(
                roms, 0, bubble, gas_fractions, mb0, 2, de, Tj)
            disp_phases.append(bent_plume_model.Particle(0, 0, 0, bubble, m0, T0,
                                        nb0, lambda_1, P, Sa, Ta, K=1., K_T=1., fdis=1.e-6, t_hyd=0.,
                                        lag_time=False))
            # Smaller free gas bubbles
            # mb0 = 222.7  # total mass flux in kg/s (mass / release time in s)
            de = 0.0005  # bubble diameter in m
            lambda_1 = 0.95
            (m0, T0, nb0, P, Sa, Ta) = dispersed_phases.initial_conditions(
                roms, z0, bubble, gas_fractions, mb0, 2, de, Tj)
            disp_phases.append(bent_plume_model.Particle(0., 0., z0, bubble, m0, T0,
                                                        nb0, lambda_1, P, Sa, Ta, K=1., K_T=1., fdis=1.e-6, t_hyd=0.,
                                                        lag_time=False))
            # Run the simulation
            bpm.simulate(np.array([0, 0, z0]), D, Vj, phi_0, theta_0, Sj, Tj, cj, tracers=chem_name, particles=disp_phases,
                        track=True, dt_max=60., sd_max=1000)

            # Added for new x_final and y_final (should be checked)
            # x_final = bpm.q_local.x[-1]
            # y_final = bpm.q_local.y[-1]

            bpm.save_sim('NaturalGas.nc', '.', 'Gas_bent_plume_model')
            bpm.save_txt('NaturalGas', '.', 'Gas_bent_plume_model')

            ### OLD CODE ###
            # sbm = single_bubble_model.Model(profile=roms)
            # composition = ['methane', 'ethane', 'propane', 'oxygen']
            # gas = dbm.FluidParticle(composition, fp_type=0.)
            # mol_frac = np.array([0.90, 0.07, 0.03, 0.0])
            # de = 0.005
            # z0 = z[-1]
            # T0 = 273.15 + 30.
            # sbm.simulate(gas, z0, de, mol_frac, T0, K=0.5, K_T=1, fdis=1e-8, delta_t=1.)
            # zz = sbm.y[:, 2]
            # xx = sbm.y[:, 0]
            # yy = sbm.y[:, 1]
            # tt = sbm.t
            # x_final = sbm.y[-1, 0]
            # y_final = sbm.y[-1, 1]

        #  2-2- Condensate Droplet
        elif release_pollutant == 'Condensate':
            Vj = 0  # Q=32 m^3/day or 4.7 m^3/s
            phi_0 = -np.pi / 2.
            theta_0 = 0.
            D = 0.8128 # 32 inch= 32*2.54= 81.28 cm
            Tj = 273.15 + 35.
            Sj = 35
            cj = 0.1
            chem_name = 'MEG'
            Ta, Sa, Pj = roms.get_values(z0, ['temperature', 'salinity','pressure'])
            # Create the stratified plume model object
            bpm = bent_plume_model.Model(roms)

            # gas_compounds = ['methane', 'ethane', 'propane', 'isobutane', 'n-butane']
            # gas_fractions = np.array([0.939, 0.042, 0.0184, 0.0003, 0.0003])
            from tamoc import dbm_utilities

            [composition, mass_frac, user_data, delta, units]=dbm_utilities.load_adios_oil('AD01728') # sharjah condensate

            drop = dbm.FluidParticle(composition,fp_type=1.,delta=delta,user_data=user_data)
            disp_phases = []
            #larger
            # mb0 = 28.86  # total mass flux in kg/s (mass / release time in s)
            de = 0.005  # bubble diameter in m
            lambda_1 = 0.85
            (m0, T0, nb0, P, Sa, Ta) = dispersed_phases.initial_conditions(
                roms, 0, drop, mass_frac, mb0*2, 2, de, Tj) # because we remove d0=0005 bubble diameter mb0 change to 2*mb0
            disp_phases.append(bent_plume_model.Particle(0, 0, 0, drop, m0, T0,
                                        nb0, lambda_1, P, Sa, Ta, K=1., K_T=1., fdis=1.e-6, t_hyd=0.,
                                        lag_time=False))

            # Smaller free gas bubbles
            # mb0 = 28.86  # total mass flux in kg/s  (mass / release time in s)
            # de = 0.0005  # bubble diameter in m
            # lambda_1 = 0.95
            # (m0, T0, nb0, P, Sa, Ta) = dispersed_phases.initial_conditions(
            #     roms, z0, drop, mass_frac, mb0, 2, de, Tj) 
            # disp_phases.append(bent_plume_model.Particle(0., 0., z0, drop, m0, T0,
            #                                             nb0, lambda_1, P, Sa, Ta, K=1., K_T=1., fdis=1.e-6, t_hyd=0.,
            #                                             lag_time=False))
            # Run the simulation
            bpm.simulate(np.array([0, 0, z0]), D, Vj, phi_0, theta_0, Sj, Tj, cj, tracers=chem_name, particles=disp_phases,
                        track=True, dt_max=60., sd_max=1000)

            bpm.save_sim('Condensate.nc', '.', 'Cond_bent_plume_model')
            bpm.save_txt('Condensate', '.', 'Cond_bent_plume_model')

            ### OLD CODE ###
            # sbm = single_bubble_model.Model(roms)
            # [composition, mass_frac, user_data, delta, units] = dbm_utilities.load_adios_oil('SHARJAH CONDENSATE') # sharjah condensate
            # drop = dbm.FluidParticle(composition, fp_type=1., delta=delta, user_data=user_data)
            # de = 0.005
            # z0 = z[-1]
            # T0 = 273.15 + 30.
            # sbm.simulate(drop, z0, de, mass_frac, T0, K_T=1, fdis=1e-8, delta_t=10.)
            # zz = sbm.y[:, 2]
            # xx = sbm.y[:, 0]
            # yy = sbm.y[:, 1]
            # tt = sbm.t

        # Plot the output
        t = bpm.t
        q = bpm.q
        LagElement = bent_plume_model.LagElement
        q_local = bpm.q_local
        p = bpm.p
        particles = bpm.particles
        q0_local = LagElement(t[0], q[0, :], q_local.D, roms, p, particles,
                            q_local.tracers, q_local.chem_names)
        n_part = q0_local.np
        #=================
        pchems = 1
        for i in range(n_part):
            if len(particles[i].composition) > pchems:
               pchems = len(particles[i].composition)
         # Store the derived variables
        M = np.zeros(t.shape)
        S = np.zeros(t.shape)
        T = np.zeros(t.shape)
        Mpf = np.zeros((len(t), n_part, pchems))
        Hp = np.zeros((len(t), n_part))
        Mp = np.zeros((len(t), n_part))
        Tp = np.zeros((len(t), n_part))
        xp = np.zeros((len(t), 3*n_part))
        u = np.zeros(t.shape)
        v = np.zeros(t.shape)
        w = np.zeros(t.shape)
        V = np.zeros(t.shape)
        h = np.zeros(t.shape)
        x = np.zeros(t.shape)
        y = np.zeros(t.shape)
        z = np.zeros(t.shape)
        s = np.zeros(t.shape)
        rho = np.zeros(t.shape)
        b = np.zeros(t.shape)
        cos_p = np.zeros(t.shape)
        sin_p = np.zeros(t.shape)
        cos_t = np.zeros(t.shape)
        sin_t = np.zeros(t.shape)
        rho_a = np.zeros(t.shape)
        Sa = np.zeros(t.shape)
        Ta = np.zeros(t.shape)
        ua = np.zeros(t.shape)
        E = np.zeros(t.shape)
                
        #=============================
        # x = np.zeros(t.shape)
        # y = np.zeros(t.shape)
        # z = np.zeros(t.shape)
        # cos_p = np.zeros(t.shape)
        # sin_p = np.zeros(t.shape)
        # cos_t = np.zeros(t.shape)
        # sin_t = np.zeros(t.shape)
        # b = np.zeros(t.shape)
        for i in range(len(t)):
            if i > 0:
                q0_local.update(t[i - 1], q[i - 1, :], roms, p, particles)
            q_local.update(t[i], q[i, :], roms, p, particles)
            #======================
            M[i] = q_local.M
            S[i] = q_local.S
            T[i] = q_local.T
            for j in range(n_part):
                Mpf[i,j,0:len(q_local.M_p[j])] = q_local.M_p[j][:]
                Mp[i,j] = np.sum(particles[j].m[:])
                Tp[i,j] = particles[j].T
                xp[i,j*3:j*3+3] = q_local.x_p[j,:]
            Hp[i,:] = q_local.H_p
            u[i] = q_local.u
            v[i] = q_local.v
            w[i] = q_local.w
            V[i] = q_local.V
            h[i] = q_local.h
            x[i] = q_local.x
            y[i] = q_local.y
            z[i] = q_local.z
            s[i] = q_local.s
            rho[i] = q_local.rho
            b[i] = q_local.b
            cos_p[i] = q_local.cos_p
            sin_p[i] = q_local.sin_p
            cos_t[i] = q_local.cos_t
            sin_t[i] = q_local.sin_t
            rho_a[i] = q_local.rho_a
            Sa[i] = q_local.Sa
            Ta[i] = q_local.Ta
            ua[i] = q_local.ua
            E[i] = lmp.entrainment(q0_local, q_local, p)
        # Compute the unit vector along the plume axis
        Sz = sin_p
        Sx = cos_p * cos_t
        Sy = cos_p * sin_t


            #========================
        # x[i] = q_local.x
        # y[i] = q_local.y
        # z[i] = q_local.z
        # b[i] = q_local.b
        # cos_p[i] = q_local.cos_p
        # sin_p[i] = q_local.sin_p
        # cos_t[i] = q_local.cos_t
        # sin_t[i] = q_local.sin_t
        # Sz = sin_p
        # Sx = cos_p * cos_t
        # Sy = cos_p * sin_t

        plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_1')
        # plt.plot(t, z * -1, 'b')
        for i in range(len(particles)):
            plt.plot(particles[i].t, -1*particles[i].z, 'ob')
            # plt.plot(xp[:,i*3], -1*xp[:,i*3+2], '.--b')
            # if particles[i].integrate is False and particles[i].z > 0.:
            #     plt.plot(particles[i].t, -1*particles[i].sbm.y[:,2],
                        # '.:b')       
        plt.ylabel('Depth (m)')
        plt.xlabel('Time (s)')
        plt.suptitle('Depth-Time Plot\n' + str(release_timestep_datetime.strftime("%Y-%m-%d %H:%M UTC")))
        # plt.xlim(0, 5)
        # plt.ylim(-z0,-57)
        plt.locator_params(tight=True, nbins=6)
        plt.grid(True)
        fig = plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_1')
        fig.savefig('DepthTime_' + str(release_timestep_datetime.strftime("%Y%m%d-%H%M")))
        # plt.show()

        plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_2')
        # plt.plot(x, z * -1, 'b')
        # [x1, z1, x2, z2] = bent_plume_model.width_projection(Sx, Sz, b)
        # plt.plot(x + x1, (z + z1) * -1, 'b--')
        # plt.plot(x + x2, (z + z2) * -1, 'b--')
        #==========
        for i in range(len(particles)):
            plt.plot(particles[i].x, -1*particles[i].z, 'ob')
            plt.plot(xp[:,i*3], -1*xp[:,i*3+2], '.--b')
            if particles[i].integrate is False and particles[i].z > 0.:
                plt.plot(particles[i].sbm.y[:,0], -1*particles[i].sbm.y[:,2],
                        '.:b')
        #===========
        plt.ylabel('Depth (m)')
        plt.xlabel('Longitude-Axis (m)')
        plt.suptitle('Depth-Longitude Plot\n' + str(release_timestep_datetime.strftime("%Y-%m-%d %H:%M UTC")))
        # plt.xlim(0, 5)
        # plt.ylim(-z0,-50)
        plt.locator_params(tight=True, nbins=6)
        plt.grid(True)
        fig = plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_2')
        fig.savefig('DepthLon_' + str(release_timestep_datetime.strftime("%Y%m%d-%H%M")))
        # plt.show()

        plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_3')
        # plt.plot(y, z * -1, 'b')
        # [y1, z1, y2, z2] = bent_plume_model.width_projection(Sy, Sz, b)
        # plt.plot(y + y1, (z + z1) * -1, 'b--')
        # plt.plot(y + y2, (z + z2) * -1, 'b--')
        #===========
        for i in range(len(particles)):
            plt.plot(particles[i].y, -1*particles[i].z, 'ob')
            plt.plot(xp[:,i*3+1], -1*xp[:,i*3+2], '.--b')
            if particles[i].integrate is False and particles[i].z > 0.:
                plt.plot(particles[i].sbm.y[:,1], -1*particles[i].sbm.y[:,2],
                        '.:b')
        #===========
        # plt.ylim(-z0,-50)
        plt.ylabel('Depth (m)')
        plt.xlabel('Latitude-Axis (m)')
        plt.suptitle('Depth-Latitude Plot\n' + str(release_timestep_datetime.strftime("%Y-%m-%d %H:%M UTC")))
        plt.locator_params(tight=True, nbins=6)
        plt.grid(True)
        fig = plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_3')
        fig.savefig('DepthLat_' + str(release_timestep_datetime.strftime("%Y%m%d-%H%M")))
        # plt.show()
        plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_1').clear()
        plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_2').clear()
        plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_3').clear()
        plt.close()
        plt.cla()
        plt.clf()

        # Added for new x_final and y_final (should be checked)
        x_final = x[-1]
        print(x_final)
        y_final = y[-1]
        print(y_final)




        ### OLD CODE ###
        # figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_1')
        # # figure(1)
        # plt.plot(tt, zz*-1)
        # plt.ylabel('Depth (m)')
        # plt.xlabel('Time (s)')
        # plt.suptitle('Depth-Time Plot\n' + str(release_timestep_datetime.strftime("%Y-%m-%d %H:%M UTC")))
        # plt.locator_params(tight = True, nbins = 6)
        # plt.grid(True)
        # fig = plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_1')
        # # fig = plt.figure(1)
        # fig.savefig('DepthTime_' + str(release_timestep_datetime.strftime("%Y%m%d-%H%M")))
        # # plt.show()

        # figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_2')
        # # figure(2)
        # plt.plot(xx, zz*-1)
        # plt.ylabel('Depth (m)')
        # plt.xlabel('Longitude-Axis (m)')
        # plt.suptitle('Depth-Longitude Plot\n' + str(release_timestep_datetime.strftime("%Y-%m-%d %H:%M UTC")))
        # plt.locator_params(tight = True, nbins = 6)
        # plt.grid(True)
        # fig = plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_2')
        # # fig = plt.figure(2)
        # fig.savefig('DepthLon_' + str(release_timestep_datetime.strftime("%Y%m%d-%H%M")))
        # # plt.show()

        # figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_3')
        # # figure(3)
        # plt.plot(yy, zz*-1)
        # plt.ylabel('Depth (m)')
        # plt.xlabel('Latitude-Axis (m)')
        # plt.suptitle('Depth-Latitude Plot\n' + str(release_timestep_datetime.strftime("%Y-%m-%d %H:%M UTC")))
        # plt.locator_params(tight = True, nbins = 6)
        # plt.grid(True)
        # fig = plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_3')
        # # fig = plt.figure(3)
        # fig.savefig('DepthLat_' + str(release_timestep_datetime.strftime("%Y%m%d-%H%M")))
        # # plt.show()

    #  2-3- MEG
    else:
        # Jet initial conditions
        vol = pollutant_amount / 1110  # MEG's rho = 1110 kg/m3
        discharge = vol / release_time # Discharge of MEG (m3/s)
        D = 0.1143 # Diameter of pipe
        pipe_cs_area = np.pi * D**2 / 4 # cross-sectional area of pipe
        Vj = discharge / pipe_cs_area # Velocity of release assuming full rupture of pipe
        phi_0 = -np.pi / 2.
        theta_0 = 0.
        Tj = 273.15 + 35.
        Sj = 150.46
        cj = 0.25
        chem_name = 'water'
        Ta, Sa, Pj = roms.get_values(z0, ['temperature', 'salinity','pressure'])
        # Create the stratified plume model object
        bpm = MEG_bent_plume_model.Model(roms)
        # Run the simulation
        bpm.simulate(np.array([0, 0, z0]), D, Vj, phi_0, theta_0, Sj, Tj, cj, tracers=chem_name,dt_max=60., sd_max=5)

        bpm.save_sim('MEG.nc', '.', 'MEG_bent_plume_model')
        bpm.save_txt('MEG', '.', 'MEG_bent_plume_model')

        ### OLD CODE ###
        # Vj = 4.7  # Q = 32 m^3/day or 4.7 m^3/s
        # phi_0 = -np.pi / 2.
        # theta_0 = 0.
        # D = 0.01
        # Tj = 273.15 + 35.
        # Sj = 150.46
        # cj = 0.
        # chem_name = 'tracer'
        # Ta, Sa, Pj = roms.get_values(z[-1], ['temperature', 'salinity', 'pressure'])
        # rho_j = seawater.density(Ta, Sa, Pj)
        # ww = dbm.seawater.density(Tj,Sj,Pj)

        # # Create the stratified plume model object
        # bpm = bent_plume_model.Model(roms)

        # composition = ['methane', 'ethane', 'propane', 'oxygen']
        # yk = np.array([0.93, 0.05, 0.02, 0.0])
        # gas = dbm.FluidParticle(composition)
        # disp_phases = []

        # # Larger free gas bubbles
        # mb0 = 0.005  # total mass flux in kg/s
        # de = 0.005  # bubble diameter in m
        # lambda_1 = 0.85
        # (m0, T0, nb0, P, Sa, Ta) = dispersed_phases.initial_conditions(roms, z0, gas, yk, mb0, 2, de, Tj)
        # disp_phases.append(bent_plume_model.Particle(0., 0., z0, gas, m0, T0, nb0, lambda_1, P, Sa, Ta,
        #                                             K=1., K_T=1., fdis=1.e-6, t_hyd=0., lag_time=False))

        # # Run the simulation
        # bpm.simulate(np.array([0., 0., z0]), D, Vj, phi_0, theta_0, Sj, Tj, cj, tracers=chem_name,particles=disp_phases, track=True, dt_max=60.,sd_max=2000.)

        # # Plot the full suite of model variables
        # t = bpm.t
        # q = bpm.q
        # LagElement = bent_plume_model.LagElement
        # q_local = bpm.q_local
        # p = bpm.p
        # particles = bpm.particles
        # q0_local = LagElement(t[0], q[0, :], q_local.D, roms, p, particles,
        #                     q_local.tracers, q_local.chem_names)
        # n_part = q0_local.np
        # x = np.zeros(t.shape)
        # y = np.zeros(t.shape)
        # z = np.zeros(t.shape)
        # cos_p = np.zeros(t.shape)
        # sin_p = np.zeros(t.shape)
        # cos_t = np.zeros(t.shape)
        # sin_t = np.zeros(t.shape)
        # b = np.zeros(t.shape)
        # for i in range(len(t)):
        #     if i > 0:
        #         q0_local.update(t[i - 1], q[i - 1, :], roms, p, particles)
        #     q_local.update(t[i], q[i, :], roms, p, particles)
        #     x[i] = q_local.x
        #     y[i] = q_local.y
        #     z[i] = q_local.z
        #     b[i] = q_local.b
        #     cos_p[i] = q_local.cos_p
        #     sin_p[i] = q_local.sin_p
        #     cos_t[i] = q_local.cos_t
        #     sin_t[i] = q_local.sin_t
        # Sz = sin_p
        # Sx = cos_p * cos_t
        # Sy = cos_p * sin_t

        # PLOT
        t = bpm.t
        q = bpm.q

        # Create a second Lagrangian element in order to compute entrainment
        q0_local = MEG_bent_plume_model.LagElement(t[0], q[0, :], bpm.q_local.D, bpm.profile, bpm.p,
                            bpm.q_local.tracers)
        # Store the derived variables
        M = np.zeros(t.shape)
        S = np.zeros(t.shape)
        T = np.zeros(t.shape)
        u = np.zeros(t.shape)
        v = np.zeros(t.shape)
        w = np.zeros(t.shape)
        V = np.zeros(t.shape)
        h = np.zeros(t.shape)
        x = np.zeros(t.shape)
        y = np.zeros(t.shape)
        z = np.zeros(t.shape)
        s = np.zeros(t.shape)
        rho = np.zeros(t.shape)
        b = np.zeros(t.shape)
        cos_p = np.zeros(t.shape)
        sin_p = np.zeros(t.shape)
        cos_t = np.zeros(t.shape)
        sin_t = np.zeros(t.shape)
        rho_a = np.zeros(t.shape)
        Sa = np.zeros(t.shape)
        Ta = np.zeros(t.shape)
        ua = np.zeros(t.shape)
        E = np.zeros(t.shape)

        for i in range(len(t)):
            if i > 0:
                q0_local.update(t[i - 1], q[i - 1, :], bpm.profile, bpm.p)
            bpm.q_local.update(t[i], q[i, :], bpm.profile, bpm.p)
            M[i] = bpm.q_local.M
            S[i] = bpm.q_local.S
            T[i] = bpm.q_local.T
            u[i] = bpm.q_local.u
            v[i] = bpm.q_local.v
            w[i] = bpm.q_local.w
            V[i] = bpm.q_local.V
            h[i] = bpm.q_local.h
            x[i] = bpm.q_local.x
            y[i] = bpm.q_local.y
            z[i] = bpm.q_local.z
            s[i] = bpm.q_local.s
            rho[i] = bpm.q_local.rho
            b[i] = bpm.q_local.b
            cos_p[i] = bpm.q_local.cos_p
            sin_p[i] = bpm.q_local.sin_p
            cos_t[i] = bpm.q_local.cos_t
            sin_t[i] = bpm.q_local.sin_t
            rho_a[i] = bpm.q_local.rho_a
            Sa[i] = bpm.q_local.Sa
            Ta[i] = bpm.q_local.Ta
            ua[i] = bpm.q_local.ua
            E[i] = MEG_lmp.entrainment(q0_local, bpm.q_local, bpm.p)

        # Compute the unit vector along the plume axis
        Sz = sin_p
        Sx = cos_p * cos_t
        Sy = cos_p * sin_t

        # x = q[:, 7]
        # y = q[:, 8]
        # z = q[:, 9]

        plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_1')
        plt.plot(t, z * -1, 'b')
        plt.ylabel('Depth (m)')
        plt.xlabel('Time (s)')
        plt.suptitle('Depth-Time Plot\n' + str(release_timestep_datetime.strftime("%Y-%m-%d %H:%M UTC")))
        # plt.xlim(0, 5)
        # plt.ylim(-z0,-57)
        plt.locator_params(tight=True, nbins=6)
        plt.grid(True)
        fig = plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_1')
        fig.savefig('DepthTime_' + str(release_timestep_datetime.strftime("%Y%m%d-%H%M")))
        # plt.show()

        plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_2')
        plt.plot(x, z * -1, 'b')
        [x1, z1, x2, z2] = MEG_bent_plume_model.width_projection(Sx, Sz, b)
        plt.plot(x + x1, (z + z1) * -1, 'b--')
        plt.plot(x + x2, (z + z2) * -1, 'b--')
        plt.ylabel('Depth (m)')
        plt.xlabel('Longitude-Axis (m)')
        plt.suptitle('Depth-Longitude Plot\n' + str(release_timestep_datetime.strftime("%Y-%m-%d %H:%M UTC")))
        # plt.xlim(0, 5)
        # plt.ylim(-z0,-50)
        plt.locator_params(tight=True, nbins=6)
        plt.grid(True)
        fig = plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_2')
        fig.savefig('DepthLon_' + str(release_timestep_datetime.strftime("%Y%m%d-%H%M")))
        # plt.show()

        plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_3')
        plt.plot(y, z * -1, 'b')
        [y1, z1, y2, z2] = MEG_bent_plume_model.width_projection(Sy, Sz, b)
        plt.plot(y + y1, (z + z1) * -1, 'b--')
        plt.plot(y + y2, (z + z2) * -1, 'b--')
        # plt.ylim(-z0,-50)
        plt.ylabel('Depth (m)')
        plt.xlabel('Latitude-Axis (m)')
        plt.suptitle('Depth-Latitude Plot\n' + str(release_timestep_datetime.strftime("%Y-%m-%d %H:%M UTC")))
        plt.locator_params(tight=True, nbins=6)
        plt.grid(True)
        fig = plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_3')
        fig.savefig('DepthLat_' + str(release_timestep_datetime.strftime("%Y%m%d-%H%M")))
        # plt.show()
        fig.clear()
        plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_1').clear()
        plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_2').clear()
        plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_3').clear()
        plt.close()
        plt.cla()
        plt.clf()
        
        # ### OLD CODE ###
        # plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_1')
        # plt.plot(t, z * -1,'b')
        # plt.ylabel('Depth (m)')
        # plt.xlabel('Time (s)')
        # plt.suptitle('Depth-Time Plot\n' + str(release_timestep_datetime.strftime("%Y-%m-%d %H:%M UTC")))
        # # plt.xlim(0, 5)
        # plt.ylim(-50.8,-38)
        # plt.locator_params(tight = True, nbins = 6)
        # plt.grid(True)
        # fig = plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_1')
        # fig.savefig('DepthTime_' + str(release_timestep_datetime.strftime("%Y%m%d-%H%M")))
        # # plt.show()

        # plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_2')
        # plt.plot(x, z * -1,'b')
        # [x1, z1, x2, z2] = bent_plume_model.width_projection(Sx, Sz, b)
        # plt.plot(x + x1, (z + z1)*-1, 'b--')
        # plt.plot(x + x2, (z + z2)*-1, 'b--')
        # plt.ylabel('Depth (m)')
        # plt.xlabel('Longitude-Axis (m)')
        # plt.suptitle('Depth-Longitude Plot\n' + str(release_timestep_datetime.strftime("%Y-%m-%d %H:%M UTC")))
        # # plt.xlim(0, 5)
        # # plt.ylim(-50.8,-38)
        # plt.locator_params(tight = True, nbins = 6)
        # plt.grid(True)
        # fig = plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_2')
        # fig.savefig('DepthLon_' + str(release_timestep_datetime.strftime("%Y%m%d-%H%M")))
        # # plt.show()

        # plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_3')
        # plt.plot(y, z*-1,'b')
        # [y1, z1, y2, z2] = bent_plume_model.width_projection(Sy, Sz, b)
        # plt.plot(y + y1, (z + z1)*-1, 'b--')
        # plt.plot(y + y2, (z + z2)*-1, 'b--')
        # # plt.ylim(-50.8,-40)
        # plt.ylabel('Depth (m)')
        # plt.xlabel('Latitude-Axis (m)')
        # plt.suptitle('Depth-Latitude Plot\n' + str(release_timestep_datetime.strftime("%Y-%m-%d %H:%M UTC")))
        # plt.locator_params(tight = True, nbins = 6)
        # plt.grid(True)
        # fig = plt.figure(str(release_timestep_datetime.strftime("%Y%m%d-%H%M"))+'_3')
        # fig.savefig('DepthLat_' + str(release_timestep_datetime.strftime("%Y%m%d-%H%M")))
        # # plt.show()

    remove('TAMOC.tmp')
    print('Aermod area source')
    print(x_final,y_final) # values of x_final, y_final
    return (x_final, y_final)
# pr#int('Aermod area source')
# print(x_final,y_final) # values of x_final, y_final
    # return (52.0, 27.0)000   000