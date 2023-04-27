#TODO: This will be a basic implementation of code to loop over all snapshots in directory/all sub-steps there-in to retrieve all data pertaining to the sink particle(s) in question
#This shouldn't be too hard, just make use of other fortran reading routines and see write_sink_fine.f90 in Ricarda's version

#Also, add all these quantities to plot_settings + dict_variables as sink_variables
#TODO: Make the sub-step loop an option, maybe just modify output_sink_csv to output coarse step values

import numpy as np



class Particle(object):
    """This is the parent class for which all Particle objects are derived
    """

    def __init__(self, obj):
        self.obj = obj
        self.ID = -1
        self.npart = -1
        self.units = 'code'
        self.mass = {}
        self.position = {}
        self.velocity  = {}
        self.angular_mom = {}

        #TODO: flesh this out properly



class Sink(Particle):
    """Sink class - for sink particles
    """
    #TODO: Make this work for multiple sink particles in the same Box

    particle_type = 'sink'
    def __init__(self, obj):
        """Initialise basic properties of the Sink Particle
        """
        super(Sink, self).__init__(obj)
        self.particle_type = 'sink'
        self.dMBHoverdt = {}
        self.spin = {}
        self.dotMbondi =  {}
        
        #TODO: In princple there should be a dictionary for each quantity which we will be saving

    def set_isolated_box_defaults(self, obj):
        """Shortcut to set default properties in the case of an isolated box sim
        """

        self.position = obj.array([0.5, 0.5, 0.5], 'code_length')
        self.velocity = obj.array([0., 0., 0.], 'km/s')
        self.angular_mom['total'] = obj.array([0., 0., 0.], 'Msun * km * km / s')

    def read_sink_particle_data(self, snapshots=[], substeps=False, n_substeps=None):
        """Routine to read all sink particle data of type sink_00000.out00000, with aim to produce time series of data
           Contains options to loop over desired snapshots, along with substeps
        """

        #TODO: Come up with a smarter way to assign the header names, maybe have some parameter file written (e.g. in the same way as for RAMSES ANR outputs)
        from scipy.io import FortranFile
        import pandas as pd

        self.path = self.obj.simulation.fullpath[:-13]# retrieves path
        for snap in snapshots:
            if substeps:
                #header_names = ['identity', 'mass', 'position_x', 'position_y', 'position_z', 'velocity_x', 'velocity_y', 'velocity_z', 'birth_time', 'dMsmbh', 'dMBH_coarse', 'dMEd_coarse', 'Esave', 'jsink_x', 'jsink_y', 'jsink_z', 'spin_x', 'spin_y', 'spin_z', 'spin_magnitude', 'epi_sink', 'sink_stat_09_01', 'sink_stat_09_02', 'sink_stat_09_03', 'sink_stat_09_04', 'sink_stat_09_05', 'sink_stat_09_06', 'sink_stat_09_07']
                #TODO
                return 
            else:
            #TODO: include 'header' array, which is looped over to initialise all dictionaries e.g. self.bondi using regex
        
                header_names = ['identity', 'mass', 'position_x', 'position_y', 'position_z', 'velocity_x', 'velocity_y', 'velocity_z', 'birth_time', 'dMBHoverdt']
                filename =  f'{self.path}/output_{str(snap).zfill(5)}/sink_{str(snap).zfill(5)}.csv'
                data = pd.read_csv(filename, names=header_names)
                
                self.position[f'{snap}'] = self.obj.array([data['position_x'][0], data['position_y'][0], data['position_z'][0]], 'code_length')
                self.velocity[f'{snap}'] = self.obj.array([data['velocity_x'][0], data['velocity_y'][0], data['velocity_z'][0]], 'km/s')
                self.dMBHoverdt[f'{snap}'] = self.obj.quantity(data['dMBHoverdt'][0], 'Msun/s')
                self.mass[f'{snap}'] = self.obj.quantity(data['mass'][0], 'Msun')

    def save_sink_particle_data(self):
        #TODO
        return



    

def create_new_particle(obj, particle_type):
    """Simple function to create a new Particle instance
    """

    if particle_type == 'sink':
        return Sink(obj)

