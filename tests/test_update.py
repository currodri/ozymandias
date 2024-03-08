import ozy
import numpy as np
obj = ozy.load('test_00010.hdf5')

gal = obj.galaxies[0]
orig_ndm = gal.ndm
orig_pos = np.array(gal.position.d)
orig_mass = float(gal.mass['gas'].d)
print(gal.ndm,gal.mass['gas'])
gal.update_attribute('ndm',100)
gal.update_attribute('position',obj.array([-1.,-1.,-1.],'code_length'))
gal.update_attribute('mass.gas',obj.quantity(0.0,'code_mass'))
print(gal.ndm,gal.mass['gas'])
topology = {'Bx':obj.quantity(1.0,'code_magnetic'),
            'By':obj.quantity(1.0,'code_magnetic'),
            'Bz':obj.quantity(1.0,'code_magnetic')}
gal.update_attribute('topology',topology)
print(gal.topology, obj._galaxy_dicts['topology']['Bx'][:])
extra_info = {'color':1,
              'hey':obj.quantity(1.0,'code_mass'),
              'acc':np.array([1.0,2.0,3.0]),
              'ne':obj.array([1.0,1.0,1.0],'code_mass')}
gal.update_attribute('extra_info',extra_info)
print(gal.extra_info)
gal.update_attribute('noldstars',300)
print(gal.noldstars)
del obj,gal

obj = ozy.load('test_00010.hdf5')
gal = obj.galaxies[0]
print(gal.ndm,gal.mass['gas'])
print(gal.topology, obj._galaxy_dicts['topology']['Bx'][:])
print(gal.noldstars)
gal.update_attribute('ndm',orig_ndm)
gal.update_attribute('position',obj.array(orig_pos,'code_length'))
gal.update_attribute('mass.gas',obj.quantity(orig_mass,'code_mass'))
topology = {'Bx':obj.quantity(1.0,'code_magnetic'),
            'By':obj.quantity(2.0,'code_magnetic'),
            'Bz':obj.quantity(3.0,'code_magnetic')}
gal.update_attribute('topology',topology)
print(gal.topology,obj._galaxy_dicts['topology']['Bx'][:])
print(gal.ndm,gal.mass['gas'])

new_bx = obj.quantity(20.0, 'code_magnetic')
gal.update_attribute('topology.Bx',new_bx)
print(gal.topology,obj._galaxy_dicts['topology']['Bx'][:])