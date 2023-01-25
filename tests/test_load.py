import ozy
obj = ozy.load('test_00035.hdf5')
print(obj._halo_dicts.keys())
print(obj._halo_dicts['radius'].keys())
print(obj._halo_dicts['radius']['BT87_simple'][:])
print(obj._halo_data.keys())
print(obj._galaxy_dicts.keys())
print(obj._galaxy_dicts['mass'].keys())
print(obj._galaxy_data.keys())