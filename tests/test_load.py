import ozy
obj = ozy.load('test_00012.hdf5')
print(obj.galaxies[0].velocity)
print(obj._halo_dicts['radius']['BT87_simple'][:])