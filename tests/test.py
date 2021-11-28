import ozy  
obj = ozy.OZY('/mnt/extraspace/currodri/NUT/cosmoNUThd/output_00010')
obj.build_HaloMaker()

obj.save('test_00010.hdf5')

