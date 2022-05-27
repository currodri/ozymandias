import ozy  
obj = ozy.OZY('/mnt/extraspace/currodri/NUT/cosmoNUThd/output_00012')
obj.build_HaloMaker()

obj.save('test_00012.hdf5')

