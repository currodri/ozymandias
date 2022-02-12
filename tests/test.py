import ozy  
obj = ozy.OZY('/mnt/extraspace/currodri/NUT/cosmoNUTcrmhd/output_00035')
obj.build_HaloMaker()

obj.save('test_00035.hdf5')

