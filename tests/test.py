import ozy  
obj = ozy.OZY('/mnt/extraspace/currodri/NUT/cosmoNUTcrmhd/output_00010')
obj.build_HaloMaker(main_gal=True)

obj.save('test_00010.hdf5')

