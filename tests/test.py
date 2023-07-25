import ozy
from ozy.read_HaloMaker import build_HaloMaker
obj = ozy.OZY('/mnt/extraspace/currodri/NUT/cosmoNUTcrmhd/output_00010',catalogue_reader=build_HaloMaker)
obj.run_halocatalogue(main_gal=True)

obj.save('test_00010.hdf5')

