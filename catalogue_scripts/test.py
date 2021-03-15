import yt
import ozy  
ds = yt.load('/mnt/extraspace/currodri/NUT/cosmoNUThd/output_00007/info_00007.txt')
obj = ozy.OZY(ds)
obj.build_HaloMaker()

obj.save('test_00007.hdf5')

