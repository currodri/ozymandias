import yt
import ozy  
ds = yt.load('/mnt/extraspace/currodri/NUT/cosmoNUThd/output_00018/info_00018.txt')
obj = ozy.OZY(ds)
obj.build_HaloMaker()

obj.save('test_00018.hdf5')

