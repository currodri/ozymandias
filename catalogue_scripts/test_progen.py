import ozy
from ozy import progen

obj1 = ozy.load('/mnt/extraspace/currodri/NUT/cosmoNUThd/Groups/ozy_00017.hdf5')
obj2 = ozy.load('/mnt/extraspace/currodri/NUT/cosmoNUThd/Groups/ozy_00016.hdf5')

my_progens = progen.progen_build(obj1, obj2,'/mnt/extraspace/currodri/NUT/cosmoNUThd/Groups/ozy_00017.hdf5')

print(my_progens)
