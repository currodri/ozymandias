import ozy
from amr2 import dictionary_commons

myvars = {'density':1,
          'vx':2,
          'vy':3,
          'vz':4,
          'thermal_pressure':5}

mydict = dictionary_commons.dictf90()

mydict.init(5)

for key,value in myvars.items():
    mydict.add(key,value)
    print(mydict.get(key.ljust(128)))
print(mydict.values)
