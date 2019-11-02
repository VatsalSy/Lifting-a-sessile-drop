# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# This script run the Basilisk code [getEnergy.c](https://www.dropbox.com/s/nc7oof3rj5oy51d/getEnergy.c?dl=0). First it looks for a file called getEnergy.dat. If the file exists, then the data is appended to this file, else a new file is created. Then it checks for existence of an intermediate snapshot file. This helps if the computer crashes (or more commonly shuts down unexpectedly) or if the snapshot file is faulty.

import numpy as np
import sys
import os

nGFS = 2001
Oh = sys.argv[1]
print("Oh = %s" % Oh)

name = "getEnergy.dat"
if os.path.exists(name):
    print("File %s found! New data will be appended to the file" % name)

for ti in range(nGFS):
    t = 0.01 * ti
    place = "intermediate/snapshot-%5.4f" % t
    if not os.path.exists(place):
        print("File %s not found!" % place)
    else:
        exe = "./getEnergy %s %s %s" % (place, name, Oh)
        os.system(exe)
        print(("Done %d of %d" % (ti, nGFS)))
