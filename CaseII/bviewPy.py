# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# This is a python script that runs [getVideo.c](https://www.dropbox.com/s/yk71hlrio9pnywz/getVideo.c?dl=0). The idea is to check for the intermediate snapshot file, then check if an image at this time-step already exists or not. And, finally generate image if and only if there is a snapshot file and there is no image. This helps if the computer crashes (or more commonly shuts down unexpectedly) or if the snapshot file is faulty.

import numpy as np
import os
nGFS = 751 # I still call it nGFS because of historical reasons :D

folder = 'bview'  # output folder
if not os.path.isdir(folder):
    os.makedirs(folder)

for ti in range(nGFS):
    t = 0.01 * ti
    place = "intermediate/snapshot-%5.4f" % t
    name = "%s/%4.4d.ppm" %(folder, ti)
    if not os.path.exists(place):
        print("File %s not found!" % place)
    else:
        if os.path.exists(name):
            print("Image %s found!" % name)
        else:
            exe = "./getVideo %s %s" % (place, name)
            os.system(exe)
            print(("Done %d of %d" % (ti, nGFS)))
