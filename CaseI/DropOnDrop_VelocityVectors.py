## Python script to generate images with velocity vectors. End result should be similar to [this](https://www.dropbox.com/s/d1wjik0emrg3yw9/CaseIII_NumericsVelocityVectors.mp4?dl=0). It runs the Basilisk script, getDataSlice.c
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

import numpy as np
import os
import subprocess as sp
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

def gettingfield(filename):
    print('Field values')
    exe = ["./getDataSlice", filename, str(xmin), str(xmax), str(ySlice), str(zmin), str(zmax), str(nx), str(nz)]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    Xtemp, Ztemp, f1temp, f2temp, Utemp, Wtemp = [] , [], [], [], [], []

    if (len(temp2) > 1e2):
        for n1 in range(len(temp2)):
            temp3 = temp2[n1].split(" ")
            if temp3 == ['']:
                pass
            else:
                Xtemp.append(float(temp3[0]))
                Ztemp.append(float(temp3[2]))
                f1temp.append(float(temp3[3]))
                f2temp.append(float(temp3[4]))
                Utemp.append(float(temp3[5]))
                Wtemp.append(float(temp3[6]))

        X = np.asarray(Xtemp)
        Z = np.asarray(Ztemp)
        f1 = np.asarray(f1temp)
        f2 = np.asarray(f2temp)
        U = np.asarray(Utemp)
        W = np.asarray(Wtemp)

        X.resize((nz+1, nx+1))
        Z.resize((nz+1, nx+1))
        f1.resize((nz+1, nx+1))
        f2.resize((nz+1, nx+1))
        U.resize((nz+1, nx+1))
        W.resize((nz+1, nx+1))

        print('Got Field values')
        return X, Z, f1, f2, U, W
    else:
        return Xtemp, Ztemp, f1temp, f2temp, Utemp, Wtemp
# ------------------------------------------------------------------------------

nGFS = 751

folder = 'TracerVelocityVideo'  # output folder
if not os.path.isdir(folder):
    os.makedirs(folder)


xmin = -2.0
xmax = 2.0
ySlice = 0.0
zmin = 0.0
zmax = 4.0
zmaxp = 4.0

LEVEL = 8
nx = 2**(LEVEL)
nz = 2**(LEVEL)

Orange = [0.9100, 0.4100, 0.1700];

for ti in range(nGFS):
    t = ti*0.01
    place = "intermediate/snapshot-%5.4f" % t
    name = "%s/%4.4d.png" %(folder, ti)
    if not os.path.exists(place):
        print("File %s not found!" % place)
    else:
        if os.path.exists(name):
            print("Image %s found!" % name)
        else:
            X, Z, f1, f2, U, W = gettingfield(place)
            if (len(X)):
                fig, ax = plt.subplots()
                fig.set_size_inches(19.20, 10.80)
                rc('axes', linewidth=2)
                plt.xticks(fontsize=30)
                plt.yticks(fontsize=30)

                ## V
                speed = np.sqrt(U**2 + W**2)
                Vmax = 2.0
                cntrl2 = ax.imshow(speed.transpose(), interpolation='bilinear', cmap="Blues", origin='lower', extent=[xmin, xmax, zmin, zmax], vmax = Vmax, vmin = 0.0)

                maxs = speed.max()
                ndx = 10
                ndz = 10
                vScale = 20.0
                if maxs > 0:
                    ax.quiver(X[::ndx,::ndz], Z[::ndx,::ndz], U[::ndx,::ndz]/Vmax, W[::ndx,::ndz]/Vmax, scale=vScale, color='black',linewidth=3)

                ax.contour(X, Z, f1, levels=[0.5], colors='Orange',linewidths=4)
                ax.contour(X, Z, f2, levels=[0.5], colors='Orange',linewidths=4)

                ax.set_xlabel(r'$X$', fontsize=50)
                ax.set_ylabel(r'$Z$', fontsize=50)
                ax.set_aspect('equal')
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(zmin, zmaxp)
                t2 = t - 0.15
                ax.set_title('t = %5.4f' % t2, fontsize=40)

                for axis in ['top','bottom','left','right']:
                    ax.spines[axis].set_linewidth(2)

                l, b, w, h = ax.get_position().bounds

                cb2 = fig.add_axes([l+w+0.01, b, 0.03, h])
                c2 = plt.colorbar(cntrl2,cax=cb2,orientation='vertical')
                c2.ax.tick_params(labelsize=30)
                c2.set_label(r'$\|U_i\|$',fontsize=40)
                c2.ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
                # plt.show()
                plt.savefig(name, bbox_inches="tight",dpi=150)
                plt.close()
            else:
                print("Problem in the available file %s" % place)

    print(("Done %d of %d" % (ti+1, nGFS)))
