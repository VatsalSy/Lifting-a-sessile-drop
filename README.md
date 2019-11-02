# Lifting-a-sessile-drop
This repository contains the codes used for simulating the cases discussed in the manuscript: Lifting a sessile oil drop with an impacting one.
We investigate the dynamics of an oil drop impacting an identical sessile drop sitting on a superamphiphobic surface.
On this page, I am presenting the code that we used to simulate the process shown in the above video. The results presented here are currently under review in Science Advances. For a detailed documentation on the code included in the manuscript, please visit [my Basilisk sandbox](http://basilisk.fr/sandbox/vatsal/DropOnDropImpact/dropOnDropImpact.c).

## Offset parameter
In the manuscript, offset parameter <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/c91091e68f0e0113ff161179172813ac.svg?invert_in_darkmode" align=middle width=10.28535419999999pt height=14.15524440000002pt/> is defined as:
<p align="center"><img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/8be2c7b35dc31b2c1be3422b699bb84a.svg?invert_in_darkmode" align=middle width=55.0032648pt height=33.81208709999999pt/></p>
Here, <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/2103f85b8b1477f430fc407cad462224.svg?invert_in_darkmode" align=middle width=8.55596444999999pt height=22.831056599999986pt/> is the distance between the axes of two drops, and <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/1e438235ef9ec72fc51ac5025516017c.svg?invert_in_darkmode" align=middle width=12.60847334999999pt height=22.465723500000017pt/> is the equivalent radius of the drop. In the simulations, we input <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/f8ae52fcf4026682a03f086114e491eb.svg?invert_in_darkmode" align=middle width=18.50456354999999pt height=21.18721440000001pt/> (as the length scale is <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/1e438235ef9ec72fc51ac5025516017c.svg?invert_in_darkmode" align=middle width=12.60847334999999pt height=22.465723500000017pt/>).

Weber number is based on the impact velocity, <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/1049ded9e1be03670fc3963966339893.svg?invert_in_darkmode" align=middle width=17.77628489999999pt height=22.465723500000017pt/>.
<p align="center"><img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/ab6f3e030eebe6eb5d41ccd935f29613.svg?invert_in_darkmode" align=middle width=95.89606455pt height=38.973783749999996pt/></p>
Note that the impact Weber number we input is slightly less than the actual Weber number that we need. As the impacting drop falls, it also gains some kinetic energy. So, while analyzing the results, the Weber number should be taken at the instant, which one decides to choose as <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/477a717e18587a5e8605780ca167c322.svg?invert_in_darkmode" align=middle width=36.07293689999999pt height=21.18721440000001pt/>. These calculations get tricky as one goes to a higher <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/c91091e68f0e0113ff161179172813ac.svg?invert_in_darkmode" align=middle width=10.28535419999999pt height=14.15524440000002pt/>. However, the results do not change much as long as <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/75fa9dc7f0d8c6fb3f37e4310669f397.svg?invert_in_darkmode" align=middle width=81.92810834999999pt height=24.65753399999998pt/>.

Navier Stokes equation for this case:
<p align="center"><img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/258c6f95c2dd30f6a4ec0dc9824e8cc1.svg?invert_in_darkmode" align=middle width=446.44984065pt height=36.1865163pt/></p>
The <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/2ece8916c80529609c5cc5d5b4e259f4.svg?invert_in_darkmode" align=middle width=9.728951099999989pt height=22.831056599999986pt/> and <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/fb2c407771af04095047a75aab1127e2.svg?invert_in_darkmode" align=middle width=9.90492359999999pt height=22.831056599999986pt/> are the VoF equations to calculate properties, given by:
<p align="center"><img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/9b8ab16b267f7c08d0729415c87acbe7.svg?invert_in_darkmode" align=middle width=232.0403052pt height=36.09514755pt/></p>

Ohnesorge number <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/2788af332492730c5b3e9e0ced0b5eb3.svg?invert_in_darkmode" align=middle width=22.46654024999999pt height=22.831056599999986pt/>: measure between surface tension and viscous forces.
<p align="center"><img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/a28e62e456a1a050953094068c748b83.svg?invert_in_darkmode" align=middle width=95.6323434pt height=33.4857765pt/></p>

Bond number <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/114800de7d48de0b11aa97b6bc03bf31.svg?invert_in_darkmode" align=middle width=21.261454499999992pt height=22.465723500000017pt/>: measure between Gravity and surface tension.
<p align="center"><img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/247bc2ba9dd1f55c5556f5b27b6cfff2.svg?invert_in_darkmode" align=middle width=87.10952414999998pt height=38.973783749999996pt/></p>

**Note:** The subscript <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/2f2322dff5bde89c37bcae4116fe20a8.svg?invert_in_darkmode" align=middle width=5.2283516999999895pt height=22.831056599999986pt/> denotes liquid. Also, the radius used in the dimensionless numbers is the equivalent radius of the drops <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/9d03528dfb8f202c7f882f6e865eef2d.svg?invert_in_darkmode" align=middle width=136.86086039999998pt height=37.80850590000001pt/>. <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/a2120f1cc6443e4cbd3a1013f8b4510e.svg?invert_in_darkmode" align=middle width=13.81283969999999pt height=22.465723500000017pt/> is the volume of the two drops.

Velocity scale as the intertial-capillary velocity,
<p align="center"><img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/5e30d2a96a3517ae8472d2244e90f985.svg?invert_in_darkmode" align=middle width=88.11261524999999pt height=39.452455349999994pt/></p>

## Initial Condition

For Bond numbers <img src="https://rawgit.com/in	git@github.com:VatsalSy/Lifting-a-sessile-drop/master/svgs/185b2578ddb033601bd1df2165be3210.svg?invert_in_darkmode" align=middle width=38.35617554999999pt height=21.18721440000001pt/>, assuming the sessile drop to be spherical is inaccurate. One can also see this in the experimental video above. So, we used the code [(here)](http://basilisk.fr/sandbox/vatsal/DropDeposition/DropDeposition.c) to get the initial shape of the sessile drop. We import an [STL file](https://www.dropbox.com/s/uenhig7lfhvss66/Sessile-Bo0.3080.stl?dl=0).

## Adaptive Mesh Refinement
We refine based on curvatures of the two drops along with the generally used VoF and velocity fields. This ensures that the refinement level along the interface is MAXlevel.

## Running the code

Use the following procedure:

**Step 1:** Importing the stl file and generating the first dump file

~~~bash
#!/bin/bash
mkdir intermediate
qcc -fopenmp -O2 -Wall dropOnDropImpact.c -o dropOnDropImpact -lm
export OMP_NUM_THREADS=8
./dropOnDropImpact
~~~

**Step 2:** Follow the method described [(here)](http://basilisk.fr/src/Tips#running-on-supercomputers). Do not forget to use the dump file generated in the previous step.
