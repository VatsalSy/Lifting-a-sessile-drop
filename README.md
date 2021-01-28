# Lifting-a-sessile-drop

This repository contains the codes used for simulating the cases discussed in the manuscript: Lifting a sessile oil drop from a superamphiphobic surface with an impacting one [(link)](https://advances.sciencemag.org/content/6/34/eaba4330.abstract).


You can use the following article citation:

```
@article{ramirez2020lifting,
  title={Lifting a sessile oil drop from a superamphiphobic surface with an impacting one},
  author={Ram{\'\i}rez-Soto, O. and Sanjay, V. and Lohse, D. and Pham, J. T. and Vollmer, D.},
  journal={Science advances},
  volume={6},
  number={34},
  pages={eaba4330},
  year={2020},
  publisher={American Association for the Advancement of Science}
} 
```

We investigate the dynamics of an oil drop impacting an identical sessile drop sitting on a superamphiphobic surface.
On this page, I am presenting the code that we used to simulate the process shown in the above video. For a more detailed documentation on the code included in the manuscript, please visit [my Basilisk sandbox](http://basilisk.fr/sandbox/vatsal/DropOnDropImpact/dropOnDropImpact.c).

## Offset parameter
In the manuscript, offset parameter $\chi$ is defined as:
$$
\chi = \frac{d}{2R}
$$
Here, $d$ is the distance between the axes of two drops, and $R$ is the equivalent radius of the drop. In the simulations, we input $2\chi$ (because the length scale is $R$).

Weber number is based on the impact velocity, $U_0$.
$$ We = \frac{\rho_lU_0^2R}{\gamma} $$
Note that the impact Weber number we input is slightly less than the actual Weber number that we need. As the impacting drop falls, it also gains some kinetic energy. So, while analyzing the results, the Weber number should be taken at the instant, which one decides to choose as $t = 0$. These calculations get tricky as one goes to a higher $\chi$. However, the results do not change much as long as $We \sim \mathcal{O}(1)$.

Navier Stokes equation for this case:
$$
\partial_tU_i+\nabla\cdot(U_iU_j) =
\frac{1}{\hat{\rho}}\left(-\nabla p + Oh\nabla\cdot(2\hat{\mu}D_{ij}) + \kappa\delta_sn_i\right) + Bog_i
$$
The $\hat{\rho}$ and $\hat{\mu}$ are the VoF equations to calculate properties, given by:
$$
\hat{A} = (f_1+f_2) + (1-f_1-f_2)\frac{A_g}{A_l}
$$

Ohnesorge number $Oh$ is a measure of viscous forces as compared to inertia and surface tension forces.
$$ Oh = \frac{\mu_l}{\sqrt{\rho_l\gamma R}} $$

Bond number $Bo$ is a measure between Gravity and surface tension.
$$ Bo = \frac{\rho_lgR^2}{\gamma} $$

**Note:** The subscript $l$ denotes liquid. Also, the radius used in the dimensionless numbers is the equivalent radius of the drops $\left(R = \left(3\pi V_l/4\right)^{1/3}\right)$. $V_l$ is the volume of the two drops.

Velocity scale as the intertial-capillary velocity,
$$ U_\gamma = \sqrt{\frac{\gamma}{\rho_l R}} $$

## Initial Condition

For Bond numbers $> 0.1$, assuming the sessile drop to be spherical is inaccurate. One can also see this in the experimental video above. So, we used the code [(here)](http://basilisk.fr/sandbox/vatsal/DropDeposition/DropDeposition.c) to get the initial shape of the sessile drop. We import an [STL file](https://www.dropbox.com/s/uenhig7lfhvss66/Sessile-Bo0.3080.stl?dl=0).

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
