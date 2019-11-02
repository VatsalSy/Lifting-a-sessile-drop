/* Title: Getting Energy
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

scalar f1[], f2[];
double ke1, ucm1, vcm1, wcm1, se1, gpe1, ke2, ucm2, vcm2, wcm2, se2, gpe2, ke3, gpe3, rho1, rho2, mu1, mu2, eps;
char nameEnergy[80], nameOut[80];

#define RHO(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#define MU(f)  (clamp(f,0.,1.)*(mu1 - mu2) + mu2)

#define Mu21 (6e-3)
#define Rho21 (1./770)
#define Ldomain 16

int main(int a, char const *arguments[]) {
  sprintf(nameOut, "%s", arguments[1]);
  sprintf(nameEnergy, "%s", arguments[2]);
  double Oh = 0.0216; // <0.003/sqrt(770*0.025*0.001)>
  double Bo = 0.308; // <770*10*0.001^2/0.025>

  FILE *fp;
  fp = fopen (nameEnergy, "a");
  restore (file = nameOut);

  rho1 = 1.0; mu1 = Oh;
  rho2 = Rho21; mu2 = Mu21*Oh;

  // boundary conditions
  u.t[back] = dirichlet(0.);
  u.r[back] = dirichlet(0.);
  f1[back] = dirichlet(0.);
  f2[back] = dirichlet(0.);

  #if TREE
    f1.prolongation = fraction_refine;
    f2.prolongation = fraction_refine;
  #endif
  boundary((scalar *){f1, f2, u.x, u.y, u.z});

  ke1 = 0., ucm1 = 0., vcm1 = 0., wcm1 = 0., se1 = 0., gpe1 = 0., ke2 = 0., ucm2 = 0., vcm2 = 0., wcm2 = 0., se2 = 0., gpe2 = 0., ke3 = 0., gpe3 = 0., eps = 0.;

  double sumU1 = 0., sumU2 = 0.;
  double sumV1 = 0., sumV2 = 0.;
  double sumW1 = 0., sumW2 = 0.;
  double wt1 = 0., wt2 = 0.;

  foreach (){
    ke1 += (0.5*clamp(f1[], 0., 1.)*rho1*(sq(u.x[]) + sq(u.y[]) + sq(u.z[])))*cube(Delta);
    ke2 += (0.5*clamp(f2[], 0., 1.)*rho1*(sq(u.x[]) + sq(u.y[]) + sq(u.z[])))*cube(Delta);
    ke3 += (0.5*clamp(1.-f1[]-f2[], 0., 1.)*rho2*(sq(u.x[]) + sq(u.y[]) + sq(u.z[])))*cube(Delta);

    sumU1 += clamp(f1[], 0., 1.)*u.x[]*cube(Delta); sumU2 += clamp(f2[], 0., 1.)*u.x[]*cube(Delta);
    sumV1 += clamp(f1[], 0., 1.)*u.y[]*cube(Delta); sumV2 += clamp(f2[], 0., 1.)*u.y[]*cube(Delta);
    sumW1 += clamp(f1[], 0., 1.)*u.z[]*cube(Delta); sumW2 += clamp(f2[], 0., 1.)*u.z[]*cube(Delta);
    wt1 += clamp(f1[], 0., 1.)*cube(Delta); wt2 += clamp(f2[], 0., 1.)*cube(Delta);

    gpe1 += clamp(f1[], 0., 1.)*rho1*z*cube(Delta)*Bo;
    gpe2 += clamp(f2[], 0., 1.)*rho1*z*cube(Delta)*Bo;
    gpe3 += clamp(1.-f1[]-f2[], 0., 1.)*rho2*z*cube(Delta)*Bo;

    double D2 = 0.;
    foreach_dimension(){
      double DII = (u.x[1,0,0]-u.x[-1,0,0])/(2*Delta);
      double DIJ = 0.5*((u.x[0,1,0]-u.x[0,-1,0] + u.y[1,0,0] - u.y[-1,0,0])/(2*Delta));
      double DIK = 0.5*((u.x[0,0,1]-u.x[0,0,-1] + u.z[1,0,0] - u.z[-1,0,0])/(2*Delta));
      D2 += sq(DII) + sq(DIJ) + sq(DIK);
    }
    eps += ( 2*MU(f1[]+f2[])*D2 )*cube(Delta);
  }

  sumU1 /= wt1; sumU2 /= wt2;
  sumV1 /= wt1; sumV2 /= wt2;
  sumW1 /= wt1; sumW2 /= wt2;

  ucm1 = sumU1; vcm1 = sumV1; wcm1 = sumW1;
  ucm2 = sumU2; vcm2 = sumV2; wcm2 = sumW2;

  se1 = (interface_area (f1)-4*pi);
  se2 = (interface_area (f2)-4*pi);
  boundary((scalar *){f1, f2, u.x, u.y, u.z});

  fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", t, ke1, ucm1, vcm1, wcm1, se1, gpe1, ke2, ucm2, vcm2, wcm2, se2, gpe2, ke3, gpe3, eps);
  fclose(fp);
}
