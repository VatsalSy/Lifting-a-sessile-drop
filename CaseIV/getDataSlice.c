/* Title: Interpolating data from dump files: gfs2oogl style
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"

char filename[80];
int nx, nz;
double xmin, zmin, ySlice, xmax, zmax;
scalar f1[], f2[], * interfaces = {f1, f2};
scalar * list = NULL;

int main(int a, char const *arguments[]){

  L0 = 16;
  origin(-L0/2, - L0/2, 0.);

  // boundary conditions
  u.t[back] = dirichlet(0.);
  u.r[back] = dirichlet(0.);
  f1[back] = dirichlet(0.);
  f2[back] = dirichlet(0.);

  sprintf (filename, "%s", arguments[1]);
  xmin = atof(arguments[2]);
  xmax = atof(arguments[3]);
  ySlice = atof(arguments[4]);
  zmin = atof(arguments[5]);
  zmax = atof(arguments[6]);

  nx = atoi(arguments[7]);
  nz = atoi(arguments[8]);

  list = list_add (list, f1);
  list = list_add (list, f2);
  list = list_add (list, u.x);
  list = list_add (list, u.z);

  restore (file = filename);
  boundary((scalar *){f1, f2, u.x, u.y, u.z});

  FILE * fp = ferr;
  nx++;
  nz++;
  double Deltax = 0.999999*(xmax-xmin)/(nx - 1);
  double Deltaz = 0.999999*(zmax-zmin)/(nz - 1);
  int len = list_len(list);
  double ** field = (double **) matrix_new (nx, nz, len*sizeof(double));
  for (int i = 0; i < nx; i++) {
    double x = Deltax*i + xmin;
    for (int j = 0; j < nz; j++) {
      double z = Deltaz*j + zmin;
      int k = 0;
      for (scalar s in list){
        field[i][len*j + k++] = interpolate (s, x, ySlice, z);
      }
    }
  }

  for (int i = 0; i < nx; i++) {
    double x = Deltax*i + xmin;
    for (int j = 0; j < nz; j++) {
      double z = Deltaz*j + zmin;
      fprintf (fp, "%g %g %g", x, ySlice, z);
      int k = 0;
      for (scalar s in list){
        fprintf (fp, " %g", field[i][len*j + k++]);
      }
      fputc ('\n', fp);
    }
  }
  fflush (fp);
  matrix_free (field);

}
