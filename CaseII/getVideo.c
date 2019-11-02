/*
## Title: Saving images with bview
# Authors:
Vatsal & Youssef<br/>
vatsalsanjay@gmail.com<br/>
Physics of Fluids<br/>
This file can be used with a [python script](https://www.dropbox.com/s/el7z0mvy68o08ct/bviewPy.py?dl=0). Python is just a choice. Any scriting language should work.
*/

#include "grid/octree.h"
#include "view.h"
#define tsnap (0.01)

scalar f1[], f2[];
char filename[80], Imagename[80];

int main(int a, char const *arguments[]) {
  sprintf (filename, "%s", arguments[1]);
  sprintf (Imagename, "%s",arguments[2]);
  restore (file = filename);
  f1[back] = dirichlet(0.);
  f2[back] = dirichlet(0.);
  boundary((scalar *){f1, f2});
  scalar m[];
  foreach(){
    m[] = 0.;
  }
  boundary((scalar *){m});
  view (fov = 4.65024, quat = {0.606693,0.0775163,0.104137,0.784265}, tx = -0.00341448, ty = -0.0938018, bg = {1,1,1}, width = 1920, height = 1136, samples = 4);
  draw_vof("f1", edges=false, fc = {0.91, 0.41, 0.17});
  draw_vof("f2", edges=false, fc = {0.91, 0.41, 0.17});
  squares("m", min=-0.5, max=0.5, n = {0, 0, 1}, map = cool_warm);
  cells(n = {0, 0, 1}, lc = {0., 0., 0.}, lw = 4);
  save (Imagename);
}
