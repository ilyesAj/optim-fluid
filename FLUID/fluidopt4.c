 
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cblas.h"
#include <omp.h>
/*
 * function used to compute the linear position in a vector express as coordinate in a two-D structure
 */
int build_index(int i, int j, int grid_size)
//entree vect sortie un autre vect
  {
	  return (i + (grid_size + 2) * j);
  }

/*
 * Utility function to push new value in the array of a preceding time step
 *
 */

void swap(float *d, float *dOld, int vector_size)
  { 
  int i;
  float tmp;

  for (i=0 ; i < vector_size ; i ++)
    {
    tmp = d[i]; 
    d[i] = dOld[i]; 
    dOld[i] = tmp; 
    }
  }

/*
 * addSource: TO COMMENT
 */

void addSource(float *x, float *x0, int vector_size, float factor)
  {
  int i;
  //( n,  alpha, x, incx, y, incy)
  //y=y+alpha*x
  cblas_saxpy(vector_size, factor, x0, 1, x, 1);
  }

/*
 * setBoundry: TO COMMENT
 * specifies simple boundry conditions.
 */
  void build_vect2(float* x,int vector_size,int deb,int step,int coef,float* res,int deb_res,int step_res)
//entree vect sortie un autre vect
  {
    //( n,x, incx, y, incy)
    //y=x
    //step=grid_size+2
    //deb+(step)*
    cblas_scopy(vector_size,(x+deb),step,(res+deb_res),step_res);
    if (coef != 0)
      //( n,  alpha, x, incx) 
      cblas_sscal(vector_size,coef,(res+deb_res),step_res);
  }
void setBoundry(int b, float* x, int grid_size)
  {
  //omp_set_num_threads(4);
  //#pragma omp parallel 
  // build_index(int i, int j, int grid_size)
  //(i + (grid_size + 2) * j)
  int step=grid_size + 2;
  if (b==1)
  {
  	build_vect2(x,grid_size,step+1,step,-1,x,step,step);
  	build_vect2(x,grid_size,1+step,step,-1,x,(grid_size+1)*step,step);
  }
  else{
  	build_vect2(x,grid_size,step+1,step,1,x,step,step);
  	build_vect2(x,grid_size,1+step,step,1,x,(grid_size+1)*step,step);
  }
  if (b==2)
  {
  	build_vect2(x,grid_size,1+step,1,-1,x,1,1);
 	build_vect2(x,grid_size,((grid_size)*step)+1,1,-1,x,((grid_size+1)*step)+1,1);
  }
  else
  {
  	build_vect2(x,grid_size,1+step,1,1,x,1,1);
 	build_vect2(x,grid_size,((grid_size)*step)+1,1,1,x,((grid_size+1)*step)+1,1);

  }

   x[build_index(0, 0, grid_size)] = 0.5f * (x[build_index(1, 0, grid_size)] + x[build_index(0, 1, grid_size)]);
   x[build_index(0, grid_size+1, grid_size)] = 0.5f * (x[build_index(1, grid_size+1, grid_size)] + x[build_index(0, grid_size, grid_size)]);
   x[build_index(grid_size+1, 0, grid_size)] = 0.5f * (x[build_index(grid_size, 0, grid_size)] + x[build_index(grid_size+1, 1, grid_size)]);
   x[build_index(grid_size+1, grid_size+1, grid_size)] = 0.5f * (x[build_index(grid_size, grid_size+1, grid_size)] + x[build_index(grid_size+1, grid_size, grid_size)]);

   }


/*
 * Iterative linear system solver using the Gauss-sidel
 * relaxation technique. Room for much improvement here...
 *
 */

void linearSolver(int b, float* x, float* x0, float a, float c, float dt, int grid_size)
  {
  int i,j,k;

  for (k = 0; k < 20; k++)
    {
    for (i = 1; i <= grid_size; i++)
      {
      for (j = 1; j <= grid_size; j++)
        {
        x[build_index(i, j, grid_size)] = (a * ( x[build_index(i-1, j, grid_size)] + x[build_index(i+1, j, grid_size)] +   x[build_index(i, j-1, grid_size)] + x[build_index(i, j+1, grid_size)]) +  x0[build_index(i, j, grid_size)]) / c;
        }
      }
    setBoundry(b, x, grid_size);
    }
  }

/*
 * Recalculate the input array with diffusion effects.
 * Here we consider a stable method of diffusion by
 * finding the densities, which when diffused backward
 * in time yield the same densities we started with.
 * This is achieved through use of a linear solver to
 * solve the sparse matrix built from this linear system.
 * TO COMMENT
 */

void diffuse(int b, float* c, float* c0, float diff, float dt, int grid_size)
  {
  float a = dt * diff * grid_size * grid_size;

  linearSolver(b, c, c0, a, 1 + 4 * a, dt, grid_size);
  }


/*
 * Calculate the curl at position (i, j) in the fluid grid.
 * Physically this represents the vortex strength at the
 * cell. Computed as follows: w = (del x U) where U is the
 * velocity vector at (i, j).
 *
 */
float calculate_curl(int i, int j, int grid_size, float *u, float *v)
  {
  float du_dy = (u[build_index(i, j + 1, grid_size)] - u[build_index(i, j - 1, grid_size)]) * 0.5f;
  float dv_dx = (v[build_index(i + 1, j, grid_size)] - v[build_index(i - 1, j, grid_size)]) * 0.5f;

  return du_dy - dv_dx;
  }

/*
 * buoyancy: TO COMMENT
 */
float buoyancy(float *dst, float *src, int grid_size)
  {
  float Tamb = 0;
  float a = 0.000625f;
  float b = 0.025f;
  int i,j;

  float current_max=0;

  for (i = 1; i <= grid_size; i++)
    {
    for (j = 1; j <= grid_size; j++)
      {
      Tamb += src[build_index(i, j, grid_size)];
      }
    }
  // get average temperature
  Tamb /= (grid_size * grid_size);

  // for each cell compute buoyancy force
  for (i = 1; i <= grid_size; i++)
    {
    for (j = 1; j <= grid_size; j++)
      {
      dst[build_index(i, j, grid_size)] = a * src[build_index(i, j, grid_size)] + -b * (src[build_index(i, j, grid_size)] - Tamb);
      }
    }
  }


/*
 * Calculate the input array after advection. We start with an
 * input array from the previous timestep and an and output array.
 * For all grid cells we need to calculate for the next timestep,
 * we trace the cell's center position backwards through the
 * velocity field. Then we interpolate from the grid of the previous
 * timestep and assign this value to the current grid cell.
 * MORE COMMENT
 */

void advect(int b, float* d, float* d0, float* du, float* dv, int grid_size, float dt)
  {
  int i, j;
  int i0, j0, i1, j1;
  float x, y, s0, t0, s1, t1, dt0;

  dt0 = dt * grid_size;

  for (i = 1; i <= grid_size; i++)
    {
    for (j = 1; j <= grid_size; j++)
      {
      // go backwards through velocity field
      x = i - dt0 * du[build_index(i, j, grid_size)];
      y = j - dt0 * dv[build_index(i, j, grid_size)];

      // interpolate results
      if (x > grid_size + 0.5) x = grid_size + 0.5f;
      if (x < 0.5)     x = 0.5f;

      i0 = (int) x;
      i1 = i0 + 1;

      if (y > grid_size + 0.5) y = grid_size + 0.5f;
      if (y < 0.5)     y = 0.5f;

      j0 = (int) y;
      j1 = j0 + 1;

      s1 = x - i0;
      s0 = 1 - s1;
      t1 = y - j0;
      t0 = 1 - t1;

      d[build_index(i, j,grid_size)] = s0 * (t0 * d0[build_index(i0, j0, grid_size)] + t1 * d0[build_index(i0, j1, grid_size)])
                           + s1 * (t0 * d0[build_index(i1, j0, grid_size)] + t1 * d0[build_index(i1, j1, grid_size)]);
      }
    }
  setBoundry(b, d, grid_size);
  }

/*
 * Calculate the vorticity confinement force for each cell
 * in the fluid grid. At a point (i,j), Fvc = N x w where
 * w is the curl at (i,j) and N = del |w| / |del |w||.
 * N is the vector pointing to the vortex center, hence we
 * add force perpendicular to N.
 *
 */
void vorticityConfinement(float* Fvc_x, float* Fvc_y, float *curl, float *u, float *v, int grid_size)
  {
  int i,j;
  float dw_dx, dw_dy;
  float length;
  float vorticity;

  // Calculate magnitude of curl(u,v) for each cell. (|w|)
  for (i = 1; i <= grid_size; i++)
    {
    for (j = 1; j <= grid_size; j++)
      {
      curl[build_index(i, j, grid_size)] = abs(calculate_curl(i, j, grid_size, u, v));
      }
    }

  for (i = 2; i < grid_size; i++)
    {
    for (j = 2; j < grid_size; j++)
      {
      // Find derivative of the magnitude (n = del |w|)
      dw_dx = (curl[build_index(i + 1, j, grid_size)] - curl[build_index(i - 1, j, grid_size)]) * 0.5f;
      dw_dy = (curl[build_index(i, j + 1, grid_size)] - curl[build_index(i, j - 1, grid_size)]) * 0.5f;

      // Calculate vector length. (|n|)
      // Add small factor to prevent divide by zeros.
      length = (float) sqrt(dw_dx * dw_dx + dw_dy * dw_dy) + 0.000001f;

      // N = ( n/|n| )
      dw_dx /= length;
      dw_dy /= length;

      vorticity = calculate_curl(i, j, grid_size, u,v);

      // N x w
      Fvc_x[build_index(i, j, grid_size)] = dw_dy * -vorticity;
      Fvc_y[build_index(i, j, grid_size)] = dw_dx *  vorticity;
      }
    }
  }


/*
 * Use project() to make the velocity a mass conserving,
 * incompressible field. Achieved through a Hodge
 * decomposition. First we calculate the divergence field
 * of our velocity using the mean finite differnce approach,
 * and apply the linear solver to compute the Poisson
 * equation and obtain a "height" field. Now we subtract
 * the gradient of this field to obtain our mass conserving
 * velocity field.
 *
 * TO COMMENT
 *
 */

void project(float* x, float* y, float* p, float* div, float dt, int grid_size)
 {
 int i, j;
    
 for (i = 1; i <= grid_size; i++)
   {
   for (j = 1; j <= grid_size; j++)
     {
     div[build_index(i, j, grid_size)] = (x[build_index(i+1, j, grid_size)] - x[build_index(i-1, j, grid_size)] + y[build_index(i, j+1, grid_size)] - y[build_index(i, j-1, grid_size)]) * - 0.5f / grid_size;
     p[build_index(i, j, grid_size)] = 0;
     }
   }
 setBoundry(0, div, grid_size);
 setBoundry(0, p, grid_size);
 linearSolver(0, p, div, 1, 4, dt, grid_size);

 for (i = 1; i <= grid_size; i++)
   {
   for (j = 1; j <= grid_size; j++)
     {
     x[build_index(i, j, grid_size)] -= 0.5f * grid_size * (p[build_index(i+1, j, grid_size)] - p[build_index(i-1, j, grid_size)]);
     y[build_index(i, j, grid_size)] -= 0.5f * grid_size * (p[build_index(i, j+1, grid_size)] - p[build_index(i, j-1, grid_size)]);
     }
   }
  setBoundry(1, x, grid_size);
  setBoundry(2, y, grid_size);
  }

/*
 * The basic density solving routine.
 */
void c_densitySolver(float *d, float *dOld, float diff, float *u, float *v , float dt, int grid_size, int vector_size)
  {
  int i;
  // add density inputted by mouse
  addSource(d, dOld, vector_size, dt);
  swap(d, dOld, vector_size);

  diffuse(0, d, dOld, diff, dt, grid_size);
  swap(d, dOld, vector_size);

  advect(0, d, dOld, u, v, grid_size, dt);
  // clear input density array for next frame
  for (i = 0; i < vector_size; i++) 
    {
    dOld[i] = 0;
    }
  }

/*
 * The basic velocity solving routine as described by Stam.
 */
void c_velocitySolver( float *u, float *uOld, float *v, float *vOld, float *curl, float *d, float visc, float dt, int grid_size, int vector_size)
  {
  int i;

  // add velocity that was input by mouse
  addSource(u, uOld, vector_size, dt);
  addSource(v, vOld, vector_size, dt);

  // add in vorticity confinement force
  vorticityConfinement(uOld, vOld, curl, u, v, grid_size);

  addSource(u, uOld, vector_size, dt);
  addSource(v, vOld, vector_size, dt);

  buoyancy(vOld, d, grid_size);
  addSource(v, vOld, vector_size, dt);

  // swapping arrays for economical mem use
  // and calculating diffusion in velocity.
  swap(u, uOld, vector_size);

  diffuse(0, u, uOld, visc, dt, grid_size);

  swap(v, vOld, vector_size);
  diffuse(0, v, vOld, visc, dt, grid_size);

  // we create an incompressible field
  // for more effective advection.
  project(u, v, uOld, vOld, dt, grid_size);

  swap(u, uOld, vector_size);
  swap(v, vOld, vector_size);

  // self advect velocities
  advect(1, u, uOld, uOld, vOld, grid_size, dt);
  advect(2, v, vOld, uOld, vOld, grid_size, dt);

  // make an incompressible field
  project(u, v, uOld, vOld, dt, grid_size);
  // clear all input velocities for next frame
  for (i = 0; i < vector_size; i++)
      { 
      uOld[i] = 0; 
      vOld[i] = 0; 
      }
  }



