#ifndef _EllipsoidFit_H
#define _EllipsoidFit_H

#include "math.h"
#include "arm_math.h"
#include <stdio.h>
#include <stdlib.h>


void * ellipsoid_fit(float32_t *rawData, int nSample);
void refine_3D_fit(double *gain, double *rotM);

void jacobi_eigenvalue ( int n, double a[], int it_max, double v[],
  double d[], int *it_num, int *rot_num );
void r8mat_diag_get_vector ( int n, double a[], double v[] );
void r8mat_identity ( int n, double a[] );




#endif
