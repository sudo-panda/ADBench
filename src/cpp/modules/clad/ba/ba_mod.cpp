// The lines in ba_clad.cpp are produced using this file. 
//
// This file can be compiled to generate the code by using:
// 
//   clang -v -Xclang -add-plugin -Xclang clad -Xclang -load -Xclang " path to clad library " 
//   -I" path to clad include directory " -x c++ -std=c++11 -lstdc++ -lm -Xclang -plugin-arg-clad
//   -Xclang -fgenerate-source-file


#include "clad/Differentiator//Differentiator.h"


double sqsum(int n, double const* x)
{
  int i;
  double res = 0;
  for (i = 0; i < n; i++)
  {
    res = res + x[i] * x[i];
  }

  return res;
}

double compute_reproj_error0(
    double const* cam,
    double const* X,
    double const w,
    double const* feat
)
{
  double proj[2];

  double Xo[3], Xcam[3];

  Xo[0] = X[0] - cam[3];
  Xo[1] = X[1] - cam[4];
  Xo[2] = X[2] - cam[5];

  int i;
  double sqtheta = sqsum(3, cam);
  if (sqtheta != 0)
  {
    double theta, costheta, sintheta, theta_inverse, tmp;
    double w_local[3], w_cross_pt[3];

    theta = sqrt(sqtheta);
    costheta = cos(theta);
    sintheta = sin(theta);
    theta_inverse = 1.0 / theta;

    for (i = 0; i < 3; i++)
    {
      w_local[i] = cam[i] * theta_inverse;
    }

    w_cross_pt[0] = w_local[1] * Xo[2] - w_local[2] * Xo[1];
    w_cross_pt[1] = w_local[2] * Xo[0] - w_local[0] * Xo[2];
    w_cross_pt[2] = w_local[0] * Xo[1] - w_local[1] * Xo[0];

    tmp = (w_local[0] * Xo[0] + w_local[1] * Xo[1] + w_local[2] * Xo[2]) *
          (1. - costheta);

    for (i = 0; i < 3; i++)
    {
      Xcam[i] = Xo[i] * costheta + w_cross_pt[i] * sintheta + w_local[i] * tmp;
    }
  }
  else
  {
    double rot_cross_pt[3];
    rot_cross_pt[0] = cam[1] * Xo[2] - cam[2] * Xo[1];
    rot_cross_pt[1] = cam[2] * Xo[0] - cam[0] * Xo[2];
    rot_cross_pt[2] = cam[0] * Xo[1] - cam[1] * Xo[0];

    for (i = 0; i < 3; i++)
    {
      Xcam[i] = Xo[i] + rot_cross_pt[i];
    }
  }

  proj[0] = Xcam[0] / Xcam[2];
  proj[1] = Xcam[1] / Xcam[2];

  double rsq, L;
  rsq = sqsum(2, proj);
  L = 1. + cam[9] * rsq + cam[10] * rsq * rsq;
  proj[0] = proj[0] * L;

  proj[0] = proj[0] * cam[6] + cam[7];

  return w * (proj[0] - feat[0]);
}

double compute_reproj_error1(
    double const* cam,
    double const* X,
    double const w,
    double const* feat
)
{
  double proj[2];

  double Xo[3], Xcam[3];

  Xo[0] = X[0] - cam[3];
  Xo[1] = X[1] - cam[4];
  Xo[2] = X[2] - cam[5];

  int i;
  double sqtheta = sqsum(3, cam);
  if (sqtheta != 0)
  {
    double theta, costheta, sintheta, theta_inverse, tmp;
    double w_local[3], w_cross_pt[3];

    theta = sqrt(sqtheta);
    costheta = cos(theta);
    sintheta = sin(theta);
    theta_inverse = 1.0 / theta;

    for (i = 0; i < 3; i++)
    {
      w_local[i] = cam[i] * theta_inverse;
    }

    w_cross_pt[0] = w_local[1] * Xo[2] - w_local[2] * Xo[1];
    w_cross_pt[1] = w_local[2] * Xo[0] - w_local[0] * Xo[2];
    w_cross_pt[2] = w_local[0] * Xo[1] - w_local[1] * Xo[0];

    tmp = (w_local[0] * Xo[0] + w_local[1] * Xo[1] + w_local[2] * Xo[2]) *
          (1. - costheta);

    for (i = 0; i < 3; i++)
    {
      Xcam[i] = Xo[i] * costheta + w_cross_pt[i] * sintheta + w_local[i] * tmp;
    }
  }
  else
  {
    double rot_cross_pt[3];
    rot_cross_pt[0] = cam[1] * Xo[2] - cam[2] * Xo[1];
    rot_cross_pt[1] = cam[2] * Xo[0] - cam[0] * Xo[2];
    rot_cross_pt[2] = cam[0] * Xo[1] - cam[1] * Xo[0];

    for (i = 0; i < 3; i++)
    {
      Xcam[i] = Xo[i] + rot_cross_pt[i];
    }
  }

  proj[0] = Xcam[0] / Xcam[2];
  proj[1] = Xcam[1] / Xcam[2];

  double rsq, L;
  rsq = sqsum(2, proj);
  L = 1. + cam[9] * rsq + cam[10] * rsq * rsq;
  proj[1] = proj[1] * L;

  proj[1] = proj[1] * cam[6] + cam[8];

  return w * (proj[1] - feat[1]);
}

double compute_zach_weight_error(double const w)
{
  return 1 - w * w;
}

int main() {
  clad::gradient(compute_reproj_error0);
  clad::gradient(compute_reproj_error1);
  clad::gradient(compute_zach_weight_error);
}