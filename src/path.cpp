////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    OpenParEM2D - A fullwave 2D electromagnetic simulator.                  //
//    Copyright (C) 2022 Brian Young                                          //
//                                                                            //
//    This program is free software: you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation, either version 3 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "path.hpp"

// with common point (x0,y0)
double angle_between_two_lines (double x0, double y0, double x1, double y1, double x2, double y2)
{
   // cross product
   double cpz=(x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);

   // dot product
   double dp=(x1-x0)*(x2-x0)+(y1-y0)*(y2-y0);

   double theta=atan2(cpz,dp);

   while (theta > M_PI) theta-=2*M_PI;
   while (theta < -M_PI) theta+=2*M_PI;

   return theta;
}

// with common point (x0,y0,z0)
double angle_between_two_lines (double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2)
{
   // cross product
   double cpx=(y1-y0)*(z2-z0)-(y2-y0)*(z1-z0);
   double cpy=(x2-x0)*(z1-z0)-(x1-x0)*(z2-z0);
   double cpz=(x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);

   // dot product
   double dp=(x1-x0)*(x2-x0)+(y1-y0)*(y2-y0)+(z1-z0)*(z2-z0);

   double theta=atan2(sqrt(cpx*cpx+cpy*cpy+cpz*cpz),dp);

   while (theta > M_PI) theta-=2*M_PI;
   while (theta < -M_PI) theta+=2*M_PI;

   return theta;
}

bool are_parallel (double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double tolerance)
{
   double theta=angle_between_two_lines(0,0,x1-x0,y1-y0,x3-x2,y3-y2);
   if (abs(theta) < tolerance) return true;
   return false;
}

bool are_parallel (double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double tolerance)
{
   double theta=angle_between_two_lines(0,0,0,x1-x0,y1-y0,z1-z0,x3-x2,y3-y2,z3-z2);
   if (abs(theta) < tolerance) return true;
   return false;
}

bool compare_xy (double x1, double y1, double x2, double y2, double tolerance)
{
   if (double_compare(x1,x2,tolerance) && double_compare(y1,y2,tolerance)) return true;
   return false;
}

bool compare_xy (double *x1, double *y1, double *x2, double *y2, double *tolerance)
{
   if (double_compare(*x1,*x2,*tolerance) && double_compare(*y1,*y2,*tolerance)) return true;
   return false;
}

bool compare_xyz (double x1, double y1, double z1, double x2, double y2, double z2, double tolerance)
{
   if (double_compare(x1,x2,tolerance) && double_compare(y1,y2,tolerance) && double_compare(z1,z2,tolerance)) return true;
   return false;
}

bool compare_xyz (double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, double *tolerance)
{
   if (double_compare(*x1,*x2,*tolerance) && double_compare(*y1,*y2,*tolerance) && double_compare(*z1,*z2,*tolerance)) return true;
   return false;
}

// checks to see if test point (xt,yt) falls on the line given by (x1,y1) to (x2,y2)
bool is_point_on_line (double xt, double yt, double x1, double y1, double x2, double y2, double tolerance)
{
   double length=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
   double lengtht=sqrt((xt-x1)*(xt-x1)+(yt-y1)*(yt-y1));

   // check ends
   if (compare_xy(xt,yt,x1,y1,tolerance)) return true;
   if (compare_xy(xt,yt,x2,y2,tolerance)) return true;

   // shift to common origin and find the angle between the vectors
   double theta=angle_between_two_lines(0,0,x2-x1,y2-y1,xt-x1,yt-y1);

   // must have the same angle
   if (abs(theta) > tolerance) return false;

   // projection must be small
   if (abs(sin(theta)*lengtht) > tolerance*length) return false;

   // test vector cannot be longer than the line
   if (lengtht > length+tolerance*length) return false;

   return true;
}

// true if t1 <= t <= t2
bool is_bound_by (double t, double t1, double t2, double tolerance)
{
   if (t2 > t1) {
      if (t >= t1-tolerance && t <= t2+tolerance) return true;
   } else {
      if (t >= t2-tolerance && t <= t1+tolerance) return true;
   }
   return false;
}

// checks to see if test point (xt,yt,zt) falls on the line given by (x1,y1,z1) to (x2,y2,z2)
bool is_point_on_line (double xt, double yt, double zt, double x1, double y1, double z1, double x2, double y2, double z2, double tolerance)
{
   double length=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
   double lengtht=sqrt((xt-x1)*(xt-x1)+(yt-y1)*(yt-y1)+(zt-z1)*(zt-z1));

   // check ends
   if (compare_xyz(xt,yt,zt,x1,y1,z1,tolerance)) return true;
   if (compare_xyz(xt,yt,zt,x2,y2,z2,tolerance)) return true;

   // shift to common origin and find the angle between the vectors
   double theta=angle_between_two_lines(0,0,0,x2-x1,y2-y1,z2-z1,xt-x1,yt-y1,zt-z1);

   // must have the same angle
   if (abs(theta) > tolerance) return false;

   // projection must be small
   if (abs(sin(theta)*lengtht) > tolerance*length) return false;

   // test vector cannot be longer than the line
   if (lengtht > length+tolerance*length) return false;

   return true;
}

bool is_point_on_line_not_ends (double xt, double yt, double x1, double y1, double x2, double y2, double tolerance)
{
   if (compare_xy(xt,yt,x1,y1,tolerance)) return false;
   if (compare_xy(xt,yt,x2,y2,tolerance)) return false;
   if (is_point_on_line (xt,yt,x1,y1,x2,y2,tolerance)) return true;
   return false;
}

bool is_point_on_line_not_ends (double xt, double yt, double zt, double x1, double y1, double z1, double x2, double y2, double z2, double tolerance)
{
   if (compare_xyz(xt,yt,zt,x1,y1,z1,tolerance)) return false;
   if (compare_xyz(xt,yt,zt,x2,y2,z2,tolerance)) return false;
   if (is_point_on_line (xt,yt,zt,x1,y1,z1,x2,y2,z2,tolerance)) return true;
   return false;
}

// check if two lines intersect, not counting the end points
// overlaid lines do not count as intersecting, either
bool do_intersect (double x1, double y1, double x2, double y2, double xt1, double yt1, double xt2, double yt2, double tolerance)
{
   // check for identical lines
   if (compare_xy(x1,y1,xt1,yt1,tolerance) && compare_xy(x2,y2,xt2,yt2,tolerance)) return false;
   if (compare_xy(x1,y1,xt2,yt2,tolerance) && compare_xy(x2,y2,xt1,yt1,tolerance)) return false;

   // check for the lines being far apart
   if (max(x1,x2) < min(xt1,xt2)+tolerance) return false;
   if (min(x1,x2) > max(xt1,xt2)-tolerance) return false;
   if (max(y1,y2) < min(yt1,yt2)+tolerance) return false;
   if (min(y1,y2) > max(yt1,yt2)-tolerance) return false;

   // general calculation with traps for infinite slopes
   if (x1 == x2) {
      if (xt1 == xt2) {
         // both vertical - can only overlap
         return false;
      } else {
         // bounding box tests passed, so there must be an intersection, but touching at an end does not count
         if (double_compare(x1,xt1,tolerance)) return false;
         if (double_compare(x1,xt2,tolerance)) return false;
         return true;
      }
   } else {
      if (xt1 == xt2) {
         // bounding box tests passed, so there must be an intersection, but touching at an end does not count
         if (double_compare(xt1,x1,tolerance)) return false;
         if (double_compare(xt1,x2,tolerance)) return false;
         return true;
      } else {
         // the general case

         double m=(y2-y1)/(x2-x1);
         double b=y1-m*x1;

         double mt=(yt2-yt1)/(xt2-xt1);
         double bt=yt1-mt*xt1;

         // check for parallel lines (can't intersect); overlapping parallel lines do not count
         if (double_compare(m,mt,tolerance)) return false;

         // find the intersecting x value by setting y values equal
         double xint=(bt-b)/(m-mt);

         // check intersections at ends, which do not count
         if (double_compare(xint,x1,tolerance)) return false;
         if (double_compare(xint,x2,tolerance)) return false;
         if (double_compare(xint,xt1,tolerance)) return false;
         if (double_compare(xint,xt2,tolerance)) return false;

         // check that x falls within each line for there to be an intersection
         if (xint > x1 && xint < x2 && xint > xt1 && xint < xt2) return true;
         if (xint > x2 && xint < x1 && xint > xt1 && xint < xt2) return true;
         if (xint > x1 && xint < x2 && xint > xt2 && xint < xt1) return true;
         if (xint > x2 && xint < x1 && xint > xt2 && xint < xt2) return true;
      }
   }

   return false;
}

void test_is_point_on_line ()
{
  int i=1;
  double x1,y1,z1,x2,y2,z2,tolerance;

  // 2D

  x1=1; y1=1; x2=10; y2=10; tolerance=1e-8/9;
  if (is_point_on_line (2, 2, x1, y1, x2, y2, tolerance))               PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (1-1e-10, 1-1e-9, x1, y1, x2, y2, tolerance))    PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (10+1e-10, 10+1e-9, x1, y1, x2, y2, tolerance))  PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, 2+1e-10, x1, y1, x2, y2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, 2-1e-10, x1, y1, x2, y2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, -2, x1, y1, x2, y2, tolerance))            PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (11, 11, x1, y1, x2, y2, tolerance))            PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 3, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (3, 2, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2+1e-7, 2, x1, y1, x2, y2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2-1e-7, 2, x1, y1, x2, y2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=1; y1=-1; x2=10; y2=-10;
  if (is_point_on_line (2, -2, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (1-1e-10, -1-1e-9, x1, y1, x2, y2, tolerance))   PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (10+1e-10, -10+1e-9, x1, y1, x2, y2, tolerance)) PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, -2+1e-10, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, -2-1e-10, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 2, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (11, -11, x1, y1, x2, y2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, -3, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (3, -2, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2+1e-7, -2, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2-1e-7, -2, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=-1; y1=1; x2=-10; y2=10;
  if (is_point_on_line (-2, 2, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-1-1e-10, 1-1e-9, x1, y1, x2, y2, tolerance))   PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-10+1e-10, 10+1e-9, x1, y1, x2, y2, tolerance)) PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, 2+1e-10, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, 2-1e-10, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, -2, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-11, 11, x1, y1, x2, y2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 3, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-3, 2, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2+1e-7, 2, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2-1e-7, 2, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=-1; y1=-1; x2=-10; y2=-10;
  if (is_point_on_line (-2, -2, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-1-1e-10, -1-1e-9, x1, y1, x2, y2, tolerance))  PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-10+1e-10, -10+1e-9, x1, y1, x2, y2, tolerance))PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, -2+1e-10, x1, y1, x2, y2, tolerance))       PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, -2-1e-10, x1, y1, x2, y2, tolerance))       PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 2, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-11, -11, x1, y1, x2, y2, tolerance))          PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, -3, x1, y1, x2, y2, tolerance))            PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-3, -2, x1, y1, x2, y2, tolerance))            PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2+1e-7, -2, x1, y1, x2, y2, tolerance))       PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2-1e-7, -2, x1, y1, x2, y2, tolerance))       PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=-1; y1=0; x2=-10; y2=0;
  if (is_point_on_line (-2, 0, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-1+1e-10, 0, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-10-1e-10, 0, x1, y1, x2, y2, tolerance))       PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, 1e-10, x1, y1, x2, y2, tolerance))          PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, -1e-10, x1, y1, x2, y2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 2, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-11, -11, x1, y1, x2, y2, tolerance))          PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, -3, x1, y1, x2, y2, tolerance))            PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 3, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 1e-7, x1, y1, x2, y2, tolerance))          PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, -1e7, x1, y1, x2, y2, tolerance))          PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=1; y1=0; x2=10; y2=0;
  if (is_point_on_line (2, 0, x1, y1, x2, y2, tolerance))               PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (1-1e-10, 0, x1, y1, x2, y2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (10+1e-10, 0, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, 1e-10, x1, y1, x2, y2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, -1e-10, x1, y1, x2, y2, tolerance))          PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 2, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (11, 0, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, -3, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 3, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 1e-7, x1, y1, x2, y2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, -1e7, x1, y1, x2, y2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=0; y1=1; x2=0; y2=10;
  if (is_point_on_line (0, 2, x1, y1, x2, y2, tolerance))               PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (0, 1-1e-10, x1, y1, x2, y2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (0, 10+1e-10, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (1e-10, 2, x1, y1, x2, y2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-1e-10, 2, x1, y1, x2, y2, tolerance))          PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 2, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (0, 11, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 2, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (3, 2, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (1e-7, 2, x1, y1, x2, y2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-1e7, 2, x1, y1, x2, y2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  // 3D

  x1=1; y1=1; z1=1; x2=10; y2=10; z2=10; 
  if (is_point_on_line (2, 2, 2, x1, y1, z1, x2, y2, z2, tolerance))                       PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (1-1e-10, 1-1e-9, 1-1e-10, x1, y1, z1, x2, y2, z2, tolerance))      PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (10+1e-10, 10+1e-9, 10+1e-9, x1, y1, z1, x2, y2, z2, tolerance))    PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, 2+1e-10, 2, x1, y1, z1, x2, y2, z2, tolerance))                 PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, 2-1e-10, 2, x1, y1, z1, x2, y2, z2, tolerance))                 PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, -2, -2, x1, y1, z1, x2, y2, z2, tolerance))                   PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (11, 11, 11, x1, y1, z1, x2, y2, z2, tolerance))                   PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 3, 2, x1, y1, z1, x2, y2, z2, tolerance))                      PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (3, 2, 2, x1, y1, z1, x2, y2, z2, tolerance))                      PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2+1e-7, 2, 2, x1, y1, z1, x2, y2, z2, tolerance))                 PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2-1e-7, 2, 2, x1, y1, z1, x2, y2, z2, tolerance))                 PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=1; y1=-1; z1=1; x2=10; y2=-10; z1=1;
  if (is_point_on_line (2, -2, 2, x1, y1, z1, x2, y2, z2, tolerance))                      PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (1-1e-10, -1-1e-9, 1-1e-10, x1, y1, z1, x2, y2, z2, tolerance))     PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (10+1e-10, -10+1e-9, 10+1e-10, x1, y1, z1, x2, y2, z2, tolerance))  PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, -2+1e-10, 2, x1, y1, z1, x2, y2, z2, tolerance))                PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, -2-1e-10, 2, x1, y1, z1, x2, y2, z2, tolerance))                PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 2, -2, x1, y1, z1, x2, y2, z2, tolerance))                    PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (11, -11, 11, x1, y1, z1, x2, y2, z1, tolerance))                  PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, -3, 2, x1, y1, z1, x2, y2, z2, tolerance))                     PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (3, -2, 3, x1, y1, z1, x2, y2, z2, tolerance))                     PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2+1e-7, -2, 2+1e-7, x1, y1, z1, x2, y2, z2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2-1e-7, -2, 2-1e-7, x1, y1, z1, x2, y2, z2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=-1; y1=1; z1=-1; x2=-10; y2=10; z2=-10;
  if (is_point_on_line (-2, 2, -2, x1, y1, z1, x2, y2, z2, tolerance))                     PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-1-1e-10, 1-1e-9, -1-1e-10, x1, y1, z1, x2, y2, z2, tolerance))    PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-10+1e-10, 10+1e-9, -10+1e-10, x1, y1, z1, x2, y2, z2, tolerance)) PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, 2+1e-10, -2, x1, y1, z1, x2, y2, z2, tolerance))               PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, 2-1e-10, -2, x1, y1, z1, x2, y2, z2, tolerance))               PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, -2, 2, x1, y1, z1, x2, y2, z2, tolerance))                     PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-11, 11, -11, x1, y1, z1, x2, y2, z2, tolerance))                 PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 3, -2, x1, y1, z1, x2, y2, z2, tolerance))                    PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-3, 2, -3, x1, y1, z1, x2, y2, z2, tolerance))                    PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2+1e-7, 2, -2+1e-7, x1, y1, z1, x2, y2, z2, tolerance))          PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2-1e-7, 2, -2-1e-7, x1, y1, z1, x2, y2, z2, tolerance))          PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;


  x1=-1; y1=-1; z1=-1; x2=-10; y2=-10; z2=-10;
  if (is_point_on_line (-2, -2, -2, x1, y1, z1, x2, y2, z2, tolerance))                    PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-1-1e-10, -1-1e-9, -1-1e-10, x1, y1, z1, x2, y2, z2, tolerance))   PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-10+1e-10, -10+1e-9, -10+1e-10, x1, y1, z1, x2, y2, z2, tolerance))PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, -2+1e-10, -2, x1, y1, z1, x2, y2, z2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, -2-1e-10, -2, x1, y1, z1, x2, y2, z2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 2, 2, x1, y1, z1, x2, y2, z2, tolerance))                      PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-11, -11, -11, x1, y1, z1, x2, y2, z2, tolerance))                PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, -3, -2, x1, y1, z1, x2, y2, z2, tolerance))                   PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-3, -2, -3, x1, y1, z1, x2, y2, z2, tolerance))                   PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2+1e-7, -2, -2+1e-7, x1, y1, z1, x2, y2, z2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2-1e-7, -2, -2-1e-7, x1, y1, z1, x2, y2, z2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=-1; y1=0; z1=-1; x2=-10; y2=0; z2=-10;
  if (is_point_on_line (-2, 0, -2, x1, y1, z1, x2, y2, z2, tolerance))                     PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-1+1e-10, 0, -1+1e-10, x1, y1, z1, x2, y2, z2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-10-1e-10, 0, -10-1e-10, x1, y1, z1, x2, y2, z2, tolerance))       PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, 1e-10, -2, x1, y1, z1, x2, y2, z2, tolerance))                 PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, -1e-10, -2, x1, y1, z1, x2, y2, z2, tolerance))                PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 2, -2, x1, y1, z1, x2, y2, z2, tolerance))                    PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-11, -11, -11, x1, y1, z1, x2, y2, z2, tolerance))                PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, -3, -2, x1, y1, z1, x2, y2, z2, tolerance))                   PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 3, -2, x1, y1, z1, x2, y2, z2, tolerance))                    PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 1e-7, -2, x1, y1, z1, x2, y2, z2, tolerance))                 PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, -1e7, -2, x1, y1, z1, x2, y2, z2, tolerance))                 PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=1; y1=0; z1=1; x2=10; y2=0; z2=10;
  if (is_point_on_line (2, 0, 2, x1, y1, z1, x2, y2, z2, tolerance))                       PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (1-1e-10, 0, 1-1e-10, x1, y1, z1, x2, y2, z2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (10+1e-10, 0, 10+1e-10, x1, y1, z1, x2, y2, z2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, 1e-10, 2, x1, y1, z1, x2, y2, z2, tolerance))                   PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, -1e-10, 2, x1, y1, z1, x2, y2, z2, tolerance))                  PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 2, 2, x1, y1, z1, x2, y2, z2, tolerance))                      PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (11, 0, 11, x1, y1, z1, x2, y2, z2, tolerance))                    PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, -3, 2, x1, y1, z1, x2, y2, z2, tolerance))                     PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 3, 2, x1, y1, z1, x2, y2, z2, tolerance))                      PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 1e-7, 2, x1, y1, z1, x2, y2, z2, tolerance))                   PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, -1e7, 2, x1, y1, z1, x2, y2, z2, tolerance))                   PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=0; y1=1; z1=0; x2=0; y2=10; z2=0;
  if (is_point_on_line (0, 2, 0, x1, y1, z1, x2, y2, z2, tolerance))                       PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (0, 1-1e-10, 0, x1, y1, z1, x2, y2, z2, tolerance))                 PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (0, 10+1e-10, 0, x1, y1, z1, x2, y2, z2, tolerance))                PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (1e-10, 2, 1e-10, x1, y1, z1, x2, y2, z2, tolerance))               PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-1e-10, 2, -1e-10, x1, y1, z1, x2, y2, z2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 2, 2, x1, y1, z1, x2, y2, z2, tolerance))                      PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (0, 11, 0, x1, y1, z1, x2, y2, z2, tolerance))                     PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 2, -2, x1, y1, z1, x2, y2, z2, tolerance))                    PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (3, 2, 3, x1, y1, z1, x2, y2, z2, tolerance))                      PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (1e-7, 2, 1e-7, x1, y1, z1, x2, y2, z2, tolerance))                PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-1e7, 2, -1e-7, x1, y1, z1, x2, y2, z2, tolerance))               PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
}

bool mergePaths (vector<Path *> *pathList, vector<long unsigned int> *pathIndexList, vector<bool> *reverseList, string boundaryType, string boundaryName, Path **mergedPath)
{
   bool fail=false;
   long unsigned int i;

   // see if any merging is required
   if (pathIndexList->size() <= 1) return fail;

   // merge required

   // none of the paths can be closed: ambiguous as to user's intent
   i=0;
   while (i < pathIndexList->size()) {
      Path *path=(*pathList)[(*pathIndexList)[i]];
      if (path->is_closed()) {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR1095: Path must be open with \"closed\" set to \"false\" at line number %d.\n",path->get_closed_lineNumber());
         fail=true;
      }
      i++;
   }
   if (fail) return fail;

   // set up a new path
   *mergedPath=new Path(-1,-1);

   // merge

   i=0;
   while (i < pathIndexList->size()) {

      Path *path=(*pathList)[(*pathIndexList)[i]];

      if ((*reverseList)[i]) {
         long unsigned int j=path->get_points_size()-1;
         while (j >= 0) {
            // always add the point to a new path
            if ((*mergedPath)->get_points_size() == 0) (*mergedPath)->push_point(path->get_point(j)->clone());           
            else {
                // add if the point does not duplicate the last one
                if (! path->get_point(j)->point_compare((*mergedPath)->get_point((*mergedPath)->get_points_size()-1))) {
                    (*mergedPath)->push_point(path->get_point(j)->clone());
                }
            }
            if (j == 0) break;
            j--;
         }
      } else {
         long unsigned int j=0;
         while (j < path->get_points_size()) {
            // always add the point to a new path
            if ((*mergedPath)->get_points_size() == 0) (*mergedPath)->push_point(path->get_point(j)->clone());
            else {
                // add if the point does not duplicate the last one
                if (! path->get_point(j)->point_compare((*mergedPath)->get_point((*mergedPath)->get_points_size()-1))) {
                    (*mergedPath)->push_point(path->get_point(j)->clone());
                }
            }
            j++;
         }
      }

      i++;
   }
   // check for closed polygon
   if ((*mergedPath)->get_point(0)->point_compare((*mergedPath)->get_point((*mergedPath)->get_points_size()-1))) {
      (*mergedPath)->pop_point();
   }

   (*mergedPath)->set_closed(true);

   // check for duplicate points along the path
   i=0;
   while (i < (*mergedPath)->get_points_size()-1) {
      long unsigned int j=i+1;
      while (j < (*mergedPath)->get_points_size()) {
         if ((*mergedPath)->get_point(i)->point_compare((*mergedPath)->get_point(j))) {
            PetscPrintf(PETSC_COMM_WORLD,"ERROR1096: Merged path for %s %s has duplicate point at (%g,%g,%g)\n",
               boundaryType.c_str(),boundaryName.c_str(),(*mergedPath)->get_point(i)->get_point_value_x(),(*mergedPath)->get_point(i)->get_point_value_y(),(*mergedPath)->get_point(i)->get_point_value_z());
            fail=true;
         }
         j++;
      }
      i++;
   }

   // ToDo - Add a check for crossing lines

   if (fail) {
      delete *mergedPath;
      *mergedPath=nullptr;
      return fail;
   }

   (*mergedPath)->calculateBoundingBox();

   return fail;
}


///////////////////////////////////////////////////////////////////////////////////////////
// Path
///////////////////////////////////////////////////////////////////////////////////////////

Path::Path(int startLine_, int endLine_)
{
   startLine=startLine_;
   endLine=endLine_;

   // name
   name.push_alias("name");
   name.set_loaded(false);
   name.set_positive_required(false);
   name.set_non_negative_required(false);
   name.set_lowerLimit(0);
   name.set_upperLimit(0);
   name.set_checkLimits(false);

   // closed
   closed.push_alias("closed");
   closed.set_loaded(false);
   closed.set_positive_required(false);
   closed.set_non_negative_required(false);
   closed.set_lowerLimit(0);
   closed.set_upperLimit(0);
   closed.set_checkLimits(false);

   rotated=false;

   theta=0;
   cos_theta_=cos(theta);
   sin_theta_=sin(theta);

   phi=0;
   cos_phi_=cos(phi);
   sin_phi_=sin(phi);

   hasNormal=false;
   nx=-1; ny=-1; nz=-1;
}


// return true if the point is close
bool Path::compare (long unsigned int i, keywordPair test_point)
{
   if (points[i]->get_point_value_x() == 0) {
      if (fabs(test_point.get_point_value_x()) > tol) return false;
   }
   if (fabs((test_point.get_point_value_x()-points[i]->get_point_value_x())/points[i]->get_point_value_x()) > tol) return false;

   if (points[i]->get_point_value_y() == 0) {
      if (fabs(test_point.get_point_value_y()) > tol) return false;
   }
   if (fabs((test_point.get_point_value_y()-points[i]->get_point_value_y())/points[i]->get_point_value_y()) > tol) return false;
   return true;
}

void Path::print (string indent)
{
   int dim=0;
   PetscPrintf(PETSC_COMM_WORLD,"%sPath %p\n",indent.c_str(),this);
   PetscPrintf(PETSC_COMM_WORLD,"%s   name=%s\n",indent.c_str(),get_name().c_str());
   long unsigned int i=0;
   while (i < points.size()) {
      if (get_point_dim(i) == 2) {dim=2; PetscPrintf(PETSC_COMM_WORLD,"%s   point=(%g,%g)\n",indent.c_str(),get_point_x(i),get_point_y(i));}
      if (get_point_dim(i) == 3) {dim=3; PetscPrintf(PETSC_COMM_WORLD,"%s   point=(%g,%g,%g)\n",indent.c_str(),get_point_x(i),get_point_y(i),get_point_z(i));}
      i++;
   }
   if (closed.get_bool_value()) PetscPrintf(PETSC_COMM_WORLD,"%s   closed=true\n",indent.c_str());
   else PetscPrintf(PETSC_COMM_WORLD,"%s   closed=false\n",indent.c_str());
   PetscPrintf(PETSC_COMM_WORLD,"%s   tol=%g\n",indent.c_str(),tol);
   if (rotated) PetscPrintf(PETSC_COMM_WORLD,"%s   rotated=true\n",indent.c_str());
   else PetscPrintf(PETSC_COMM_WORLD,"%s   rotated=false\n",indent.c_str());
   PetscPrintf(PETSC_COMM_WORLD,"%s   theta=%g\n",indent.c_str(),theta);
   PetscPrintf(PETSC_COMM_WORLD,"%s   phi=%g\n",indent.c_str(),phi);
   PetscPrintf(PETSC_COMM_WORLD,"%s   sin_theta_=%g\n",indent.c_str(),sin_theta_);
   PetscPrintf(PETSC_COMM_WORLD,"%s   cos_theta_=%g\n",indent.c_str(),cos_theta_);
   PetscPrintf(PETSC_COMM_WORLD,"%s   sin_phi_=%g\n",indent.c_str(),sin_phi_);
   PetscPrintf(PETSC_COMM_WORLD,"%s   cos_phi_=%g\n",indent.c_str(),cos_phi_);
   PetscPrintf(PETSC_COMM_WORLD,"%s   xmax=%g\n",indent.c_str(),xmax);
   PetscPrintf(PETSC_COMM_WORLD,"%s   xmin=%g\n",indent.c_str(),xmin);
   PetscPrintf(PETSC_COMM_WORLD,"%s   ymax=%g\n",indent.c_str(),ymax);
   PetscPrintf(PETSC_COMM_WORLD,"%s   ymin=%g\n",indent.c_str(),ymin);
   if (dim == 3) PetscPrintf(PETSC_COMM_WORLD,"%s   zmax=%g\n",indent.c_str(),zmax);
   if (dim == 3) PetscPrintf(PETSC_COMM_WORLD,"%s   zmin=%g\n",indent.c_str(),zmin);
   if (hasNormal) PetscPrintf(PETSC_COMM_WORLD,"%s   hasNormal=true\n",indent.c_str());
   else PetscPrintf(PETSC_COMM_WORLD,"%s   hasNormal=false\n",indent.c_str());
   PetscPrintf(PETSC_COMM_WORLD,"%s   nx=%g\n",indent.c_str(),nx);
   PetscPrintf(PETSC_COMM_WORLD,"%s   ny=%g\n",indent.c_str(),ny);
   PetscPrintf(PETSC_COMM_WORLD,"%s   nz=%g\n",indent.c_str(),nz);
   PetscPrintf(PETSC_COMM_WORLD,"%sEndPath\n",indent.c_str());
}

void Path::output (ofstream *out, int force_dim)
{
   *out << "Path" << endl;
   *out << "   name=" << get_name() << endl;
   long unsigned int i=0;
   while (i < points.size()) {
      if (force_dim == 2) {*out << setprecision(16) << "   point=(" << get_point_x(i) << "," << get_point_y(i) << ")" << endl;}
      if (force_dim == 3) {*out << setprecision(16) << "   point=(" << get_point_x(i) << "," << get_point_y(i) << "," << get_point_z(i) << ")" << endl;}
      i++;
   }
   if (closed.get_bool_value()) *out << "   closed=true" << endl;
   else *out << "   closed=false" << endl;

   *out << "EndPath" << endl;
}

bool Path::load(int dim, string *indent, inputFile *inputs)
{
   bool fail=false;

   int lineNumber=inputs->get_next_lineNumber(startLine);
   int stopLineNumber=inputs->get_previous_lineNumber(endLine);
   while (lineNumber <= stopLineNumber) {

      string token,value,line;
      line=inputs->get_line(lineNumber);
      get_token_pair(&line,&token,&value,&lineNumber,*indent);

      int recognized=0;

      if (name.match_alias(&token)) {
         if (name.is_loaded()) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1097: Duplicate entry at line %d for previous entry at line %d.\n",
                                         indent->c_str(),indent->c_str(),lineNumber,name.get_lineNumber());
            fail=true;
         } else {
            name.set_keyword(token);
            name.set_value(value);
            name.set_lineNumber(lineNumber);
            name.set_loaded(true);
         }
         recognized++;
      }

      if (token.compare("point") == 0) {
         keywordPair *point=new keywordPair;
         point->push_alias("point");
         point->set_keyword(token);
         point->set_value(value);
         point->set_lineNumber(lineNumber);
         point->set_positive_required(false);
         point->set_non_negative_required(false);
         point->set_lowerLimit(-100);
         point->set_upperLimit(100);
         point->set_checkLimits(true);
         point->set_loaded(false);
         point->set_point_value(-DBL_MAX,-DBL_MAX,-DBL_MAX);

         if (point->loadPoint(dim,&token,&value,lineNumber)) delete point;
         else points.push_back(point);

         recognized++;
      }

      if (closed.match_alias(&token)) {
         recognized++;
         if (closed.loadBool(&token, &value, lineNumber)) fail=true;
      }

      // should recognize one keyword
      if (recognized != 1) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1098: Unrecognized keyword at line %d.\n",indent->c_str(),indent->c_str(),lineNumber);
         fail=true;
      }
      lineNumber=inputs->get_next_lineNumber(lineNumber);
   }

   calculateBoundingBox();
   calculateNormal();

   return fail;
}

bool Path::checkBoundingBox(Vector *lowerLeft, Vector *upperRight, string *indent, double tol)
{
   bool fail=false;

   long unsigned int i=0;
   while (i < points.size()) {

      bool point_fail=false;

      if (points[i]->get_point_value_x() < lowerLeft->Elem(0)-tol) point_fail=true;
      if (points[i]->get_point_value_x() > upperRight->Elem(0)+tol) point_fail=true;

      if (points[i]->get_point_value_y() < lowerLeft->Elem(1)-tol) point_fail=true;
      if (points[i]->get_point_value_y() > upperRight->Elem(1)+tol) point_fail=true;

      if (points[i]->get_point_value_dim() == 3) {
         if (points[i]->get_point_value_z() < lowerLeft->Elem(2)-tol) point_fail=true;
         if (points[i]->get_point_value_z() > upperRight->Elem(2)+tol) point_fail=true;
      }

      if (point_fail) {
         if (points[i]->get_point_value_dim() == 2) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1099: Path block at line %d has point (%g,%g) outside of the mesh bounding box.\n",
                                         indent->c_str(),indent->c_str(),startLine,points[i]->get_point_value_x(),points[i]->get_point_value_y());
         }
         if (points[i]->get_point_value_dim() == 3) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1100: Path block at line %d has point (%g,%g,%g) outside of the mesh bounding box.\n",
                                         indent->c_str(),indent->c_str(),startLine,
                                         points[i]->get_point_value_x(),points[i]->get_point_value_y(),points[i]->get_point_value_z());
         }
         fail=true;
      }

      i++;
   }

   return fail;
}

bool Path::inBlock (int lineNumber)
{
   if (lineNumber >= startLine && lineNumber <= endLine) return true;
   return false;
}

bool Path::check(string *indent)
{
   bool fail=false;

   if (!name.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1101: Path block at line %d must specify a name.\n",indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   if (!closed.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1102: Path block at line %d must specify \"closed\".\n",indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   if (points.size() == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1103: Path block at line %d must specify points.\n",indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   if (points.size() == 1) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1104: Path block at line %d must specify more than one point.\n",indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   if (points.size() == 2 && closed.get_bool_value()) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1105: Path block at line %d cannot be closed with just two points.\n",indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   return fail;
}

// returns the segment index on which the segment falls
long unsigned int Path::is_segmentOnLine (double x1, double y1, double x2, double y2)
{
   long unsigned int max=-1;

   long unsigned int i=0;
   while (points.size() > 0 && i < points.size()-1) {
      if (is_point_on_line (x1,y1,points[i]->get_point_value_x(),points[i]->get_point_value_y(),
                                  points[i+1]->get_point_value_x(),points[i+1]->get_point_value_y(),1e-8) &&
          is_point_on_line (x2,y2,points[i]->get_point_value_x(),points[i]->get_point_value_y(),
                                  points[i+1]->get_point_value_x(),points[i+1]->get_point_value_y(),1e-8)) return i;
      i++;
   }

   if (is_closed()) {
      if (is_point_on_line (x1,y1,points[points.size()-1]->get_point_value_x(),points[points.size()-1]->get_point_value_y(),
                                  points[0]->get_point_value_x(),points[0]->get_point_value_y(),1e-8) &&
          is_point_on_line (x2,y2,points[points.size()-1]->get_point_value_x(),points[points.size()-1]->get_point_value_y(),
                                  points[0]->get_point_value_x(),points[0]->get_point_value_y(),1e-8)) return points.size()-1;
   }

   return max;
}

// returns the segment index on which the segment falls
long unsigned int Path::is_segmentOnLine (double x1, double y1, double z1, double x2, double y2, double z2)
{
   long unsigned int max=-1;

   long unsigned int i=0;
   while (points.size() > 0 && i < points.size()-1) {
      if (is_point_on_line (x1,y1,z1,points[i]->get_point_value_x(),points[i]->get_point_value_y(),points[i]->get_point_value_z(),
                                     points[i+1]->get_point_value_x(),points[i+1]->get_point_value_y(),points[i+1]->get_point_value_z(),1e-8) &&
          is_point_on_line (x2,y2,z2,points[i]->get_point_value_x(),points[i]->get_point_value_y(),points[i]->get_point_value_z(),
                                     points[i+1]->get_point_value_x(),points[i+1]->get_point_value_y(),points[i+1]->get_point_value_z(),1e-8)) return i;
      i++;
   }

   if (is_closed()) {
      if (is_point_on_line (x1,y1,z1,points[points.size()-1]->get_point_value_x(),points[points.size()-1]->get_point_value_y(),points[points.size()-1]->get_point_value_z(),
                                     points[0]->get_point_value_x(),points[0]->get_point_value_y(),points[0]->get_point_value_z(),1e-8) &&
          is_point_on_line (x2,y2,z2,points[points.size()-1]->get_point_value_x(),points[points.size()-1]->get_point_value_y(),points[points.size()-1]->get_point_value_z(),
                                     points[0]->get_point_value_x(),points[0]->get_point_value_y(),points[0]->get_point_value_z(),1e-8)) return points.size()-1;
   }

   return max;
}

double Path::sum_of_angles (double x, double y)
{
   double theta=0;
   long unsigned int i=0;
   while (i < points.size()-1) {
      theta+=angle_between_two_lines (x,y,points[i]->get_point_value_x(),points[i]->get_point_value_y(),
                                          points[i+1]->get_point_value_x(),points[i+1]->get_point_value_y());
      i++;
   }
   if (is_closed()) {
      theta+=angle_between_two_lines (x,y,points[i]->get_point_value_x(),points[i]->get_point_value_y(),
                                          points[0]->get_point_value_x(),points[0]->get_point_value_y());
   }
   return theta;
}

// eliminate partial overlaps of paths by subdividing
// crossing paths are ok
void Path::subdivide2D(Path *test)
{
   bool modified=false;
   vector<keywordPair *> newPoints;

   if (points.size() == 0) return;
   if (test->points.size() == 0) return;

   newPoints.push_back(test->points[0]->clone());

   long unsigned int i=0;
   while (i < points.size()-1) {

      double x1=points[i]->get_point_value_x();
      double y1=points[i]->get_point_value_y();
      double x2=points[i+1]->get_point_value_x();
      double y2=points[i+1]->get_point_value_y();

      long unsigned int j=0;
      while (j < test->points.size()-1) {

         double xt1=test->points[j]->get_point_value_x();
         double yt1=test->points[j]->get_point_value_y();
         double xt2=test->points[j+1]->get_point_value_x();
         double yt2=test->points[j+1]->get_point_value_y();

         bool break_on_1=false;
         bool break_on_2=false;

         bool parallel=are_parallel (xt1,yt1,xt2,yt2,x1,y1,x2,y2,1e-12);

         if (parallel && is_point_on_line_not_ends(xt1,yt1,x1,y1,x2,y2,1e-8)) break_on_1=true;
         if (parallel && is_point_on_line_not_ends(xt2,yt2,x1,y1,x2,y2,1e-8)) break_on_2=true;

         if (break_on_1) {
            if (break_on_2) {
               // segment is fully enclosed

               // maintain ordering along the line
               if (points[i]->get_point_distance(test->points[j]) < points[i]->get_point_distance(test->points[j+1])) {
                  if (! points[points.size()-1]->is_close_point (test->points[j])) {
                     newPoints.push_back(test->points[j]->clone());
                  }
                  newPoints.push_back(test->points[j+1]->clone());
               } else {
                  if (! points[points.size()-1]->is_close_point (test->points[j+1])) {
                     newPoints.push_back(test->points[j+1]->clone());
                  }
                  newPoints.push_back(test->points[j]->clone());
               }
               modified=true;
            } else {
               // partial overlap - break the segment at test point 1
               newPoints.push_back(test->points[j]->clone());
               modified=true;
            }
         } else {
            if (break_on_2) {
               // partial overlap - break the segment at test point 2
               newPoints.push_back(test->points[j+1]->clone());
               modified=true;
            } else {
               // nothing to do
            }
         }

         j++;
      }

      // finish the segment
      newPoints.push_back(points[i+1]->clone());

      i++;
   }

   if (modified) {
      long unsigned int k=0;
      while (k < points.size()) {
         delete points[k];
         k++;
      }
      points.clear();

      k=0;
      while (k < newPoints.size()) {
         points.push_back(newPoints[k]);
         k++;
      }
   } else {
      long unsigned int k=0;
      while (k < newPoints.size()) {
         delete newPoints[k];
         k++;
      }
   }
}

// eliminate partial overlaps of paths by subdividing
// crossing paths are ok
void Path::subdivide3D(Path *test)
{
   bool modified=false;
   vector<keywordPair *> newPoints;

   if (points.size() == 0) return;
   if (test->points.size() == 0) return;

   newPoints.push_back(test->points[0]->clone());

   long unsigned int i=0;
   while (i < points.size()-1) {

      double x1=points[i]->get_point_value_x();
      double y1=points[i]->get_point_value_y();
      double z1=points[i]->get_point_value_z();
      double x2=points[i+1]->get_point_value_x();
      double y2=points[i+1]->get_point_value_y();
      double z2=points[i+1]->get_point_value_z();

      long unsigned int j=0;
      while (j < test->points.size()-1) {

         double xt1=test->points[j]->get_point_value_x();
         double yt1=test->points[j]->get_point_value_y();
         double zt1=test->points[j]->get_point_value_z();
         double xt2=test->points[j+1]->get_point_value_x();
         double yt2=test->points[j+1]->get_point_value_y();
         double zt2=test->points[j+1]->get_point_value_z();

         bool break_on_1=false;
         bool break_on_2=false;

         bool parallel=are_parallel (xt1,yt1,zt1,xt2,yt2,zt2,x1,y1,z1,x2,y2,z2,1e-12);

         if (parallel && is_point_on_line_not_ends(xt1,yt1,zt1,x1,y1,z1,x2,y2,z2,1e-8)) break_on_1=true;
         if (parallel && is_point_on_line_not_ends(xt2,yt2,zt2,x1,y1,z1,x2,y2,z2,1e-8)) break_on_2=true;

         if (break_on_1) {
            if (break_on_2) {
               // segment is fully enclosed

               // maintain ordering along the line
               if (points[i]->get_point_distance(test->points[j]) < points[i]->get_point_distance(test->points[j+1])) {
                  if (! points[points.size()-1]->is_close_point (test->points[j])) {
                     newPoints.push_back(test->points[j]->clone());
                  }
                  newPoints.push_back(test->points[j+1]->clone());
               } else {
                  if (! points[points.size()-1]->is_close_point (test->points[j+1])) {
                     newPoints.push_back(test->points[j+1]->clone());
                  }
                  newPoints.push_back(test->points[j]->clone());
               }
               modified=true;
            } else {
               // partial overlap - break the segment at test point 1
               newPoints.push_back(test->points[j]->clone());
               modified=true;
            }
         } else {
            if (break_on_2) {
               // partial overlap - break the segment at test point 2
               newPoints.push_back(test->points[j+1]->clone());
               modified=true;
            } else {
               // nothing to do
            }
         }

         j++;
      }

      // finish the segment
      newPoints.push_back(points[i+1]->clone());

      i++;
   }

   if (modified) {
      long unsigned int k=0;
      while (k < points.size()) {
         delete points[k];
         k++;
      }
      points.clear();

      k=0;
      while (k < newPoints.size()) {
         points.push_back(newPoints[k]);
         k++;
      }
   } else {
      long unsigned int k=0;
      while (k < newPoints.size()) {
         delete newPoints[k];
         k++;
      }
   }
}

Path* Path::clone()
{
   Path *newPath=new Path(startLine,endLine);
   newPath->name=name;
   newPath->closed=closed;
   newPath->tol=tol;
   newPath->rotated=rotated;
   newPath->theta=theta;
   newPath->phi=phi;
   newPath->sin_theta_=sin_theta_;
   newPath->cos_theta_=cos_theta_;
   newPath->sin_phi_=sin_phi_;
   newPath->cos_phi_=cos_phi_;
   newPath->xmax=xmax;
   newPath->xmin=xmin;
   newPath->ymax=ymax;
   newPath->ymin=ymin;
   newPath->zmax=zmax;
   newPath->zmin=zmin;
   newPath->hasNormal=hasNormal;
   newPath->nx=nx;
   newPath->ny=ny;
   newPath->nz=nz;

   long unsigned int i=0;
   while (i < points.size()) {
      newPath->points.push_back(points[i]->clone());
      i++;
   }
   return newPath;
}

void Path::calculateBoundingBox()
{
   xmax=-DBL_MAX;
   xmin=DBL_MAX;
   ymax=-DBL_MAX;
   ymin=DBL_MAX;
   zmax=-DBL_MAX;
   zmin=DBL_MAX;
   long unsigned int i=0;
   while (i < points.size()) {
      if (points[i]->get_point_value_x() > xmax) xmax=points[i]->get_point_value_x();
      if (points[i]->get_point_value_x() < xmin) xmin=points[i]->get_point_value_x();
      if (points[i]->get_point_value_y() > ymax) ymax=points[i]->get_point_value_y();
      if (points[i]->get_point_value_y() < ymin) ymin=points[i]->get_point_value_y();
      if (points[i]->get_point_value_z() > zmax) zmax=points[i]->get_point_value_z();
      if (points[i]->get_point_value_z() < zmin) zmin=points[i]->get_point_value_z();
      i++;
   }
}

// calculate a normal to the path
// assumes that the path is planar
bool Path::calculateNormal ()
{
   // this is a 3D operation,
   // just need to check one point
   if (points[0]->get_point_value_dim() == 2) {
      return true;
   }

   // exit if 0, 1, or 2 points - cannot determine a normal vector
   if (points.size() < 3) {
      //PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Path::calculateNormal passed path with fewer than 3 points.\n");
      return true;
   }

   // find the two consecutive segments with the angle between them closest to 90 degrees
   double smallest_difference=1e30;
   long unsigned int index=0;
   long unsigned int i=1;
   while (i < points.size()-1) {
      double theta=angle_between_two_lines (points[i]->get_point_value_x(), points[i]->get_point_value_y(), points[i]->get_point_value_z(),
                                            points[i-1]->get_point_value_x(), points[i-1]->get_point_value_y(), points[i-1]->get_point_value_z(),
                                            points[i+1]->get_point_value_x(), points[i+1]->get_point_value_y(), points[i+1]->get_point_value_z());
      if (abs(abs(theta)-M_PI/2) < smallest_difference) {smallest_difference=abs(abs(theta)-M_PI/2); index=i;}
      i++;
   }

   if (is_closed()) {
      i=0;
      double theta=angle_between_two_lines (points[i]->get_point_value_x(), points[i]->get_point_value_y(), points[i]->get_point_value_z(),
                                            points[points.size()-1]->get_point_value_x(), points[points.size()-1]->get_point_value_y(), points[points.size()-1]->get_point_value_z(),
                                            points[i+1]->get_point_value_x(), points[i+1]->get_point_value_y(), points[i+1]->get_point_value_z());
      if (abs(abs(theta)-M_PI/2) < smallest_difference) {smallest_difference=abs(abs(theta)-M_PI/2); index=i;}

      i=points.size()-1;
      theta=angle_between_two_lines (points[i]->get_point_value_x(), points[i]->get_point_value_y(), points[i]->get_point_value_z(),
                                     points[i-1]->get_point_value_x(), points[i-1]->get_point_value_y(), points[i-1]->get_point_value_z(),
                                     points[0]->get_point_value_x(), points[0]->get_point_value_y(), points[0]->get_point_value_z());
      if (abs(abs(theta)-M_PI/2) < smallest_difference) {smallest_difference=abs(abs(theta)-M_PI/2); index=i;}
   }

   // cannot continue if the points are colinear - cannot determine a normal vector
   if (double_compare(abs(theta),M_PI,tol)) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Path::calculateNormal passed a colinear path.\n");
      return true;
   }

   // find the normal to these two segments
   double vx,vy,vz;
   double ux,uy,uz;

   if (index == 0) {
      vx=points[points.size()-1]->get_point_value_x()-points[0]->get_point_value_x();
      vy=points[points.size()-1]->get_point_value_y()-points[0]->get_point_value_y();
      vz=points[points.size()-1]->get_point_value_z()-points[0]->get_point_value_z();
      ux=points[1]->get_point_value_x()-points[0]->get_point_value_x();
      uy=points[1]->get_point_value_y()-points[0]->get_point_value_y();
      uz=points[1]->get_point_value_z()-points[0]->get_point_value_z();
   } else if (index == points.size()-1) {
      vx=points[points.size()-2]->get_point_value_x()-points[points.size()-1]->get_point_value_x();
      vy=points[points.size()-2]->get_point_value_y()-points[points.size()-1]->get_point_value_y();
      vz=points[points.size()-2]->get_point_value_z()-points[points.size()-1]->get_point_value_z();
      ux=points[0]->get_point_value_x()-points[points.size()-1]->get_point_value_x();
      uy=points[0]->get_point_value_y()-points[points.size()-1]->get_point_value_y();
      uz=points[0]->get_point_value_z()-points[points.size()-1]->get_point_value_z();
   } else {
      vx=points[index-1]->get_point_value_x()-points[index]->get_point_value_x();
      vy=points[index-1]->get_point_value_y()-points[index]->get_point_value_y();
      vz=points[index-1]->get_point_value_z()-points[index]->get_point_value_z();
      ux=points[index+1]->get_point_value_x()-points[index]->get_point_value_x();
      uy=points[index+1]->get_point_value_y()-points[index]->get_point_value_y();
      uz=points[index+1]->get_point_value_z()-points[index]->get_point_value_z();
   }

   // perpendicular vector: n=uxv
   nx=uy*vz-vy*uz;
   ny=-ux*vz+vx*uz;
   nz=ux*vy-vx*uy;

   // normalize
   double mag=sqrt(nx*nx+ny*ny+nz*nz);
   nx/=mag;
   ny/=mag;
   nz/=mag;

   hasNormal=true;

   return false;
}

Path* Path::rotateToXYplane ()
{
   if (!hasNormal) {
      return nullptr;
   }

   // see if already rotated
   if (rotated) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Path::rotateToXYplane passed a rotated path.\n");
      return nullptr;
   }

   // rotated structure
   Path *rotated=this->clone();
   rotated->rotated=true;

   // make copies
   double cnx=rotated->nx;
   double cny=rotated->ny;
   double cnz=rotated->nz;

   // rotation angles

   // rotation about the z-axis to put n in the x-z plane
   rotated->theta=-atan2(cny,cnx);
   cnx=sqrt(cnx*cnx+cny*cny);
   cny=0;

   // rotation about the y-axis to put n in the z direction
   rotated->phi=-angle_between_two_lines(0,0,0,cnx,cny,cnz,0,0,1);

   if (rotated->theta > M_PI/2) {
      rotated->theta-=M_PI;
      rotated->phi=-rotated->phi;
   }

   if (rotated->theta < -M_PI/2) {
      rotated->theta+=M_PI;
      rotated->phi=-rotated->phi;
   }

   rotated->cos_theta_=cos(rotated->theta);
   rotated->sin_theta_=sin(rotated->theta);

   rotated->cos_phi_=cos(rotated->phi);
   rotated->sin_phi_=sin(rotated->phi);

   long unsigned int i=0;
   while (i < rotated->points.size()) {

      //  theta - rotation about z-axis
      //  phi - rotation about the y-axis
      //  (px,py,pz) - point to be rotated
      //  (px',py',pz') - rotated point
      //  rotation matrix for z x rotation matrix for y x point:
      //  [px']   [  cos(phi)  0  sin(phi) ] [ cos(theta)  -sin(theta)  0 ] [px]
      //  [py'] = [      0     1      0    ] [ sin(theta)   cos(theta)  0 ] [py]
      //  [pz']   [ -sin(phi)  0  cos(phi) ] [     0            0       1 ] [pz]
      //
      //  multiplied out to avoid the hassle of setting up matrix*matrix then matrix*vector multiplications
      //  not too bad with all the zeros

      rotated->points[i]->set_point_value(rotated->cos_phi_*rotated->cos_theta_*rotated->points[i]->get_point_value_x()
                                             -rotated->cos_phi_*rotated->sin_theta_*rotated->points[i]->get_point_value_y()
                                             +rotated->sin_phi_*rotated->points[i]->get_point_value_z(),
                                          rotated->sin_theta_*rotated->points[i]->get_point_value_x()
                                             +rotated->cos_theta_*rotated->points[i]->get_point_value_y(),
                                          -rotated->sin_phi_*rotated->cos_theta_*rotated->points[i]->get_point_value_x()
                                             +rotated->sin_phi_*rotated->sin_theta_*rotated->points[i]->get_point_value_y()
                                             +rotated->cos_phi_*rotated->points[i]->get_point_value_z());
      i++;
   }

   // check that the path is planar - all z components should be close
   i=1;
   while (i < rotated->points.size()) {
      if (! double_compare(rotated->points[i]->get_point_value_z(),rotated->points[0]->get_point_value_z(),rotated->tol)) {
         PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Path::rotateToXYplane passed a 3D path that is not planar.\n");
      }
      i++;
   }

   rotated->calculateBoundingBox();
   rotated->calculateNormal();

   return rotated;
}

void Path::rotatePoint (double *x, double *y, double *z)
{
   // nothing to do
   if (! rotated) return;

   // rotate the given test point by the plane's theta and phi angles
   // see Path::rotateToXYplane for notes on the calculation
   double xr=cos_phi_*cos_theta_*(*x)-cos_phi_*sin_theta_*(*y)+sin_phi_*(*z);
   double yr=sin_theta_*(*x)+cos_theta_*(*y);
   double zr=-sin_phi_*cos_theta_*(*x)+sin_phi_*sin_theta_*(*y)+cos_phi_*(*z);

   *x=xr;
   *y=yr;
   *z=zr;
}

// rotate an extra 180 degrees around the y axis
void Path::rotatePoint (double *x, double *y, double *z, bool spin180degrees)
{
   // nothing to do
   if (! rotated) return;

   double scale=1;
   if (spin180degrees) scale=-1;

   // rotate the given test point by the plane's theta and phi angles
   // see Path::rotateToXYplane for notes on the calculation
   double xr=scale*cos_phi_*cos_theta_*(*x)-scale*cos_phi_*sin_theta_*(*y)+scale*sin_phi_*(*z);
   double yr=sin_theta_*(*x)+cos_theta_*(*y);
   double zr=-scale*sin_phi_*cos_theta_*(*x)+scale*sin_phi_*sin_theta_*(*y)+scale*cos_phi_*(*z);

   *x=xr;
   *y=yr;
   *z=zr;
}

// rotate to a given path's rotation
void Path::rotateToPath (Path *rotatedPath)
{
   double x,y,z;

   if (!rotatedPath->rotated) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Path::rotateToPath passed a path that is not rotated.\n");
      return;
   }

   if (rotated) return;

   long unsigned int i=0;
   while (i < points.size()) {
      x=get_point_x(i); y=get_point_y(i); z=get_point_z(i);
      rotatedPath->rotatePoint(&x,&y,&z);
      points[i]->set_point_value(x,y,z);
      i++;
   }

   calculateBoundingBox();
   rotated=true;

   return;
}

// rotate to a given path's rotation with spin
void Path::rotateToPath (Path *rotatedPath, bool spin180degrees)
{
   double x,y,z;

   if (!rotatedPath->rotated) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Path::rotateToPath passed a path that is not rotated.\n");
      return;
   }

   if (rotated) return;

   long unsigned int i=0;
   while (i < points.size()) {
      x=get_point_x(i); y=get_point_y(i); z=get_point_z(i);
      rotatedPath->rotatePoint(&x,&y,&z,spin180degrees);
      points[i]->set_point_value(x,y,z);
      i++;
   }

   calculateBoundingBox();
   rotated=true;

   return;
}

// calculate the perpendicular distance from the given point to the plane
double Path::distanceFromPoint (double x, double y, double z)
{
  double x0,y0,z0,d,distance;

  if (!hasNormal) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Path::distanceFromPoint does not have a normal defined.\n");
      return -1;
  }

  if (points.size() == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Path::distanceFromPoint does not have any points.\n");
      return -1;
  }

  // get a point on the plane
  x0=points[0]->get_point_value_x();
  y0=points[0]->get_point_value_y();
  z0=points[0]->get_point_value_z();

  // calculate d for the equation of the plane
  d=-(nx*x0+ny*y0+nz*z0);

  // calculate the distance
  distance=abs(nx*x+ny*y+nz*z+d)/sqrt(nx*nx+ny*ny+nz*nz);

  return distance;
}

// define inside as including the lines themselves
// assumes that the path has been rotated to an x-y plane parallel to the z-plane
bool Path::is_point_inside (double xtr, double ytr, double ztr)
{
   long unsigned int i;
   if (points.size() == 0) return false;

   if (! rotated) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Path::is_point_inside passed an unrotated path.\n");
      return false;
   }

   // rotate the given test point by the plane's theta and phi angles
   rotatePoint(&xtr,&ytr,&ztr);

//   // cannot be inside if the z's do not align - just need to check one point
//   if (! double_compare(points[0]->get_point_value_z(),ztr,tol)) return false;

   // check all the points and accept if there is a match
   bool found=false;
   i=0;
   while (i < points.size()) {
      if (double_compare(points[i]->get_point_value_z(),ztr,tol)) {found=true; break;}
      i++;
   }
   if (!found) return false;

   // simple checks for points outside the path's bounding box
   if (!is_bound_by (xtr,xmin,xmax,tol)) return false;
   if (!is_bound_by (ytr,ymin,ymax,tol)) return false;

   // check the lines first - inside if on the lines
   i=0;
   while (i < points.size()-1) {
      if (is_point_on_line (xtr,ytr,points[i]->get_point_value_x(),points[i]->get_point_value_y(),
                                    points[i+1]->get_point_value_x(),points[i+1]->get_point_value_y(),tol)) return true;
      i++;
   }

   if (is_closed()) {
      if (is_point_on_line (xtr,ytr,points[i]->get_point_value_x(),points[i]->get_point_value_y(), 
                                    points[0]->get_point_value_x(),points[0]->get_point_value_y(),tol)) return true;
   }

   // check for interior points

   double theta=sum_of_angles(xtr,ytr);
   if (double_compare(theta,M_PI*2,tol)) return true;
   if (double_compare(theta,-M_PI*2,tol)) return true;

   return false;
}


// define interior as NOT including the lines themselves
// assumes that the path has been rotated to an x-y plane parallel to the z-plane
bool Path::is_point_interior (double xtr, double ytr, double ztr)
{
   long unsigned int i;
   if (points.size() == 0) return false;

   if (! rotated) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Path::is_point_interior passed an unrotated path.\n");
      return false;
   }

   // rotate the given test point by the plane's theta and phi angles
   rotatePoint(&xtr,&ytr,&ztr);

   // cannot be inside if the z's do not align - just need to check one point
   if (! double_compare(points[0]->get_point_value_z(),ztr,tol)) return false;

   // simple checks for points outside the path's bounding box
   if (!is_bound_by (xtr,xmin,xmax,tol)) return false;
   if (!is_bound_by (ytr,ymin,ymax,tol)) return false;

   // check the lines first - not interior if on the lines
   i=0;
   while (i < points.size()-1) {
      if (is_point_on_line (xtr,ytr,points[i]->get_point_value_x(),points[i]->get_point_value_y(),
                                    points[i+1]->get_point_value_x(),points[i+1]->get_point_value_y(),tol)) return false;
      i++;
   }

   if (is_closed()) {
      if (is_point_on_line (xtr,ytr,points[i]->get_point_value_x(),points[i]->get_point_value_y(),
                                    points[0]->get_point_value_x(),points[0]->get_point_value_y(),tol)) return false;
   }

   // check for interior points

   double theta=sum_of_angles(xtr,ytr);
   if (double_compare(theta,M_PI*2,tol)) return true;
   if (double_compare(theta,-M_PI*2,tol)) return true;

   return false;
}

bool Path::does_line_intersect (double x1, double y1, double z1, double x2, double y2, double z2)
{
   long unsigned int i;
   if (points.size() == 0) return false;

   if (! rotated) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Path::does_line_intersect passed an unrotated path.\n");
      return false;
   }

   // rotate the given end points by the plane's theta and phi angles
   rotatePoint(&x1,&y1,&z1);
   rotatePoint(&x2,&y2,&z2);

   // Only work with lines that are in the plane of the path.
   // It is assumed that the model is trying to be well-constructed and that there is no need
   // to handle the general case of a line not in the plane of the path.

   if (! double_compare(points[0]->get_point_value_z(),z1,tol)) return false;
   if (! double_compare(points[0]->get_point_value_z(),z2,tol)) return false;

   // check each segment
   i=0;
   while (i < points.size()-1) {
      if (do_intersect (points[i]->get_point_value_x(),points[i]->get_point_value_y(),points[i+1]->get_point_value_x(),points[i+1]->get_point_value_y(),
                        x1,y1,x2,y2,tol)) return true;
      i++;
   }

   if (is_closed()) {
      if (do_intersect (points[i]->get_point_value_x(),points[i]->get_point_value_y(),points[0]->get_point_value_x(),points[0]->get_point_value_y(),
                        x1,y1,x2,y2,tol)) return true;
   }

   return false;
}

bool Path::is_path_overlap (Path *test)
{
   long unsigned int i=0;
   while (i < test->points.size()-1) {
      if (does_line_intersect (test->get_point_x(i),test->get_point_y(i),test->get_point_z(i),test->get_point_x(i+1),test->get_point_y(i+1),test->get_point_z(i+1))) {
         return true;
      }
      i++;
   }

   if (is_closed()) {
      if (does_line_intersect (test->get_point_x(i),test->get_point_y(i),test->get_point_z(i),test->get_point_x(0),test->get_point_y(0),test->get_point_z(0))) {
         return true;
      }
   }

   return false;
}

bool Path::is_path_inside (Path *test)
{
   // check for points inside the path
   long unsigned int i=0;
   while (i < test->points.size()) {
      if (! is_point_inside(test->get_point_x(i),test->get_point_y(i),test->get_point_z(i))) return false;
      i++;
   }

   // check for crossing lines for the case where overlapping paths do not place a point inside another
   if (is_path_overlap(test)) return false;

   return true;
}


// Assumes the path:
// Path
//    name=m
//    point=(0,0,3)
//    point=(1,2,3)
//    point=(2,1,3)
//    point=(3,2,3)
//    point=(4,0,3)
//    point=(4,-1,3)
//    point=(0,-1,3)
//    closed=true
// EndPath
//
// the is value and expected values should be the same for all cases
//
void Path::test_is_point_inside_m ()
{
   double xt,yt,zt;
   string expected="";

   if (points.size() < 3) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Path::test_is_point_inside was passed a point or a line.\n");
      return;
   }

   Path *rotated=this->rotateToXYplane();

   xt=1; yt=0.5; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1; yt=1; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1; yt=1.5; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=3; yt=1.5; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=2; yt=1.001; zt=3; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=0.5; yt=1.5; zt=3; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=2; yt=2; zt=3; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=0.5; yt=2; zt=3; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1; yt=2; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=-1; yt=0; zt=3; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1; yt=0; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=2; yt=0.9999; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=3.00001; yt=2; zt=3; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=3; yt=2.0001; zt=3; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=3; yt=1.9999; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=2.9999; yt=2; zt=3; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1.5; yt=1.5; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1.5; yt=1.4999; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1.5; yt=1.50001; zt=3; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=2; yt=-1; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=2; yt=-1.0001; zt=3; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=3; yt=0.0001; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=3; yt=-0.0001; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=2.5; yt=1.5; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=2.4999; yt=1.5; zt=3; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=2.50001; yt=1.5; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=3.75; yt=1; zt=3; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   delete rotated;
}

// Assumes the path:
// Path
//    name=mr
//    point=(0,3,0)
//    point=(2,3,1)
//    point=(1,3,2)
//    point=(2,3,3)
//    point=(0,3,4)
//    point=(-1,3,4)
//    point=(-1,3,0)
//    closed=true
// EndPath
//
// the is value and expected values should be the same for all cases
//
void Path::test_is_point_inside_mr ()
{
   double xt,yt,zt;
   string expected="";

   if (points.size() < 3) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Path::test_is_point_inside was passed a point or a line.\n");
      return;
   }

   Path *rotated=this->rotateToXYplane();

   xt=0.5; yt=3; zt=1; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1; yt=3; zt=1; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1.5; yt=3; zt=1; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1.5; yt=3; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1.001; yt=3; zt=2; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1.5; yt=3; zt=0.5; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=2; yt=3; zt=2; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=2; yt=3; zt=0.5; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=2; yt=3; zt=1; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=0; yt=3; zt=-1; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=0; yt=3; zt=1; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=0.9999; yt=3; zt=2; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=2; yt=3; zt=3.00001; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=2.0001; yt=3; zt=3; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1.9999; yt=3; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=2; yt=3; zt=2.99999; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1.5; yt=3; zt=1.5; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1.4999; yt=3; zt=1.5; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1.50001; yt=3; zt=1.5; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=-1; yt=3; zt=2; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=-1.0001; yt=3; zt=2; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=0.0001; yt=3; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=-0.0001; yt=3; zt=3; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1.5; yt=3; zt=2.5; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1.5; yt=3; zt=2.49999; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1.5; yt=3; zt=2.50001; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=1; yt=3; zt=3.75; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   delete rotated;
}

// Assumes the path:
//Path
//   name=sqr2
//   point=(0.183012701892219,1.36602540378444,1.04903810567666)
//   point=(0.616025403784439,2.23205080756888,0.799038105676658)
//   point=(-0.133974596215561,2.73205080756888,1.23205080756888)
//   point=(-0.566987298107781,1.86602540378444,1.48205080756888)
//   closed=true
//EndPath
//
// the is value and expected values should be the same for all cases
//
void Path::test_is_point_inside_sqr2 ()
{
   double xt,yt,zt;
   string expected="";

   if (points.size() < 3) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Path::test_is_point_inside was passed a point or a line.\n");
      return;
   }

   Path *rotated=this->rotateToXYplane();

   xt=0.0245190528383291; yt=2.04903810567666; zt=1.14054445662277; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=0.182695714594112; yt=1.36739142918822; zt=1.04922111837855; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=-0.130804723234483; yt=2.71839055353103; zt=1.23022068054996; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=0.399519052838329; yt=1.79903810567666; zt=0.924038105676658; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=-0.350480947161671; yt=2.29903810567666; zt=1.35705080756888; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=-0.191987298107781; yt=1.61602540378444; zt=1.26554445662277; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=0.241025403784439; yt=2.48205080756888; zt=1.01554445662277; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=-0.555157171088858; yt=1.86968565782228; zt=1.47522068054995; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=0.604195276765517; yt=2.22839055353103; zt=0.80586823269558; expected="inside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=0.182579689190327; yt=1.36515937838065; zt=1.04928810567666; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=0.620355530803361; yt=2.24071106160672; zt=0.796538105676658; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=-0.425480947161671; yt=2.34903810567666; zt=1.4003520777581; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=0.245355530803361; yt=2.49071106160672; zt=1.01304445662277; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   xt=-0.0334936490538903; yt=0.933012701892219; zt=1.17403810567666; expected="outside";
   if (rotated->is_point_inside (xt,yt,zt)) cout << "(" << xt << "," << yt << "," << zt << ") is inside, expected=" << expected << endl;
   else cout << "(" << xt << "," << yt << "," << zt << ") is outside, expected=" << expected << endl;

   delete rotated;
}

Path::~Path ()
{
   long unsigned int i=0;
   while (i < points.size()) {
      delete points[i];
      i++;
   }
}

