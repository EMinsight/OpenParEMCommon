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

#ifndef PATH_H
#define PATH_H

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <limits>
#include "petscsys.h"
#include "keywordPair.hpp"
#include "misc.hpp"

using namespace std;
using namespace mfem;

bool compare_xy (double, double, double, double, double);
bool compare_xy (double*, double*, double*, double*, double*);
bool compare_xyz (double, double, double, double, double, double, double);
bool compare_xyz (double*, double*, double*, double*, double*, double*, double*);

class Path {
   private:
      int startLine;
      int endLine;
      keywordPair name;
      vector<keywordPair *> points;
      keywordPair closed;
      double tol=1e-11; // 1e-12
      bool hasNormal;
      double nx,ny,nz;
      bool rotated;
      double theta;        // rotation about z-axis
      double phi;          // rotation about y-axis
      double sin_theta_;
      double cos_theta_;
      double sin_phi_;
      double cos_phi_;
      double xmax,xmin;
      double ymax,ymin;
      double zmax,zmin;
   public:
      Path (int, int);
      ~Path();
      bool load(int, string *, inputFile *);
      bool inBlock (int);
      bool check(string *);
      bool checkBoundingBox(Vector *, Vector *, string *, double);
      string get_name() {return name.get_value();}
      bool name_is_loaded() {return name.is_loaded();}
      int get_name_lineNumber() {return name.get_lineNumber();}
      bool get_closed() {return closed.get_bool_value();}
      int get_closed_lineNumber() {return closed.get_lineNumber();}
      int get_startLine() {return startLine;}
      int get_endLine() {return endLine;}
      long unsigned int get_points_size() {return points.size();}
      keywordPair* get_point (long unsigned int i) {return points[i];}
      double get_point_x (long unsigned int i) {return points[i]->get_point_value_x();}
      double get_point_y (long unsigned int i) {return points[i]->get_point_value_y();}
      double get_point_z (long unsigned int i) {return points[i]->get_point_value_z();}
      int get_point_dim (long unsigned int i) {return points[i]->get_point_value_dim();}
      void push_point (keywordPair *point) {points.push_back(point);}
      void pop_point () {points.pop_back();}
      bool compare (long unsigned int i, keywordPair test_point);
      void set_closed(bool value) {closed.set_bool_value(value); closed.set_loaded(true);}
      bool is_closed() {return closed.get_bool_value();}
      void set_name (string name_) {name.set_value(name_); name.set_loaded(true);}
      keywordPair* get_startPoint() {return points[0];}
      keywordPair* get_endPoint() {if (closed.get_bool_value()) return points[0]; return points[points.size()-1];}
      long unsigned int is_segmentOnLine (double, double, double, double);
      long unsigned int is_segmentOnLine (double, double, double, double, double, double);
      bool does_line_intersect (double, double, double, double, double, double);
      bool is_path_overlap (Path *);
      void subdivide2D (Path *);
      void subdivide3D (Path *);
      double sum_of_angles (double, double);
      bool calculateNormal ();
      bool is_rotated() {return rotated;}
      Path* rotateToXYplane ();
      void rotatePoint (double *, double *, double *);
      void rotatePoint (double *, double *, double *, bool);
      void rotateToPath (Path *);
      void rotateToPath (Path *, bool);
      bool is_point_inside (double, double, double);
      double distanceFromPoint (double, double, double);
      bool is_point_interior (double, double, double);
      bool is_path_inside (Path *);
      void test_is_point_inside_m ();
      void test_is_point_inside_mr ();
      void test_is_point_inside_sqr2 ();
      Path* clone();
      void calculateBoundingBox();
      void get_normal(double *nx_, double *ny_, double *nz_) {*nx_=nx; *ny_=ny; *nz_=nz;}
      void print(string);
      void output (ofstream *, int);
};

bool mergePaths (vector<Path *> *, vector<long unsigned int> *, vector<bool> *, string, string, Path **);

#endif

