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

#include "OpenParEMmaterials.hpp"

///////////////////////////////////////////////////////////////////////////////////////////
// keywordPair
///////////////////////////////////////////////////////////////////////////////////////////

bool keywordPair::int_limit_checks (string *keyword, int lineNumber)
{
   if (positive_required && int_value <= 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1017: %s at line %d is required to be positive.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber);
      return false;
   }

   if (non_negative_required && int_value < 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1018: %s at line %d is required to be non-negative.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber);
      return false;
   }

   if (int_value < lowerLimit) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1019: %s at line %d is required to be >= %g.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber,lowerLimit);
      return false;
   }

   if (int_value > upperLimit) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1020: %s at line %d is required to be <= %g.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber,upperLimit);
      return false;
   }

   return true;
}

bool keywordPair::dbl_limit_checks (string *keyword, int lineNumber)
{
   if (positive_required && dbl_value <= 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1021: %s at line %d is required to be positive.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber);
      return false; 
   }

   if (non_negative_required && dbl_value < 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1022: %s at line %d is required to be non-negative.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber);
      return false;
   }

   if (dbl_value < lowerLimit*(1-dbl_tolerance)) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1023: %s at line %d is required to be >= %g.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber,lowerLimit);
      return false;
   }

   if (dbl_value > upperLimit*(1+dbl_tolerance)) { 
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1024: %s at line %d is required to be <= %g.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber,upperLimit);
      return false;
   }

   return true;
}

bool keywordPair::point_limit_checks (string *keyword, int lineNumber)
{
   if (positive_required && (point_value.x <= 0 || point_value.y <= 0 || (point_value.dim == 3 && point_value.z <= 0))) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1025: %s at line %d is required to be positive.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber);
      return false;
   }

   if (non_negative_required && (point_value.x < 0 || point_value.y < 0 || (point_value.dim == 3 && point_value.z < 0))) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1026: %s at line %d is required to be non-negative.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber);
      return false;
   }

   if (point_value.x < lowerLimit*(1-dbl_tolerance) || point_value.y < lowerLimit*(1-dbl_tolerance) || (point_value.dim == 3 && point_value.z < lowerLimit*(1-dbl_tolerance))) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1027: %s at line %d is required to be >= %g.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber,lowerLimit);

      return false;
   }

   if (point_value.x > upperLimit*(1+dbl_tolerance) || point_value.y > upperLimit*(1+dbl_tolerance) || (point_value.dim == 3 && point_value.z > upperLimit*(1+dbl_tolerance))) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1028: %s at line %d is required to be <= %g.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber,upperLimit);
      return false;
   }

   return true;
}

bool keywordPair::limit_check (string type)
{
   bool fail=false;

   if (type.compare("int") == 0) {
      if (! int_limit_checks (&keyword, lineNumber)) fail=true;
   } else if (type.compare("double") == 0) {
      if (! dbl_limit_checks (&keyword, lineNumber)) fail=true;
   } else if (type.compare("point") == 0) {
      if (! point_limit_checks (&keyword, lineNumber)) fail=true;
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: bad selection in keywordPair::limit_check\n");
   }

   return fail;
}

bool keywordPair::match_alias (string *token)
{
   long unsigned int i=0;
   while (i < aliases.size()) {
      if (aliases[i].compare(*token) == 0) return true;
      i++;
   }
   return false;
}

bool keywordPair::loadBool (string *token, string *value_, int lineNumber_)
{
   // check for duplicate
   if (loaded) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1029: Duplicate entry at line %d for previous entry at line %d.\n",
                                   indent.c_str(),indent.c_str(),lineNumber_,lineNumber);
      return true;
   }

   // check for a boolean
   if (value_->compare("true") != 0 && value_->compare("false") != 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1030: %s value at line %d is invalid.\n",
                                   indent.c_str(),indent.c_str(),token->c_str(),lineNumber_);
      return true;
   }

   // get the value
   if (value_->compare("true") == 0) bool_value=true;
   else bool_value=false;

   // save it
   keyword=*token;
   value=*value_;
   lineNumber=lineNumber_;
   loaded=true;

   return false;
}

bool keywordPair::loadInt (string *token, string *value_, int lineNumber_)
{
   // check for duplicate
   if (loaded) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1031: Duplicate entry at line %d for previous entry at line %d.\n",
                                   indent.c_str(),indent.c_str(),lineNumber_,lineNumber);
      return true;
   }

   // check for a pure number
   if (!is_int(value_)) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1032: %s value at line %d is invalid.\n",
                                   indent.c_str(),indent.c_str(),token->c_str(),lineNumber_);
      return true;
   }

   // get the value
   try {int_value=stoi(*value_);}
   catch (const std::invalid_argument& ia) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1033: %s value at line %d is invalid.\n",
                                   indent.c_str(),indent.c_str(),token->c_str(),lineNumber_);
      return true;
   }

   // check the limits
   if (checkLimits && ! int_limit_checks (token, lineNumber_)) return true;

   // save it
   keyword=*token;
   value=*value_;
   lineNumber=lineNumber_;
   loaded=true;

   return false;
}

bool keywordPair::loadDouble (string *token, string *value_, int lineNumber_)
{
   // check for duplicate
   if (loaded) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1112: Duplicate entry at line %d for previous entry at line %d.\n",
                                   indent.c_str(),indent.c_str(),lineNumber_,lineNumber);
      return true;
   }

   // check for a pure number
   if (!is_double(value_)) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1111: %s value at line %d is invalid.\n",
                                   indent.c_str(),indent.c_str(),token->c_str(),lineNumber_);
      return true;
   }

   // get the value
   try {dbl_value=stod(*value_);}
   catch (const std::invalid_argument& ia) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1113: %s value at line %d is invalid.\n",
                                   indent.c_str(),indent.c_str(),token->c_str(),lineNumber_);
      return true;
   }

   // check the limits
   if (checkLimits && ! dbl_limit_checks (token, lineNumber_)) return true;

   // save it
   keyword=*token;
   value=*value_;
   lineNumber=lineNumber_;
   loaded=true;

   return false;
}

bool keywordPair::loadPoint (int dim, string *token, string *value_, int lineNumber_)
{
   // check for duplicate
   if (loaded) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1114: Duplicate entry at line %d for previous entry at line %d.\n",
                                   indent.c_str(),indent.c_str(),lineNumber_,lineNumber);
      return true;
   }

   // check
   if (!is_point(value_,dim)) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR1115: %s value at line %d is invalid.\n",
                                   indent.c_str(),indent.c_str(),token->c_str(),lineNumber_);
      return true;
   }

   // get values
   if (point_get (value_, &point_value.x, &point_value.y, &point_value.z, dim, indent, lineNumber_)) return true;
   point_value.dim=dim;

   // check the limits
   if (checkLimits && ! point_limit_checks (token, lineNumber_)) return true;

   point_value.dim=dim;
   loaded=true;

   return false;
}

bool keywordPair::dbl_compare (keywordPair *test)
{
   if (dbl_value == test->dbl_value) return true;
   if (dbl_value == 0 && fabs(test->dbl_value) < dbl_tolerance) return true;
   if (test->dbl_value == 0 && fabs(dbl_value) < dbl_tolerance) return true;
   if (fabs((dbl_value-test->dbl_value)/dbl_value) < dbl_tolerance) return true;
   return false;
}

bool keywordPair::value_compare (keywordPair *test)
{
   if (value.compare(test->value) == 0) return true;
   return false;
}

bool keywordPair::point_compare (keywordPair *a)
{
   if (point_value.dim != a->point_value.dim) return false;
   if (! double_compare(point_value.x,a->point_value.x,dbl_tolerance)) return false;
   if (! double_compare(point_value.y,a->point_value.y,dbl_tolerance)) return false;
   if (point_value.dim == 3 && ! double_compare(point_value.z,a->point_value.z,dbl_tolerance)) return false;
   return true;
}

double keywordPair::get_point_distance (keywordPair *a)
{
   if (get_point_value_dim() == 2) return sqrt(pow(get_point_value_x()-a->get_point_value_x(),2)+pow(get_point_value_y()-a->get_point_value_y(),2));
   // dim == 3
   return sqrt(pow(get_point_value_x()-a->get_point_value_x(),2)+pow(get_point_value_y()-a->get_point_value_y(),2)+pow(get_point_value_z()-a->get_point_value_z(),2));
}

bool keywordPair::is_close_point (keywordPair *a)
{
   if (! double_compare(get_point_value_x(),a->get_point_value_x(),1e-12)) return false;
   if (! double_compare(get_point_value_y(),a->get_point_value_y(),1e-12)) return false;
   if (get_point_value_dim() == 3 && ! double_compare(get_point_value_z(),a->get_point_value_z(),1e-12)) return false;
   return true;
}

void keywordPair::copy (keywordPair a)
{
   aliases.clear();
   long unsigned int i=0;
   while (i < a.aliases.size()) {
      aliases.push_back(a.aliases[i]);
      i++;
   }

   keyword=a.keyword;
   value=a.value;
   lineNumber=a.lineNumber;
   int_value=a.int_value;
   dbl_value=a.dbl_value;
   bool_value=a.bool_value;
   point_value=a.point_value;
   loaded=a.loaded;
   lowerLimit=a.lowerLimit;
   upperLimit=a.upperLimit;
   positive_required=a.positive_required;
   non_negative_required=a.non_negative_required;
   indent=a.indent;
   dbl_tolerance=a.dbl_tolerance;
   checkLimits=a.checkLimits;
}

keywordPair* keywordPair::clone ()
{
   keywordPair *b=new keywordPair();
   b->copy(*this);
   return b;
}

void keywordPair::print()
{
  long unsigned int i=0;
  while (i < aliases.size()) {
     PetscPrintf(PETSC_COMM_WORLD,"alias: %s\n",aliases[i].c_str());
     i++;
  }
  PetscPrintf(PETSC_COMM_WORLD,"keyword: %s\n",keyword.c_str());
  PetscPrintf(PETSC_COMM_WORLD,"value: %s\n",value.c_str());
  PetscPrintf(PETSC_COMM_WORLD,"lineNumber: %d\n",lineNumber);
  PetscPrintf(PETSC_COMM_WORLD,"int_value: %d\n",int_value);
  PetscPrintf(PETSC_COMM_WORLD,"dbl_value: %g\n",dbl_value);
  PetscPrintf(PETSC_COMM_WORLD,"loaded: %d\n",loaded);
  PetscPrintf(PETSC_COMM_WORLD,"lowerLimit: %g\n",lowerLimit);
  PetscPrintf(PETSC_COMM_WORLD,"upperLimit: %g\n",upperLimit);
  PetscPrintf(PETSC_COMM_WORLD,"postive_required: %d\n",positive_required);
  PetscPrintf(PETSC_COMM_WORLD,"non_negative_required: %d\n",non_negative_required);
  PetscPrintf(PETSC_COMM_WORLD,"indent: [%s]\n",indent.c_str());
  PetscPrintf(PETSC_COMM_WORLD,"dbl_tolerance: %g\n",dbl_tolerance);
  PetscPrintf(PETSC_COMM_WORLD,"checkLimits: %d\n",checkLimits);
}

