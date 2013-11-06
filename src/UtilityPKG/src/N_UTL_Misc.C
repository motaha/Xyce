//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2013  Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename      : $RCSfile: N_UTL_Misc.C,v $
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Lon Waters, SNL
//
// Creation Date : 4/30/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.30.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:52 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>
#include <iostream>

#ifdef HAVE_CSTDLIB
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#ifdef Xyce_PARALLEL_MPI
 #include <mpi.h>
#endif

// ---------  Other Includes  -----------

//-----------------------------------------------------------------------------
// Function      : Xyce_exit
// Purpose       : 
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 07/26/04
//-----------------------------------------------------------------------------
void Xyce_exit( int code )
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Finalize();
#endif
  exit(code);
}


// Class string_param

//-----------------------------------------------------------------------------
// Function      : string_param
// Purpose       : Copy constructor
// Special Notes :
// Scope         : Public
// Creator       : Lon Waters, SNL
// Creation Date : 04/30/01
//-----------------------------------------------------------------------------

string_param::string_param(const string_param &right)
  :
  tag(right.tag),
  val(right.val)

{

}

//-----------------------------------------------------------------------------
// Function      : ~string_param
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Lon Waters, SNL
// Creation Date : 04/30/01
//-----------------------------------------------------------------------------

string_param::~string_param()

{
}

//-----------------------------------------------------------------------------
// Function      : getString
// Purpose       : Return the string value of the parameter.
// Special Notes :
// Scope         : Public
// Creator       : Lon Waters, SNL
// Creation Date : 04/30/01
//-----------------------------------------------------------------------------
string string_param::getString()
{
  return val;
}

//-----------------------------------------------------------------------------
// Function      : getReal
// Purpose       : Return the double value of the parameter.
// Special Notes :
// Scope         : Public
// Creator       : Lon Waters, SNL
// Creation Date : 04/30/01
//-----------------------------------------------------------------------------
double string_param::getReal()
{
  return ExtendedString(val).Value();
}

//-----------------------------------------------------------------------------
// Function      : getInteger
// Purpose       : Return the integer value of the parameter.
// Special Notes :
// Scope         : Public
// Creator       : Lon Waters, SNL
// Creation Date : 04/30/01
//-----------------------------------------------------------------------------
int string_param::getInteger()
{
  return atoi( val.c_str() );
}

//-----------------------------------------------------------------------------
// Description   : Turn string into integer value
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 12/01/05
//-----------------------------------------------------------------------------
int ExtendedString::Ival()
{
  int value;

  if (isInt())
    value = atoi((*this).c_str());
  else
    return 0;

  return value;
}

//-----------------------------------------------------------------------------
// Description   : Adjustment of data values when recognizable scale factor is
//                 present (e.g. 2.5n ==> 2.5E-9).
// Special Notes :
// Creator       : Alan Lundin
// Creation Date :
//-----------------------------------------------------------------------------
double ExtendedString::Value()
{
   char tmp[3];
   double value = atof((*this).c_str());
   int j = (*this).find_first_not_of("0123456789.-+eE", 0);
   if (j == (*this).npos) return value;
   switch ((*this)[j])
   {
      case 'T' : case 't' :
         return value*1.0e12;
      case 'G' : case 'g' :
         return value*1.0e9;
      case 'K' : case 'k' :
         return value*1.0e3;
      case 'M' :
      case 'm' :
         tmp[0] = tolower ((*this)[j]);
         if ((*this).size() > j+2)
         {
           tmp[1] = tolower ((*this)[j+1]);
           tmp[2] = tolower ((*this)[j+2]);
           if (tmp[1]  == 'i' && tmp[2] == 'l')
             return value*25.4e-6;
           else if (tmp[1] == 'e' && tmp[2] == 'g')
             return value*1.0e6;
         }
         return value*1.0e-3;
      case 'u' : case 'U' :
         return value*1.0e-6;
      case 'n' : case 'N' :
         return value*1.0e-9;
      case 'p' : case 'P' :
         return value*1.0e-12;
      case 'f' : case 'F' :
         return value*1.0e-15;
      default :
         return value;
   }
   return 0;
}

//-----------------------------------------------------------------------------
// Function      : isValue
// Purpose       : Test is string is a valid Value
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 10/28/05
//-----------------------------------------------------------------------------
bool ExtendedString::isValue()
{ 
  int stringPos ( 0 );
  int i = ( 0 );
  int stringSize ( (*this).size() );
  static const char *units[] = {"V", "VOLT", "VDC", "A", "AMP", "AMPERE", "F", 
                          "FARAD", "HENRY", "HY", "IL", "EG", "H", "HZ", 
                          "HERTZ", "OHM", "SECOND", "S", "METER", "M",
                          "MEG", "MIL", NULL};

  // Check for leading + or - sign.
  char ch ( (*this)[stringPos] );
  if ( ch == '+' || ch == '-' )
    ++stringPos;

  if ( stringPos == stringSize )
    return false; // string ended to soon to be a numeric value.

  ch = (*this)[stringPos];
  if ( (!isdigit(ch)) && ch != '.' )
    return false;

  while (isdigit(ch))
  {
    ++stringPos;
    if ( stringPos == stringSize ) return true; // all is well.
    ch = (*this)[stringPos];
  }

  if ( ch == '.' )
  {
    // Must have a digit before or after the decimal.
    if ( stringPos + 1 >= stringSize )
    {
      // Decimal is at end of string, must have digit before decimal.
      if ( stringPos == 0 ) 
        return false;
      else if ( !isdigit((*this)[stringPos-1]) )
        return false;
    }
    else if ( stringPos == 0 )
    {
      // Decimal is at beginning of string, must have digit after decimal.
      if ( stringPos + 1 >= stringSize )
        return false;
      else if ( !isdigit((*this)[stringPos+1]) )
        return false;
    }
    else if ( !isdigit((*this)[stringPos-1]) &&
              !isdigit((*this)[stringPos+1]) )
      return false;

    ++stringPos;
    if ( stringPos == stringSize ) return true; // all is well.
    ch = (*this)[stringPos];
  }

  while (isdigit(ch))
  {
    ++stringPos;
    if ( stringPos == stringSize ) return true; // all is well.
    ch = (*this)[stringPos];
  }

  // Check for exponent.
  if ( ch == 'E' || ch == 'e' )
  {
    ++stringPos;
    if ( stringPos == stringSize ) return true; // all is well.
    ch = (*this)[stringPos];

    // Look for exponent sign.
    if ( ch == '+' || ch == '-' )
    {
      ++stringPos;
      if ( stringPos == stringSize ) return true; // all is well.
      ch = (*this)[stringPos];
    }

    while (isdigit(ch))
    {
      ++stringPos;
      if ( stringPos == stringSize ) return true; // all is well.
      ch = (*this)[stringPos];
    }
  }

  if ((*this)[stringPos] == 'T' || (*this)[stringPos] == 't' ||
      (*this)[stringPos] == 'G' || (*this)[stringPos] == 'g' ||
      (*this)[stringPos] == 'K' || (*this)[stringPos] == 'k' ||
      (*this)[stringPos] == 'U' || (*this)[stringPos] == 'u' ||
      (*this)[stringPos] == 'N' || (*this)[stringPos] == 'n' ||
      (*this)[stringPos] == 'P' || (*this)[stringPos] == 'p' ||
      (*this)[stringPos] == 'F' || (*this)[stringPos] == 'f' )
    ++stringPos; 

  if ((*this)[stringPos] == 'M' || (*this)[stringPos] == 'm') 
  {
     char tmp[3];
     tmp[0] = tolower ((*this)[stringPos]);
     ++stringPos;

    if (stringSize >= stringPos +2)
    {
           tmp[1] = tolower ((*this)[stringPos]);
           tmp[2] = tolower ((*this)[stringPos + 1]);
           if ((tmp[1]  == 'i' && tmp[2] == 'l') || (tmp[1] == 'e' && tmp[2] == 'g'))
             stringPos = stringPos + 2;
    }

  }

  if (stringPos == stringSize)
    return true;

  ExtendedString u((*this).substr(stringPos, stringSize-stringPos));
  u = (*this).substr(stringPos, stringSize-stringPos);
  u.toUpper();
  i = 0;
  while (units[i] != NULL)
  {
    if (u == units[i++])
      return true;
  }

  return false;
}

//-----------------------------------------------------------------------------
// Function      : Xycepow10
// Purpose       : Compute powers of ten without using the "pow" function
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
double Xycepow10(int power)
{

  double retval=1.0;
  double fact;

  if (power < 0)
    {
      fact=0.1;
      power *= -1;
    }
  else
    {
      fact=10.0;
    }

  // handle odd powers
  if (power%2 == 1)
    {
      retval *= fact;
      power--;
    }

  // we now know power is even (or zero)
  // factors of 100 or .01
  for (; power > 0 ; power -= 2)
    {
      retval *= fact*fact;
    }


  return (retval);
}

namespace N_UTL {

  

int Ival(const string & tmpStr)
{
  int value;

  if (isInt(tmpStr))
   value = atoi(tmpStr.c_str());
  else
   return 0;

  return value;
}

double Value(const string & tmpStr)
{
  char tmp[3];
  double value = atof(tmpStr.c_str());
  int j = tmpStr.find_first_not_of("0123456789.-+eE", 0);
  if (j == tmpStr.npos) return value;
  switch (tmpStr[j])
  {
    case 'T' : case 't' :
       return value*1.0e12;
    case 'G' : case 'g' :
       return value*1.0e9;
    case 'K' : case 'k' :
       return value*1.0e3;
    case 'M' :
    case 'm' :
       tmp[0] = tolower (tmpStr[j]);
       if (tmpStr.size() > j+2)
       {
	 tmp[1] = tolower (tmpStr[j+1]);
	 tmp[2] = tolower (tmpStr[j+2]);
	 if (tmp[1]  == 'i' && tmp[2] == 'l')
	   return value*25.4e-6;
	 else if (tmp[1] == 'e' && tmp[2] == 'g')
	   return value*1.0e6;
       }
       return value*1.0e-3;
    case 'u' : case 'U' :
       return value*1.0e-6;
    case 'n' : case 'N' :
       return value*1.0e-9;
    case 'p' : case 'P' :
       return value*1.0e-12;
    case 'f' : case 'F' :
       return value*1.0e-15;
    default :
       return value;
  }
  return 0;
}

bool isValue(const string & tmpStr)
{
  int stringPos ( 0 );
  int i ( 0 );
  int stringSize ( tmpStr.size() );
  static const char *units[] = {"V", "VOLT", "VDC", "A", "AMP", "AMPERE", "F", 
                          "FARAD", "HENRY", "HY", "IL", "EG", "H", "HZ", 
                          "HERTZ", "OHM", "SECOND", "S", "METER", "M",
                          "MEG", "MIL", NULL};

  // Check for leading + or - sign.
  char ch ( tmpStr[stringPos] );
  if ( ch == '+' || ch == '-' )
    ++stringPos;

  if ( stringPos == stringSize )
    return false; // string ended to soon to be a numeric value.

  ch = tmpStr[stringPos];
  if ( (!isdigit(ch)) && ch != '.' )
    return false;

  while (isdigit(ch))
  {
    ++stringPos;
    if ( stringPos == stringSize ) return true; // all is well.
    ch = tmpStr[stringPos];
  }

  if ( ch == '.' )
  {
    // Must have a digit before or after the decimal.
    if ( stringPos + 1 >= stringSize )
    {
      // Decimal is at end of string, must have digit before decimal.
      if ( stringPos == 0 ) 
        return false;
      else if ( !isdigit(tmpStr[stringPos-1]) )
        return false;
    }
    else if ( stringPos == 0 )
    {
      // Decimal is at beginning of string, must have digit after decimal.
      if ( stringPos + 1 >= stringSize )
        return false;
      else if ( !isdigit(tmpStr[stringPos+1]) )
        return false;
    }
    else if ( !isdigit(tmpStr[stringPos-1]) &&
              !isdigit(tmpStr[stringPos+1]) )
      return false;

    ++stringPos;
    if ( stringPos == stringSize ) return true; // all is well.
    ch = tmpStr[stringPos];
  }

  while (isdigit(ch))
  {
    ++stringPos;
    if ( stringPos == stringSize ) return true; // all is well.
    ch = tmpStr[stringPos];
  }

  // Check for exponent.
  if ( ch == 'E' || ch == 'e' )
  {
    ++stringPos;
    if ( stringPos == stringSize ) return true; // all is well.
    ch = tmpStr[stringPos];

    // Look for exponent sign.
    if ( ch == '+' || ch == '-' )
    {
      ++stringPos;
      if ( stringPos == stringSize ) return true; // all is well.
      ch = tmpStr[stringPos];
    }

    while (isdigit(ch))
    {
      ++stringPos;
      if ( stringPos == stringSize ) return true; // all is well.
      ch = tmpStr[stringPos];
    }
  }

  if (tmpStr[stringPos] == 'T' || tmpStr[stringPos] == 't' ||
      tmpStr[stringPos] == 'G' || tmpStr[stringPos] == 'g' ||
      tmpStr[stringPos] == 'K' || tmpStr[stringPos] == 'k' ||
      tmpStr[stringPos] == 'U' || tmpStr[stringPos] == 'u' ||
      tmpStr[stringPos] == 'N' || tmpStr[stringPos] == 'n' ||
      tmpStr[stringPos] == 'P' || tmpStr[stringPos] == 'p' ||
      tmpStr[stringPos] == 'F' || tmpStr[stringPos] == 'f')
    ++stringPos;

  if (tmpStr[stringPos] == 'M' || tmpStr[stringPos] == 'm')
  {
     char tmp[3];
     tmp[0] = tolower (tmpStr[stringPos]);
     ++stringPos;

    if (stringSize >= stringPos +2)
    {
           tmp[1] = tolower (tmpStr[stringPos]);
           tmp[2] = tolower (tmpStr[stringPos + 1]);
           if ((tmp[1]  == 'i' && tmp[2] == 'l') || (tmp[1] == 'e' && tmp[2] == 'g'))
             stringPos = stringPos + 2;
    }

  }

  if (stringPos == stringSize)
    return true;

  ExtendedString u(tmpStr.substr(stringPos, stringSize-stringPos));
  u = tmpStr.substr(stringPos, stringSize-stringPos);
  u.toUpper();
  i = 0;
  while (units[i] != NULL)
  {
    if (u == units[i++])
      return true;
  }

  return false;
}

}

