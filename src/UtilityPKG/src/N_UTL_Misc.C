//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
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
// Revision Number: $Revision: 1.39 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

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
  : tag(right.tag),
    val(right.val)
{}

//-----------------------------------------------------------------------------
// Function      : ~string_param
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Lon Waters, SNL
// Creation Date : 04/30/01
//-----------------------------------------------------------------------------

string_param::~string_param()
{}

//-----------------------------------------------------------------------------
// Function      : getString
// Purpose       : Return the string value of the parameter.
// Special Notes :
// Scope         : Public
// Creator       : Lon Waters, SNL
// Creation Date : 04/30/01
//-----------------------------------------------------------------------------
std::string string_param::getString()
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

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Description   : Adjustment of data values when recognizable scale factor is
//                 present (e.g. 2.5n ==> 2.5E-9).
// Special Notes :
// Creator       : Alan Lundin
// Creation Date :
//-----------------------------------------------------------------------------
double Value(const std::string &tmpStr)
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


bool isValue(const std::string & tmpStr)
{
  int stringPos = 0;
  int i = 0;
  int stringSize = tmpStr.size();
  
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
      tmpStr[stringPos] == 'F' || tmpStr[stringPos] == 'f' )
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

//-----------------------------------------------------------------------------
// Function      : isInt
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool isInt(const std::string & tmpStr)
{
  int j;

  if (tmpStr.empty())
    return false;

  if (tmpStr[0] == '-' || tmpStr[0] == '+')
    j = tmpStr.find_first_not_of("0123456789", 1);
  else
    j = tmpStr.find_first_not_of("0123456789");

  if (j == (int)std::string::npos)
    return true;

  // But there's one case where we *could* still be an int.  That would be
  // if we've got .[0]* at the current point.

  if (tmpStr[j] == '.')
  {
    std::string::size_type i = tmpStr.find_first_not_of("0",j+1);

    // If we find nothing but 0 after the ., we are still an int.
    if (i == std::string::npos)
     return true;
  }

  return false;
}

//-----------------------------------------------------------------------------
// Function      : possibleParam
// Purpose       : Test is string is a valid Value
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/02/05
//-----------------------------------------------------------------------------
bool possibleParam(const std::string &tmpStr)
{
  std::string first("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_$");
  std::string legal("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_0123456789.");
  std::string::iterator iterS, iterEnd;
  int i;
  bool ok = false;

  for (i=0 ; i<(int)(tmpStr.size()) ; ++i)
  {
    if (i == 0)
    {
      iterS=first.begin();
      iterEnd=first.end();
    }
    else
    {
      iterS=legal.begin();
      iterEnd=legal.end();
    }
    ok = false;
    while (iterS!=iterEnd)
    {
      if (*(iterS++) == tmpStr[i])
      {
        ok = true;
        break;
      }
    }
    if (!ok)
      break;
  }
  if (ok && isBool(tmpStr))
    ok = false;

  return ok;
}

namespace {

std::string::const_iterator
find_next_char(
  std::string::const_iterator	p,
  std::string::const_iterator	end,
  char				c)
{
  while (p != end && *p != c)
    p++;
  return p;
}

std::string::const_iterator
find_next_not_char(
  std::string::const_iterator	p,
  std::string::const_iterator	end,
  char				c)
{
  while (p != end && *p == c)
    p++;
  return p;
}

inline std::string::const_iterator find_next_space(std::string::const_iterator p, std::string::const_iterator end) {
  return find_next_char(p, end, ' ');
}

inline std::string::const_iterator find_next_endl(std::string::const_iterator p, std::string::const_iterator end) {
  return find_next_char(p, end, '\n');
}

inline std::string::const_iterator find_next_nonspace(std::string::const_iterator p, std::string::const_iterator end) {
  return find_next_not_char(p, end, ' ');
}

} // namespace <null>

std::ostream &
word_wrap(
  std::ostream &                os,
  const std::string &	        s,
  std::string::size_type        line_length,
  const std::string &	        prefix,
  const std::string &	        prefix_first_line)
{
  const std::string *u = &prefix_first_line;

  std::string::const_iterator p0, p1, p2, p3;
  p0 = p1 = p2 = s.begin();

  while (p2 != s.end() ) {

    // skip preceeding whitespace
    p1 = find_next_nonspace(p0, s.end());
    p3 = find_next_endl(p0, s.end());
    p2 = p1 = find_next_space(p1, s.end());
    do { // skip words
      p1 = find_next_nonspace(p1, s.end());
      p1 = find_next_space(p1, s.end());
      if (p3 < p1) {
	p2 = p3;
	break;
      }
      if (p1 - p0 + u->size() > line_length) // hit word end past line_length
	break;
      p2 = p1;
    } while (p2 != s.end());

    os << *u << std::string(p0, p2) << "\n";

    // if (p2 == p3) // If you want an embedded newline to mean
    //   u = &prefix_first_line;
    // else
    u = &prefix;

    p0 = p2 + 1;
  }

  return os;
}

} // namespace Util
} // namespace Xyce

