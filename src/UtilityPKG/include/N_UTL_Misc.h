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
// Filename       : $RCSfile: N_UTL_Misc.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 5/22/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.45.4.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:52 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef  _MISC_H
#define  _MISC_H

typedef unsigned int    UINT;

#define TRUE  1
#define FALSE 0

#define STATUS_SUCCESS 1
#define STATUS_FAILURE 0

#define INIT_SUPER_FLAG 1

// ---------- Standard Includes ----------
#include <string>

// Some implementations seem not to get this included by the time we use "toupper"
#ifdef HAVE_CCTYPE
#include <cctype>
#endif

// ----------   Xyce Includes   ----------
#include <N_UTL_Xyce.h>
#include <N_UTL_MachDepParams.h>

#ifndef M_PI
#ifdef MATH_H_HAS_M_PI
#include <math.h>
#else
  #define M_PI (2.0*asin(1.0))
#endif
#endif

inline double Xycemax ( double f1, double f2) { return f1 > f2 ? f1 : f2; }
inline double Xycemin ( double f1, double f2) { return f1 < f2 ? f1 : f2; }

// integer version:
inline int Xycemax ( int f1, int f2) { return f1 > f2 ? f1 : f2; }
inline int Xycemin ( int f1, int f2) { return f1 < f2 ? f1 : f2; }

double Xycepow10(int power);

//-----------------------------------------------------------------------------
// Function      : Xyce_exit
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/26/04
//-----------------------------------------------------------------------------
void Xyce_exit( int code );

//-----------------------------------------------------------------------------
// Class         : tagged_param
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
struct tagged_param
{

  tagged_param() : tag(string("")), param(0.0), given(false) {}

  tagged_param(const string & t, const double & p, bool g = false ) :
	tag(t), param(p), given(g) {}

  tagged_param( const tagged_param & right ) :
	tag(right.tag), param(right.param), given(right.given) {}

  tagged_param & operator=( const tagged_param & right ) 
  { tag = right.tag; param = right.param; given = right.given; return *this; }

  bool operator<( const tagged_param & rhs ) const
  { return (tag < rhs.tag); }

  friend bool operator==(const tagged_param & tp1, const tagged_param & tp2);

  friend bool operator!=(const tagged_param & tp1, const tagged_param & tp2);

  string tag;
  double      param;
  bool        given;
};

//-----------------------------------------------------------------------------
// Function      : operator==
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline bool operator==(const tagged_param & tp1, const tagged_param & tp2)
{
  return tp1.tag == tp2.tag;
}

//-----------------------------------------------------------------------------
// Function      : operator!=
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline bool operator!=(const tagged_param & tp1, const tagged_param & tp2)
{
  return !(tp1 == tp2);
}


//-----------------------------------------------------------------------------
// Class         : string_param
// Purpose       :
//                 Similar to "tagged_param" except designed where both the
//                 name and the value are stored in strings.
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
// Last Modified : 4/30/01
//-----------------------------------------------------------------------------
class string_param
{
  public:
    string_param(const string & t, const string & v) :
	  tag(t), val(v) {}
    // Constructor

    string_param( string_param const& right );
    // Copy constructor.

    ~string_param();
    // Destructor.

    string tag, val;

    string getString();
    double getReal();
    int getInteger();
    // Methods to get the "val" in the desired format.

    friend bool operator==(const string_param & sp1, const string_param & sp2);
    friend bool operator!=(const string_param & sp1, const string_param & sp2);

  private:
    double Value(const string& str);
};


//-----------------------------------------------------------------------------
// Purpose       : string_param "==" operator
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline bool operator == (const string_param & sp1, const string_param & sp2)
{
  return sp1.tag == sp2.tag;
}

//-----------------------------------------------------------------------------
// Purpose       : string_param "!=" operator
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline bool operator != (const string_param & sp1, const string_param & sp2)
{
  return !(sp1 == sp2);
}

//-----------------------------------------------------------------------------
// Class         : index_pair
// Purpose       : contains a pair of indices.  Mostly used for matrix
//                 access.
//
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/23/00
//-----------------------------------------------------------------------------
struct index_pair
{
  index_pair()
     : row(0), col(0), superFlag(0) {}

  index_pair (const int & rowTmp,  const int & colTmp,
	const int & superTmp = -1)
     : row(rowTmp), col(colTmp), superFlag(superTmp) {}

  friend bool operator==(const index_pair & ip1, const index_pair & ip2);

  int row, col, superFlag;
};

//-----------------------------------------------------------------------------
// Function      : operator==
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline bool operator==(const index_pair & ip1, const index_pair & ip2)
{
  return ((ip1.row == ip2.row) && (ip1.col == ip2.col));
}

//-----------------------------------------------------------------------------
// Class         : ExtendedString
// Purpose       :
//                  This  class  extends the C++ string class by adding
//                  functions toUpper, toLower, and removeWhiteSpace.
//                  Also includes spice style conversion to numeric Value.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class ExtendedString : public string
{
  public:
    ExtendedString        (const char *X) : string(X) {};
    ExtendedString        (const string & X) : string(X) {};
    ExtendedString & toUpper();
    ExtendedString & toLower();
    ExtendedString & removeWhiteSpace();
    double Value();
    bool isValue();
    int Ival();
    bool isInt();
    bool Bval();
    bool isBool();
    bool possibleParam();
};

//-----------------------------------------------------------------------------
// Function      : ExtendedString::toUpper
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline ExtendedString & ExtendedString::toUpper(void)
{
  //string::size_type size (length());
  //for( string::size_type i = 0; i < size; ++i )
  int size (length());
  for( int i = 0; i < size; ++i )
  {
    if (islower ( (*this)[i] ) )
    {
      (*this)[i] = toupper( (*this)[i] );
    }
  }

  return (*this);
}

//-----------------------------------------------------------------------------
// Function      : ExtendedString::toLower
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline ExtendedString & ExtendedString::toLower(void)
{
  //string::size_type size (length());
  //for( string::size_type i =0; i < size; ++i )
  int size (length());
  for( int i =0; i < size; ++i )
  {
    if (isupper ( (*this)[i] ) )
    {
      (*this)[i] = tolower( (*this)[i] );
    }
  }

  return (*this);
}

//-----------------------------------------------------------------------------
// Function      : ExtendedString::removeWhiteSpace
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline ExtendedString & ExtendedString::removeWhiteSpace(void)
{
  string::size_type N = string::npos - 1;

  while (N != string::npos)
  {
    N = find(" ");
    if (N != string::npos) erase (N, 1);
  }

  return (*this);
}

//-----------------------------------------------------------------------------
// Description   : Test if string is a bool value
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
inline bool ExtendedString::isBool()
{
  if (isValue())
    return true;
  ExtendedString tmp(*this);
  tmp.toUpper();
  if ( tmp == "TRUE" )
    return true;
  if ( tmp == "FALSE" )
    return true;
  return false;
}



//-----------------------------------------------------------------------------
// Description   : Test if string is an integer
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 12/01/05
//-----------------------------------------------------------------------------
inline bool ExtendedString::isInt()
{
  int j;

  //if ((*this).size() == 0)
  if ((*this).empty())
    return false;

  if ((*this)[0] == '-' || (*this)[0] == '+')
    j = (*this).find_first_not_of("0123456789", 1);
  else
    j = (*this).find_first_not_of("0123456789");

  if (j == (int)string::npos)
    return true;

  // But there's one case where we *could* still be an int.  That would be
  // if we've got .[0]* at the current point.

  if ((*this)[j] == '.')
  {
    int i;
    i=(*this).find_first_not_of("0",j+1);

    // If we find nothing but 0 after the ., we are still an int.
    if (i == (int)string::npos)
     return true;
  }

  return false;
}

//-----------------------------------------------------------------------------
// Description   : Turn string into bool value
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
inline bool ExtendedString::Bval()
{
  if (isValue())
  {
    if (Value() == 0)
      return false;
    else
      return true;
  }
  ExtendedString tmp(*this);
  tmp.toUpper();
  if ( tmp == "TRUE" )
    return true;
  if ( tmp == "FALSE" )
    return false;
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
inline bool ExtendedString::possibleParam()
{ 
  string first("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_$");
  string legal("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_0123456789.");
  string::iterator iterS, iterEnd;
  int i;
  bool ok = false;

  for (i=0 ; i<(int)((*this).size()) ; ++i)
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
      if (*(iterS++) == (*this)[i])
      {
        ok = true;
        break;
      }
    }
    if (!ok)
      break;
  }
  if (ok && isBool())
    ok = false;

  return ok;
}


namespace N_UTL {

void toUpper(string & tmp);
void toLower(string & tmp);
void removeWhiteSpace(string & tmp);
int Ival(const string & tmpStr);
bool isInt(const string & tmpStr);
bool Bval(const string & tmpStr);
bool isBool(const string & tmpStr);
bool possibleParam(const string & tmpStr);

double Value(const string & tmp);
bool isValue(const string & tmp);

//-----------------------------------------------------------------------------
// Function      : toUpper
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
inline void toUpper(string & tmp)
{
  int size (tmp.length());
  for( int i=0; i < size; ++i )
  {
    if (islower ( tmp[i] ) )
    {
      tmp[i] = toupper( tmp[i] );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : toLower
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
inline void toLower(string & tmp)
{
  int size (tmp.length());
  for( int i=0; i < size; ++i )
  {
    if (isupper ( tmp[i] ) )
    {
      tmp[i] = tolower( tmp[i] );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : removeWhiteSpace
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
inline void removeWhiteSpace(string & tmp)
{
  string::size_type N = string::npos - 1;
  while (N != string::npos)
  {
    N = tmp.find(" ");
    if (N != string::npos) tmp.erase (N, 1);
  }
}

//-----------------------------------------------------------------------------
// Function      : isInt
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
inline bool isInt(const string & tmpStr)
{
  int j(0);
  if (tmpStr.empty())
   return false;

  if (tmpStr[0] == '-' || tmpStr[0] == '+')
    j = tmpStr.find_first_not_of("0123456789", 1);
  else
    j = tmpStr.find_first_not_of("0123456789");

  if (j == (int)string::npos)
   return true;

  // But there's one case where we *could* still be an int.  That would be
  // if we've got .[0]* at the current point.
  if (tmpStr[j] == '.')
  {
    int i=tmpStr.find_first_not_of("0",j+1);

    // If we find nothing but 0 after the ., we are still an int.
    if (i == (int)string::npos)
     return true;
  }
  return false;
}

//-----------------------------------------------------------------------------
// Function      : Bval
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
inline bool Bval(const string & tmpStr)
{
  if (isValue(tmpStr))
  {
    if (Value(tmpStr) == 0)
      return false;
    else
      return true;
  }

  ExtendedString tmp(tmpStr);
  tmp.toUpper();
  if ( tmp == "TRUE" )
    return true;
  if ( tmp == "FALSE" )
    return false;
  return false;
}

//-----------------------------------------------------------------------------
// Function      : isBool
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
inline bool isBool(const string & tmpStr)
{
  if (isValue(tmpStr))
    return true;
  ExtendedString tmp(tmpStr);
  tmp.toUpper();
  if ( tmp == "TRUE" )
    return true;
  if ( tmp == "FALSE" )
    return true;
  return false;
}

//-----------------------------------------------------------------------------
// Function      : possibleParam
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
inline bool possibleParam(const string & tmpStr)
{
  string first("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_");
  string legal("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_0123456789");
  string::iterator iterS, iterEnd;
  int i;
  bool ok(false);

  int size ( tmpStr.size() );
  for (i=0 ; i< size; ++i)
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

}

#ifndef HAVE_IOTA
// iota is not part of the C++ standard.  It is an extension.

template <class _ForwardIterator, class _Tp>
void
iota(_ForwardIterator __first, _ForwardIterator __last, _Tp __value)
{
  while (__first != __last)
    *__first++ = __value++;
}
#endif

#endif

