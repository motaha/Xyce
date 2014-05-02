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
// Revision Number: $Revision: 1.56 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_Misc_h
#define Xyce_N_UTL_Misc_h

typedef unsigned int    UINT;

#define INIT_SUPER_FLAG 1

#include <string>
#include <algorithm>

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

#include <N_UTL_NoCase.h>

namespace Xyce {
namespace Util {

void toUpper(std::string & tmp);
void toLower(std::string & tmp);
void removeWhiteSpace(std::string & tmp);
int Ival(const std::string & tmpStr);
bool isInt(const std::string & tmpStr);
bool Bval(const std::string & tmpStr);
bool isBool(const std::string & tmpStr);
bool possibleParam(const std::string & tmpStr);

double Value(const std::string & tmp);
bool isValue(const std::string & tmp);

std::ostream &word_wrap(std::ostream &os, const std::string &s, std::string::size_type line_length, const std::string &prefix, const std::string &prefix_first_line);


} // namespace Util
} // namespace Xyce

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

  tagged_param() : tag(std::string("")), param(0.0), given(false) {}

  tagged_param(const std::string & t, const double & p, bool g = false ) :
	tag(t), param(p), given(g) {}

  tagged_param( const tagged_param & right ) :
	tag(right.tag), param(right.param), given(right.given) {}

  tagged_param & operator=( const tagged_param & right )
  { tag = right.tag; param = right.param; given = right.given; return *this; }

  bool operator<( const tagged_param & rhs ) const
  { return (tag < rhs.tag); }

  friend bool operator==(const tagged_param & tp1, const tagged_param & tp2);

  friend bool operator!=(const tagged_param & tp1, const tagged_param & tp2);

  std::string tag;
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
    string_param(const std::string & t, const std::string & v) :
	  tag(t), val(v) {}
    // Constructor

    string_param( string_param const& right );
    // Copy constructor.

    ~string_param();
    // Destructor.

    std::string tag, val;

    std::string getString();
    double getReal();
    int getInteger();
    // Methods to get the "val" in the desired format.

    friend bool operator==(const string_param & sp1, const string_param & sp2);
    friend bool operator!=(const string_param & sp1, const string_param & sp2);

  private:
    double Value(const std::string& str);
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
class ExtendedString : public std::string
{
  public:
    ExtendedString(const char *X)
      : std::string(X)
    {}

    ExtendedString(const std::string & X)
      : std::string(X)
    {}

    ~ExtendedString()
    {}

    ExtendedString &toUpper() {
      Xyce::Util::toUpper(*this);

      return *this;
    }

    ExtendedString &toLower() {
      Xyce::Util::toLower(*this);

      return *this;
    }

    ExtendedString &removeWhiteSpace() {
      Xyce::Util::removeWhiteSpace(*this);

      return *this;
    }

    double Value() const {
      return Xyce::Util::Value(*this);
    }

    bool isValue() const {
      return Xyce::Util::isValue(*this);
    }

    int Ival() const {
      return Xyce::Util::Ival(*this);
    }

    bool isInt() const {
      return Xyce::Util::isInt(*this);
    }

    bool Bval() const {
      return Xyce::Util::Bval(*this);
    }

    bool isBool() const {
      return Xyce::Util::isBool(*this);
    }

    bool possibleParam() {
      return Xyce::Util::possibleParam(*this);
    }
};

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Function      : toUpper
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline void toUpper(std::string &s)
{
  for (std::string::iterator it = s.begin(); it != s.end(); ++it)
  {
    *it = toupper(*it);
  }
}

//-----------------------------------------------------------------------------
// Function      : toLower
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline void toLower(std::string &s)
{
  for(std::string::iterator it = s.begin(); it != s.end(); ++it)
  {
    *it = tolower(*it);
  }
}

//-----------------------------------------------------------------------------
// Function      : removeWhiteSpace
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline void removeWhiteSpace(std::string & tmp)
{
  std::string::iterator x = std::remove(tmp.begin(), tmp.end(), ' ');
  tmp.erase(x, tmp.end());
}

//-----------------------------------------------------------------------------
// Function      : Bval
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline bool Bval(const std::string & tmpStr)
{
  if (isValue(tmpStr))
  {
    return Value(tmpStr) != 0;
  }

  return equal_nocase(tmpStr, "TRUE");
}

//-----------------------------------------------------------------------------
// Function      : isBool
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline bool isBool(const std::string & s)
{
  return isValue(s) || equal_nocase(s, "TRUE") || equal_nocase(s, "FALSE");
}

inline int Ival(const std::string & tmpStr)
{
  if (isInt(tmpStr))
    return atoi(tmpStr.c_str());

  return 0;
}

} // namespace Util
} // namespace Xyce

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

#endif // Xyce_N_UTL_Misc_h
