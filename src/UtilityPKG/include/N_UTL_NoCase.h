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
// Filename       : $RCSfile: N_UTL_NoCase.h,v $
//
// Purpose        : Case insensitive functions and functors
//
// Special Notes  : 
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11.2.1 $
//
// Revision Date  : $Date: 2014/03/03 18:29:29 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_NoCase_h
#define Xyce_N_UTL_NoCase_h

#include <string>
#include <functional>

namespace Xyce {

/**
 * Compare two strings, case insensitively.
 *
 * Return negative if the first differing character of s1 is less than
 * s2.  Return 0 is they are equal, otherize positive.
 *
 * @param s1
 * @param s2
 *
 * @return Negative if s1 is less than s2, positive if s1 is greater than s2, zero otherwise
 */
int compare_nocase(const char *s1, const char *s2);

/**
 * Test if first string is less than second string, case insensitively.
 *
 * Return true if s0 is less than s1, case insensitively.
 *
 * @param s0
 * @param s1
 *
 * @return true if s0 is less than s1, case insensitively
 */
inline bool less_nocase(const std::string &s0, const std::string &s1) 
{
  return compare_nocase(s0.c_str(), s1.c_str()) < 0;
}

/**
 * Functor to test if first string is less than second string, case insensitively.
 *
 * Return true if s0 is less than s1, case insensitively.
 *
 * @param s0
 * @param s1
 *
 * @return true if s0 is less than s1, case insensitively
 */
struct LessNoCase : public std::binary_function<std::string,std::string,bool> 
 {
  bool operator()(const std::string &s0, const std::string &s1 ) const 
  {
    return less_nocase(s0, s1);
  }
};

 /**
 * Compare two strings for equality, case insensitively.
 *
 * Return true if s0 is equal to s1, case insensitively.
 *
 * @param s0
 * @param s1
 *
 * @return true if s0 is equal to s1, case insensitively
 */
inline bool equal_nocase(const std::string &s0, const std::string &s1) 
{
  return compare_nocase(s0.c_str(), s1.c_str()) == 0;
}

/**
 * Functor to compare two strings for equality, case insensitively.
 *
 * Return true if s0 is equal to s1, case insensitively.
 *
 * @param s0
 * @param s1
 *
 * @return true if s0 is equal to s1, case insensitively
 */
struct EqualNoCase : public std::binary_function<std::string,std::string,bool> 
 {
  bool operator()(const std::string &s0, const std::string &s1 ) const 
  {
    return equal_nocase(s0, s1);
  }
};

/**
 * Predicate to compare string for equality, case insensitively.
 *
 * Return true if predicate string is equal to s, case insensitively.
 *
 * @param s0
 * @param s1
 *
 * @return true if predicate string is equal to s, case insensitively
 */
struct EqualNoCasePred : public std::binary_function<std::string,std::string,bool> 
{
  EqualNoCasePred(const std::string &s)
    : s_(s)
  {}
  
  bool operator()(const std::string &s) const 
  {
    return equal_nocase(s_, s);
  }

  private:
    const std::string s_;
};

} // namepace Xyce

#endif // Xyce_N_UTL_NoCase_h
