//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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
// Filename       : $RCSfile: N_UTL_NameLevelKey.h,v $
//
// Purpose        : Map key for Device name and level
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
// Revision Number: $Revision: 1.1.2.2 $
//
// Revision Date  : $Date: 2014/03/03 18:29:29 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_NameLevelKey_h
#define Xyce_N_UTL_NameLevelKey_h

#include <utility>
#include <string>
#include <functional>

#include <N_UTL_NoCase.h>

namespace Xyce {

struct NameLevelKey : public std::pair<std::string, int>
{
  NameLevelKey()
    : std::pair<std::string, int>()
  {}

  NameLevelKey(const std::string &s, int i)
    : std::pair<std::string, int>(s, i)
  {}
};

struct NameLevelLess : public std::binary_function<NameLevelKey, NameLevelKey, bool> 
{
  bool operator()( const NameLevelKey &lhs , const NameLevelKey &rhs ) const 
  {
    //    int i = strcasecmp( lhs.first.c_str() , rhs.first.c_str() );
    int i = compare_nocase(lhs.first.c_str(), rhs.first.c_str());

    if (i == 0)
      return lhs.second < rhs.second;
    else
      return i < 0;
  }
};

std::ostream &operator<<(std::ostream &os, const NameLevelKey &device_level_key);

} // namespace Xyce

#endif // Xyce_N_UTL_NameLevelKey_h
