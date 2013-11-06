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
// Filename       : $RCSfile: N_DEV_DeviceLevelKey.h,v $
//
// Purpose        : Map key for Device name and level
//
// Special Notes  :
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
//
// Creation Date  : 2013/04/18 18:01:27
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.4.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DeviceLevelKey_h
#define Xyce_N_DEV_DeviceLevelKey_h

#include <utility>
#include <string>
#include <string>
#include <functional>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#else
#ifdef HAVE_STRING_H
// On windows, we don't have this and cannot use strcasecmp, and need to use
// _stricmp instead.  Groan.
#include <string.h>
#define strcasecmp _stricmp
#else
error FIXME: Neither strings.h (POSIX) nor string.h (Windows) found.  Cannot find a case-insensitive string compare to use.
#endif  // HAVE_STRING_H
#endif  // HAVE_STRINGS_H

namespace Xyce {
namespace Device {

struct DeviceLevelKey : public std::pair<std::string, int>
{
    DeviceLevelKey()
      : std::pair<std::string, int>()
    {}

    DeviceLevelKey(const std::string &s, int i)
      : std::pair<std::string, int>(s, i)
    {}
};

struct DeviceLevelLess : public std::binary_function<DeviceLevelKey, DeviceLevelKey, bool> {
  bool operator()( const DeviceLevelKey &lhs , const DeviceLevelKey &rhs ) const {
    int i = strcasecmp( lhs.first.c_str() , rhs.first.c_str() );

    if (i == 0)
      return lhs.second < rhs.second;
    else
      return i < 0;
  }
};

std::ostream &operator<<(std::ostream &os, const Xyce::Device::DeviceLevelKey &device_level_key);

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_DeviceLevelKey_h
