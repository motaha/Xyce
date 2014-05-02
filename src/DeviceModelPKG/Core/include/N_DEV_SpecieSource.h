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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_SpecieSource.h,v $
//
// Purpose        : This file contains the PDE device instance base class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/15/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7 $
//
// Revision Date  : $Date: 2014/02/24 23:49:14 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_SpecieSource_h
#define Xyce_N_DEV_SpecieSource_h

// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_Const.h>
#include <N_DEV_CompositeParam.h>

namespace Xyce {
namespace Device {


//-----------------------------------------------------------------------------
// Class         : SpecieSource
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
class SpecieSource : public CompositeParam
{
  friend class ParametricData<SpecieSource>;

public:
  static ParametricData<SpecieSource> &getParametricData();

  SpecieSource ();
  bool processParam (Param & ndParam, std::string & param, DevicePDEInstance & di);
  void processParams ();

public:
  std::string name;         // this is also the index into the map.
};

// inline functions
//-----------------------------------------------------------------------------
// Function      : SpecieSource::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline std::ostream & operator<<(std::ostream & os, const SpecieSource & ds)
{
  os << ds.name << ":\n";
  os << std::endl;
  return os;
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::SpecieSource N_DEV_SpecieSource;

#endif

