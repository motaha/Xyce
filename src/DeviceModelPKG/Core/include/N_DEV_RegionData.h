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
// Filename       : $RCSfile: N_DEV_RegionData.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 07/19/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8 $
//
// Revision Date  : $Date: 2014/02/24 23:49:14 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_RegionData_h
#define Xyce_N_DEV_RegionData_h

// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_CompositeParam.h>

// ---------- Forward Declarations -------

//-----------------------------------------------------------------------------
// Class         : N_DEV_RegionData
// Purpose       :
// Special Notes : This class is intended to be a class for passing data into
//                 the Rxn Region constructor.
//
// Creator       : Eric Keiter, SNL
// Creation Date : 4/17/11
//-----------------------------------------------------------------------------

namespace Xyce {
namespace Device {

class RegionData : public CompositeParam
{
  friend class ParametricData<RegionData>;

public:
  static ParametricData<RegionData> &getParametricData();

  RegionData ();

  void processParams ();

#ifdef Xyce_DEBUG_DEVICE
  friend std::ostream & operator<<(std::ostream & os, const RegionData & rd);
#endif

private:

public:
  std::string name;
  std::string outName;
  std::string type;
  std::string reactionFile;
  double area;
  double xloc;

  bool doNothing; // set to true if reaction set should be ignored.
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::RegionData N_DEV_RegionData;

#endif

