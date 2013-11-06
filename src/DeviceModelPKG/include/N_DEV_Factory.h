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

/**
 * @file   N_DEV_Factory.h
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Mon Aug  5 09:38:00 2013
 * 
 * @brief  Device registry factory singleton declarations
 */
// Purpose        :
//
// Special Notes  :
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.6.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $


#ifndef Xyce_N_DEV_Factory_h
#define Xyce_N_DEV_Factory_h

#include <iosfwd>
#include <string>
#include <functional>

#include <N_DEV_fwd.h>
#include <N_DEV_DeviceLevelKey.h>
#include <N_DEV_Pars.h>
#include <N_UTL_Registry.h>

namespace Xyce {
namespace Device {

typedef Xyce::Plugin::Registry<int> Registry;
typedef Xyce::Plugin::UserPlugin<Registry, Device, Device *(*)(SolverState &, DeviceOptions &) > Fred;

typedef Xyce::Plugin::Registry<DeviceLevelKey, DeviceLevelLess> Registry2;
typedef Xyce::Plugin::UserPlugin<Registry2, Device, Device *(*)(SolverState &, DeviceOptions &) > Bob;

typedef Xyce::Plugin::UserSubroutine<Registry2, ParametricData<void> &() > Tom;

Registry &getXyceRegistry();
Registry2 &getXyceRegistry2();
Registry2 &getXyceInstanceRegistry();
Registry2 &getXyceModelRegistry();

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_Factory_h
