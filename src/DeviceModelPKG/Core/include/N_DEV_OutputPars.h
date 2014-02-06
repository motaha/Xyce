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
// Filename       : $RCSfile: N_DEV_OutputPars.h,v $
//
// Purpose        : Functions to output parameter definitions
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
// Revision Number: $Revision: 1.3.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_OutputPars_h
#define Xyce_N_DEV_OutputPars_h

#include <iosfwd>

#include <N_DEV_fwd.h>
#include <N_DEV_Pars.h>
#include <N_UTL_Param.h>

namespace Xyce {
namespace Device {

namespace OutputMode {

enum Mode {
  DEFAULT, PARAM, INFO, DOC, DOC_CAT
};

} // namespace OutputMode

void outputParams(const DeviceEntity &entity, OutputMode::Mode mode);
void outputParams(const DeviceEntity &entity, OutputMode::Mode mode, map<int, DeviceEntity *> & base);
void getUnitDescription(const DeviceEntity &entity, std::string &, ParameterUnit &, ParameterCategory &, std::string &, map<int, DeviceEntity *> & base);

std::ostream &outputParams(std::ostream &os, const ParametricData<void> &parametric_data, OutputMode::Mode mode);
std::ostream &outputParameterMap(std::ostream &os, const ParametricData<void>::ParameterMap &parameter_map);
std::ostream &outputDescriptor(std::ostream &os, const Descriptor &descriptor);

std::ostream &laTexDevice(std::ostream &os, const std::string &name, const int level, const int type, const std::string &device_description, const ParametricData<void> &parameters, const OutputMode::Mode format);

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_OutputPars_h
