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
// Filename       : $RCSfile: N_DEV_LaTexDoc.h,v $
//
// Purpose        : Functions to output parameter definitions
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
// Revision Date  : $Date: 2014/03/17 21:28:42 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_OutputPars_h
#define Xyce_N_DEV_OutputPars_h

#include <iosfwd>
#include <string>
#include <map>

#include <N_DEV_fwd.h>
#include <N_DEV_Pars.h>

namespace Xyce {
namespace Device {

namespace OutputMode {

enum Mode {
  DEFAULT, PARAM, INFO, DOC, DOC_CAT
};

} // namespace OutputMode

std::ostream &laTexDevice(std::ostream &os, const std::string &name, const int level, const int type, const std::string &device_description, const ParametricData<void> &parameters, const OutputMode::Mode format);

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_OutputPars_h
