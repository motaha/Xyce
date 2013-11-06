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
// Filename       : $RCSfile: N_IO_fwd.h,v $
//
// Purpose        : Forward declarations
//
// Special Notes  : Forward declaring everything as a class breaks if the implementation of the type changes (like during
//                  templatization)
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
//
// Creation Date  : 2013/04/18 18:01:27
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.5.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_IO_fwd_h
#define Xyce_N_IO_fwd_h

namespace Xyce {
namespace IO {

class CircuitMetadata;
class DeviceMetadata;
class OutputMgr;
struct PrintParameters;
struct Table;

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::CircuitMetadata N_IO_CircuitMetadata;
typedef Xyce::IO::DeviceMetadata N_IO_DeviceMetadata;
typedef Xyce::IO::OutputMgr N_IO_OutputMgr;

#endif // Xyce_N_DEV_fwd_h

