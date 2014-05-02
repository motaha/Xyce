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
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10.2.1 $
//
// Revision Date  : $Date: 2014/03/03 18:29:28 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_IO_fwd_h
#define Xyce_N_IO_fwd_h

namespace Xyce {
namespace IO {

class CircuitBlock;
class CircuitContext;
class CircuitMetadata;
class CmdParse;
class DeviceBlock;
class DeviceMetadata;
class DistributionTool;
class FourierMgr;
class FunctionBlock;
class NetlistImportTool;
class Objective;
class OptionBlock;
class OutputMgr;
class ParameterBlock;
class PkgOptionsMgr;
class PkgOptionsReg;
class RestartMgr;
class RestartNode;
class SpiceSeparatedFieldTool;

class OutputFileBase;

struct PrintParameters;
struct Table;

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::CircuitBlock N_IO_CircuitBlock;
typedef Xyce::IO::CircuitContext N_IO_CircuitContext;
typedef Xyce::IO::CircuitMetadata N_IO_CircuitMetadata;
typedef Xyce::IO::CmdParse N_IO_CmdParse;
typedef Xyce::IO::DeviceBlock N_IO_DeviceBlock;
typedef Xyce::IO::DeviceMetadata N_IO_DeviceMetadata;
typedef Xyce::IO::DistributionTool N_IO_DistributionTool;
typedef Xyce::IO::FourierMgr N_IO_FourierMgr;
typedef Xyce::IO::FunctionBlock N_IO_FunctionBlock;
typedef Xyce::IO::Objective N_IO_Objective;
typedef Xyce::IO::OptionBlock N_IO_OptionBlock;
typedef Xyce::IO::OutputMgr N_IO_OutputMgr;
typedef Xyce::IO::ParameterBlock N_IO_ParameterBlock;
typedef Xyce::IO::PkgOptionsMgr N_IO_PkgOptionsMgr;
typedef Xyce::IO::PkgOptionsReg N_IO_PkgOptionsReg;
typedef Xyce::IO::RestartMgr N_IO_RestartMgr;
typedef Xyce::IO::RestartNode N_IO_RestartNode;
typedef Xyce::IO::SpiceSeparatedFieldTool N_IO_SpiceSeparatedFieldTool;

typedef Xyce::IO::OutputFileBase N_IO_OutputFileBase;

#endif // Xyce_N_IO_fwd_h

