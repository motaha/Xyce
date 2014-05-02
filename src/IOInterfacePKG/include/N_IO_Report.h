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
// Filename      : $RCSfile: N_IO_Report.h,v $
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : David Baur
//
// Creation Date : 10/25/2013
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2.2.1 $
//
// Revision Date  : $Date: 2014/02/27 00:52:17 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_Report_H
#define Xyce_N_IO_Report_H

#include <N_IO_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_UTL_LogStream.h>

namespace Xyce {
namespace IO {

} // namespace IO
} // namespace Xyce

#define AssertLIDs(cmp) Assert(cmp, *this, #cmp)

#endif // Xyce_N_IO_Report_H

