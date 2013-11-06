//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2013  Sandia Corporation
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
// Filename       : $RCSfile: N_NLS_NOX.h,v $
//
// Purpose        : Error-handling for N_NLS_NOX
//
// Special Notes  :
//
// Creator        : Tammy Kolda, NLS, 8950
//
// Creation Date  : 01/31/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.5.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:47 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_h
#define Xyce_N_NLS_NOX_h

// ---------- Standard Includes ----------

#include <string>

// ----------   Xyce Includes   ----------
#include <N_UTL_Xyce.h>
// ----------   NOX Includes   ----------

// ---------- Forward Declarations ----------

// ---------- Namespace Declarations ----------

// ---------- Function Declarations ----------

namespace N_NLS_NOX {

  //---------------------------------------------------------------------------
  // Function      : error
  // Purpose       : This is a wrapper for
  //                 N_ERH_ErrorMsg::report(N_ERH_ErrorMsg::USR_ERROR_0, msg).
  //                 It throws "N_NLS_NOX Error" after calling report.
  //---------------------------------------------------------------------------
  void error(const string msg);

  //---------------------------------------------------------------------------
  // Function      : error
  // Purpose       : This is simply a wrapper for
  //                 N_ERH_ErrorMsg::report(N_ERH_ErrorMsg::USR_WARNING_0, msg)
  //---------------------------------------------------------------------------
  void warning(const string msg);

  //---------------------------------------------------------------------------
  // Function      : error
  // Purpose       : This is simply a wrapper for
  //                 N_ERH_ErrorMsg::report(N_ERH_ErrorMsg::USR_INFO_0, msg)
  //---------------------------------------------------------------------------
  void info(const string msg);

} // namespace N_NLS_NOX

#endif

