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
// Filename       : $RCSfile: N_NLS_NOX.C,v $
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
// Revision Number: $Revision: 1.6.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:48 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include "N_NLS_NOX.h"
#include "N_ERH_ErrorMgr.h"

// ----------   Namespaces   ----------

// ----------   Code   ----------

void N_NLS_NOX::error(const string msg)
{
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
}

void N_NLS_NOX::warning(const string msg)
{
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
}

void N_NLS_NOX::info(const string msg)
{
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, msg);
}

