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
// Filename      : $RCSfile: N_ANP_Dakota.C,v $
// Purpose       : Class for handling Dakota optimization analysis.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.13 $
// Revision Date  : $Date: 2014/02/24 23:49:12 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_Dakota.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Analysis {

#ifndef Xyce_Dakota

//-----------------------------------------------------------------------------
// Function      : Dakota::run()
// Purpose       : provide stub function here for linking and
//                 generate an error if they're called
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool Dakota::run()
{
  std::string msg;
  msg = "Dakota::run() - Dakota analysis requested in a non-Dakota enabled build of Xyce";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return false;
}

#else

//-----------------------------------------------------------------------------
// Function      : Dakota::run()
// Purpose       : This is the main controlling loop for Dakota analysis.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/04/00
//-----------------------------------------------------------------------------
bool Dakota::run()
{
  bool integration_status = false;
  mainAnalysisRCPtr_->resetForStepAnalysis();
  integration_status = mainAnalysisRCPtr_->run();
  return integration_status;
}

#endif // Xyce_Dakota

} // namespace Analysis
} // namespace Xyce
