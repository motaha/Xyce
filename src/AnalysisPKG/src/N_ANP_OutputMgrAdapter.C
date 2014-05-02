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
// Filename       : $RCSfile: N_ANP_OutputMgrAdapter.C,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11 $
//
// Revision Date  : $Date: 2014/02/24 23:49:12 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_ANP_OutputMgrAdapter.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::OutputMgrAdapter( )
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
OutputMgrAdapter::OutputMgrAdapter( ):
    stepAnalysisStepNumber_(0), 
    stepAnalysisMaxSteps_(0),
    dcAnalysisStepNumber_(0),
    dcAnalysisMaxSteps_(0)
{
  stepParamVecRCPtr_ = Teuchos::rcp( new std::vector<SweepParam> );
  dcParamVecRCPtr_ = Teuchos::rcp( new std::vector<SweepParam> );
}

} // namespace Analysis
} // namespace Xyce
