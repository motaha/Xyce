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

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: N_IO_MeasureEquationEvaluation.C,v $
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.4.2.2 $
// Revision Date  : $Date: 2013/10/03 17:23:42 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------


// ----------   Xyce Includes   ----------
#include <N_IO_MeasureEquationEvaluation.h>
#include <N_ERH_ErrorMgr.h>


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureEquationEvaluation::N_IO_MeasureEquationEvaluation()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
N_IO_MeasureEquationEvaluation::N_IO_MeasureEquationEvaluation( const N_UTL_OptionBlock & measureBlock, N_IO_OutputMgr &outputMgr ):
  N_IO_MeasureBase( measureBlock, outputMgr)
{
  string msg = "MEASURE equation evaluation is not supported.  The measure item " + name_ +
    " will be ignored.";
  N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_WARNING, msg);
}


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureEquationEvaluation::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void N_IO_MeasureEquationEvaluation::updateTran( const double circuitTime, RCP< N_LAS_Vector > solnVecRCP)
{

}

//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureEquationEvaluation::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void N_IO_MeasureEquationEvaluation::updateDC( const vector<N_ANP_SweepParam> & dcParamsVec, RCP< N_LAS_Vector > solnVecRCP)
{

}
