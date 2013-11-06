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
// Filename      : $RCSfile: N_IO_MeasureMin.C,v $
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.3.2.2 $
// Revision Date  : $Date: 2013/10/03 17:23:42 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------


// ----------   Xyce Includes   ----------
#include <N_IO_MeasureMin.h>
#include <N_ERH_ErrorMgr.h>


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureMin::N_IO_MeasureMin()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
N_IO_MeasureMin::N_IO_MeasureMin( const N_UTL_OptionBlock & measureBlock, N_IO_OutputMgr &outputMgr ):
  N_IO_MeasureBase(measureBlock, outputMgr),
  minimumValue_(0.0),
  initialized_(false)

{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // this measurement should have only one dependent variable.
  // Error for now if it doesn't
  numOutVars_ = depSolVarIterVector_.size();

  if ( numOutVars_ > 1 )
  {
    string msg = "Too many dependent variables for statistical measure, \"" + name_ + "\" Exiting.";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL, msg);
  }

  outVarValues_.resize( numOutVars_ );
  for( int i=0; i< numOutVars_; i++ )
  {
    outVarValues_[i] = 0.0;
  }

}


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureMin::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void N_IO_MeasureMin::updateTran( const double circuitTime, RCP< N_LAS_Vector > solnVecRCP)
{
  if( !calculationDone_ && withinTransientWindow( circuitTime ) )
  {
    // we're in the transient window, now we need to calculate the value of this
    // measure and see if it triggers any specified rise, fall, cross windowing.
    double tempResult = 0.0;

    // update our outVarValues_ vector
    for( int i=0; i< numOutVars_; i++ )
    {
      outVarValues_[i] = getOutputValue(depSolVarIterVector_[i], solnVecRCP);
    }

    if( !initialized_  )
    {
      minimumValue_ = outVarValues_[0];
      initialized_ = true;
    }

    if( minimumValue_ > outVarValues_[0] )
    {
      minimumValue_ = outVarValues_[0];
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureMin::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void N_IO_MeasureMin::updateDC( const vector<N_ANP_SweepParam> & dcParamsVec, RCP< N_LAS_Vector > solnVecRCP)
{

}


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureMin::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double N_IO_MeasureMin::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ =  minimumValue_;
  }
  return calculationResult_;
};
