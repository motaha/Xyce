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
// Filename      : $RCSfile: N_IO_MeasureFindWhen.C,v $
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.5.2.2 $
// Revision Date  : $Date: 2013/10/03 17:23:42 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------


// ----------   Xyce Includes   ----------
#include <N_IO_MeasureFindWhen.h>
#include <N_ERH_ErrorMgr.h>


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureFindWhen::N_IO_MeasureFindWhen()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
N_IO_MeasureFindWhen::N_IO_MeasureFindWhen( const N_UTL_OptionBlock & measureBlock, N_IO_OutputMgr &outputMgr ):
  N_IO_MeasureBase(measureBlock, outputMgr)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // this measurement can involve up to three solution variables
  numOutVars_ = depSolVarIterVector_.size();
  outVarValues_.resize( numOutVars_ );
  for( int i=0; i< numOutVars_; i++ )
  {
    outVarValues_[i] = 0.0;
  }
}


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureFindWhen::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void N_IO_MeasureFindWhen::updateTran( const double circuitTime, RCP< N_LAS_Vector > solnVecRCP)
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

    if( withinRiseFallCrossWindow( outVarValues_[0] ) )
    {
      if( (type_ == "WHEN") && outputValueTargetGiven_ )
      {
        // This is the form WHEN v(a)=value
        if( fabs(outVarValues_[0] - outputValueTarget_) < minval_ )
        {
          calculationResult_ = circuitTime;
          calculationDone_ = true;
        }
      }
      else if( (type_ == "WHEN") && (numOutVars_==2) )
      {
        // This is the form WHEN v(a)=value
        if( fabs(outVarValues_[0] - outVarValues_[1]) < minval_ )
        {
          calculationResult_ = circuitTime;
          calculationDone_ = true;
        }
      }
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureFindWhen::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void N_IO_MeasureFindWhen::updateDC( const vector<N_ANP_SweepParam> & dcParamsVec, RCP< N_LAS_Vector > solnVecRCP)
{

}
