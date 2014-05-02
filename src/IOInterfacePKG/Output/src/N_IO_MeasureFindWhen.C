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
// Filename      : $RCSfile: N_IO_MeasureFindWhen.C,v $
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.13.2.1 $
// Revision Date  : $Date: 2014/03/10 19:28:04 $
// Current Owner  : $Author: rlschie $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_MeasureFindWhen.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : FindWhen::FindWhen()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
FindWhen::FindWhen( const Util::OptionBlock & measureBlock, N_IO_OutputMgr &outputMgr ):
  Base(measureBlock, outputMgr),
  initialized_(false),
  lastIndepVarValue_(0.0),
  lastDepVarValue_(0.0)

{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

}

void FindWhen::prepareOutputVariables() 
{
  // this measurement can involve up to three solution variables
  numOutVars_ = outputVars_.size();
  outVarValues_.resize( numOutVars_, 0.0 );
}


//-----------------------------------------------------------------------------
// Function      : FindWhen::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void FindWhen::updateTran( const double circuitTime, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{
  

  if( !calculationDone_ && withinTransientWindow( circuitTime ) )
  {
    // we're in the transient window, now we need to calculate the value of this
    // measure and see if it triggers any specified rise, fall, cross windowing.
    double tempResult = 0.0;

    // update our outVarValues_ vector
    updateOutputVars( outVarValues_, circuitTime, solnVec, stateVec, storeVec, 0);
    if( !initialized_ )
    {
      // assigned last dependent and independent var to current time and outVarValue_[0] 
      // while we can't interpolate on this step, it ensurs that the initial history is
      // something realistic
      lastIndepVarValue_=circuitTime;
      lastDepVarValue_=outVarValues_[0];
      initialized_=true;
    }  

    double targVal=0.0;
    bool doneIfFound=false;
    if( outputValueTargetGiven_ )
    {
      // This is the form WHEN v(a)=fixed value
      targVal = outputValueTarget_;
      doneIfFound = true;
    }
    else
    {
      // This is the form WHEN v(a)= potentially changing value
      targVal = outVarValues_[1];
      // since we can't determine if the calculation is done at this piont
        // we don't set calculationDone_ = true;
      doneIfFound = false;
    }
    
    // I think type_ is always "WHEN" here need to check 
    if( (type_ == "WHEN") &&  withinRiseFallCrossWindow( outVarValues_[0], targVal ) )
    {
      // this is the simple case where Xyce output a value within tolerace 
      // to the target value 
      if( fabs(outVarValues_[0] - targVal) < minval_ )
      {
        calculationResult_ = circuitTime;
        calculationDone_ = doneIfFound;
      }
      else
      {
        // check an see if last point and this point bound the target point 
        double backDiff    = lastDepVarValue_ - targVal;
        double forwardDiff = outVarValues_[0] - targVal;

        // if we bound the target then either
        //  (backDiff < 0) && (forwardDiff > 0)  
        //   OR
        //  (backDiff > 0) && (forwardDiff < 0) 
        // or more simpley sgn( backDiff ) = - sgn( forwardDiff ) 
        if( ((backDiff < 0.0) && (forwardDiff > 0.0)) || ((backDiff > 0.0) && (forwardDiff < 0.0)) )
        {
          // bound the solution so interpolate to find the target time (or frequency etc)
          calculationResult_ = circuitTime - ( ((circuitTime - lastIndepVarValue_)/(outVarValues_[0]-lastDepVarValue_)) * (outVarValues_[0]-targVal) );
          calculationDone_ = doneIfFound;
        }
      }
    }

  }

  // remember the last point in case we need to interpolate
  // to the time when v(a)=x
  lastIndepVarValue_=circuitTime;
  lastDepVarValue_=outVarValues_[0]; 
}


//-----------------------------------------------------------------------------
// Function      : FindWhen::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void FindWhen::updateDC( const std::vector<N_ANP_SweepParam> & dcParamsVec, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{}

} // namespace Measure
} // namespace IO
} // namespace Xyce
