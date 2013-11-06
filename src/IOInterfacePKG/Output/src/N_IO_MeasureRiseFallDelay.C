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
// Filename      : $RCSfile: N_IO_MeasureRiseFallDelay.C,v $
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.7.2.2 $
// Revision Date  : $Date: 2013/10/03 17:23:42 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------


// ----------   Xyce Includes   ----------
#include <N_IO_MeasureRiseFallDelay.h>
#include <N_ERH_ErrorMgr.h>


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureRiseFallDelay::N_IO_MeasureRiseFallDelay()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
N_IO_MeasureRiseFallDelay::N_IO_MeasureRiseFallDelay( const N_UTL_OptionBlock & measureBlock, N_IO_OutputMgr &outputMgr):
  N_IO_MeasureBase(measureBlock, outputMgr),
  trigHistoryNeeded_( false ),
  targHistoryNeeded_( false ),
  trigMax_(0.0),
  targMax_(0.0),
  trigResultIndex_(0),
  targResultIndex_(0),
  timeForTrig_(0.0),
  timeForTarg_(0.0)
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

  // check if this measure will need to record any history.  It will need history
  // info if trig_frac_max or targ_frac_max are given.  (These are cases where the
  // trigger / target are a fraction of the TRIG/TARG maximum.  Thus a history will
  // be needed to find the extrema and than use it find the fraction of that extrema
  if( trigFracMaxGiven_ )
  {
    trigHistoryNeeded_ = true;
  }
  if( targFracMaxGiven_ )
  {
    targHistoryNeeded_ = true;
  }
}


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureRiseFallDelay::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void N_IO_MeasureRiseFallDelay::updateTran( const double circuitTime, RCP< N_LAS_Vector > solnVecRCP)
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

    // outVarValues has our TRIG and TARG values so, first store them in
    // our history array if needed.

    // one memory savings.  We won't store the value if it hasn't changed
    bool recordData = false;

    if( trigHistoryNeeded_ )
    {
      int trigSize = trigVarHistory_.size();
      if( numOutVars_ > 0 )
      {
        if( trigSize == 0 )
        {
          // must record the first point
          recordData = true;
        }
        else if(fabs(trigVarHistory_[trigSize-1] - outVarValues_[0]) > minval_ )
        {
          //  trig changed enough that we need to record it.
          recordData = true;
        }
      }
    }

    if( targHistoryNeeded_ )
    {
      int targSize = targetVarHistory_.size();
      if( numOutVars_ > 0 )
      {
        if( targSize == 0 )
        {
          // must record the first point
          recordData = true;
        }
        else if( (numOutVars_ > 1) && (fabs(targetVarHistory_[targSize-1] - outVarValues_[1]) > minval_ ) )
        {
          // either trig or targ changed enough that we need to record it.
          recordData = true;
        }
      }
    }


    if ( recordData )
    {
      // earlier checks indicated that we needed to record history
      // so do that for both trig and targ here.  We record both at each independent var
      // point so that we don't have to have two independent var histories.
      if( numOutVars_ >0 )
      {
        trigIndepVarHistory_.push_back( circuitTime );
        trigVarHistory_.push_back( outVarValues_[0] );
        // update trigMax_ if needed.
        if( outVarValues_[0] > trigMax_ )
        {
          trigMax_ =  outVarValues_[0];
          //std::cout << "trigMax = " << trigMax_;
        }
      }
      if( numOutVars_ >1 )
      {
        targIndepVarHistory_.push_back( circuitTime );
        targetVarHistory_.push_back( outVarValues_[1] );
        // update trigMax_ if needed.
        if( outVarValues_[1] > targMax_ )
        {
          targMax_ =  outVarValues_[1];
          //std::cout << "  targMax = " << targMax_ << std::endl;
        }
      }

    }

    if( trigHistoryNeeded_ )
    {
      // find time for TRIG
      int trigSize = trigVarHistory_.size();
      for( int i=trigResultIndex_;i<(trigSize-1);i++ )
      {
        // goal is to find trigFracMax_ * trigMax_
        double difference = trigVarHistory_[i] - trigFracMax_ * trigMax_;
        double nextDifference = trigVarHistory_[i+1] - trigFracMax_ * trigMax_;

        bool diffNeg = ( difference < 0 );
        bool nextDiffNeg = ( nextDifference < 0 );

        if( (fabs(difference) < minval_) || ( diffNeg != nextDiffNeg ) )
        {
          timeForTrig_ = trigIndepVarHistory_[i];
          trigResultIndex_=i;
          //std::cout << "Setting timeForTrig_ = " << timeForTrig_ << " and trigResultIndex_ = " << trigResultIndex_ << std::endl;
          break;
        }
      }
    }

    if( targHistoryNeeded_ )
    {
      // find time for TARG
      int targSize = targetVarHistory_.size();
      for( int i=targResultIndex_;i<(targSize-1); i++ )
      {
        double difference = targetVarHistory_[i] - targFracMax_ * targMax_;
        double nextDifference = targetVarHistory_[i+1] - targFracMax_ * targMax_;
        bool diffNeg = ( difference < 0 );
        bool nextDiffNeg = ( nextDifference < 0 );

        // goal is to find targFracMax_ * trigMax_
        if( (fabs(difference) < minval_) || ( diffNeg != nextDiffNeg ) )
        {
          timeForTarg_ = targIndepVarHistory_[i];
          targResultIndex_=i;
          //std::cout << "Setting timeForTarg_ = " << timeForTarg_ << " and targResultIndex_ = " << targResultIndex_ << std::endl;
          break;
        }
      }
    }

    // prune old history that isn't needed
    // don't prune history every timestep.  Just when the trig/tarResultIndex_ is over
    // some max (1000 here)
    const int pruningHistoryThreshold=1000;
    if( trigResultIndex_> pruningHistoryThreshold )
    {
      vector<double>::iterator pruneIttr = trigVarHistory_.begin();
      pruneIttr += trigResultIndex_;
      trigVarHistory_.erase( trigVarHistory_.begin(), pruneIttr );
      pruneIttr = trigIndepVarHistory_.begin();
      pruneIttr += trigResultIndex_;
      trigIndepVarHistory_.erase( trigIndepVarHistory_.begin(), pruneIttr );
      trigResultIndex_ = 0;
    }

    if( targResultIndex_ > pruningHistoryThreshold )
    {
      vector<double>::iterator pruneIttr = targetVarHistory_.begin();
      pruneIttr += targResultIndex_;
      targetVarHistory_.erase( targetVarHistory_.begin(), pruneIttr );
      pruneIttr = targIndepVarHistory_.begin();
      pruneIttr += targResultIndex_;
      targIndepVarHistory_.erase( targIndepVarHistory_.begin(), pruneIttr );
      targResultIndex_ = 0;
    }


  }
}


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureRiseFallDelay::updateMeasures()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void N_IO_MeasureRiseFallDelay::updateDC( const vector<N_ANP_SweepParam> & dcParmsVec, RCP< N_LAS_Vector > solnVecRCP)
{

}


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureRiseFallDelay::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/27/2012
//-----------------------------------------------------------------------------
double N_IO_MeasureRiseFallDelay::getMeasureResult()
{
  calculationResult_=timeForTarg_ - timeForTrig_;
  return calculationResult_;
}
