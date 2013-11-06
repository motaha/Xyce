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
// Filename       : $RCSfile: N_IO_MeasureRiseFallDelay.h,v $
//
// Purpose        : Measure rise/fall delay times
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 03/10/2009
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.6.2.2 $
// Revision Date  : $Date: 2013/10/03 17:23:42 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureRiseFallDelay_h
#define Xyce_N_IO_MeasureRiseFallDelay_h

// ----------   Xyce Includes   ----------
#include <N_IO_MeasureBase.h>

// ---------- Forward Declarations ----------


//-------------------------------------------------------------------------
// Class         : N_IO_MeasureRiseFallDelay
// Purpose       : Measure Rise/fall/delay times
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-------------------------------------------------------------------------
class N_IO_MeasureRiseFallDelay : public N_IO_MeasureBase
{
public:
  N_IO_MeasureRiseFallDelay( const N_UTL_OptionBlock & measureBlock, N_IO_OutputMgr &outputMgr );
  ~N_IO_MeasureRiseFallDelay() {};

  void updateTran( const double circuitTime, RCP< N_LAS_Vector > solnVecRCP);
  void updateDC( const vector<N_ANP_SweepParam> & dcParamsVec, RCP< N_LAS_Vector > solnVecRCP);
  double getMeasureResult();

private:
  bool trigHistoryNeeded_;
  bool targHistoryNeeded_;
  double trigMax_;
  double targMax_;
  int trigResultIndex_;
  int targResultIndex_;
  double timeForTrig_;
  double timeForTarg_;
  int numOutVars_;
  vector<double> outVarValues_;
  // these are vectors to store history information.
  // we need two independent var vectors because we will
  // trim the vectors dynamically to keep down on memory use
  vector<double> trigIndepVarHistory_; // usually time
  vector<double> trigVarHistory_;      // the trigger history
  vector<double> targIndepVarHistory_; // usually time
  vector<double> targetVarHistory_;    // the target history

};

#endif
