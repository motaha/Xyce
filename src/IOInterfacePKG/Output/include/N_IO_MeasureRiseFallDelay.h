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
// Revision Number: $Revision: 1.15 $
// Revision Date  : $Date: 2014/02/24 23:49:20 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureRiseFallDelay_h
#define Xyce_N_IO_MeasureRiseFallDelay_h

#include <N_IO_MeasureBase.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : RiseFallDelay
// Purpose       : Measure Rise/fall/delay times
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-------------------------------------------------------------------------
class RiseFallDelay : public Base
{
public:
  RiseFallDelay( const Util::OptionBlock & measureBlock, N_IO_OutputMgr &outputMgr );
  ~RiseFallDelay() {};

    void prepareOutputVariables();
  void updateTran( const double circuitTime, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec);
  void updateDC( const std::vector<N_ANP_SweepParam> & dcParamsVec, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec);
  double getMeasureResult();

private:
  bool trigVariableLengthHistoryNeeded_;
  bool targVariableLengthHistoryNeeded_;
  double trigMax_;
  double targMax_;
  int trigResultIndex_;
  int targResultIndex_;
  double timeForTrig_;
  double timeForTarg_;
  bool trigMaxChanged_;
  bool targMaxChanged_;
  bool timeForTrigFound_;
  bool timeForTargFound_;
  bool trigOutputValueTargetChanged_;
  bool targOutputValueTargetChanged_;
  int numOutVars_;
  std::vector<double> outVarValues_;
  // these are vectors to store history information.
  // we need two independent var vectors because we will
  // trim the vectors dynamically to keep down on memory use
  std::vector<double> trigIndepVarHistory_; // usually time
  std::vector<double> trigVarHistory_;      // the trigger history
  std::vector<double> targIndepVarHistory_; // usually time
  std::vector<double> targetVarHistory_;    // the target history

};

} // namespace Measure
} // namespace IO
} // namespace Xyce

typedef Xyce::IO::Measure::RelativeError N_IO_MeasureRelativeError;

#endif
