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
// Filename       : $RCSfile: N_IO_MeasureFindWhen.h,v $
//
// Purpose        : Find time when a variable hits a target value
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
// Revision Number: $Revision: 1.11.2.1 $
// Revision Date  : $Date: 2014/03/10 19:28:03 $
// Current Owner  : $Author: rlschie $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureFindWhen_h
#define Xyce_N_IO_MeasureFindWhen_h

#include <N_IO_MeasureBase.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : FindWhen
// Purpose       : Find time when a variable hits a target value
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-------------------------------------------------------------------------
class FindWhen : public Base
{
public:
  FindWhen( const Util::OptionBlock & measureBlock, N_IO_OutputMgr &outputMgr );
  ~FindWhen() {};

    void prepareOutputVariables();
  void updateTran( const double circuitTime, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec);
  void updateDC( const std::vector<N_ANP_SweepParam> & dcParamsVec, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec);

private:
  int numOutVars_;
  std::vector<double> outVarValues_;

  // these are used to interpolate to the target independent var (time, frequency etc)  when 
  // the simulation has reported values that bound the target value
  bool initialized_;
  double lastIndepVarValue_;
  double lastDepVarValue_;

};

} // namespace Measure
} // namespace IO
} // namespace Xyce

typedef Xyce::IO::Measure::FindWhen N_IO_MeasureFindWhen;

#endif
