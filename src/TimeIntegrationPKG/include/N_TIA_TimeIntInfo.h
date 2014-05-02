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
// Filename      : $RCSfile: N_TIA_TimeIntInfo.h,v $
// Purpose       : 
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 1/29/07
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.16 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_TIME_INT_INFO_H
#define Xyce_N_TIA_TIME_INT_INFO_H

// ---------- Standard Declarations ----------
#include <iosfwd>

// ---------- Forward Declarations ----------


//-----------------------------------------------------------------------------
// Class         : N_TIA_TimeIntInfo
//
// Purpose       : This class contains time integration info that needs to
//                 be passed from the time integration package to other
//                 parts of Xyce.
//
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 1/29/07
//-----------------------------------------------------------------------------
class N_TIA_TimeIntInfo
{
public:
  N_TIA_TimeIntInfo():
  currentOrder(1),
  numberOfSteps(0),
  usedOrder(1),
  nscsco(0),
  pdt(0.0),
  nextTimeStep(0.0),
  currTimeStep(0.0),
  currentTime(0.0),
  nextTime(0.0),
  finalTime(0.0),
  startingTimeStep(0.0),
  bpTol(0.0),
  dcopFlag(false),
  acopFlag(false),
  inputOPFlag(false),
  tranopFlag(false),
  transientFlag(false),
  dcsweepFlag(false),
  timeStepNumber(0),
  initTranFlag(false),
  beginIntegrationFlag(false),
  doubleDCOPStep(0),
  doubleDCOPEnabled(false),
  stepLoopIter(0),
  timeIntMode(0),
  sweepSourceResetFlag(false)   
  {};

  virtual ~N_TIA_TimeIntInfo() {};

  int currentOrder;
  int numberOfSteps;
  int usedOrder;
  int nscsco;

  double pdt;
  double nextTimeStep;
  double currTimeStep;
  double currentTime;
  double nextTime;
  double finalTime;
  double startingTimeStep;
  double bpTol;
  bool dcopFlag;
  bool acopFlag;
  bool inputOPFlag;
  bool tranopFlag;
  bool transientFlag;
  bool dcsweepFlag;
  int timeStepNumber;
  bool initTranFlag;
  bool beginIntegrationFlag;
  int doubleDCOPStep;
  bool doubleDCOPEnabled;
  int stepLoopIter;
  int timeIntMode;
  bool sweepSourceResetFlag;
};

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : "<<" operator for two-level info class.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 1/29/07
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const N_TIA_TimeIntInfo & tiInfo);

#endif // Xyce_N_TIA_TIME_INT_INFO_H


