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
// Filename       : $RCSfile: N_DEV_SolverState.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/25/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.42.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_UTL_Expression.h>
#include <N_DEV_SolverState.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : SolverState::SolverState
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/25/03
//-----------------------------------------------------------------------------
SolverState::SolverState () :
  pdt              (0.0),
  currTimeStep     (0.0),
  lastTimeStep     (0.0),
  currTime         (0.0),
  finalTime        (0.0),
  startingTimeStep (0.0),
  bpTol            (0.0),
  acceptedTime     (0.0),
  currentOrder     (0.0),
  usedOrder        (0.0), // BNB, integration order for 2-level stamp
  currFastTime     (0.0),
  finalFastTime    (0.0),
  currentHoldTime  (0.0),
  mpdeOnFlag       (false),
  blockAnalysisFlag(false),
  forceFinalOutput (false),
  doubleDCOPEnabled(false),
  doubleDCOPStep   (0),
  timeStepNumber   (0),
  ltraTimeIndex    (0),
  ltraTimeHistorySize (0),
  ltraDoCompact    (false),
  initTranFlag     (false),
  sweepSourceResetFlag(false),
  beginIntegrationFlag(true),
  dcopFlag         (true),
  inputOPFlag      (false),
  transientFlag    (true),
  dcsweepFlag      (false),
  tranopFlag       (true),
  acopFlag         (false),
  PDESystemFlag    (false),
  newtonIter       (0),
  stepLoopIter     (0),
  locaEnabledFlag  (false),
  continuationStepNumber (0),
  firstContinuationParam (true),
  firstSolveComplete (false),
  initJctFlag      (false),
  initFixFlag      (false),
  debugTimeFlag    (false),
  twoLevelNewtonCouplingMode (FULL_PROBLEM),
  PDEcontinuationFlag(false),
  pdeAlpha (1.0),
  chargeAlpha (1.0),
  chargeHomotopy (false),
  maxPDEContinuationSteps(10),
  currPDEContinuationStep(0),
  prevPDEContinuationStep(0),
  artParameterFlag(false),
  sizeParameterFlag(false),
  gainScale   (1, 1.0),
  nltermScale (1.0),
  sizeScale   (1.0),
  previousSizeScale   (1.0),
  bjtArtParameterFlag(false),
  ACspecified(false),
  TRANspecified(false),
  DCspecified(false),
  STEPspecified(false),
  OPspecified(false),
  MPDEspecified(false),
  HBspecified(false)
{
}


//-----------------------------------------------------------------------------
// Function      : SolverState::InitializeHomotopyBlockSize
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/25/03
//-----------------------------------------------------------------------------
void SolverState::InitializeHomotopyBlockSize(int numBlocks)
{
  gainScale.resize(numBlocks, 1.0);
}


//-----------------------------------------------------------------------------
// Function      : SolverState::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/05/2005
//-----------------------------------------------------------------------------
ostream& operator<<(ostream & os, const SolverState & ss)
{

  os << "-----------------------------------------------------" << endl;
  os << "  Device Package Solver State:" << endl;

  os << "  pdt = " << ss.pdt << endl;
  os << "  currTimeStep = " << ss.currTimeStep << endl;
  os << "  lastTimeStep = " << ss.lastTimeStep << endl;
  os << "  currTime = " << ss.currTime << endl;
  os << "  finalTime = " << ss.finalTime << endl;
  os << "  startingTimeStep = " << ss.startingTimeStep << endl;
  os << "  bpTol = " << ss.bpTol << endl;

  os << "  acceptedTime = " << ss.acceptedTime << endl;
  os << "  currentOrder = " << ss.currentOrder << endl;
  os << "  usedOrder =  " << ss.usedOrder << endl;

  os << "  mpdeOnFlag = ";
  if (ss.mpdeOnFlag)
  {
    os << "yes" << endl;
    os << "  currFastTime = " << ss.currFastTime << endl;
    os << "  finalFastTime = " << ss.finalFastTime << endl;
    os << "  blockAnalysisFlag  = " << ss.blockAnalysisFlag << endl;
  }
  else
  {
    os << "no" << endl;
  }

  os << "  timeStepNumber = " << ss.timeStepNumber << endl;
  os << "  ltraTimeIndex = " << ss.ltraTimeIndex << endl;
  os << "  ltraTimeStepHistorySize = " << ss.ltraTimeHistorySize << endl;
  os << "  ltraDoCompact = " << ss.ltraDoCompact << endl;
  os << "  newtonIter = " << ss.newtonIter << endl;
  os << "  stepLoopIter = "  << ss.stepLoopIter << endl;
  os << "  continuationStepNumber = " << ss.continuationStepNumber << endl;
  os << "  firstContinuationParam = ";
  if (ss.firstContinuationParam) os << "yes" << endl;
  else                           os << "no" << endl;

  os << "  firstSolveComplete = ";
  if (ss.firstSolveComplete) os << "yes" << endl;
  else                       os << "no" << endl;

  os << "  initTranFlag = ";
  if (ss.initTranFlag) os << "yes" << endl;
  else                 os << "no" << endl;

  os << "  beginIntegrationFlag = ";
  if (ss.beginIntegrationFlag) os << "yes" << endl;
  else                         os << "no" << endl;

  os << "  dcopFlag = ";
  if (ss.dcopFlag) os << "yes" << endl;
  else             os << "no" << endl;

  os << "  inputOPFlag = ";
  if (ss.inputOPFlag) os << "yes" << endl;
  else             os << "no" << endl;

  os << "  transientFlag = ";
  if (ss.transientFlag) os << "yes" << endl;
  else                  os << "no" << endl;

  os << "  dcsweepFlag = ";
  if (ss.dcsweepFlag) os << "yes" << endl;
  else                os << "no" << endl;

  os << "  tranopFlag = ";
  if (ss.tranopFlag) os << "yes" << endl;
  else               os << "no" << endl;

  os << "  acopFlag = ";
  if (ss.acopFlag) os << "yes" << endl;
  else             os << "no" << endl;

  os << "  PDESystemFlag = ";
  if (ss.PDESystemFlag) os << "yes" << endl;
  else                  os << "no" << endl;

  os << "  locaEnabledFlag = ";
  if (ss.locaEnabledFlag) os << "yes" << endl;
  else                    os << "no" << endl;

  os << "  initJctFlag = ";
  if (ss.initJctFlag) os << "yes" << endl;
  else                os << "no" << endl;

  os << "  initFixFlag = ";
  if (ss.initFixFlag) os << "yes" << endl;
  else                os << "no" << endl;

  os << "  sweepSourceResetFlag = ";
  if (ss.sweepSourceResetFlag) os << "yes" << endl;
  else                         os << "no" << endl;

  os << "  debugTimeFlag = ";
  if (ss.debugTimeFlag) os << "yes" << endl;
  else                  os << "no" << endl;

  os << "-----------------------------------------------------" << endl;
  os << endl;

  return os;
}

} // namespace Device
} // namespace Xyce

