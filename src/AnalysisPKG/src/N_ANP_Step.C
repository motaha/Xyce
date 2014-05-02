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
// Filename      : $RCSfile: N_ANP_Step.C,v $
// Purpose       : .STEP Sweep class analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.26 $
// Revision Date  : $Date: 2014/02/24 23:49:12 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_Step.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : Step::setAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 6/22/10
//-----------------------------------------------------------------------------
bool Step::setAnalysisParams(const N_UTL_OptionBlock & paramsBlock)
{
#ifdef Xyce_DEBUG_ANALYSIS
  if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
  {
    Xyce::dout() << "In Step::setAnalysisParams" << std::endl;
  }
#endif

  std::list<N_UTL_Param>::const_iterator it_tp;
  std::list<N_UTL_Param>::const_iterator it_param;
  std::list<N_UTL_Param>::const_iterator it_type;
  std::list<N_UTL_Param>::const_iterator first = paramsBlock.getParams().begin();
  std::list<N_UTL_Param>::const_iterator last = paramsBlock.getParams().end();

  std::string msg;

  // first check to see that there is only 1 PARAM set.  They need to be in
  // separate lines for this to work.
  int countPar = 0;
  for (it_tp = first; it_tp != last; ++it_tp)
  {
    if (it_tp->uTag() == "PARAM")
    {
      ++countPar;
    }
  }
  if (countPar > 1)
  {
    msg =  "Step::setSTEPAnalysisParams\n";
    msg += "You have more than one step parameter on a single line.\n";
    msg += "Each parameter needs its own line.\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

  SweepParam sp;
#ifdef Xyce_DEBUG_ANALYSIS
  if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
  {
    for (it_tp = first; it_tp != last; ++it_tp)
    {
      Xyce::dout() << it_tp->uTag() ;
      Xyce::dout() << "\t";
      if (it_tp->uTag() == "PARAM" || it_tp->uTag() == "TYPE")
      {
        Xyce::dout() << it_tp->stringValue();
      }
      else
      {
        Xyce::dout() << it_tp->getImmutableValue<double>();
      }
      Xyce::dout() << std::endl;
    }
  }
#endif

  for (it_tp = first; it_tp != last; ++it_tp)
  {
    if (it_tp->uTag() == "TYPE")
    {
      it_type = it_tp;
      sp.type = it_tp->stringValue();
    }

    if (it_tp->uTag() == "PARAM")
    {
      it_param = it_tp;
      sp.name = it_tp->stringValue();
    }
  }

  it_tp = it_param;
  ++it_tp;
  if (sp.type == "LIN") // default
  {
    sp.startVal = it_tp->getImmutableValue<double>(); ++it_tp;
    sp.stopVal  = it_tp->getImmutableValue<double>(); ++it_tp;
    sp.stepVal  = it_tp->getImmutableValue<double>(); ++it_tp;
  }
  else if (sp.type == "DEC" || sp.type == "OCT")
  {
    sp.startVal = it_tp->getImmutableValue<double>(); ++it_tp;
    sp.stopVal  = it_tp->getImmutableValue<double>(); ++it_tp;
    sp.numSteps = it_tp->getImmutableValue<int>(); ++it_tp;
  }
  else if (sp.type == "LIST")
  {
    for (;it_tp!=last;++it_tp)
    {
      sp.valList.push_back(it_tp->getImmutableValue<double>());
    }
  }
  else
  {
    msg =  "Step::setSTEPAnalysisParams: ";
    msg += " unsupported STEP type\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

  // need to do a bunch of stuff to initialize the step loop.
  (*stepParamVec_).push_back(sp);

  // good idea to put this here, but this code isn't touched in
  // non .step runs and check_outputs will expcect the param vec
  // to exist (but be empty) even if this isn't a .step run
  //outputMgrAdapterRCPtr_->setStepParamVec( & stepParamVec_ );

  //stepLoopFlag_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Step::run()
// Purpose       : This is the main controlling loop for Step analysis.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/04/00
//-----------------------------------------------------------------------------
bool Step::run()
{
  bool bsuccess = true;
  bsuccess = init();
  bsuccess &= loopProcess();
  bsuccess &= finish();
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Step::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool Step::init()
{
#ifdef Xyce_DEBUG_ANALYSIS
  if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
  {
    Xyce::dout() << std::endl << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "Step::init" << std::endl;
  }
#endif

  stepLoopSize_ = setupSweepLoop_(*stepParamVec_);

  outputMgrAdapterRCPtr_->setStepAnalysisMaxSteps( stepLoopSize_ );

  anaManagerRCPtr_->stepLoopInitialized_ = true;

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Step::loopProcess()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool Step::loopProcess()
{
  std::string msg;
  bool integration_status = true;

  for (stepLoopIter_=0; stepLoopIter_< stepLoopSize_; ++stepLoopIter_)
  {
    updateSweepParams_(stepLoopIter_, *stepParamVec_);

#ifdef Xyce_DEBUG_ANALYSIS
    if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
    {
      // output parameter(s)
      std::vector <SweepParam>::iterator iterParam;
      std::vector <SweepParam>::iterator firstParam = stepParamVec_->begin();
      std::vector <SweepParam>::iterator lastParam = stepParamVec_->end ();
      for (iterParam=firstParam; iterParam != lastParam;++iterParam)
      {
        Xyce::dout() << "Step Analysis # " << stepLoopIter_<<"\t";
        Xyce::dout() << (*iterParam);
      }
    }
#endif

    secRCPtr_->resetAll ();
    anaManagerRCPtr_->getTIADataStore()->setZeroHistory();

    // solve the loop.
    outputMgrAdapterRCPtr_->setStepAnalysisStepNumber( stepLoopIter_);
    mainAnalysisRCPtr_->resetForStepAnalysis();
    integration_status &= mainAnalysisRCPtr_->run();

    outputMgrAdapterRCPtr_->outputRESULT( *(anaManagerRCPtr_->getTIADataStore()->currSolutionPtr), *(anaManagerRCPtr_->getTIADataStore()->currStatePtr), *(anaManagerRCPtr_->getTIADataStore()->currStorePtr) );

  } // end of for loop, and end of step analysis.

  outputMgrAdapterRCPtr_->finishOutputSTEP ();

  return integration_status;
}

//-----------------------------------------------------------------------------
// Function      : Step::processSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool Step::processSuccessfulStep()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Step::processFailedStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool Step::processFailedStep()
{
  return true;
}


//-----------------------------------------------------------------------------
// Function      : Step::finish()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool Step::finish()
{
  return true;
}

} // namespace Analysis
} // namespace Xyce

