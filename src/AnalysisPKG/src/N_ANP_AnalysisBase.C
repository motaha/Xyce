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
// Filename      : $RCSfile: N_ANP_AnalysisBase.C,v $
// Purpose       : Base class analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.35 $
// Revision Date  : $Date: 2014/02/24 23:49:12 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_ANP_AnalysisBase.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ERH_Message.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::AnalysisBase
// Purpose       : Constructor
// Special Notes : 
// Scope         : public
// Creator       : Todd S. Coffey, SNL.
// Creation Date : 01/29/08
//-----------------------------------------------------------------------------
AnalysisBase::AnalysisBase( AnalysisManager * anaManagerPtr )
  : tiaParams(anaManagerPtr->tiaParams),
    beginningIntegration(true),
    integrationMethod_(TIAMethod_NONE),
    stepNumber(0),
    tranStepNumber(0),
    totalNumberSuccessfulStepsTaken_(0),
    totalNumberSuccessStepsThisParameter_(0),
    totalNumberFailedStepsAttempted_(0),
    totalNumberJacobiansEvaluated_(0),
    totalNumberIterationMatrixFactorizations_(0),
    totalNumberLinearSolves_(0),
    totalNumberFailedLinearSolves_(0),
    totalNumberLinearIters_(0),
    totalNumberResidualEvaluations_(0),
    totalNonlinearConvergenceFailures_(0),
    totalLinearSolutionTime_(0.0),
    totalResidualLoadTime_(0.0),
    totalJacobianLoadTime_(0.0),
    doubleDCOPFlag_(false),
    doubleDCOPStep_(0),
    sensFlag_(false),
    inputOPFlag_(false),
    commandLine_(anaManagerPtr->getCommandLine())
{
  anaManagerRCPtr_ = rcp(anaManagerPtr, false );
  assemblerRCPtr_ = anaManagerPtr->assemblerPtr;
  lasSystemRCPtr_ = anaManagerPtr->lasSysPtr;
  loaderRCPtr_ = anaManagerPtr->loaderPtr;
  nlsMgrRCPtr_ = anaManagerPtr->nlsMgrPtr;
  outputMgrAdapterRCPtr_ = anaManagerPtr->outputMgrAdapterRCPtr_;
  secRCPtr_ = anaManagerPtr->secPtr_;
  wimRCPtr_ = anaManagerPtr->wimPtr;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::~AnalysisBase()
// Purpose       : destructor
// Special Notes : 
// Scope         : 
// Creator       : Todd S. Coffey, SNL.
// Creation Date : 01/29/08
//-----------------------------------------------------------------------------
AnalysisBase::~AnalysisBase()
{
}


//-----------------------------------------------------------------------------
// Function      : AnalysisBase::printStepHeader()
// Purpose       : Prints out step information.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
void AnalysisBase::printStepHeader(std::ostream &os) 
{}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::printProgress()
// Purpose       : Outputs run completion percentage and estimated
//                 time-to-completion.
// Special Notes : 
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 06/07/2002
//-----------------------------------------------------------------------------
void AnalysisBase::printProgress(std::ostream &os) 
{}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::resetForStepAnalysis() 
// Purpose       : When doing a .STEP sweep, some data must be reset to its 
//                 initial state.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 8/26/04
//-----------------------------------------------------------------------------
bool AnalysisBase::resetForStepAnalysis()
{
  totalNumberSuccessStepsThisParameter_ = 0;
  stepNumber = 0;
  beginningIntegration = true;

  return true;
}


//-----------------------------------------------------------------------------
// Function      : AnalysisBase::setupSweepLoop
// Purpose       : Processes sweep parameters.
// Special Notes : Used for DC and STEP analysis classes.
// Scope         : public
// Creator       : Eric R. Keiter, SNL.
// Creation Date : 08/21/04
//-----------------------------------------------------------------------------
int AnalysisBase::setupSweepLoop_( std::vector<SweepParam> & sweepParamVec )
{
  double pinterval = 1.0;
  double pcount = 0.0, pstart, pstop, pstep;

  // loop over the param containers, and check that all the params exist.
  // (the device package will complain if it can't find the param)
  double tmp = 0.0;
  for (int iparam=0;iparam<(int)sweepParamVec.size();++iparam)
  {
    tmp = loaderRCPtr_->getParamAndReduce(sweepParamVec[iparam].name);
  }


#ifdef Xyce_DEBUG_ANALYSIS
  if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
  {
    Xyce::dout() << std::endl << std::endl
                 << Xyce::subsection_divider << std::endl
                 << "AnalysisManager::setupSweepLoop" << std::endl;
  }
#endif

  // loop over the param containers:
  for (int iparam=0;iparam<(int)sweepParamVec.size();++iparam)
  {
#ifdef Xyce_DEBUG_ANALYSIS
    if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
    {
      Xyce::dout() << "name = " << sweepParamVec[iparam].name << std::endl;
    }
#endif
    // set interval:
    sweepParamVec[iparam].interval = static_cast<int> (pinterval);

    // This stuff should probably be moved up into the SweepParam class.

    // obtain next pinterval:
    if (sweepParamVec[iparam].type=="LIN")
    {
      pstart = sweepParamVec[iparam].startVal;
      pstop  = sweepParamVec[iparam].stopVal;
      pstep  = sweepParamVec[iparam].stepVal;
      // ----------
      // pcount = floor(((pstop - pstart)/pstep) + 1.0);
      // The computation of "pcount" above is notoriously prone to roundoff
      // error, especially if the "floor" function is an in-line function
      // and is subject to high levels of optimization on x86 processors.
      // The symptom is that the last step of a DC sweep (or other sweep)
      // gets lost for very specific combinations of start/stop/step.
      // The next few lines are an attempt to mitigate this roundoff issue,
      // which was present in Xyce for years, and was inherited from SPICE3F5,
      // from which the above expression was taken verbatim.

      // Compute the number of steps of size pstep between pstart and pstop
      pcount = floor(((pstop - pstart)/pstep));
      // Here we're checking that adding one more step doesn't pass pstop
      // by more than machine precision.  If we're within machine precision
      // of pstop by taking one more step, that must mean we undercounted 
      // due to roundoff in the division --- if we hadn't undercounted, we'd 
      // exceed pstop by (nearly) a full pstep.
      if ( fabs(pstop-(pstart+(pcount+1.0)*pstep)) < 2.0*N_UTL_MachineDependentParams::MachinePrecision())
      {
        pcount += 1.0;
      }
      
      // Pcount is now the exact number of steps of size pstep between pstart
      // and pstop, with roundoff handled mostly cleanly.

      // finally, because our actual loop does a loop from zero to maxStep-1,
      // we have to pad maxStep (as is done in the original pcount expression
      // above) to get the full range.
      pcount += 1.0;

      // done this way, we should no longer miss final steps of DC sweeps.
      // Fixed 31 Jul 2012.  This was bug 695 in Bugzilla, and had plagued
      // us since Xyce was first ported to Linux with GCC.
      // ----------

      sweepParamVec[iparam].maxStep = static_cast<int>(pcount);
#ifdef Xyce_DEBUG_ANALYSIS
      if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
      {
        Xyce::dout() << "pstart  = " << pstart << std::endl;
        Xyce::dout() << "pstop   = " << pstop  << std::endl;
        Xyce::dout() << "pstep   = " << pstep  << std::endl;
        Xyce::dout() << "pstop-pstart/pstep = " << ((pstop - pstart)/pstep) << std::endl;
        Xyce::dout() << "floor ()= " << floor(((pstop - pstart)/pstep)+1.0) << std::endl;
        Xyce::dout() << "pcount  = " << pcount << std::endl;
        Xyce::dout() << "maxStep = " << sweepParamVec[iparam].maxStep << std::endl;
      }
#endif
    }
    else if(sweepParamVec[iparam].type=="DEC")
    {
      double numSteps = static_cast<double>(sweepParamVec[iparam].numSteps);
      // stepMult could also be calculated as pow(10,(1/numSteps))
      double stepMult = exp(log(10.0)/numSteps);
      sweepParamVec[iparam].stepMult = stepMult;

      pstart   = sweepParamVec[iparam].startVal;
      pstop    = sweepParamVec[iparam].stopVal;
      pcount   = floor(fabs(log10(pstart) - log10(pstop)) * numSteps + 1);
      sweepParamVec[iparam].maxStep = static_cast<int>(pcount);
    }
    else if(sweepParamVec[iparam].type=="OCT")
    {
      double numSteps = static_cast<double>(sweepParamVec[iparam].numSteps);
      // stepMult could also be calculated as pow(2,1/(numSteps))
      double stepMult = exp(log(2.0)/numSteps);

      // changed to remove dependence on "log2" function, which apparently
      // doesn't exist in the math libraries of FreeBSD or the mingw 
      // cross-compilation suite.   Log_2(x)=log_e(x)/log_e(2.0)
      double ln2=log(2.0);

      sweepParamVec[iparam].stepMult = stepMult;
      pstart   = sweepParamVec[iparam].startVal;
      pstop    = sweepParamVec[iparam].stopVal;
      pcount   = floor(fabs(log(pstart) - log(pstop))/ln2 * numSteps + 1);
      sweepParamVec[iparam].maxStep = static_cast<int>(pcount);
    }
    else if(sweepParamVec[iparam].type=="LIST")
    {
      pcount = sweepParamVec[iparam].valList.size();
      sweepParamVec[iparam].maxStep = sweepParamVec[iparam].valList.size();
    }
    else
    {
      Report::UserError0() << " Unsupported STEP type";
    }
    pinterval *= pcount;

#ifdef Xyce_DEBUG_ANALYSIS
    if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
    {
      Xyce::dout() << "parameter = " << sweepParamVec[iparam].name << std::endl;
      Xyce::dout() << "pcount    = " << pcount << std::endl;
      Xyce::dout() << "pinterval = " << pinterval << std::endl;
    }
#endif
  }

  // At this point, pinterval equals the total number of steps 
  // for the step loop.
  return static_cast<int> (pinterval);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::updateSweepParams_ 
// Purpose       : Update parameters either for DC or STEP sweeps
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 8/26/04
//-----------------------------------------------------------------------------
bool AnalysisBase::updateSweepParams_ 
  (int loopIter, std::vector<SweepParam> & sweepParamVec)
{
  std::vector<SweepParam>::iterator iterParam;
  std::vector<SweepParam>::iterator firstParam = sweepParamVec.begin();
  std::vector<SweepParam>::iterator lastParam = sweepParamVec.end ();
  bool resetFlag=false;

  // set parameter(s)
  for (iterParam=firstParam; iterParam != lastParam;++iterParam)
  {
    iterParam->updateCurrentVal (loopIter);
    resetFlag = resetFlag || iterParam->getSweepResetFlag();
    loaderRCPtr_->setParam (iterParam->name, iterParam->currentVal);
  }

  // Tell the manager if any of our sweeps are being reset in this loop 
  // iteration.
  anaManagerRCPtr_->setSweepSourceResetFlag(resetFlag);
  
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::resetAll
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 12/16/10
//-----------------------------------------------------------------------------
void AnalysisBase::resetAll()
{
  stepNumber = 0;
  tranStepNumber = 0;
  totalNumberSuccessStepsThisParameter_ = 0;

  totalNumberSuccessfulStepsTaken_ = 0;
  totalNumberFailedStepsAttempted_ = 0;
  totalNumberJacobiansEvaluated_ = 0;
  totalNumberIterationMatrixFactorizations_ = 0;
  totalNumberLinearSolves_ = 0;
  totalNumberFailedLinearSolves_ = 0;
  totalNumberLinearIters_ = 0;
  totalNumberResidualEvaluations_ = 0;
  totalNonlinearConvergenceFailures_ = 0;
  totalResidualLoadTime_ = 0.0;
  totalJacobianLoadTime_ = 0.0;
  totalLinearSolutionTime_ = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::saveLoopInfo
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 12/16/10
//-----------------------------------------------------------------------------
int AnalysisBase::saveLoopInfo ()
{
  int pos = saveTimeI.size();
  saveTimeI.resize(pos+1);
  saveTimeD.resize(pos+1);

  saveTimeI[pos].push_back(totalNumberSuccessfulStepsTaken_);
  saveTimeI[pos].push_back(totalNumberFailedStepsAttempted_);
  saveTimeI[pos].push_back(totalNumberJacobiansEvaluated_);
  saveTimeI[pos].push_back(totalNumberIterationMatrixFactorizations_);
  saveTimeI[pos].push_back(totalNumberLinearSolves_);
  saveTimeI[pos].push_back(totalNumberFailedLinearSolves_);
  saveTimeI[pos].push_back(totalNumberLinearIters_);
  saveTimeI[pos].push_back(totalNumberResidualEvaluations_);
  saveTimeI[pos].push_back(totalNonlinearConvergenceFailures_);
  saveTimeD[pos].push_back(totalResidualLoadTime_);
  saveTimeD[pos].push_back(totalJacobianLoadTime_);
  saveTimeD[pos].push_back(totalLinearSolutionTime_);
	
  if (pos == 0)
  {
    int i,j;

    pos++;
    saveTimeI.resize(pos+1);
    saveTimeD.resize(pos+1);

    j = saveTimeI[0].size();
    saveTimeI[1].resize(j);
    for (i=0 ; i<j ; ++i)
    {
      saveTimeI[1][i] = saveTimeI[0][i];
      saveTimeI[0][i] = 0;
    }
    j = saveTimeD[0].size();
    saveTimeD[1].resize(j);
    for (i=0 ; i<j ; ++i)
    {
      saveTimeD[1][i] = saveTimeD[0][i];
      saveTimeD[0][i] = 0;
    }
  }

  return pos;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::printLoopInfo
// Purpose       : Prints out time loop information.
// Special Notes : Prints stats from save point start to save point finish.
//                 Special case 0,0 is entire run to this point
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
bool AnalysisBase::printLoopInfo(int start, int finish)
{
  bool bsuccess = true;
  int t1, t2, indI = 0, indD = 0;
  if (start == 0 && finish == 0)
  {
    t2 = saveLoopInfo();
    t1 = 0;
  }
  else
  {
    t1 = start;
    t2 = finish;
  }

  lout() << "\tNumber Successful Steps Taken:\t\t" << saveTimeI[t2][indI]-saveTimeI[t1][indI] << std::endl;
  indI++;

  lout() << "\tNumber Failed Steps Attempted:\t\t" << saveTimeI[t2][indI]-saveTimeI[t1][indI] << std::endl;
  indI++;

  lout() << "\tNumber Jacobians Evaluated:\t\t" << saveTimeI[t2][indI]-saveTimeI[t1][indI] << std::endl;
  indI++;

  lout() << "\tNumber Iteration Matrix Factorizations:\t" << saveTimeI[t2][indI]-saveTimeI[t1][indI] << std::endl;
  indI++;

  lout() << "\tNumber Linear Solves:\t\t\t" << saveTimeI[t2][indI]-saveTimeI[t1][indI] << std::endl;
  indI++;

  lout() << "\tNumber Failed Linear Solves:\t\t" << saveTimeI[t2][indI]-saveTimeI[t1][indI] << std::endl;
  indI++;

  if (saveTimeI[t2][indI] > saveTimeI[t1][indI])
  {
    lout() << "\tNumber Linear Solver Iterations:\t" << saveTimeI[t2][indI]-saveTimeI[t1][indI] << std::endl;
  }
  indI++;  // still needs to be incremented, even if info is not printed 

  lout() << "\tNumber Residual Evaluations:\t\t" << saveTimeI[t2][indI]-saveTimeI[t1][indI] << std::endl;
  indI++;

  lout() << "\tNumber Nonlinear Convergence Failures:\t" << saveTimeI[t2][indI]-saveTimeI[t1][indI] << std::endl;
  indI++;

  lout() << std::endl;

  lout() << "\tTotal Residual Load Time:\t\t" << saveTimeD[t2][indD]-saveTimeD[t1][indD] << " seconds" << std::endl;
  indD++;

  lout() << "\tTotal Jacobian Load Time:\t\t" << saveTimeD[t2][indD]-saveTimeD[t1][indD] << " seconds" << std::endl;
  indD++;

  lout() << "\tTotal Linear Solution Time:\t\t" << saveTimeD[t2][indD]-saveTimeD[t1][indD] << " seconds" << std::endl;
  indD++;

  return bsuccess;
}

void AnalysisBase::gatherStepStatistics_ ()
{
  if (secRCPtr_->newtonConvergenceStatus <= 0)
  {
    ++totalNonlinearConvergenceFailures_;
  }

  totalNumberJacobiansEvaluated_      += nlsMgrRCPtr_->getNumJacobianLoads();
  totalNumberLinearSolves_            += nlsMgrRCPtr_->getNumLinearSolves();
  totalNumberFailedLinearSolves_      += nlsMgrRCPtr_->getNumFailedLinearSolves();
  totalNumberLinearIters_             += nlsMgrRCPtr_->getTotalNumLinearIters();
  totalNumberResidualEvaluations_     += nlsMgrRCPtr_->getNumResidualLoads();
  totalNumberIterationMatrixFactorizations_ += nlsMgrRCPtr_->getNumJacobianFactorizations();
  totalLinearSolutionTime_            += nlsMgrRCPtr_->getTotalLinearSolveTime();
  totalResidualLoadTime_              += nlsMgrRCPtr_->getTotalResidualLoadTime();
  totalJacobianLoadTime_              += nlsMgrRCPtr_->getTotalJacobianLoadTime();
}

} // namespace Analysis
} // namespace Xyce

