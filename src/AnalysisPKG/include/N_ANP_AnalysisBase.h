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
// Filename       : $RCSfile: N_ANP_AnalysisBase.h,v $
//
// Purpose        : Base class for Analysis types
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.21 $
//
// Revision Date  : $Date: 2014/02/24 23:49:11 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_AnalysisBase_h
#define Xyce_N_ANP_AnalysisBase_h

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

#include <N_ANP_fwd.h>

#include <N_UTL_OptionBlock.h>
#include <N_TIA_TIAParams.h>

class N_TIA_Assembler;
class N_LAS_System;
class N_LOA_Loader;
class N_NLS_Manager;
class N_TIA_StepErrorControl;
class N_TIA_WorkingIntegrationMethod;

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : AnalysisBase
// Purpose       : Base class for common analysis functions
// Special Notes : 
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class AnalysisBase
{
  public: 

    AnalysisBase( AnalysisManager * anaManagerPtr );
    virtual ~AnalysisBase();

    virtual bool setAnalysisParams(const N_UTL_OptionBlock & paramsBlock) {return true;};
    virtual bool outputFailureStats () {return true;};

    virtual void setParamsWithOutputMgrAdapter 
      (RefCountPtr< OutputMgrAdapter > & outputMgrAdapterRCPtr) {};

    virtual int getStepIter () { return 0; }
    virtual int getStepNumber () { return stepNumber; }
    virtual void setStepNumber (int step) { stepNumber=step; }

    virtual void setTranStepNumber (int step) { tranStepNumber=step; }
    virtual int getTranStepNumber () { return tranStepNumber; }
    virtual void setSensFlag() {sensFlag_=true;}

    virtual bool run() = 0;
    virtual bool init() = 0;
    virtual bool loopProcess() = 0;
    virtual bool processSuccessfulStep() = 0;
    virtual bool processFailedStep() = 0;
    virtual bool finish() = 0;
    virtual bool handlePredictor() = 0;

    virtual void printStepHeader(std::ostream &os);
    virtual void printProgress(std::ostream &os);

    // mixed-signal
    virtual void preStepDetails(double maxTimeStepFromHabanero)
    {}
    
    virtual bool mixedSignalStep() {
      return true;
    }

    virtual bool finalizeStep () {
      return true;
    }

    // Two Level specific
    virtual bool twoLevelStep() {
      return true;
    }

    virtual bool resetForStepAnalysis();

    // utility functions common to more than one analysis
    int setupSweepLoop_( std::vector<SweepParam> & sweepParamVec );
    bool updateSweepParams_ (int loopIter, std::vector<SweepParam> & sweepParamVec);

    // Utility function for AnalysisManager to determine if a complex analysis type (like HB)
    // is performing another type of analysis under the hood.  This is necessary for populating
    // the N_TIA_TimeIntInfo struct.  For straightforward analysis types, this method is not
    // needed because the AnalysisManager already knows which analysis is being performed.
    virtual bool isAnalysis( int analysis_type ) { return false; }

    void resetAll();
    int saveLoopInfo ();
    virtual bool printLoopInfo(int start, int finish);

    virtual void setBeginningIntegrationFlag(bool bif) {beginningIntegration = bif;}
    virtual bool getBeginningIntegrationFlag() {return beginningIntegration;}

    virtual void setIntegrationMethod (int im) {integrationMethod_= im;}
    virtual unsigned int getIntegrationMethod () {return integrationMethod_;}

    virtual bool getInputOPFlag(){return inputOPFlag_;}

    // step statistic functions
    void gatherStepStatistics_ ();
    double getTotalLinearSolutionTime() const;
    double getTotalResidualLoadTime() const;
    double getTotalJacobianLoadTime() const;
    bool getDoubleDCOPEnabled ();
    int getDoubleDCOPStep();
    bool firstDoubleDCOPStep_ ();

  public:
    RefCountPtr< AnalysisManager  > anaManagerRCPtr_;
    RefCountPtr< N_TIA_Assembler > assemblerRCPtr_;
    RefCountPtr< N_LAS_System > lasSystemRCPtr_;
    RefCountPtr< N_LOA_Loader > loaderRCPtr_;
    RefCountPtr< N_NLS_Manager > nlsMgrRCPtr_;
    RefCountPtr< OutputMgrAdapter > outputMgrAdapterRCPtr_;
    RefCountPtr< N_TIA_StepErrorControl > secRCPtr_;
    RefCountPtr< N_TIA_WorkingIntegrationMethod > wimRCPtr_;

    N_TIA_TIAParams & tiaParams;

  //protected:
    bool beginningIntegration;

    // Current time-integration method flag.
    unsigned int integrationMethod_;

    // Time-integration step number counter.
    unsigned int stepNumber;
    // NOTE: For now tranStepNumber is the same as stepNumber, but later I will
    //       change it.  stepNumber will later include both dcop and tran steps.
    //       I haven't changed it yet b/c I need to check what devices call
    //       getStepNumber, and what they expect to get.
    unsigned int tranStepNumber;

    // Counters used in Monitoring the Cost of the Process:

    // Number of consecutive successful time-integration steps.
    unsigned int totalNumberSuccessfulStepsTaken_;
    unsigned int totalNumberSuccessStepsThisParameter_;

    // Total number of failed time-integration steps.
    unsigned int totalNumberFailedStepsAttempted_;

    unsigned int totalNumberJacobiansEvaluated_;
    unsigned int totalNumberIterationMatrixFactorizations_;
    unsigned int totalNumberLinearSolves_;
    unsigned int totalNumberFailedLinearSolves_;
    unsigned int totalNumberLinearIters_;
    unsigned int totalNumberResidualEvaluations_;
    unsigned int totalNonlinearConvergenceFailures_;

    double totalLinearSolutionTime_;
    double totalResidualLoadTime_;
    double totalJacobianLoadTime_;

    bool doubleDCOPFlag_;         // true if doing a double-DCOP is possible.
    int doubleDCOPStep_;          // current step in the DCOP loop.

    bool sensFlag_;
    bool inputOPFlag_;            // true if starting from an initial condition.


    std::vector<std::vector<int> > saveTimeI;
    std::vector<std::vector<double> > saveTimeD;

  protected:
    // command line object
    const N_IO_CmdParse & commandLine_;

  private:

};

//-----------------------------------------------------------------------------
inline double AnalysisBase::getTotalLinearSolutionTime() const
{
  return totalLinearSolutionTime_;
}

//-----------------------------------------------------------------------------
inline double AnalysisBase::getTotalResidualLoadTime() const
{
  return totalResidualLoadTime_;
}

//-----------------------------------------------------------------------------
inline double AnalysisBase::getTotalJacobianLoadTime() const
{
  return totalJacobianLoadTime_;
}

//-----------------------------------------------------------------------------
inline bool AnalysisBase::getDoubleDCOPEnabled () 
{
  return doubleDCOPFlag_;
}

//-----------------------------------------------------------------------------
inline int AnalysisBase::getDoubleDCOPStep () 
{
  return doubleDCOPStep_ ;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::firstDoubleDCOPStep_
// Purpose       : If the current step is the first step of
//                  a "doubleDCOP", then return "true".
//
//  Explanation:
//
//  If there are PDE semiconductor devices as part of this problem,
//  there may need to be a "double-pass"
//
//  first pass  = nonlinear poisson solution
//  second pass = drift diffusion solution
//
// Special Notes : Only PDE problems can ever return true.
// Scope         : 
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/21/04
//-----------------------------------------------------------------------------
inline bool AnalysisBase::firstDoubleDCOPStep_ ()
{
  return (getDoubleDCOPEnabled() && getDoubleDCOPStep() != tiaParams.lastDCOPStep);
}

} // namespace Analysis
} // namespace Xyce

typedef Xyce::Analysis::AnalysisBase N_ANP_AnalysisBase;

#endif

