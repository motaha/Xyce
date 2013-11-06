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
// Filename      : $RCSfile: N_TIA_StepErrorControl.h,v $
//
// Purpose       : This file defines the class for the time integration
//                 stepsize control algorithms.
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.77.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:49 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_STEP_ERROR_CONTROL_H
#define Xyce_N_TIA_STEP_ERROR_CONTROL_H

// ---------- Standard Declarations ----------
#include <iostream>
#include <set>
#include <N_UTL_BreakPoint.h>

// ---------- Forward Declarations ----------

class N_ANP_AnalysisManager;
class N_TIA_TimeIntegrationMethod;
class N_TIA_WorkingIntegrationMethod;
class N_TIA_TIAParams;
class N_TIA_BackwardDifferentiation15;
class N_TIA_OneStep;
class N_TIA_TimeIntInfo;

class N_PDS_Comm;
class N_IO_CmdParse;
class N_LOA_Loader;

//-----------------------------------------------------------------------------
// Class         : N_TIA_StepErrorControl
// Purpose       : This is the time integration step & error control class.
//
// Special Notes :  ERK. 6/19/2010  This class was originally designed to 
//                  handle much of the time step selection logic.  This 
//                  design intent was associated with the original 
//                  "old-DAE(ODE)" version of the integrator, which 
//                  is no longer in the code.  However numerous artifacts
//                  of this design are still here.  For example, the step size
//                  variables (such as currentTimeStep, currentTime) 
//                  are still here.  So at this point, this is mostly a data
//                  storage class for step variables that have a global 
//                  scope.
//
//                  The implementation of the "new-DAE" form 
//                  moved away from this design idea, and placed most of the
//                  step-size selection directly inside of specific algorithms.
//                  So, BDF15 handles stepsize selection specific to it in the
//                  BDF15 class, and the same goes for OneStep.
//
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
class N_TIA_StepErrorControl
{
  // member functions:
  public:
    // constructor  
    N_TIA_StepErrorControl(N_IO_CmdParse & cp,
                          N_ANP_AnalysisManager & anaManager,
                          N_TIA_TIAParams & tiaP);

    // Destructor:
    virtual ~N_TIA_StepErrorControl();

    void initializeStepSizeVariables();

    // The "stop time" is either the next discontinuity point, or the final time,
    // whichever comes first.
    void updateStopTime();

    double findNextStopTime();

    virtual void updateTwoLevelTimeInfo (const N_TIA_TimeIntInfo & tiInfo);

    // Print out the time-step information.
#ifdef Xyce_VERBOSE_TIME
    virtual void outputTimeInfo();
#endif

    bool initializeBreakPoints();

    // Requests dynamic breakpoint information from the loader.  Adds, subtracts
    // from the breakpoints array.
    bool updateBreakPoints();
    // add individual breakpoints
    void setBreakPoint(N_UTL_BreakPoint bp);
    void setBreakPoint(double bp);
    // signal that pause breakpoint has been reached
    void simulationPaused();
    // is the current time a pause time?
    bool isPauseTime();

    // Requests dynamic time step information from the loader.
    bool updateMaxTimeStep(double suggestedMaxTimeStep=0.0);

    // Sets the minimum time step based on machine precision.
    bool updateMinTimeStep();

    bool finished();

    void evaluateStepError ();

    void integrationStepReport_
      (bool step_attempt_status, bool sAStatus, bool testedError);

    void terseIntegrationStepReport_ 
      (bool step_attempt_status, bool sAStatus, bool testedError);

    // Gets the size of the restart data.
    virtual int restartDataSize( bool pack );

    // Output restart data.
    virtual bool dumpRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack );

    // Load restart data.
    virtual bool restoreRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack );

    // called if you want to start the time integration over from scratch.
    bool resetAll ();

    // Method which sets the TIA parameters object pointer.
    virtual bool setTIAParams();

    // Method which registers the TIA integration object pointer.
    bool registerWIMPtr(N_TIA_WorkingIntegrationMethod * wimPtr);

    void printBreakPoints (ostream & os) const;
    
    double getEstOverTol() const;

    void setTimeStep(double newTimeStep);
    
    bool getTranOPFlag() const;

  protected:

    void updatePauseTime(N_UTL_BreakPoint bp);

  private :
    // make default constructor private:
    N_TIA_StepErrorControl();

  // member data:
  public:

    // StepSize Variables 
    double startingTimeStep;
    double currentTimeStep;
    double lastAttemptedTimeStep;
    double lastTimeStep;
    double minTimeStep;
    double maxTimeStep;
    double maxTimeStepUser;
    double maxTimeStepBP;
    
    double savedTimeStep;
    
    double lastTime;
    double currentTime;
    double nextTime;
    double stopTime;
    double initialTime;
    double finalTime;

    double currentTimeStepRatio;
    double currentTimeStepSum;

    double lastTimeStepRatio;
    double lastTimeStepSum;

    int newtonConvergenceStatus;
  // number of newton iters
    int nIterations;
    
    int numberSuccessiveFailures;
    bool stepAttemptStatus;

    // Needed for 2-level:
    bool previousCallStepSuccessful;
    double estOverTol_;

  protected:

    // This initialization flag is true if tiaParams have been
    // registered and used.  Otherwise false.
    bool initializeFlag_;

    double minStepPrecisionFac_;

    double newtonStepReduction_;

    double restartTimeStepScale_;

    double tolAimFac_;

    // Pointer to the TIA parameters object.
    N_TIA_TIAParams & tiaParams_;

    // Pointer to the TIA control algorithm.
    N_ANP_AnalysisManager & anaManager_;

    // Pointer to the current integration method.
    N_TIA_WorkingIntegrationMethod * wimPtr_;

    N_LOA_Loader & loader_;

    set < N_UTL_BreakPoint > breakPoints_;
    set < N_UTL_BreakPoint >::iterator currentPauseBP;

  private :
    // command line object
    N_IO_CmdParse & commandLine_;

    // 03/08/04 tscoffe:  Local data for BDF 1-5 method
    int currentOrder_;      // Current order of integration
    int oldOrder_;          // previous order of integration
    int minOrder_;          // minimum order = max(1,user option minord)
    int maxOrder_;          // maximum order = min(5,user option maxord)
    int usedOrder_;         // order used in current step (used after currentOrder is updated)
    double alphas_;         // $\alpha_s$ fixed-leading coefficient of this BDF method
    vector<double> alpha_;  // $\alpha_j(n)=h_n/\psi_j(n)$ coefficient used in local error test
                            // note:   $h_n$ = current step size, n = current time step
    double alpha0_;         // $-\sum_{j=1}^k \alpha_j(n)$ coefficient used in local error test
    double cj_;             // $-\alpha_s/h_n$ coefficient used in local error test
    double ck_;             // local error coefficient 
    vector<double> sigma_;  // $\sigma_j(n) = \frac{h_n^j(j-1)!}{\psi_1(n)*\cdots *\psi_j(n)}$
    vector<double> gamma_;  // $\gamma_j(n)=\sum_{l=1}^{j-1}1/\psi_l(n)$ coefficient used to
                            // calculate time derivative of history array for predictor 
    vector<double> beta_;   // coefficients used to evaluate predictor from history array
    vector<double> psi_;    // $\psi_j(n) = t_n-t_{n-j}$ intermediary variable used to 
                            // compute $\beta_j(n)$
    int numberOfSteps_;     // number of total time integration steps taken
    int nef_;               // number of successive error failures (yes I know this is duplicated above)
    double  usedStep_;      // step-size used in current step
    int nscsco_;            // Number of Steps taken with Constant Step-size and Constant Order 
    double Ek_,Ekm1_,Ekm2_,Ekp1_; // Estimated LTE at currentOrder, currentOrder-1, and currentOrder-2
    double Est_;            // Represents which of Ek,Ekm1,Ekm2 will be used for step-size selection
    double Tk_,Tkm1_,Tkm2_,Tkp1_; // Estimated $\|h_n^{k+1} y_n^{(k+1)}\|$,  basically, order * E
    int  newOrder_;         // order of next step
    bool initialPhase_;     // determines if we're in the initial phase of integration where we
                            // double the step-size & increase the order. 
    double h0_safety_;      // safety factor in picking initial step-size
    double h0_max_factor_;  // h0 <= h0_max_factor * length of integration time
    double h_phase0_incr_;  // amount to increase step-sizes during initial phase
    double h_max_inv_;      // inverse of maximum step size
    double Tkm1_Tk_safety_; // magic number for determining order reduction
    double Tkp1_Tk_safety_; // magic number for determining order reduction
    double r_factor_;       // basic reduction factor for rr
    double r_safety_;       // safety factor in computing rr
    double r_fudge_;        // fudge factor in computing rr
    double r_min_;          // minimum reduction factor
    double r_max_;          // maximum reduction factor
    double r_hincr_test_;   // threshold for increasing the step-size
    double r_hincr_;        // factor used for increasing the step-size
    int max_LET_fail_;      // max number of error test failures before quitting.

    friend class N_TIA_BackwardDifferentiation15;
    friend class N_TIA_Gear12;
    friend class N_TIA_OneStep;
};

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::finished
// Purpose       : This function is called by the time integraion algorithm
//                 class to determine if the main time integration while
//                 loop should end.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel ComputationalSciences.
// Creation Date : 10/11/00
//-----------------------------------------------------------------------------
inline bool N_TIA_StepErrorControl::finished()
{
  double delta = currentTime - finalTime;

  delta = fabs(delta);
  return (delta < 1.0e-10 * (finalTime - initialTime));
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::registerWIMPtr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel ComputationalSciences.
// Creation Date : 09/06/01
//-----------------------------------------------------------------------------
inline bool N_TIA_StepErrorControl::registerWIMPtr(
  N_TIA_WorkingIntegrationMethod * wimPtr)
{
  wimPtr_ = wimPtr;
  return (wimPtr != NULL);
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::isPauseTime
// Purpose       : public method to flag whether the current time is a pause
//                 time
// Special Notes : Need this method because the prior values get in the way 
//                 when we resume.
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/30/04
//-----------------------------------------------------------------------------
inline bool N_TIA_StepErrorControl::isPauseTime()
{
  // Check whether pauseTime == currentTime (to within bpTol).
  // If the pause time is equal to the final time, we are at the end, not
  // at a pause time.  
#ifdef Xyce_DEBUG_TIME
  if (tiaParams_.debugLevel >0)
  {
    cout << "N_TIA_StepErrorControl::isPauseTime\n";
    cout << "currentPauseBP.value = " << currentPauseBP->value () << endl;
    cout << "final Time           = " << finalTime  << endl;
    cout << "current Time         = " << currentTime  << endl;
  }
#endif 
  return ( !(*currentPauseBP == finalTime) &&
           *currentPauseBP == currentTime );
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : "<<" operator for step error control class.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/17/05
//-----------------------------------------------------------------------------
inline ostream & operator<<(ostream & os, const N_TIA_StepErrorControl & sec)
{
  os << "\n\n-----------------------------------------\n";
  os << "\tStepErrorControl:\n";
  os << "\t\tstartingTimeStep      = " << sec.startingTimeStep << endl; 
  os << "\t\tcurrentTimeStep       = " << sec.currentTimeStep << endl; 
  os << "\t\tlastAttemptedTimeStep = " << sec.lastAttemptedTimeStep << endl; 
  os << "\t\tlastTimeStep          = " << sec.lastTimeStep << endl; 
  os << "\t\tminTimeStep           = " << sec.minTimeStep << endl; 
  os << "\t\tmaxTimeStep           = " << sec.maxTimeStep << endl; 
  os << "\t\tmaxTimeStepUser       = " << sec.maxTimeStepUser << endl; 
  os << "\t\tlastTime              = " << sec.lastTime << endl; 
  os << "\t\tcurrentTime           = " << sec.currentTime << endl; 
  os << "\t\tnextTime              = " << sec.nextTime << endl; 
  os << "\t\tstopTime              = " << sec.stopTime << endl; 
  os << "\t\tinitialTime           = " << sec.initialTime << endl; 
  os << "\t\tfinalTime             = " << sec.finalTime << endl; 
  os << endl;
  os << "\t\tBreak Points:" << endl;

  sec.printBreakPoints (os);


  os << "-----------------------------------------\n";
  os << endl;

  return os;
}

#endif // Xyce_N_TIA_STEP_ERROR_CONTROL_H


