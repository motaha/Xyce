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
// Filename      : $RCSfile: N_TIA_BackwardDifferentiation15.C,v $
//
// Purpose       : This file contains the functions which define the
//		             backward differentiation, order 1-5, class.
//
// Special Notes :
//
// Creator       : Todd Coffey
//
// Creation Date : 2/16/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.134.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:49 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

#ifdef HAVE_IOSTREAM
#include <iostream>
#else
#include <iostream.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_BackwardDifferentiation15.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_TIAParams.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_System.h>

#include <N_PDS_Comm.h>
#include <N_PDS_Manager.h>

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::
//                                               N_TIA_BackwardDifferentiation15
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 2/16/04
//-----------------------------------------------------------------------------

N_TIA_BackwardDifferentiation15::N_TIA_BackwardDifferentiation15
  (N_TIA_TIAParams & tiaP,
   N_TIA_StepErrorControl & secTmp,
   N_TIA_DataStore & dsTmp)
: N_TIA_TimeIntegrationMethod(tiaP,secTmp,dsTmp)
{
  leadingCoeff = 1;
  sec.maxOrder_=(Xycemin(5,tiaParams.maxOrder));
  sec.minOrder_=(Xycemax(1,tiaParams.minOrder));
  if (sec.minOrder_ > sec.maxOrder_) 
  {
    sec.minOrder_ = sec.maxOrder_;
  }
  timept_ = -1.0;
  return ;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::obtainPredictor
// Purpose       : Calculate predictor 
// Special Notes : stored in ds.xn0Ptr,qn0Ptr,qpn0Ptr
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 2/16/04
//-----------------------------------------------------------------------------

void N_TIA_BackwardDifferentiation15::obtainPredictor()
{
  // prepare history array for prediction
  for (int i=sec.nscsco_;i<=sec.currentOrder_;++i)
  {
    (ds.xHistory[i])->scale(sec.beta_[i]);
    (ds.qHistory[i])->scale(sec.beta_[i]);
    (ds.sHistory[i])->scale(sec.beta_[i]);
    (ds.stoHistory[i])->scale(sec.beta_[i]);
    (ds.stoLeadCurrQCompHistory[i])->scale(sec.beta_[i]);
  }
  
  // evaluate predictor
  *ds.xn0Ptr = *(ds.xHistory[0]);
  *ds.qn0Ptr = *(ds.qHistory[0]);
  *ds.sn0Ptr = *(ds.sHistory[0]);
  *ds.ston0Ptr = *(ds.stoHistory[0]);
  *ds.stoQCn0Ptr = *(ds.stoLeadCurrQCompHistory[0]);
  
  ds.qpn0Ptr->putScalar(0.0);
  ds.spn0Ptr->putScalar(0.0);
  ds.stopn0Ptr->putScalar(0.0);
  ds.stoQCpn0Ptr->putScalar(0.0);
  for (int i=1;i<=sec.currentOrder_;++i)
  {
    ds.xn0Ptr->linearCombo(1.0,*(ds.xHistory[i]),1.0,*ds.xn0Ptr);
    ds.qn0Ptr->linearCombo(1.0,*(ds.qHistory[i]),1.0,*ds.qn0Ptr);
    ds.sn0Ptr->linearCombo(1.0,*(ds.sHistory[i]),1.0,*ds.sn0Ptr);
    ds.ston0Ptr->linearCombo(1.0,*(ds.stoHistory[i]),1.0,*ds.ston0Ptr);
    ds.stoQCn0Ptr->linearCombo(1.0,*(ds.stoLeadCurrQCompHistory[i]),1.0,*ds.stoQCn0Ptr);
    ds.qpn0Ptr->linearCombo(sec.gamma_[i],*(ds.qHistory[i]),1.0,*ds.qpn0Ptr);
    ds.spn0Ptr->linearCombo(sec.gamma_[i],*(ds.sHistory[i]),1.0,*ds.spn0Ptr);
    ds.stopn0Ptr->linearCombo(sec.gamma_[i],*(ds.stoHistory[i]),1.0,*ds.stopn0Ptr);
    ds.stoQCpn0Ptr->linearCombo(sec.gamma_[i],*(ds.stoLeadCurrQCompHistory[i]),1.0,*ds.stoQCpn0Ptr);
    
  }

#ifdef Xyce_DEBUG_TIME
  cout.width(21); cout.precision(13); cout.setf(ios::scientific);

  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 1)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::obtainPredictor");
    cout << "\n currentOrder = " << sec.currentOrder_ << endl;
    cout << "\n sec.nscsco_: " << sec.nscsco_ << endl;
    for (int i=0; i<=sec.currentOrder_ ; ++i)
      cout << "\n sec.beta_[" << i << "] = " << sec.beta_[i] << "\n" << endl;
    for (int i=0; i<=sec.currentOrder_ ; ++i)
    {
      cout << "\n xHistory["<< i << "]: \n" << endl;
      (ds.xHistory[i])->printPetraObject();
      cout << endl;
    }
    for (int i=0; i<=sec.currentOrder_ ; ++i)
    {
      cout << "\n qHistory["<< i << "]: \n" << endl;
      (ds.qHistory[i])->printPetraObject();
      cout << endl;
    }
    for (int i=0; i<=sec.currentOrder_ ; ++i)
    {
      cout << "\n sHistory["<< i << "]: \n" << endl;
      (ds.sHistory[i])->printPetraObject();
      cout << endl;
    }
    cout << "\n xn0: \n" << endl;
    ds.xn0Ptr->printPetraObject();
    cout << endl;
    cout << "\n qn0: \n" << endl;
    ds.qn0Ptr->printPetraObject();
    cout << endl;
    cout << "\n qpn0: \n" << endl;
    ds.qpn0Ptr->printPetraObject();
    cout << endl;
    cout << "\n sn0: \n" << endl;
    ds.sn0Ptr->printPetraObject();
    cout << endl;
    cout << "\n spn0: \n" << endl;
    ds.spn0Ptr->printPetraObject();
    cout << endl;
    cout << "\n stoQCn0Ptr: \n" << endl;
    ds.stoQCn0Ptr->printPetraObject();
    cout << endl;
    cout << "\n stoQCpn0Ptr: \n" << endl;
    ds.stoQCpn0Ptr->printPetraObject();
    cout << endl;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME
  // copy the prediction into the next solution:
  *(ds.nextSolutionPtr) = *(ds.xn0Ptr);

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::obtainResidual
// Purpose       : Calculate Residual
// Special Notes : 
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 3/8/04
//-----------------------------------------------------------------------------

void N_TIA_BackwardDifferentiation15::obtainResidual()
{
  // output: ds.RHSVectorPtr

  // This function returns the following residual:
  // $qpn0 - (sec.alphas_/hn)(Q(x)-qn0)+F(x)-B(t)$

  // Note:  ds.nextSolutionPtr is used to get Q,F,B in N_ANP_AnalysisManager::loadRHS.
  ds.RHSVectorPtr->linearCombo(1.0,*ds.daeQVectorPtr,-1.0,*ds.qn0Ptr);
  ds.RHSVectorPtr->linearCombo(1.0,*ds.qpn0Ptr,-sec.alphas_/sec.currentTimeStep,*ds.RHSVectorPtr);

#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 1)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::obtainResidual");
    cout << "\n t = " << sec.nextTime << "\n" << endl;
    cout << "\n solution: \n" << endl;
    ds.nextSolutionPtr->printPetraObject();
    cout << "\n daeQVector: \n" << endl;
    ds.daeQVectorPtr->printPetraObject();
    cout << "\n qn0: \n" << endl;
    ds.qn0Ptr->printPetraObject();
    cout << "\n qpn0: \n" << endl;
    ds.qpn0Ptr->printPetraObject();
    cout << "\n sec.alphas_/hn: " << sec.alphas_/sec.currentTimeStep << "\n" << endl;
    cout << "\n daeFVector: \n" << endl;
    ds.daeFVectorPtr->printPetraObject();

    cout << "\n dQdt-vector: \n" << endl;
    ds.RHSVectorPtr->printPetraObject();
    cout << endl;
  }
#endif

  ds.RHSVectorPtr->linearCombo(1.0,*ds.RHSVectorPtr,+1.0,*ds.daeFVectorPtr);

  // since the nonlinear solver is expecting a -f, scale by -1.0:
  ds.RHSVectorPtr->scale(-1.0);

  // if voltage limiting is on, add it in:
  if (ds.limiterFlag)
  {
    (ds.dQdxdVpVectorPtr)->scale( -sec.alphas_/sec.currentTimeStep );

    (ds.RHSVectorPtr)->daxpy(
      *(ds.RHSVectorPtr), +1.0, *(ds.dQdxdVpVectorPtr));

    (ds.RHSVectorPtr)->daxpy(
      *(ds.RHSVectorPtr), +1.0, *(ds.dFdxdVpVectorPtr));
  }

#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 1)
  {
    cout << "\n Residual-vector: \n" << endl;
    cout << "-(qpn0-(sec.alpha_s/h)*(Q-qn0)+F-B) \n" << endl;
    ds.RHSVectorPtr->printPetraObject();
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    cout << endl;
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::obtainJacobian
// Purpose       : Calculate Jacobian
// Special Notes : 
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 3/8/04
//-----------------------------------------------------------------------------

void N_TIA_BackwardDifferentiation15::obtainJacobian()
{

#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 1)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::obtainJacobian");
  }
#endif
  // output: ds.JMatrixPtr

  // This function returns the following matrix:
  // $-(sec.alphas_/hn)dQdx(x)+dFdx$

  // Note:  ds.nextSolutionPtr is used to get dQdx,dFdx in N_ANP_AnalysisManager::loadJacobian.
//  ds.JmatrixPtr->linearCombo(-sec.alphas_/sec.currentTimeStep,*ds.dQdxMatrixPtr,+1.0,*ds.dFdxMatrixPtr);

  N_LAS_Matrix & dQdx = *(ds.dQdxMatrixPtr);
  N_LAS_Matrix & dFdx = *(ds.dFdxMatrixPtr);
  N_LAS_Matrix & Jac = *(ds.JMatrixPtr);

  double qscalar(-sec.alphas_/sec.currentTimeStep);

  Jac.linearCombo( qscalar, dQdx, 1.0, dFdx );

#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 1)
  {
    cout << "\n dFdx:" <<endl;
    dFdx.printPetraObject();
    cout << "\n Total Jacobian:" <<endl;
    Jac.printPetraObject();
//    for (int i=0;i<3;++i)
//    {
//      printf("[ %25.20g\t%25.20g\t%25.20g ]\n",Jac[i][0],Jac[i][1],Jac[i][2]);
//    }

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    cout << endl;
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::interpolateSolution
// Purpose       : Interpolate solution approximation at prescribed time point.
// Special Notes : This routine computes the solution at the output 
//               : timepoint by intepolation of the history using the order
//               : used for the most recent completed step, orderUsed.
//               : The output is put into ds.tmpSolVectorPtr.
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 2/16/04
//-----------------------------------------------------------------------------

bool N_TIA_BackwardDifferentiation15::interpolateSolution(double timepoint, 
    	N_LAS_Vector * tmpSolVectorPtr, vector <N_LAS_Vector*> & historyVec)
{
  // 03/23/04 tscoffe:  Currently this code is nearly identical to the IDA code
  // for interpolating to an output time.  Either we acknowledge the copyright,
  // the list of conditions in the license and the disclaimer or we rewrite this
  // function.  The IDA license is included after this routine.
  double tfuzz;   // fuzz factor to check for valid output time
  double tp;      // approximately t_{n-1}
  double delt;    // distance between timepoint and currentTime
  double c = 1.0; // coefficient for interpolation
  double gam;     // coefficient for interpolation
  int kord;       // order of interpolation
  double tn = sec.currentTime;
  double hh = sec.currentTimeStep;
  double hused = sec.usedStep_;
  int kused = sec.usedOrder_;
  // TVR: changed line below in addressing bug 1258
  // This is what was here prior to my change:
  //  double uround = 0.0;  // unit round-off (set to zero for now)
  // The value below is a WAG
  double uround = N_UTL_MachineDependentParams::MachinePrecision();  

  tfuzz = 100 * uround * (tn + hh);
  tp = tn - hused - tfuzz;

  if ( (timepoint - tp)*hh < 0.0 ) 
    return false;

  *tmpSolVectorPtr = *(historyVec[0]);
  kord = kused;
  if ( (kused == 0) || (timepoint == tn) ) 
    kord = 1;

  delt = timepoint - tn;
  gam = delt/sec.psi_[0];
  for (int j=1 ; j <= kord ; ++j)
  {
    c = c*gam;
    gam = (delt + sec.psi_[j-1])/sec.psi_[j];
    tmpSolVectorPtr->linearCombo(1.0,*tmpSolVectorPtr,c,*(historyVec[j]));
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::interpolateMPDESolution
// Purpose       : Interpolate solution approximation at prescribed time points.
// Special Notes : This routine computes the solution at the output 
//               : timepoints by intepolation of the history using the order
//               : used for the most recent completed step, orderUsed.
//               : The output is put into provided N_LAS_Vector pointer.
//               : The interpolation is as follows:
//               : tmpSolVectorPtr->block(i) is interpolated at timepoint(i)
//               : Therefore, if you want them all interpolated at the same time, 
//               : then use timepoint(i) = timepoint(0) forall i
//               : or use interpolateSolution. 
// Scope         : public
// Creator       : Todd Coffey, Eric Keiter, SNL 
// Creation Date : 11/28/06
//-----------------------------------------------------------------------------

bool N_TIA_BackwardDifferentiation15::interpolateMPDESolution(std::vector<double>& timepoint, 
    	N_LAS_Vector * tmpSolVectorPtr)
{
#ifdef Xyce_PARALLEL_MPI
  string msg = "N_TIA_BackwardDifferentiation15::interpolateMPDESolution: ";
  msg += "Not set up for Parallel yet";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return(false);
#endif

  N_LAS_BlockVector * blockTempSolVectorPtr = 
     dynamic_cast<N_LAS_BlockVector*>(tmpSolVectorPtr);
  if (blockTempSolVectorPtr == NULL)
  {
    string msg = "N_TIA_BackwardDifferentiation15::interpolateMPDESolution: ";
    msg += "N_LAS_Vector tmpSolVectorPtr is not of type N_LAS_BlockVector";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    return(false);
  }

  double tfuzz;   // fuzz factor to check for valid output time
  double tp;      // approximately t_{n-1}
  int numblocks = timepoint.size();
  int blockCount = blockTempSolVectorPtr->blockCount();
  if (numblocks > blockCount)
  {
    string msg = "N_TIA_BackwardDifferentiation15::interpolateMPDESolution: ";
    msg += "Number of time points requested is greater than number of fast time points in MPDE block vector";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    return(false);
  }
  double delt;  
  double c = 1.0;
  double gam;  
  int kord;       // order of interpolation
  double tn = sec.currentTime;
  double hh = sec.currentTimeStep;
  double hused = sec.usedStep_;
  int kused = sec.usedOrder_;
  double uround = 0.0;  // unit round-off (set to zero for now)

  tfuzz = 100 * uround * (tn + hh);
  tp = tn - hused - tfuzz;
  for (int i=0; i<numblocks ; ++i)
  {
    if ( (timepoint[i] - tp)*hh < 0.0 ) 
      return false;
  }

  *tmpSolVectorPtr = *(ds.xHistory[0]);

  N_LAS_Vector * solVectorPtr;
  N_LAS_Vector * xHistoryVectorPtr;
  // Loop over blocks first so that maximal order can be maintained
  for (int i=0; i < numblocks ; ++i)
  {
    if ((kused == 0) || (timepoint[i] == tn)) { kord = 1; }
    else { kord = kused; }
    solVectorPtr = &(blockTempSolVectorPtr->block(i));
    c = 1.0;
    delt = timepoint[i] - tn;
    gam = delt/sec.psi_[0];
    for (int j=1 ; j <= kord ; ++j)
    {
      c = c*gam;
      gam = (delt + sec.psi_[j-1])/sec.psi_[j];
      N_LAS_BlockVector * blockXHistoryVectorPtr = 
        dynamic_cast<N_LAS_BlockVector*>(ds.xHistory[j]);
      if (blockXHistoryVectorPtr == NULL)
      {
        string msg = "N_TIA_BackwardDifferentiation15::interpolateMPDESolution: ";
        msg += "N_LAS_Vector ds.xHistory[j] is not of type N_LAS_BlockVector\n j = ";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg, j);
        return(false);
      }
      xHistoryVectorPtr = &(blockXHistoryVectorPtr->block(i));
      solVectorPtr->linearCombo(1.0,*solVectorPtr,c,*xHistoryVectorPtr);
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::printMPDEOutputSolution()
// Purpose       : Print transient output from MPDE simulation
// Special Notes : This routine uses interpolateMPDESolution.
// Scope         : public
// Creator       : Todd Coffey, SNL, 1414
// Creation Date : 11/28/06
//-----------------------------------------------------------------------------
bool N_TIA_BackwardDifferentiation15::printMPDEOutputSolution(
        RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
        const double time,
        N_LAS_Vector * solnVecPtr,
        const std::vector<double> & fastTimes )
{
#ifdef Xyce_DEBUG_TIME
  const string dashedline =
	"---------------------------------------------------------------"
	"-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
      "  N_TIA_BackwardDifferentiation15::printMPDEOutputSolution");
  }
#endif // Xyce_DEBUG_TIME
  double timestep = sec.lastAttemptedTimeStep;
  double lasttime = sec.currentTime - timestep;
  double tn = sec.currentTime;
  // Set these values up to read output time intervals.  FIXME
  double beg_of_output_time_interval = lasttime;
  double end_of_output_time_interval = tn;
  double start_time = max(lasttime,beg_of_output_time_interval);
  double stop_time = min(tn,end_of_output_time_interval);
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 0)
  {
    cout << "timestep = " << timestep << endl;
    cout << "lasttime = " << lasttime << endl;
    cout << "tn = " << tn << endl;
    cout << "beg_of_output_time_interval = " << beg_of_output_time_interval << endl;
    cout << "end_of_output_time_interval = " << end_of_output_time_interval << endl;
    cout << "start_time = " << start_time << endl;
    cout << "stop_time = " << stop_time << endl;
  }
#endif // Xyce_DEBUG_TIME


  N_LAS_BlockVector * blockTmpSolVectorPtr = 
    dynamic_cast<N_LAS_BlockVector*>(ds.tmpSolVectorPtr);
  
  if (blockTmpSolVectorPtr == NULL)
  {
    string msg = "N_TIA_BackwardDifferentiation15::printMPDEOutputSolution: ";
    msg += "N_LAS_Vector ds.tmpSolVectorPtr is not of type N_LAS_BlockVector";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    return(false);
  }
  int blockCount = blockTmpSolVectorPtr->blockCount();

  // Create list of timepoints to interpolate (along characteristic curve)
  double T2 = fastTimes.back();
  //double charcross = start_time - floor(start_time/T2)*T2; // (start_time mod T2)
  double charcross = fmod(start_time,T2); // (start_time mod T2)
  int s_ind_0 = -1;
  // find s_ind_0 = first fast time point >= charcross.
  // This could probably be made faster: FIXME
  for (int i=0 ; i<=blockCount ; ++i)
  {
    if (fastTimes[i] >= charcross)
    {
      s_ind_0 = i;
      break;
    }
  }
  if (s_ind_0 == -1)
  {
    string msg = "N_TIA_BackwardDifferentiation15::printMPDEOutputSolution: ";
    msg += "Cannot find where characteristic curve crosses fast time slice at start_time";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    return(false);
  }
  vector<double> h2(blockCount,0);
  for (int j=0 ; j < blockCount ; ++j)
  {
    h2[j] = fastTimes[j+1] - fastTimes[j];
  }
  vector<double> ti;
  //double first_interp = floor(start_time/T2)*T2 + fastTimes[s_ind_0];
  double first_interp = start_time - charcross + fastTimes[s_ind_0];
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 0)
  {
    cout << "first_interp = " << first_interp << endl;
  }
#endif // Xyce_DEBUG_TIME
  if (s_ind_0 == blockCount) { s_ind_0 = 0; };
  // Don't interpolate the first point 
  double eps = fabs(start_time)*1.0e-6; 
  if ( fabs(first_interp-timept_) <= eps )
  {
    first_interp += h2[s_ind_0];
    s_ind_0++;
    if (s_ind_0 == blockCount) { s_ind_0 = 0; };
#ifdef Xyce_DEBUG_TIME
    if (tiaParams.debugLevel > 0)
    {
      cout << "Moving first_interp forward to avoid duplicate outputs:  " << first_interp << endl;
    }
#endif // Xyce_DEBUG_TIME
  }
  int sn = s_ind_0;
  double t = first_interp;
  while (t <= stop_time)
  {
    ti.push_back(t);
    t += h2[sn];
    sn++;
    if (sn >= blockCount) { sn = 0; }
  }
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 0)
  {
    cout << "T2 = " << T2 << endl;
    cout << "charcross = " << charcross << endl;
    cout << "s_ind_0 = " << s_ind_0 << endl;
    cout << "Expecting to interpolate the following points:" << endl;
    unsigned int numinterp = ti.size();
    for (unsigned int i=0 ; i < numinterp ; ++i)
    {
      cout << ti[i] << endl;
    }
    cout << "Total of " << numinterp << " points" << endl;
  }
#endif // Xyce_DEBUG_TIME
  timept_ = start_time;  // used later for interpolating stop_time
  unsigned int tinum = ti.size();
  int total_interp = 0;
  std::vector<double> timepoint_vec(blockCount,stop_time);
  int num_interp_this_cycle = 0;
  int s_ind = s_ind_0;
  for (unsigned int i=0; i < tinum ; ++i)
  {
    timepoint_vec[s_ind] = ti[i];
    num_interp_this_cycle++;
    s_ind++;
    if (s_ind >= blockCount) { s_ind = 0; };
    // If we're back to s_ind_0 or out of ti points, then interpolate:
    if ((s_ind == s_ind_0) || (i == tinum-1))
    {
      interpolateMPDESolution(timepoint_vec, ds.tmpSolVectorPtr);
      // Now we print these points out
      int s = s_ind_0;
      for (int j=0 ; j < num_interp_this_cycle ; ++j)
      {
        timept_ = timepoint_vec[s];
        outputMgrAdapterRCPtr->tranOutput(timept_, blockTmpSolVectorPtr->block(s), *ds.tmpStaVectorPtr, *ds.tmpStoVectorPtr );
        total_interp++;
        s++;
        if (s >= blockCount) { s = 0; }
#ifdef Xyce_DEBUG_TIME
        cout << "Interpolated to t = " << timept_ << endl;
#endif // Xyce_DEBUG_TIME
      }
      num_interp_this_cycle = 0;
    }
  }
#ifdef Xyce_DEBUG_TIME
  cout << "Total of " << total_interp << " points" << endl;
#endif // Xyce_DEBUG_TIME

  // Now we interpolate stop_time unless its too close to the last timept interpolated.
  eps = fabs(stop_time)*1.0e-8; 
  // fudge factor for printing, this should come from elsewhere FIXME
  if (fabs(timept_ - stop_time) >= eps)
  {
#ifdef Xyce_DEBUG_TIME
    if (tiaParams.debugLevel > 0)
    {
      cout << "Previous timept = " << timept_ << endl;
      cout << "Expecting to interpolate the following point: " << stop_time << endl;
    }
#endif // Xyce_DEBUG_TIME
    N_LAS_Vector* tmpSolnVecPtr = solnVecPtr;
    N_LAS_Vector* tmpVecPtr = ds.tmpXn0APtr;
    if (stop_time < tn)
    {
      interpolateSolution(stop_time,ds.tmpXn0BPtr, ds.xHistory);
      tmpSolnVecPtr = ds.tmpXn0BPtr;
    }
    N_LAS_BlockVector * blockTmpSolnVecPtr = 
      dynamic_cast<N_LAS_BlockVector*>(tmpSolnVecPtr);
    N_LAS_BlockVector * blockTmpVecPtr = 
      dynamic_cast<N_LAS_BlockVector*>(tmpVecPtr);
    if (blockTmpSolnVecPtr == NULL)
    {
      string msg = "N_TIA_BackwardDifferentiation15::printMPDEOutputSolution: ";
      msg += "N_LAS_Vector tmpSolnVecPtr is not of type N_LAS_BlockVector";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
      return(false);
    }
    if (blockTmpVecPtr == NULL)
    {
      string msg = "N_TIA_BackwardDifferentiation15::printMPDEOutputSolution: ";
      msg += "N_LAS_Vector tmpVecPtr is not of type N_LAS_BlockVector";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
      return(false);
    }
    // Interpolate where the caracteristic crosses stop_time (instead of start_time).
    //charcross = stop_time - floor(stop_time/T2)*T2;
    charcross = fmod(stop_time,T2);
    int s_ind_1 = -1;
    // Find index just before charcross:
    if( charcross < fastTimes[0] )
    {
      // by periodicity, the one before in this case is the one at the end
      s_ind_1 = blockCount-1;
    }
    else
    {
      for (int i=blockCount-1 ; i>=0 ; --i)
      { 
        if (fastTimes[i] <= charcross)
        {
          s_ind_1 = i;
          break;
        }
      }
    }
    if (s_ind_1 == -1)
    {
      string msg = "N_TIA_BackwardDifferentiation15::printMPDEOutputSolution: ";
      msg += "Cannot find where characteristic curve crosses fast time slice at stop_time";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
      return(false);
    }
    int sm = s_ind_1;
    int sp = s_ind_1+1;
    double coeff_sm = fastTimes[sp]-charcross;
    double coeff_sp = charcross-fastTimes[sm];
    if (sp == blockCount) { sp = 0; }
    double dt = h2[s_ind_1];
    timept_ = stop_time;
#ifdef Xyce_DEBUG_TIME
    cout << "charcross = " << charcross << endl;
    cout << "s_ind_1 = " << s_ind_1 << endl;
    cout << "sp = " << sp << endl;
    cout << "sm = " << sm << endl;
    cout << "dt = " << dt << endl;
    cout << "timept = " << timept_ << endl;
    cout << "coeff_sm = " << coeff_sm << endl;
    cout << "coeff_sp = " << coeff_sp << endl;
#endif // Xyce_DEBUG_TIME
    blockTmpVecPtr->block(0).linearCombo(
        coeff_sm/dt, blockTmpSolnVecPtr->block(sm),
        coeff_sp/dt, blockTmpSolnVecPtr->block(sp)  
        );
    outputMgrAdapterRCPtr->tranOutput(timept_, blockTmpVecPtr->block(0), *ds.tmpStaVectorPtr, *ds.tmpStoVectorPtr );

#ifdef Xyce_DEBUG_TIME
    cout << "Interpolated to t = " << timept_ << endl;
#endif // Xyce_DEBUG_TIME
  }
#ifdef Xyce_DEBUG_TIME
  else // no extra interpolation
  {
    if (tiaParams.debugLevel > 0)
    {
      cout << "No further interpolation required." << endl;
    }
  }
#endif // Xyce_DEBUG_TIME

#ifdef Xyce_DEBUG_TIME
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
#endif // Xyce_DEBUG_TIME
  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::printWaMPDEOutputSolution()
// Purpose       : Print transient output from WaMPDE simulation
// Special Notes : This routine uses interpolateSolution.
// Scope         : public
// Creator       : Todd Coffey, SNL, 1414
// Creation Date : 12/15/06
//-----------------------------------------------------------------------------
bool N_TIA_BackwardDifferentiation15::printWaMPDEOutputSolution(
        RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
        const double time,
        N_LAS_Vector * solnVecPtr,
        const std::vector<double> & fastTimes,
        const int phiGID )
{
#ifdef Xyce_DEBUG_TIME
  const string dashedline =
	"---------------------------------------------------------------"
	"-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
      "  N_TIA_BackwardDifferentiation15::printWaMPDEOutputSolution");
  }
#endif // Xyce_DEBUG_TIME
  double timestep = sec.lastAttemptedTimeStep;
  double lasttime = sec.currentTime - timestep;
  double tn = sec.currentTime;
  // Set these values up to read output time intervals.  FIXME
  double beg_of_output_time_interval = lasttime;
  double end_of_output_time_interval = tn;
  double start_time = max(lasttime,beg_of_output_time_interval);
  double stop_time = min(tn,end_of_output_time_interval);
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 0)
  {
    cout << "start_time = " << start_time << endl;
    cout << "stop_time = " << stop_time << endl;
  }
#endif // Xyce_DEBUG_TIME
  
  // 12/15/06 tscoffe:
  // First break the interval up as printOutputSolution does:
  // hh = timestep/(sec.usedOrder_) and interpolate in the intervals:
  // [tn+(i-1)*hh,tn+i*hh], i=1..usedOrder
  // Assume phi(t_1) is linear in these intervals and approximate with:
  // phi(t) = (1/hh)*(phi(tn)(tn+hh-t)+phi(tn+hh)(t-tn))
  // Then the t_1 values we want to interpolate are:
  // n2 = number of fast time points.
  // T2 = fast time period
  // h2 = T2/n2 = average spacing on fast time scale
  // t1_vals = [tn:h2:tn+hh]
  // And the t_2 values we want to interpolate are:
  // t2_vals = phi(t1_vals) mod T2
  // Then take the N_LAS blocks and do 2D linear interpolation on the intervals:
  // (t1,s1), (t1,s2), (t2,s1), (t2,s2)
  // x(t) = x(t,s) approx = 
  // (1/(t2-t1))(1/(s2-s1))[  x(t1,s1)(t2-t)(s2-s)
  //                         +x(t1,s2)(t2-t)(s-s1)
  //                         +x(t2,s1)(t-t1)(s2-s)
  //                         +x(t2,s2)(t-t1)(s-s1) ]
  // where t = t1_vals and s = t2_vals

  N_LAS_BlockVector * blockTmpSolVectorPtr = 
    dynamic_cast<N_LAS_BlockVector*>(ds.tmpSolVectorPtr);
  N_LAS_BlockVector * blockTmpXn0APtr = 
    dynamic_cast<N_LAS_BlockVector*>(ds.tmpXn0APtr);
  N_LAS_BlockVector * blockTmpXn0BPtr = 
    dynamic_cast<N_LAS_BlockVector*>(ds.tmpXn0BPtr);
  if (blockTmpSolVectorPtr == NULL)
  {
    string msg = "N_TIA_BackwardDifferentiation15::printWaMPDEOutputSolution: ";
    msg += "N_LAS_Vector ds.tmpSolVectorPtr is not of type N_LAS_BlockVector";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    return(false);
  }
  if (blockTmpXn0APtr == NULL)
  {
    string msg = "N_TIA_BackwardDifferentiation15::printWaMPDEOutputSolution: ";
    msg += "N_LAS_Vector ds.tmpXn0APtr is not of type N_LAS_BlockVector";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    return(false);
  }
  if (blockTmpXn0BPtr == NULL)
  {
    string msg = "N_TIA_BackwardDifferentiation15::printWaMPDEOutputSolution: ";
    msg += "N_LAS_Vector ds.tmpXn0BPtr is not of type N_LAS_BlockVector";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    return(false);
  }
  int phiLID = blockTmpSolVectorPtr->pmap()->globalToLocalIndex(phiGID);
  double hh = timestep/(sec.usedOrder_);
  double timeA = -1.0;
  double timeB = -1.0;
#ifdef Xyce_DEBUG_TIME
    if (tiaParams.debugLevel > 0)
    {
      cout << " sec.usedOrder_ = " << sec.usedOrder_ << endl;
      cout << " sec.currentTime_ = " << sec.currentTime << endl;
      cout << " lasttime = " << lasttime << endl;
    }
#endif // Xyce_DEBUG_TIME

  for (int i=0 ; i < sec.usedOrder_ ; ++i)
  {
    if (i == 0)
    {
      bool junk;
      timeA = lasttime + hh*i;
      junk = interpolateSolution(timeA,ds.tmpXn0APtr, ds.xHistory);
      if (!junk) 
      {
        string msg = "interpolateSolution returned false!";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0,msg);
      }
    }
    else
    {
      // We don't need to interpolate this again.
      *ds.tmpXn0APtr = *ds.tmpXn0BPtr;
      timeA = timeB;
    }
    timeB = lasttime + hh*(i+1);
    interpolateSolution(timeB,ds.tmpXn0BPtr, ds.xHistory);
#ifdef Xyce_DEBUG_TIME
    if (tiaParams.debugLevel > 0)
    {
      cout << "Interpolating in [ " << timeA << ", " << timeB << " ]" << endl;
      cout << "timeA = " << timeA << endl;
      cout << "timeB = " << timeB << endl;
    }
#endif // Xyce_DEBUG_TIME

    // Now we can interpolate [tmpXn0APtr,tmpXn0BPtr] in [timeA,timeB].
    vector<double> t1vals;
    double T2 = fastTimes.back(); 
    int blockCount = blockTmpSolVectorPtr->blockCount();
    double h2 = T2/blockCount; // Average mesh width in fast time-scale
    double tval = timeA+h2;
    while (tval <= timeB)
    {
      t1vals.push_back(tval);
      tval += h2;
    }
    // fudge factor for printing, this should come from elsewhere FIXME
    double eps = fabs(timeB)*1.0e-8; 
    if ( (t1vals.size() == 0) || (fabs(t1vals.back() - timeB) >= eps) )
    {
      t1vals.push_back(timeB);
    }
    vector<double> t2vals;
    double phiA, phiB;
    phiA = (*ds.tmpXn0APtr)[phiLID]; // Get from MPDE Manager
    phiB = (*ds.tmpXn0BPtr)[phiLID]; // Get from MPDE Manager
    for (unsigned int j=0 ; j<t1vals.size() ; ++j)
    {
      double phi = (1/(timeB-timeA))*(phiA*(timeB-t1vals[j])+phiB*(t1vals[j]-timeA));
      t2vals.push_back(fmod(phi,T2)); // phi(t1vals[j]) mod T2
    }
#ifdef Xyce_DEBUG_TIME
    if (tiaParams.debugLevel > 0)
    {
      cout << "t1vals = " << endl;
      for (unsigned int j=0 ; j < t1vals.size() ; ++j)
      {
        cout << t1vals[j] << endl;
      }
      cout << "phi(" << timeA << ") = " << phiA << endl;
      cout << "phi(" << timeB << ") = " << phiB << endl;
      cout << "t2vals = " << endl;
      for (unsigned int j=0 ; j< t2vals.size() ; ++j)
      {
        cout << t2vals[j] << endl;
      }
    }
#endif // Xyce_DEBUG_TIME
    // Now we can do our block 2D interpolations
    // loop through t1vals and move the fast time blocks as we go
    double t1 = timeA; // slow time indices
    double t2 = timeB;
    double t = t1vals[0]; // current time indices
    double s = t2vals[0];
    int b1,b2; // fast time block indices corresponding to s1,s2
    // Find the block surrounding s:
    b1 = -2;
    for (int j=0 ; j < blockCount ; ++j)
    {
      if ((fastTimes[j] <= s) && (s < fastTimes[j+1]))
      {
        b1 = j;
      }
    }
    b2 = b1+1;
    if (b2 == blockCount)
    {
      b2 = 0;
    }
    double s1 = fastTimes[b1];
    double s2 = fastTimes[b1+1]; // Note:  fastTimes[blockCount] = T2
    if ((s < s1) || (s > s2)) 
    {
      string msg = "N_TIA_BackwardDifferentiation15::printWaMPDEOutputSolution: ";
      msg += "  Interpolator cannot find a fast time block containing the first point  ";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    }
#ifdef Xyce_DEBUG_TIME
    cout << "Found s = " << s << " in block " << b1; 
    cout << " with boundary = [" << s1 << "," << s2 << "]" << endl;
#endif // Xyce_DEBUG_TIME
    for (unsigned int j=0 ; j < t1vals.size() ; ++j)
    {
      t = t1vals[j];
      s = t2vals[j];
      if (t > t2) break; // This should never occur
      // If s has moved outside our block, then increment block.
      if ( (s < s1) || (s > s2) ) 
      {
#ifdef Xyce_DEBUG_TIME
        cout << "Incrementing fast time block for next interpolation." << endl;
#endif // Xyce_DEBUG_TIME
        b1++;
        if (b1 == blockCount)
        {
          b1 = 0;
        }
        b2 = b1+1;
        if (b2 == blockCount)
        {
          b2 = 0;
        }
        s1 = fastTimes[b1];
        s2 = fastTimes[b1+1];
      }
      // If s isn't in the next block, then search for it.
      if ((s < s1) || (s > s2)) 
      {
#ifdef Xyce_DEBUG_TIME
        cout << "Searching for fast time block for next interpolation." << endl;
#endif // Xyce_DEBUG_TIME
        b1 = -2;
        for (int j2=0 ; j2 < blockCount ; ++j2)
        {
          if ((fastTimes[j2] <= s) && (s < fastTimes[j2+1]))
          {
            b1 = j2;
          }
        }
        b2 = b1+1;
        if (b2 == blockCount)
        {
          b2 = 0;
        }
        s1 = fastTimes[b1];
        s2 = fastTimes[b1+1]; 
      }
      // If a block surrounding s can't be found, then quit.
      if ((s < s1) || (s > s2)) 
      {
        string msg = "N_TIA_BackwardDifferentiation15::printWaMPDEOutputSolution: ";
        msg += "  Interpolator moved fast time block but new point is not in this block  ";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
      }
      if (s > T2) break; // Just double checking...
      //blockTmpXn0APtr->block(b1) // x(t1,s1)
      //blockTmpXn0APtr->block(b2) // x(t1,s2)
      //blockTmpXn0BPtr->block(b1) // x(t2,s1)
      //blockTmpXn0BPtr->block(b2) // x(t2,s2)
      // (1/(t2-t1))(1/(s2-s1))[  x(t1,s1)(t2-t)(s2-s)
      //                         +x(t1,s2)(t2-t)(s-s1)
      //                         +x(t2,s1)(t-t1)(s2-s)
      //                         +x(t2,s2)(t-t1)(s-s1) ]
#ifdef Xyce_DEBUG_TIME
      if (tiaParams.debugLevel > 0)
      {
        cout << "Interpolating in the block:" << endl;
        cout << "(t1,t2) = (" << t1 << "," << t2 << ")" << endl;
        cout << "(s1,s2) = (" << s1 << "," << s2 << ")" << endl;
      }
#endif // Xyce_DEBUG_TIME
      double denom = (t2-t1)*(s2-s1);
      double coeff0 = (t2-t)*(s2-s)/denom;
      double coeff1 = (t2-t)*(s-s1)/denom;
      double coeff2 = (t-t1)*(s2-s)/denom;
      double coeff3 = (t-t1)*(s-s1)/denom;
      (blockTmpSolVectorPtr->block(b1)).linearCombo(
        coeff0, blockTmpXn0APtr->block(b1),
        coeff1, blockTmpXn0APtr->block(b2),
        coeff2, blockTmpXn0BPtr->block(b1)  );
      (blockTmpSolVectorPtr->block(b1)).update(
        coeff3, blockTmpXn0BPtr->block(b2), 1.0  );

      // erkeite 2/24/07. This is needed because currently the interpolation goes back to t=-1.0.
      if (t >= 0.0) 
      {
        outputMgrAdapterRCPtr->tranOutput(t, blockTmpSolVectorPtr->block(b1), *ds.tmpStaVectorPtr, *ds.tmpStoVectorPtr );

#ifdef Xyce_DEBUG_TIME
        cout << "Interpolated to (t,phi(t)) = (" << t << "," << s << ")" << endl;
#endif // Xyce_DEBUG_TIME
      }
    }
  }
  
#ifdef Xyce_DEBUG_TIME
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
#endif // Xyce_DEBUG_TIME
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::printOutputSolution()
// Purpose       : Print output that is dumbed down in order.
// Special Notes : This routine picks smaller time steps to approximate first
//               : order integration from the perspective of the output.
// Scope         : public
// Creator       : Todd Coffey, SNL, 1414
// Creation Date : 11/21/05
//-----------------------------------------------------------------------------
bool N_TIA_BackwardDifferentiation15::printOutputSolution(
        RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
        const double time,
        N_LAS_Vector * solnVecPtr,
        const bool doNotInterpolate,
        const vector<double> &outputInterpolationTimes,
        bool skipPrintLineOutput)
{
#ifdef Xyce_DEBUG_TIME
  const string dashedline =
	"---------------------------------------------------------------"
	"-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
      "  N_TIA_BackwardDifferentiation15::printOutputSolution");
    cout << "usedOrder_ = " << sec.usedOrder_ << endl;
  }
#endif // Xyce_DEBUG_TIME
  double timestep = sec.lastAttemptedTimeStep;
  double lasttime = sec.currentTime - timestep;
  bool dointerp = true;
  double hh = timestep/(sec.usedOrder_);
  if (hh <= 10*sec.minTimeStep)
  {
    dointerp = false;
  }

  if (!(tiaParams.interpOutputFlag))
  {
    dointerp = false;
  }

  if (doNotInterpolate)
  {
    dointerp = false;
  }

  if (dointerp && !outputInterpolationTimes.empty())
  {
    for (unsigned int i=0;i<outputInterpolationTimes.size();++i)
    {
      interpolateSolution(outputInterpolationTimes[i], ds.tmpSolVectorPtr, ds.xHistory);    // interpolate solution vector
      interpolateSolution(outputInterpolationTimes[i], ds.tmpStaVectorPtr, ds.sHistory);    // interpolate state vector
      interpolateSolution(outputInterpolationTimes[i], ds.tmpStoVectorPtr, ds.stoHistory);  // interpolate store vector
      outputMgrAdapterRCPtr->tranOutput(outputInterpolationTimes[i], *ds.tmpSolVectorPtr, 
        *ds.tmpStaVectorPtr, *ds.tmpStoVectorPtr, skipPrintLineOutput);
    }
  }
  else if ( (sec.usedOrder_ > 2) && dointerp )
  {
    for (int i=1 ; i<sec.usedOrder_; ++i)
    {
      double timept = lasttime + hh*i;
      interpolateSolution(timept,ds.tmpSolVectorPtr, ds.xHistory);    // interpolate solution vector
      interpolateSolution(timept,ds.tmpStaVectorPtr, ds.sHistory);    // interpolate state vector
      interpolateSolution(timept,ds.tmpStoVectorPtr, ds.stoHistory);  // interpolate store vector
      outputMgrAdapterRCPtr->tranOutput(timept, *ds.tmpSolVectorPtr, *ds.tmpStaVectorPtr, *ds.tmpStoVectorPtr, skipPrintLineOutput);

#ifdef Xyce_DEBUG_TIME
      cout << "Interpolated to t = " << timept << endl;
#endif // Xyce_DEBUG_TIME
    }
  }

  // Either way, do an output on the actual computed time step, but only
  // if we weren't given a list of specific times *or* we were told not to
  // interpoloate.
  if (outputInterpolationTimes.empty() || doNotInterpolate)
    outputMgrAdapterRCPtr->tranOutput(time, *ds.currSolutionPtr, *ds.currStatePtr, *ds.currStorePtr, skipPrintLineOutput);

#ifdef Xyce_DEBUG_TIME
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
#endif // Xyce_DEBUG_TIME
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::saveOutputSolution
// Purpose       : This is similar to printOutputSolution, but is in support of
//                 the .SAVE capability, rather than .PRINT.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
bool N_TIA_BackwardDifferentiation15::saveOutputSolution  (
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  N_LAS_Vector * solnVecPtr,
                  const double saveTime,
                  const bool doNotInterpolate)
{
#ifdef Xyce_DEBUG_TIME
  const string dashedline =
	"---------------------------------------------------------------"
	"-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
      "  N_TIA_BackwardDifferentiation15::saveOutputSolution");
  }
#endif // Xyce_DEBUG_TIME

  double timestep = sec.lastAttemptedTimeStep;
  double lasttime = sec.currentTime - timestep;
  bool dointerp = true;
  double hh = timestep/(sec.usedOrder_);
  if (hh <= 10*sec.minTimeStep)
  {
    dointerp = false;
  }

  if (!(tiaParams.interpOutputFlag))
  {
    dointerp = false;
  }

  if (doNotInterpolate)
  {
    dointerp = false;
  }

  if (dointerp)
  {
    interpolateSolution(saveTime,ds.tmpSolVectorPtr, ds.xHistory);
    outputMgrAdapterRCPtr->outputDCOP( *(ds.tmpSolVectorPtr) );
  }
  else
  {
    outputMgrAdapterRCPtr->outputDCOP( *(solnVecPtr) );
  }

#ifdef Xyce_DEBUG_TIME
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
#endif // Xyce_DEBUG_TIME
  return true;
}

//-----------------------------------------------------------------------------
// IDA source license:
//-----------------------------------------------------------------------------
// Copyright (c) 2002, The Regents of the University of California.
// Produced at the Lawrence Livermore National Laboratory.
// Written by Alan Hindmarsh, Allan Taylor, Radu Serban.
// UCRL-CODE-2002-59
// All rights reserved.
//                                                                                       
// This file is part of IDA.
//  
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//  
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the disclaimer below.
//  
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the disclaimer (as noted below)
// in the documentation and/or other materials provided with the
// distribution.
//  
// 3. Neither the name of the UC/LLNL nor the names of its contributors
// may be used to endorse or promote products derived from this software
// without specific prior written permission.
//  
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// REGENTS OF THE UNIVERSITY OF CALIFORNIA, THE U.S. DEPARTMENT OF ENERGY
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  
// Additional BSD Notice
// ---------------------
// 1. This notice is required to be provided under our contract with
// the U.S. Department of Energy (DOE). This work was produced at the
// University of California, Lawrence Livermore National Laboratory
// under Contract No. W-7405-ENG-48 with the DOE.
//  
// 2. Neither the United States Government nor the University of
// California nor any of their employees, makes any warranty, express
// or implied, or assumes any liability or responsibility for the
// accuracy, completeness, or usefulness of any information, apparatus,
// product, or process disclosed, or represents that its use would not
// infringe privately-owned rights.
//  
// 3. Also, reference herein to any specific commercial products,
// process, or services by trade name, trademark, manufacturer or
// otherwise does not necessarily constitute or imply its endorsement,
// recommendation, or favoring by the United States Government or the
// University of California. The views and opinions of authors expressed
// herein do not necessarily state or reflect those of the United States
// Government or the University of California, and shall not be used for
// advertising or product endorsement purposes.
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::updateHistory
// Purpose       : Update history array after a successful step 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 03/08/04
//-----------------------------------------------------------------------------
void N_TIA_BackwardDifferentiation15::updateHistory()
{

#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 1)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::updateHistory");
    cout << "\n Before updates \n" << endl;
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      cout << "\n xHistory["<< i << "]: \n" << endl;
      (ds.xHistory[i])->printPetraObject();
      cout << endl;
    }
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      cout << "\n qHistory["<< i << "]: \n" << endl;
      (ds.qHistory[i])->printPetraObject();
      cout << endl;
    }
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      cout << "\n sHistory["<< i << "]: \n" << endl;
      (ds.sHistory[i])->printPetraObject();
      cout << endl;
    }
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME

  // Save Newton correction for potential order increase on next step.
  if (sec.usedOrder_ < sec.maxOrder_)  
  {
    *(ds.xHistory[sec.usedOrder_+1]) = *ds.newtonCorrectionPtr;
    *(ds.qHistory[sec.usedOrder_+1]) = *ds.qNewtonCorrectionPtr;
  }
  // Update history arrays
  (ds.xHistory[sec.usedOrder_])->linearCombo(1.0,*(ds.xHistory[sec.usedOrder_]),1.0,*ds.newtonCorrectionPtr);
  (ds.qHistory[sec.usedOrder_])->linearCombo(1.0,*(ds.qHistory[sec.usedOrder_]),1.0,*ds.qNewtonCorrectionPtr);
  for (int j=sec.usedOrder_-1;j>=0;j--) 
  {
    (ds.xHistory[j])->linearCombo(1.0,*(ds.xHistory[j]),1.0,*(ds.xHistory[j+1]));
    (ds.qHistory[j])->linearCombo(1.0,*(ds.qHistory[j]),1.0,*(ds.qHistory[j+1]));
  }

  // Update State History
  if (sec.usedOrder_ < sec.maxOrder_)  
  {
    *(ds.sHistory[sec.usedOrder_+1]) = *ds.sNewtonCorrectionPtr;
  }
  // Update history arrays
  (ds.sHistory[sec.usedOrder_])->linearCombo(1.0,*(ds.sHistory[sec.usedOrder_]),1.0,*ds.sNewtonCorrectionPtr);
  for (int j=sec.usedOrder_-1;j>=0;j--) 
  {
    (ds.sHistory[j])->linearCombo(1.0,*(ds.sHistory[j]),1.0,*(ds.sHistory[j+1]));
  }

  // Update Store History
  if (sec.usedOrder_ < sec.maxOrder_)  
  {
    *(ds.stoHistory[sec.usedOrder_+1]) = *ds.stoNewtonCorrectionPtr;
    *(ds.stoLeadCurrQCompHistory[sec.usedOrder_+1]) = *ds.stoLeadCurrQCompNewtonCorrectionPtr;
  }
  // Update history arrays
  (ds.stoHistory[sec.usedOrder_])->linearCombo(1.0,*(ds.stoHistory[sec.usedOrder_]),1.0,*ds.stoNewtonCorrectionPtr);
  (ds.stoLeadCurrQCompHistory[sec.usedOrder_])->
    linearCombo(1.0,*(ds.stoLeadCurrQCompHistory[sec.usedOrder_]),1.0,*ds.stoLeadCurrQCompNewtonCorrectionPtr);
  for (int j=sec.usedOrder_-1;j>=0;j--) 
  {
    (ds.stoHistory[j])->linearCombo(1.0,*(ds.stoHistory[j]),1.0,*(ds.stoHistory[j+1]));
    (ds.stoLeadCurrQCompHistory[j])->
      linearCombo(1.0,*(ds.stoLeadCurrQCompHistory[j]),1.0,*(ds.stoLeadCurrQCompHistory[j+1]));
  }
  
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 1)
  {
    cout << "\n After updates \n" << endl;
    cout << "\n newtonCorrectionPtr: " << endl;
    ds.newtonCorrectionPtr->printPetraObject();
    cout << "\n qnewtonCorrectionPtr: " << endl;
    ds.qNewtonCorrectionPtr->printPetraObject();
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      cout << "\n xHistory["<< i << "]: \n" << endl;
      (ds.xHistory[i])->printPetraObject();
      cout << endl;
    }
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      cout << "\n qHistory["<< i << "]: \n" << endl;
      (ds.qHistory[i])->printPetraObject();
      cout << endl;
    }
    cout << "\n sNewtonCorrectionPtr: " << endl;
    ds.sNewtonCorrectionPtr->printPetraObject();
    cout << endl;
    cout << "\n nextStatePtr: " << endl;
    ds.nextStatePtr->printPetraObject();
    cout << endl;
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      cout << "\n sHistory["<< i << "]: \n" << endl;
      (ds.sHistory[i])->printPetraObject();
      cout << endl;
    }
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::restoreHistory
// Purpose       : Restore history array after a failed step
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 03/08/04
//-----------------------------------------------------------------------------
void N_TIA_BackwardDifferentiation15::restoreHistory()
{

  // undo preparation of history array for prediction
  for (int i=sec.nscsco_;i<=sec.currentOrder_;++i)
  {
    (ds.xHistory[i])->scale(1/sec.beta_[i]);
    (ds.qHistory[i])->scale(1/sec.beta_[i]);
    (ds.sHistory[i])->scale(1/sec.beta_[i]);
    (ds.stoHistory[i])->scale(1/sec.beta_[i]);
  }
  for (int i=1;i<=sec.currentOrder_;++i)
  {
    sec.psi_[i-1] = sec.psi_[i] - (sec.currentTimeStep);
  }
#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 1)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::restoreHistory");
    for (int i=1;i<=sec.currentOrder_;++i)
      cout << "\n sec.psi_[i] = " << sec.psi_[i] << endl;
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      cout << "\n xHistory["<< i << "]: \n" << endl;
      (ds.xHistory[i])->printPetraObject();
      cout << endl;
    }
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      cout << "\n qHistory["<< i << "]: \n" << endl;
      (ds.qHistory[i])->printPetraObject();
      cout << endl;
    }
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      cout << "\n sHistory["<< i << "]: \n" << endl;
      (ds.sHistory[i])->printPetraObject();
      cout << endl;
    }
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME
} 

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::updateCoeffs
// Purpose       : Update method coefficients
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 03/08/04
//-----------------------------------------------------------------------------
void N_TIA_BackwardDifferentiation15::updateCoeffs()
{
  // synchronize with Step Error Control
//  sec.psi_[0] = sec.currentTimeStep;
#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::updateCoeffs");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  currentTimeStep = ", sec.currentTimeStep);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  numberOfSteps_ = ", sec.numberOfSteps_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  currentOrder_ = ", sec.currentOrder_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  nscsco_ = ", sec.nscsco_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi_[0] = ", sec.psi_[0]);
  }
#endif
  // If the number of steps taken with constant order and constant stepsize is
  // more than the current order + 1 then we don't bother to update the
  // coefficients because we've reached a constant step-size formula.  When
  // this is is not true, then we update the coefficients for the variable
  // step-sizes. 
  if ((sec.currentTimeStep != sec.usedStep_) || (sec.currentOrder_ != sec.usedOrder_))
    sec.nscsco_ = 0;
  sec.nscsco_ = min(sec.nscsco_+1,sec.usedOrder_+2);
  if (sec.currentOrder_+1 >= sec.nscsco_)
  {
    sec.beta_[0] = 1.0;
    sec.alpha_[0] = 1.0;
    double temp1 = sec.currentTimeStep;
    sec.sigma_[0] = 1.0;
    sec.gamma_[0] = 0.0;
    for (int i=1;i<=sec.currentOrder_;++i)
    {
      double temp2 = sec.psi_[i-1];
      sec.psi_[i-1] = temp1;
      sec.beta_[i] = sec.beta_[i-1]*sec.psi_[i-1]/temp2;
      temp1 = temp2 + sec.currentTimeStep;
      sec.alpha_[i] = (sec.currentTimeStep)/temp1;
      sec.sigma_[i] = (i+1)*sec.sigma_[i-1]*sec.alpha_[i];
      sec.gamma_[i] = sec.gamma_[i-1]+sec.alpha_[i-1]/(sec.currentTimeStep);
    }
    sec.psi_[sec.currentOrder_] = temp1;
    sec.alphas_ = 0.0;
    sec.alpha0_ = 0.0;
    for (int i=0;i<sec.currentOrder_;++i)
    {
      sec.alphas_ = sec.alphas_ - 1.0/(i+1.0);
      sec.alpha0_ = sec.alpha0_ - sec.alpha_[i];
    }
    sec.cj_ = -sec.alphas_/(sec.currentTimeStep);
    sec.ck_ = abs(sec.alpha_[sec.currentOrder_]+sec.alphas_-sec.alpha0_);
    sec.ck_ = max(sec.ck_,sec.alpha_[sec.currentOrder_]);
  }
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  nscsco_ = ", sec.nscsco_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  beta_[0] = ", sec.beta_[0]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  beta_[1] = ", sec.beta_[1]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  beta_[2] = ", sec.beta_[2]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  beta_[3] = ", sec.beta_[3]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  beta_[4] = ", sec.beta_[4]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha_[0] = ", sec.alpha_[0]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha_[1] = ", sec.alpha_[1]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha_[2] = ", sec.alpha_[2]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha_[3] = ", sec.alpha_[3]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha_[4] = ", sec.alpha_[4]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alphas_ = ", sec.alphas_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha0_ = ", sec.alpha0_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  gamma_[0] = ", sec.gamma_[0]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  gamma_[1] = ", sec.gamma_[1]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  gamma_[2] = ", sec.gamma_[2]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  gamma_[3] = ", sec.gamma_[3]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  gamma_[4] = ", sec.gamma_[4]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi_[0] = ", sec.psi_[0]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi_[1] = ", sec.psi_[1]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi_[2] = ", sec.psi_[2]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi_[3] = ", sec.psi_[3]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi_[4] = ", sec.psi_[4]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  sigma_[0] = ", sec.sigma_[0]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  sigma_[1] = ", sec.sigma_[1]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  sigma_[2] = ", sec.sigma_[2]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  sigma_[3] = ", sec.sigma_[3]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  sigma_[4] = ", sec.sigma_[4]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  ck_ = ", sec.ck_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::initialize
// Purpose       : Initialize method with initial solution & step-size
// Special Notes : 
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 3/09/04
//-----------------------------------------------------------------------------
void N_TIA_BackwardDifferentiation15::initialize()
{

  // we assume the solution vector is available here 
  // Note that I'm using currSolutionPtr instead of
  // nextSolutionPtr because this is the first step.

  // Update next stop time from StepErrorControl:
  // ERK.  Commenting this out, as it is already called from N_ANP_AnalysisManager,
  // right before this initialize call.  It should not be called 2x, as 
  // it is history dependent (unfortunately), so calling it 2x in a row changes 
  // the stop time to a different number.
  // sec.updateStopTime();

  // Choose initial step-size
  double time_to_stop = sec.stopTime - sec.currentTime;
  double currentTimeStep;
  if (tiaParams.constantStepSize)
  {
    currentTimeStep = 0.1 * time_to_stop;
    currentTimeStep = Xycemin(sec.startingTimeStep, currentTimeStep);
    sec.currentTimeStep = currentTimeStep;
  }
  else
  {
    // compute an initial step-size based on rate of change in the 
    // solution initially
#ifdef Xyce_INCOMPLETE_2LEVEL_NORMS
    double dnorm_q = 0.0;
    (ds.qHistory[1])->wRMSNorm(*ds.qErrWtVecPtr, &dnorm_q);
#else
    double dnorm_q = ds.delta_x_errorNorm_q1();
#endif
    if (dnorm_q > 0.0)  // time-dependent DAE
    {
      currentTimeStep = Xycemin(sec.h0_max_factor_*abs(time_to_stop),sqrt(2.0)/(sec.h0_safety_*dnorm_q));
    } 
    else  // non-time-dependent DAE
    {
      currentTimeStep = sec.h0_max_factor_*abs(time_to_stop);
    }
    // choose min of user specified value and our value:
    if (sec.startingTimeStep > 0.0)
      currentTimeStep = Xycemin(sec.startingTimeStep, currentTimeStep);
    // check for maximum step-size:
    double rh = abs(currentTimeStep)*sec.h_max_inv_; 
    if (rh>1.0) currentTimeStep = currentTimeStep/rh;


    // Apply this new stepsize only if it is smaller than the one preceding 
    // the breakpoint, but only do this if this is a non-DCOP breakpoint.
    if (sec.currentTime != sec.initialTime) // if not DCOP:
    {
      sec.currentTimeStep = Xycemin(sec.lastTimeStep, currentTimeStep);
    }
    else // else if DCOP:
    {
      sec.currentTimeStep = currentTimeStep;
    }
  }

  sec.currentTimeStepRatio = 1.0;
  sec.currentTimeStepSum   = 2.0*sec.currentTimeStep;

  sec.lastTimeStep      = sec.currentTimeStep;
  sec.lastTimeStepRatio = sec.currentTimeStepRatio;
  sec.lastTimeStepSum   = sec.currentTimeStepSum;

  sec.numberSuccessiveFailures = 0;
  sec.stepAttemptStatus        = true;

//  sec.tolAimFac_ = 0.5;

  sec.nextTime = sec.currentTime + sec.currentTimeStep;

  // x history
  *(ds.xHistory[0]) = *(ds.currSolutionPtr);
  (ds.xHistory[1])->putScalar(0.0); // no need to multiply by dt here

  // q history
  *(ds.qHistory[0]) = *(ds.daeQVectorPtr);
  *(ds.qHistory[1]) = *(ds.daeFVectorPtr);
  (ds.qHistory[1])->scale(-sec.currentTimeStep);

  // state history
  *(ds.sHistory[0]) = *(ds.currStatePtr);
  (ds.sHistory[1])->putScalar(0.0); 

  // store history
  *(ds.stoHistory[0]) = *(ds.currStorePtr);
  (ds.stoHistory[1])->putScalar(0.0); 
  
  // lead current Q compontent history
  *(ds.stoLeadCurrQCompHistory[0]) = *(ds.currStoreLeadCurrQCompPtr);
  (ds.stoLeadCurrQCompHistory[1])->putScalar(0.0);
    
  // Coefficient initialization 
  sec.numberOfSteps_ = 0;    // number of total time integration steps taken
  sec.currentOrder_ = 1;
  sec.usedOrder_ = 1;
  sec.psi_[0] = sec.currentTimeStep;
  sec.cj_ = 1/sec.psi_[0];
  sec.nscsco_ = 0;
#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::initialize");
  if (tiaParams.debugLevel > 1)
  {
    cout << "\n xHistory: \n" << endl;
    (ds.xHistory[0])->printPetraObject();
    cout << endl;
    (ds.xHistory[1])->printPetraObject();
    cout << endl;
    cout << "\n qHistory: \n" << endl;
    (ds.qHistory[0])->printPetraObject();
    cout << endl;
    (ds.qHistory[1])->printPetraObject();
    cout << endl;
    cout << "\n sHistory: \n" << endl;
    (ds.sHistory[0])->printPetraObject();
    cout << endl;
    (ds.sHistory[1])->printPetraObject();
    cout << endl;
    cout << "\n" << "currentTimeStep = " << currentTimeStep << "\n" << endl;
    cout << "\n" << "time_to_stop = " << time_to_stop << "\n" << endl;
  }
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
#endif // Xyce_DEBUG_TIME
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::setTwoLevelTimeInfo
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/01/07
//-----------------------------------------------------------------------------
void N_TIA_BackwardDifferentiation15::setTwoLevelTimeInfo
  (const N_TIA_TimeIntInfo & tiInfo)
{
  // set initial step-size
  double time_to_stop = sec.stopTime - sec.currentTime;


  // x history
  *(ds.xHistory[0]) = *(ds.currSolutionPtr);
  (ds.xHistory[1])->putScalar(0.0); // no need to multiply by dt here

  // q history
  *(ds.qHistory[0]) = *(ds.daeQVectorPtr);
  *(ds.qHistory[1]) = *(ds.daeFVectorPtr);
  (ds.qHistory[1])->scale(-sec.currentTimeStep);

  // state history
  *(ds.sHistory[0]) = *(ds.nextStatePtr);
  (ds.sHistory[1])->putScalar(0.0);
  
  // store history
  *(ds.stoHistory[0]) = *(ds.nextStatePtr);
  (ds.stoHistory[1])->putScalar(0.0);
  
  // lead current Q compontent history
  *(ds.stoLeadCurrQCompHistory[0]) = *(ds.nextStoreLeadCurrQCompPtr);
  (ds.stoLeadCurrQCompHistory[1])->putScalar(0.0);

  // Coefficient initialization 
  sec.numberOfSteps_ = 0;    // number of total time integration steps taken
  sec.usedOrder_ = 1;
  sec.psi_[0] = sec.currentTimeStep;
  sec.cj_ = 1/sec.psi_[0];
  sec.nscsco_ = 0;

}


//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::checkReduceOrder()
// Purpose       : check whether to reduce order independent of local error test
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 06/08/04
//-----------------------------------------------------------------------------
void N_TIA_BackwardDifferentiation15::checkReduceOrder()
{

// This routine puts its output in sec.newOrder_

// This routine changes the following variables:
//    sec.Ek_, sec.Tk_, sec.Est_, sec.newOrder_, ds.delta_x, ds.delta_q,
//    sec.Ekm1_, sec.Tkm1_, sec.Ekm2_, sec.Tkm2_ 

// This routine reads but does not change the following variables:
//    sec.currentOrder_, sec.sigma_, ds.newtonCorrectionPtr, ds.qNewtonCorrectionPtr,
//    ds.errWtVecPtr, ds.qErrWtVecPtr, ds.xHistory, ds.qHistory

#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::checkReduceOrder");
  }
#endif // Xyce_DEBUG_TIME

  // 03/11/04 tscoffe:  I only want to run this block after a step has been
  // attempted, but I want to do this regardless of the status from the local
  // error test.
  // 03/10/04 tscoffe:  Decide whether to reduce the order before considering
  // the local error test result.
#ifndef Xyce_USE_Q_NORM
  //double dnorm_x = 0.0;
  //ds.newtonCorrectionPtr->wRMSNorm(&dnorm_x); // delta = newtonCorrection
  //double dnorm_x *= sec.ck_;
  //double dnorm = dnorm_x;
  double dnorm_x = sec.estOverTol_;
  double dnorm = dnorm_x;
#else
  double dnorm_x = 0.0, dnorm_q = 0.0;
  ds.newtonCorrectionPtr->wRMSNorm(*ds.errWtVecPtr, &dnorm_x); // delta = newtonCorrection
  ds.qNewtonCorrectionPtr->wRMSNorm(*ds.qErrWtVecPtr, &dnorm_q); // dnorm = norm of delta
  double dnorm = sqrt(0.5*dnorm_x*dnorm_x+0.5*dnorm_q*dnorm_q);
#endif

  sec.Ek_ = sec.sigma_[sec.currentOrder_]*dnorm;
  sec.Tk_ = (sec.currentOrder_+1)*sec.Ek_;
  sec.Est_ = sec.Ek_;
  sec.newOrder_ = sec.currentOrder_;
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 0)
  {  
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  Est_= ", sec.Est_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  Ek_= ", sec.Ek_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  dnorm = ", dnorm );
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  sigma[order] = ",sec.sigma_[sec.currentOrder_] );
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  Tk_= ", sec.Tk_);
  }
#endif // Xyce_DEBUG_TIME

     
  if (sec.currentOrder_>1)
  {

#ifndef Xyce_USE_Q_NORM
    ds.delta_x->linearCombo(1.0,*(ds.xHistory[sec.currentOrder_]),1.0,*ds.newtonCorrectionPtr);
#ifdef Xyce_INCOMPLETE_2LEVEL_NORMS
    ds.delta_x->wRMSNorm(*ds.errWtVecPtr, &dnorm_x);
    dnorm_x *= sec.ck_;
#else
    dnorm_x = sec.ck_ * ds.delta_x_errorNorm_m1();
#endif
    dnorm = dnorm_x;
#else
    ds.delta_x->linearCombo(1.0,*(ds.xHistory[sec.currentOrder_]),1.0,*ds.newtonCorrectionPtr);
    ds.delta_q->linearCombo(1.0,*(ds.qHistory[sec.currentOrder_]),1.0,*ds.qNewtonCorrectionPtr);
    ds.delta_x->wRMSNorm(*ds.errWtVecPtr, &dnorm_x);
    ds.delta_q->wRMSNorm(*ds.qErrWtVecPtr, &dnorm_q);
    dnorm = sqrt(0.5*dnorm_x*dnorm_x+0.5*dnorm_q*dnorm_q);
#endif

    sec.Ekm1_ = sec.sigma_[sec.currentOrder_-1]*dnorm;
    sec.Tkm1_ = sec.currentOrder_*sec.Ekm1_;
#ifdef Xyce_DEBUG_TIME
    if (tiaParams.debugLevel > 0)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  Ekm1_= ", sec.Ekm1_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  Tkm1_= ", sec.Tkm1_);
    }
#endif // Xyce_DEBUG_TIME

    if (sec.currentOrder_>2)
    {
#ifndef Xyce_USE_Q_NORM
      ds.delta_x->linearCombo(1.0,*(ds.xHistory[sec.currentOrder_-1]),1.0,*ds.delta_x);
#ifdef Xyce_INCOMPLETE_2LEVEL_NORMS
      ds.delta_x->wRMSNorm(*ds.errWtVecPtr, &dnorm_x);
      dnorm_x *= sec.ck_;
#else
      dnorm_x = sec.ck_ * ds.delta_x_errorNorm_m2();
#endif
      dnorm = dnorm_x;
#else
      ds.delta_x->linearCombo(1.0,*(ds.xHistory[sec.currentOrder_-1]),1.0,*ds.delta_x);
      ds.delta_q->linearCombo(1.0,*(ds.qHistory[sec.currentOrder_-1]),1.0,*ds.delta_q);
      ds.delta_x->wRMSNorm(*ds.errWtVecPtr, &dnorm_x);
      ds.delta_q->wRMSNorm(*ds.qErrWtVecPtr, &dnorm_q);
      dnorm = sqrt(0.5*dnorm_x*dnorm_x+0.5*dnorm_q*dnorm_q);
#endif

      sec.Ekm2_ = sec.sigma_[sec.currentOrder_-2]*dnorm;
      sec.Tkm2_ = (sec.currentOrder_-1)*sec.Ekm2_;
      if (Xycemax(sec.Tkm1_,sec.Tkm2_)<=sec.Tk_)
      {
        sec.newOrder_--;
        sec.Est_ = sec.Ekm1_;
      }
    }
    else if (sec.Tkm1_ <= sec.Tkm1_Tk_safety_ * sec.Tk_)
    {
#ifdef Xyce_DEBUG_TIME
    if (tiaParams.debugLevel > 0)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  Tkm1_Tk_safety= ", sec.Tkm1_Tk_safety_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  Tkm1_= ", sec.Tkm1_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  Tk_= ", sec.Tk_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  Tk_*safety= ", sec.Tk_*sec.Tkm1_Tk_safety_);
    }
#endif // Xyce_DEBUG_TIME

      sec.newOrder_--;
      sec.Est_ = sec.Ekm1_;
    }
  }
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  newOrder = ", sec.newOrder_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::rejectStep()
// Purpose       : code to restore history, choose new order/step-size
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 06/07/04
//-----------------------------------------------------------------------------
void N_TIA_BackwardDifferentiation15::rejectStep()
{

// This routine puts its output in newTimeStep_ and sec.newOrder_

// This routine changes the following variables:
//    lastAttemptedTimeStep, sec.initialPhase_, sec.nef_, sec.psi_, newTimeStep_,
//    sec.newOrder_, sec.currentOrder_, currentTimeStep_, ds.xHistory,
//    ds.qHistory, nextTimePt, nextTime, currentTimeStepRatio,
//    currentTimeStepSum, nextTimePt

// This routine reades but does not change the following variables:
//    stepAttemptStatus, sec.r_factor_, sec.r_safety_, sec.Est_, sec.r_fudge_, sec.r_min_, sec.r_max_,
//    minTimeStep, maxTimeStep, currentTime, stopTime, lastTimeStep


  // First we decide if we'll reduce the order independent of the local error test:
  checkReduceOrder();

#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";

  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::rejectStep");
  }
#endif

  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = (!(tiaParams.constantStepSize) );

  sec.lastAttemptedTimeStep = sec.currentTimeStep;


  double newTimeStep_ = sec.currentTimeStep;
  double rr = 1.0; // step size ratio = new step / old step
  if ((sec.stepAttemptStatus == false) && (adjustStep))
  {
    sec.initialPhase_ = false;
    sec.nef_++;
    restoreHistory();
    // restore sec.psi_
//    for (int i=1;i<=sec.currentOrder_;++i)
//      sec.psi_[i-1] = sec.psi_[i] - sec.currentTimeStep;



    // erkeite:  11/13/2007.  Commenting this max failure error out.  The 
    // Code already checks against a minimum time step, based on machine-precision.
    // The machine-precision limit is a better test.
#if 0
    if (sec.nef_ >= sec.max_LET_fail_)  
    {
      string msg = "N_TIA_BackwardDifferentiation15::rejectStep: ";
      msg += "  Maximum number of local error test failures.  ";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    }
#endif

    if ((sec.newtonConvergenceStatus <= 0))
    {
      /// 11/11/05 erkeite:  If the Newton solver fails, don't 
      // rely on the error estimate - it may be full of Nan's.
      rr = sec.r_min_;
      newTimeStep_ = rr * sec.currentTimeStep;

      if (sec.nef_ > 2) sec.newOrder_ = 1;//consistent with block below.
#ifdef Xyce_DEBUG_TIME
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "rejectStep 1");
#endif
    }
    else
    {
      // 03/11/04 tscoffe:  Here is the block for choosing order & 
      // step-size when the local error test FAILS (but Newton 
      // succeeded). 
      if (sec.nef_ == 1) // first local error test failure
      {
#ifndef Xyce_USE_Q_NORM
        rr = sec.r_factor_*pow(sec.r_safety_*(sec.Est_+sec.r_fudge_),-1.0/(sec.newOrder_+1.0));
        rr = Xycemax(sec.r_min_,Xycemin(sec.r_max_,rr));
#else
        rr = sec.r_factor_*pow(sec.r_safety_*sec.Est_+sec.r_fudge_,-1.0/(sec.newOrder_+1.0));
        rr = Xycemax(sec.r_min_,Xycemin(sec.r_max_,rr));
#endif
        newTimeStep_ = rr * sec.currentTimeStep;
#ifdef Xyce_DEBUG_TIME
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "rejectStep 2");
#endif
      }
      else if (sec.nef_ == 2) // second Dae failure
      {
        rr = sec.r_min_;
        newTimeStep_ = rr * sec.currentTimeStep;
#ifdef Xyce_DEBUG_TIME
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "rejectStep 3");
#endif
      }
      else if (sec.nef_ > 2) // third and later failures
      {
        sec.newOrder_ = 1;
        rr = sec.r_min_;
        newTimeStep_ = rr * sec.currentTimeStep;
#ifdef Xyce_DEBUG_TIME
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "rejectStep 4");
#endif
      }
    }
    if (sec.newOrder_ >= sec.minOrder_) 
    {
      sec.currentOrder_ = sec.newOrder_;
#ifdef Xyce_DEBUG_TIME
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "rejectStep 5");
#endif
    }
    if (sec.numberOfSteps_ == 0) // still first step
    {
      sec.psi_[0] = newTimeStep_;
      (ds.xHistory[1])->scale(rr);
      (ds.qHistory[1])->scale(rr);
#ifdef Xyce_DEBUG_TIME
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "rejectStep 6");
#endif
    }
#ifdef Xyce_DEBUG_TIME
    if (tiaParams.debugLevel > 0)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  currentTimeStep = ", sec.currentTimeStep);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  numberOfSteps_ = ", sec.numberOfSteps_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  currentOrder_ = ", sec.currentOrder_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  nscsco_ = ", sec.nscsco_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  alpha_[0] = ", sec.alpha_[0]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  alpha_[1] = ", sec.alpha_[1]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  alpha_[2] = ", sec.alpha_[2]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  alpha_[3] = ", sec.alpha_[3]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  alpha_[4] = ", sec.alpha_[4]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  psi_[0] = ", sec.psi_[0]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  psi_[1] = ", sec.psi_[1]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  psi_[2] = ", sec.psi_[2]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  psi_[3] = ", sec.psi_[3]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  psi_[4] = ", sec.psi_[4]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  sigma_[0] = ", sec.sigma_[0]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  sigma_[1] = ", sec.sigma_[1]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  sigma_[2] = ", sec.sigma_[2]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  sigma_[3] = ", sec.sigma_[3]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  sigma_[4] = ", sec.sigma_[4]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  rr = ", rr);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_factor_ = ", sec.r_factor_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_safety_ = ", sec.r_safety_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  Est_ = ", sec.Est_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_fudge_ = ", sec.r_fudge_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  newOrder_ = ", sec.newOrder_);

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  currentTimeStep = ", sec.currentTimeStep);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  newTimeStep_ = ", newTimeStep_);
    }
#endif // Xyce_DEBUG_TIME
  }
  else if ((sec.stepAttemptStatus == false) & (!adjustStep))
  {
    string tmp = "  BackwardDifferentiation15:rejectStep: Warning: Local error test failed with constant step-size.\n";
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_INFO_0, tmp);
  }

  // If the step needs to be adjusted:
  if (adjustStep)
  {
    newTimeStep_ = Xycemax(newTimeStep_, sec.minTimeStep);
    newTimeStep_ = Xycemin(newTimeStep_, sec.maxTimeStep);

    double nextTimePt = sec.currentTime + newTimeStep_;

    if (nextTimePt > sec.stopTime)
    {
      nextTimePt  = sec.stopTime;
      newTimeStep_ = sec.stopTime - sec.currentTime;
    }

    sec.nextTime = nextTimePt;

    sec.currentTimeStepRatio = newTimeStep_/sec.lastTimeStep;
    sec.currentTimeStepSum   = newTimeStep_ + sec.lastTimeStep;

#ifdef Xyce_DEBUG_TIME
    if (tiaParams.debugLevel >0)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  newTimeStep_ = ", newTimeStep_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  nextTime = ", sec.nextTime);
    }
#endif // Xyce_DEBUG_TIME

    sec.currentTimeStep = newTimeStep_;
  }
  else // if time step is constant for this step:
  {
    double nextTimePt = sec.currentTime + sec.currentTimeStep;

    if (nextTimePt > sec.stopTime)
    {
      nextTimePt      = sec.stopTime;
      sec.currentTimeStep = sec.stopTime - sec.currentTime;
    }

    sec.currentTimeStepRatio = sec.currentTimeStep / sec.lastTimeStep;
    sec.currentTimeStepSum   = sec.currentTimeStep + sec.lastTimeStep;

    sec.nextTime = nextTimePt;
  }
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel >0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::rejectStepForHabanero
// Purpose       : step rejection, but from an outside program (Habanero API)
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 08/11/09
//-----------------------------------------------------------------------------
void N_TIA_BackwardDifferentiation15::rejectStepForHabanero()
{
  restoreHistory();
  sec.setTimeStep(sec.currentTimeStep);
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::completeStep()
// Purpose       : code to update history, choose new order/step-size
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 06/07/04
//-----------------------------------------------------------------------------
void N_TIA_BackwardDifferentiation15::completeStep()
{

  sec.numberOfSteps_ ++;
  sec.nef_ = 0;
  sec.lastTime    = sec.currentTime;
  sec.currentTime = sec.nextTime;
  // First we decide if we'll reduce the order independent of the local error test:
  checkReduceOrder();


#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";

  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::completeStep");
  }
#endif

  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = (!(tiaParams.constantStepSize) );

  sec.lastAttemptedTimeStep = sec.currentTimeStep;

  double newTimeStep_ = sec.currentTimeStep;
  double rr = 1.0; // step size ratio = new step / old step
  // 03/11/04 tscoffe:  Here is the block for choosing order & step-size when
  // the local error test PASSES (and Newton succeeded). 
  sec.lastTimeStep = sec.currentTimeStep;
  sec.lastTimeStepRatio = sec.currentTimeStepRatio; // copied from calcTStep1
  sec.lastTimeStepSum   = sec.currentTimeStepSum; // copied from calcTStep1
  int orderDiff = sec.currentOrder_ - sec.usedOrder_;
  sec.usedOrder_ = sec.currentOrder_;
  sec.usedStep_ = sec.currentTimeStep;
  if ((sec.newOrder_ == sec.currentOrder_-1) || (sec.currentOrder_ == sec.maxOrder_))
  {
    // If we reduced our order or reached max order then move to the next phase
    // of integration where we don't automatically double the step-size and
    // increase the order.
    sec.initialPhase_ = false;
  }
  if (sec.initialPhase_)
  {
    sec.currentOrder_++;
    newTimeStep_ = sec.h_phase0_incr_ * sec.currentTimeStep;
  }
  else // not in the initial phase of integration
  {
    int action = TIAAction_UNSET;
    if (sec.newOrder_ == sec.currentOrder_-1)
      action = TIAAction_LOWER;
    else if (sec.newOrder_ == sec.maxOrder_)
      action = TIAAction_MAINTAIN;
    else if ((sec.currentOrder_+1>=sec.nscsco_) || (orderDiff == 1))
    {
      // If we just raised the order last time then we won't raise it again
      // until we've taken sec.currentOrder_+1 steps at order sec.currentOrder_.
      action = TIAAction_MAINTAIN;
    }
    else // consider changing the order 
    {

#ifndef Xyce_USE_Q_NORM
      ds.delta_x->linearCombo(1.0,*ds.newtonCorrectionPtr,-1.0,*(ds.xHistory[sec.currentOrder_+1]));
#ifdef Xyce_INCOMPLETE_2LEVEL_NORMS
      double dnorm_x = 0.0;
      ds.delta_x->wRMSNorm(*ds.errWtVecPtr, &dnorm_x);
      dnorm_x *= sec.ck_;
#else
      double dnorm_x = sec.ck_ * ds.delta_x_errorNorm_p1();
#endif
      double dnorm = dnorm_x;
#else
      ds.delta_x->linearCombo(1.0,*ds.newtonCorrectionPtr,-1.0,*(ds.xHistory[sec.currentOrder_+1]));
      ds.delta_q->linearCombo(1.0,*ds.qNewtonCorrectionPtr,-1.0,*(ds.qHistory[sec.currentOrder_+1]));
      double dnorm_x = 0.0, dnorm_q = 0.0;
      ds.delta_x->wRMSNorm(*ds.errWtVecPtr, &dnorm_x);
      ds.delta_q->wRMSNorm(*ds.qErrWtVecPtr, &dnorm_q);
      double dnorm = sqrt(0.5*dnorm_x*dnorm_x+0.5*dnorm_q*dnorm_q);
#endif

      sec.Tkp1_ = dnorm;
      sec.Ekp1_ = sec.Tkp1_/(sec.currentOrder_+2);
      if (sec.currentOrder_ == 1)
      {
        if (sec.Tkp1_ >= sec.Tkp1_Tk_safety_ * sec.Tk_)
          action = TIAAction_MAINTAIN;
        else
          action = TIAAction_RAISE;
      }
      else
      {
        if (sec.Tkm1_ <= Xycemin(sec.Tk_,sec.Tkp1_))
          action = TIAAction_LOWER;
        else if (sec.Tkp1_ >= sec.Tk_)
          action = TIAAction_MAINTAIN;
        else
          action = TIAAction_RAISE;
      }
    }
    if (sec.currentOrder_ < sec.minOrder_) 
    {
      action = TIAAction_RAISE;
    }
    else if ((sec.currentOrder_ == sec.minOrder_) && (action == TIAAction_LOWER))
    {
      action = TIAAction_MAINTAIN;
    }
    if (action == TIAAction_RAISE)
    {
      sec.currentOrder_++;
      sec.Est_ = sec.Ekp1_;
    }
    else if (action == TIAAction_LOWER)
    {
      sec.currentOrder_--;
      sec.Est_ = sec.Ekm1_;
    }
    newTimeStep_ = sec.currentTimeStep;

    // ERK:  if errorAnalysisOption==1, that means that we're not using LTE to determine the
    // new step size.  We're only considering whether or not the Newton solve succeeded.
    // If we are in this function, then it succeeded, as otherwise we'd be in the "rejectStep"
    // function.
    if (tiaParams.errorAnalysisOption == 1)
    {
      rr = 0.4/(sec.r_min_);
      newTimeStep_ = rr*sec.currentTimeStep;
    }
    else
    {
#ifndef Xyce_USE_Q_NORM
      rr = pow(sec.r_safety_*(sec.Est_+sec.r_fudge_),-1.0/(sec.currentOrder_+1.0));
#else
      rr = pow(sec.r_safety_*sec.Est_+sec.r_fudge_,-1.0/(sec.currentOrder_+1.0));
#endif

#ifdef Xyce_DEBUG_TIME
      if (tiaParams.debugLevel > 0)
      {
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                              "  currentOrder_ = ", sec.currentOrder_);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                              "  r_safety = ", sec.r_safety_);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                              "  r_fudge_ = ", sec.r_fudge_);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                              "  r_hincr_ = ", sec.r_hincr_);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                              "  r_hincr_test_ = ", sec.r_hincr_test_);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                              "  Est = ", sec.Est_);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                              "  Ekp1_ = ", sec.Ekp1_);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                              "  Ekm1_ = ", sec.Ekm1_);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                              "  Tkp1_ = ", sec.Tkp1_);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                              "  raw rr = ", rr);
      }
#endif
      if (rr >= sec.r_hincr_test_)
      {
        rr = sec.r_hincr_;
        newTimeStep_ = rr*sec.currentTimeStep;
      }
      else if (rr <= 1)
      {
        rr = Xycemax(sec.r_min_,Xycemin(sec.r_max_,rr));
        newTimeStep_ = rr*sec.currentTimeStep;
      }
    }
  }
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  initialPhase_ = ", sec.initialPhase_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  rr = ", rr);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  currentTimeStep = ", sec.currentTimeStep);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  currentTime = ", sec.currentTime);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  nextTime = ", sec.nextTime);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  newTimeStep_ = ", newTimeStep_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  minTimeStep = ", sec.minTimeStep);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  maxTimeStep = ", sec.maxTimeStep);
  }
#endif
  // 03/22/04 tscoffe:  Note that updating the history has nothing to do with
  // the step-size and everything to do with the newton correction vectors.
  updateHistory();


  // 12/01/05 tscoffe:  This is necessary to avoid currentTimeStep == 0 right
  // before a breakpoint.  So I'm checking to see if currentTime is identically
  // equal to stopTime, in which case we are right before a breakpoint and we
  // should not adjust currentStepSize because that would result in
  // currentStepSize == 0.
  if (sec.currentTime < sec.stopTime)
  {
    // If the step needs to be adjusted:
    if (adjustStep)
    {
      newTimeStep_ = Xycemax(newTimeStep_, sec.minTimeStep);
      newTimeStep_ = Xycemin(newTimeStep_, sec.maxTimeStep);

      double nextTimePt = sec.currentTime + newTimeStep_;

      if (nextTimePt > sec.stopTime)
      {
        nextTimePt  = sec.stopTime;
        newTimeStep_ = sec.stopTime - sec.currentTime;
      }

      sec.nextTime = nextTimePt;

      sec.currentTimeStepRatio = newTimeStep_/sec.lastTimeStep;
      sec.currentTimeStepSum   = newTimeStep_ + sec.lastTimeStep;

#ifdef Xyce_DEBUG_TIME
      if (tiaParams.debugLevel >0)
      {
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                "  nextTime = ", sec.nextTime);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                "  newTimeStep_ = ", newTimeStep_);
      }
#endif

      sec.currentTimeStep = newTimeStep_;
    }
    else // if time step is constant for this step:
    {
      double nextTimePt = sec.currentTime + sec.currentTimeStep;

      if (nextTimePt > sec.stopTime)
      {
        nextTimePt      = sec.stopTime;
        sec.currentTimeStep = sec.stopTime - sec.currentTime;
      }

      sec.currentTimeStepRatio = sec.currentTimeStep / sec.lastTimeStep;
      sec.currentTimeStepSum   = sec.currentTimeStep + sec.lastTimeStep;

      sec.nextTime = nextTimePt;
    }
  }
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel >0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif

// 11/02/04 tscoffe:  This should be done at the beginning of an integration
//                    stop rather than at the end, so its been moved into
//                    ControlAlgorithm::transientLoop_
  // Update next stop time in StepErrorControl:
  // sec.updateStopTime();
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::updateStateDeriv
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 11/17/05
//-----------------------------------------------------------------------------
void N_TIA_BackwardDifferentiation15::updateStateDeriv ()
{
  // dS/dt = spn0 - (sec.alpha_/hn)(S(x)-sn0)
  ds.nextStateDerivPtr->linearCombo(1.0,*ds.nextStatePtr,-1.0,*ds.sn0Ptr);
  ds.nextStateDerivPtr->linearCombo(1.0,*ds.spn0Ptr,-sec.alphas_/sec.currentTimeStep,*ds.nextStateDerivPtr);
}


//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::updateLeadCurrent
// Purpose       : calculates lead currents in store vector with 
//                 the storeLeadCurrQCompVec.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling, SNL
// Creation Date : 03/22/2013
//-----------------------------------------------------------------------------
void N_TIA_BackwardDifferentiation15::updateLeadCurrent ()
{
  // dStoreQ/dt = spn0 - (sec.alpha_/hn)(S(x)-sn0)
  ds.nextStoreLeadCurrQCompDerivPtr->
    linearCombo(1.0,*ds.nextStoreLeadCurrQCompPtr,-1.0,*ds.stoQCn0Ptr);
  ds.nextStoreLeadCurrQCompDerivPtr->
    linearCombo(1.0,*ds.stoQCpn0Ptr,-sec.alphas_/sec.currentTimeStep,*ds.nextStoreLeadCurrQCompDerivPtr);
  ds.nextStorePtr->linearCombo(1.0,*ds.nextStorePtr,1.0,*ds.nextStoreLeadCurrQCompDerivPtr);

}



//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::getInitialQnorm
// Purpose       : Needed by 2-level solves.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/18/07
//-----------------------------------------------------------------------------
void N_TIA_BackwardDifferentiation15::getInitialQnorm (N_TIA_TwoLevelError & tle)
{
  tle.q1HistorySum = ds.partialSum_q1();
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::setupTwoLevelError
// Purpose       : Needed by 2-level solves.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
void N_TIA_BackwardDifferentiation15::setupTwoLevelError(N_TIA_TwoLevelError & tle)
{
  tle.xErrorSum    = ds.partialErrorNormSum ();
  tle.qErrorSum    = ds.partialQErrorNormSum ();
  tle.xErrorSum_m1 = ds.partialSum_m1 (sec.currentOrder_);
  tle.xErrorSum_m2 = ds.partialSum_m2 (sec.currentOrder_);
  tle.xErrorSum_p1 = ds.partialSum_p1 (sec.currentOrder_, sec.maxOrder_);
  tle.innerSize    = ds.globalLength ();

#ifdef Xyce_DEBUG_TIME
  cout << tle;
#endif
}

