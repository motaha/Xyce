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
// Filename      : $RCSfile: N_TIA_Gear12.C,v $
//
// Purpose       : This file contains the functions which define the
//		             backward differentiation, order 1-2, class.
//
// Special Notes :
//
// Creator       : Ting Mei
//
// Creation Date : 2/16/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10.2.3 $
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
#include <N_TIA_Gear12.h>
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
// Function      : N_TIA_Gear::N_TIA_Gear
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/10/12
//-----------------------------------------------------------------------------

N_TIA_Gear12::N_TIA_Gear12
  (N_TIA_TIAParams & tiaP,
   N_TIA_StepErrorControl & secTmp,
   N_TIA_DataStore & dsTmp)
: N_TIA_TimeIntegrationMethod(tiaP,secTmp,dsTmp)
{
  leadingCoeff = 1;
  sec.maxOrder_=(Xycemin(2,tiaParams.maxOrder));
  sec.minOrder_=(Xycemax(1,tiaParams.minOrder));

  if (sec.minOrder_ > sec.maxOrder_)
  {
    sec.minOrder_ = sec.maxOrder_;
  }
//  sec.maxOrder_ = 2;
  timept_ = -1.0;
  return ;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_Gear12::obtainPredictor
// Purpose       : Calculate predictor 
// Special Notes : stored in ds.xn0Ptr,qn0Ptr,qpn0Ptr
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/10/12
//-----------------------------------------------------------------------------

void N_TIA_Gear12::obtainPredictor()
{
 
  // evaluate predictor
//  *ds.xn0Ptr = *(ds.xHistory[0]);
//  *ds.qn0Ptr = *(ds.qHistory[0]);
  *ds.sn0Ptr = *(ds.sHistory[0]);
  ds.xn0Ptr->putScalar(0.0);
  ds.qn0Ptr->putScalar(0.0);
  *ds.stoQCn0Ptr = *(ds.stoLeadCurrQCompHistory[0]);

  ds.spn0Ptr->putScalar(0.0);
  for (int i=0;i<=sec.currentOrder_;++i)
  {
    ds.xn0Ptr->linearCombo(sec.beta_[i],*(ds.xHistory[i]),1.0,*ds.xn0Ptr);
    ds.qn0Ptr->linearCombo(sec.beta_[i],*(ds.qHistory[i]),1.0,*ds.qn0Ptr);
//    ds.sn0Ptr->linearCombo(sec.beta_[i],*(ds.sHistory[i]),1.0,*ds.sn0Ptr);
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
                           "  N_TIA_Gear12::obtainPredictor");
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
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME
  // copy the prediction into the next solution:
  *(ds.nextSolutionPtr) = *(ds.xn0Ptr);

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_Gear12::obtainResidual
// Purpose       : Calculate Residual
// Special Notes : 
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------

void N_TIA_Gear12::obtainResidual()
{
  // output: ds.RHSVectorPtr
  // Note:  ds.nextSolutionPtr is used to get Q,F,B in N_ANP_AnalysisManager::loadRHS.
  ds.RHSVectorPtr->linearCombo(sec.alpha_[0],*ds.daeQVectorPtr, sec.alpha_[1],*(ds.qHistory[0]) );
  
#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 1)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_Gear12::obtainResidual");
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

  if (sec.currentOrder_  == 2)
  {
    ds.RHSVectorPtr->linearCombo(1.0, *ds.RHSVectorPtr, sec.alpha_[2],*(ds.qHistory[1]));
  }

  ds.RHSVectorPtr->linearCombo(1.0/sec.currentTimeStep,*ds.RHSVectorPtr,+1.0,*ds.daeFVectorPtr);

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
// Function      : N_TIA_Gear12::obtainJacobian
// Purpose       : Calculate Jacobian
// Special Notes : 
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------

void N_TIA_Gear12::obtainJacobian()
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
                           "  N_TIA_Gear12::obtainJacobian");
  }
#endif
  // output: ds.JMatrixPtr

  // This function returns the following matrix:
  // $-(sec.alphas_/hn)dQdx(x)+dFdx$

  // Note:  ds.nextSolutionPtr is used to get dQdx,dFdx in N_ANP_AnalysisManager::loadJacobian.

  N_LAS_Matrix & dQdx = *(ds.dQdxMatrixPtr);
  N_LAS_Matrix & dFdx = *(ds.dFdxMatrixPtr);
  N_LAS_Matrix & Jac = *(ds.JMatrixPtr);

  double qscalar(sec.alpha_[0]/sec.currentTimeStep);
  double fscalar(1.0);

#ifdef Xyce_NO_MATRIX_LINEAR_COMBO

  Jac.put( 0.0 );

#ifndef Xyce_FAST_MATRIX_ADD
  Jac.add( dQdx );
#else
  Jac.addLocal( dQdx );
#endif

#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 1)
  {
    cout << "\n dQdx: \n" <<endl;
    dQdx.printPetraObject();
  }
#endif

  Jac.scale( qscalar );

#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 1)
  {
    cout << "\n scaled dQdx by " << qscalar << " :" <<endl;
    Jac.printPetraObject();
  }
#endif

#ifndef Xyce_FAST_MATRIX_ADD
  Jac.add( dFdx );
#else
  Jac.addLocal( dFdx );
#endif
#else

  Jac.linearCombo( qscalar, dQdx, fscalar, dFdx );

#endif

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
// Function      : N_TIA_Gear12::interpolateSolution
// Purpose       : Interpolate solution approximation at prescribed time point.
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------

bool N_TIA_Gear12::interpolateSolution(double timepoint, 
        N_LAS_Vector * tmpSolVectorPtr, vector<N_LAS_Vector*> & historyVec)

{
  // this is a very course approximation to determine if we are too 
  // close to the actual time step to do an interpolation.
  // it could use more work. 
  double dtr = timepoint - sec.currentTime;  // the delta time requested.
  if( -dtr < 100 * N_UTL_MachineDependentParams::MachinePrecision() )
  {
    *tmpSolVectorPtr = *(historyVec[0]);
    return false;
  }

  tmpSolVectorPtr->linearCombo(1.0, *(historyVec[0]), -1.0, *(historyVec[1]));
  
  if( sec.usedOrder_ <= 2)
  {
    // do first order interpolation
    // X_interp = X + delta_t_requested * delta_X/delta_t[last step]
    dtr = dtr / sec.lastTimeStep;
    tmpSolVectorPtr->linearCombo(1.0, *(historyVec[0]), dtr, *tmpSolVectorPtr);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_Gear12::interpolateMPDESolution
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
// Creator       : Ting Mei, Eric Keiter, SNL 
// Creation Date : 11/28/06
//-----------------------------------------------------------------------------

bool N_TIA_Gear12::interpolateMPDESolution(std::vector<double>& timepoint, 
    	N_LAS_Vector * tmpSolVectorPtr)
{
#ifdef Xyce_PARALLEL_MPI
  string msg = "N_TIA_Gear12::interpolateMPDESolution: ";
  msg += "Not set up for Parallel yet";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return(false);
#endif

  N_LAS_BlockVector * blockTempSolVectorPtr = 
     dynamic_cast<N_LAS_BlockVector*>(tmpSolVectorPtr);
  if (blockTempSolVectorPtr == NULL)
  {
    string msg = "N_TIA_Gear12::interpolateMPDESolution: ";
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
    string msg = "N_TIA_Gear12::interpolateMPDESolution: ";
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
        string msg = "N_TIA_Gear12::interpolateMPDESolution: ";
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
// Function      : N_TIA_Gear12::printMPDEOutputSolution()
// Purpose       : Print transient output from MPDE simulation
// Special Notes : This routine uses interpolateMPDESolution.
// Scope         : public
// Creator       : Ting Mei, SNL, 1414
// Creation Date : 11/28/06
//-----------------------------------------------------------------------------

bool N_TIA_Gear12::printMPDEOutputSolution(
        RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
        const double time,
        N_LAS_Vector * solnVecPtr,
        const std::vector<double> & fastTimes )
{
/*
#ifdef Xyce_DEBUG_TIME
  const string dashedline =
	"---------------------------------------------------------------"
	"-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
      "  N_TIA_Gear12::printMPDEOutputSolution");
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
    string msg = "N_TIA_Gear12::printMPDEOutputSolution: ";
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
    string msg = "N_TIA_Gear12::printMPDEOutputSolution: ";
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
        outputMgrAdapterRCPtr->tranOutput(timept_, blockTmpSolVectorPtr->block(s), *ds.tmpStaVectorPtr, *ds.tmpStoVectorPtr  );
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
      interpolateSolution(stop_time,ds.tmpXn0BPtr);
      tmpSolnVecPtr = ds.tmpXn0BPtr;
    }
    N_LAS_BlockVector * blockTmpSolnVecPtr = 
      dynamic_cast<N_LAS_BlockVector*>(tmpSolnVecPtr);
    N_LAS_BlockVector * blockTmpVecPtr = 
      dynamic_cast<N_LAS_BlockVector*>(tmpVecPtr);
    if (blockTmpSolnVecPtr == NULL)
    {
      string msg = "N_TIA_Gear12::printMPDEOutputSolution: ";
      msg += "N_LAS_Vector tmpSolnVecPtr is not of type N_LAS_BlockVector";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
      return(false);
    }
    if (blockTmpVecPtr == NULL)
    {
      string msg = "N_TIA_Gear12::printMPDEOutputSolution: ";
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
      string msg = "N_TIA_Gear12::printMPDEOutputSolution: ";
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
    outputMgrAdapterRCPtr->tranOutput(timept_, blockTmpVecPtr->block(0), *ds.tmpStaVectorPtr, *ds.tmpStoVectorPtr  );
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
*/
}


//-----------------------------------------------------------------------------
// Function      : N_TIA_Gear12::printWaMPDEOutputSolution()
// Purpose       : Print transient output from WaMPDE simulation
// Special Notes : This routine uses interpolateSolution.
// Scope         : public
// Creator       : Ting Mei, SNL, 1414
// Creation Date : 12/15/06
//-----------------------------------------------------------------------------

bool N_TIA_Gear12::printWaMPDEOutputSolution(
        RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
        const double time,
        N_LAS_Vector * solnVecPtr,
        const std::vector<double> & fastTimes,
        const int phiGID )
{
/*
#ifdef Xyce_DEBUG_TIME
  const string dashedline =
	"---------------------------------------------------------------"
	"-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
      "  N_TIA_Gear12::printWaMPDEOutputSolution");
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
    string msg = "N_TIA_Gear12::printWaMPDEOutputSolution: ";
    msg += "N_LAS_Vector ds.tmpSolVectorPtr is not of type N_LAS_BlockVector";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    return(false);
  }
  if (blockTmpXn0APtr == NULL)
  {
    string msg = "N_TIA_Gear12::printWaMPDEOutputSolution: ";
    msg += "N_LAS_Vector ds.tmpXn0APtr is not of type N_LAS_BlockVector";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    return(false);
  }
  if (blockTmpXn0BPtr == NULL)
  {
    string msg = "N_TIA_Gear12::printWaMPDEOutputSolution: ";
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
      junk = interpolateSolution(timeA,ds.tmpXn0APtr);
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
    interpolateSolution(timeB,ds.tmpXn0BPtr);
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
      string msg = "N_TIA_Gear12::printWaMPDEOutputSolution: ";
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
        string msg = "N_TIA_Gear12::printWaMPDEOutputSolution: ";
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
        outputMgrAdapterRCPtr->tranOutput(t, blockTmpSolVectorPtr->block(b1), *ds.tmpStaVectorPtr, *ds.tmpStoVectorPtr  );
#ifdef Xyce_DEBUG_TIME
        cout << "Interpolated to (t,phi(t)) = (" << t << "," << s << ")" << endl;
#endif // Xyce_DEBUG_TIME
      }
    }
  }
  
#ifdef Xyce_DEBUG_TIME
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
#endif // Xyce_DEBUG_TIME
  return true;         */
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_Gear12::printOutputSolution()
// Purpose       : Print output that is dumbed down in order.
// Special Notes : This routine picks smaller time steps to approximate first
//               : order integration from the perspective of the output.
// Scope         : public
// Creator       : Ting Mei, SNL, 1414
// Creation Date : 11/16/07 
//-----------------------------------------------------------------------------

bool N_TIA_Gear12::printOutputSolution(
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
      "  N_TIA_Gear12::printOutputSolution");
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

  // Either way, do an output on the actual computed time step, but only
  // if we weren't given a list of specific times *or* we were told not to
  // interpoloate.
  if (outputInterpolationTimes.empty() || doNotInterpolate)
  {
    outputMgrAdapterRCPtr->tranOutput(time, *ds.currSolutionPtr, *ds.currStatePtr, *ds.currStorePtr, skipPrintLineOutput);
  }

#ifdef Xyce_DEBUG_TIME
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
#endif // Xyce_DEBUG_TIME
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_Gear12::saveOutputSolution
// Purpose       : This is similar to printOutputSolution, but is in support of
//                 the .SAVE capability, rather than .PRINT.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
bool N_TIA_Gear12::saveOutputSolution  ( 
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
  
  outputMgrAdapterRCPtr->outputDCOP( *(ds.currSolutionPtr) );

#ifdef Xyce_DEBUG_TIME
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
#endif // Xyce_DEBUG_TIME
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_Gear12::updateHistory
// Purpose       : Update history array after a successful step 
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
void N_TIA_Gear12::updateHistory()
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
                           "  N_TIA_Gear12::updateHistory");
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


   if (sec.currentOrder_ == 2)
   {
     *(ds.xHistory[2]) = *(ds.xHistory[1]);
     *(ds.qHistory[1]) = *(ds.qHistory[0]);
     *(ds.stoLeadCurrQCompHistory[1]) = *(ds.stoLeadCurrQCompHistory[0]);
     *(ds.sHistory[1]) = *(ds.sHistory[0]);
   }
     
   *(ds.xHistory[1]) = *(ds.xHistory[0]);
//   (ds.qHistory[1])->linearCombo(1.0, *ds.daeQVectorPtr, -1.0,*(ds.qHistory[0]));

   *(ds.xHistory[0]) = *ds.nextSolutionPtr;
   *(ds.qHistory[0]) =  *ds.daeQVectorPtr;
   *(ds.sHistory[0]) =  *ds.nextStatePtr;    
   *(ds.stoLeadCurrQCompHistory[0]) = *ds.nextStoreLeadCurrQCompPtr;

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
// Function      : N_TIA_Gear12::restoreHistory
// Purpose       : Restore history array after a failed step
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07 
//-----------------------------------------------------------------------------
void N_TIA_Gear12::restoreHistory()
{

  for (int i=1;i<=sec.currentOrder_;++i)
  {
    sec.psi_[i-1] = sec.psi_[i];
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
                           "  N_TIA_Gear12::restoreHistory");
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
// Function      : N_TIA_Gear12::updateCoeffs
// Purpose       : Update method coefficients
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
void N_TIA_Gear12::updateCoeffs()
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
                           "  N_TIA_Gear12::updateCoeffs");
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
  
  
  double temp1 = sec.currentTimeStep;
    
  if (sec.currentOrder_ == 2)
  {
	sec.psi_[2] = sec.psi_[1];
  }
  sec.psi_[1] = sec.psi_[0];
  sec.psi_[0] = temp1;
   
//    sec.beta_[0] = 1.0;
//    sec.alpha_[0] = 1.0; 
    sec.ck_ = 1.0;
    sec.alphas_ = -1.0;
    
  if (sec.currentOrder_ == 2)
  {

// the coeffs of predictor
    sec.beta_[2] = temp1/sec.psi_[2] * (temp1 + sec.psi_[1])/(sec.psi_[1] + sec.psi_[2]);
    sec.beta_[1] = -temp1/sec.psi_[1] - sec.beta_[2] * (sec.psi_[1] + sec.psi_[2])/sec.psi_[1];
    sec.beta_[0] = 1.0 - sec.beta_[2] - sec.beta_[1];

    sec.alpha_[2] = -temp1/sec.psi_[1] * temp1/(2 * temp1 + sec.psi_[1]); 
    sec.alpha_[1] = 1 - sec.alpha_[2];
    sec.alpha_[0] = -sec.alpha_[1] - sec.alpha_[2] * (1 + sec.psi_[1]/temp1);

    sec.alpha_[2] = sec.alpha_[2]/sec.alpha_[0];
    sec.alpha_[1] = sec.alpha_[1]/sec.alpha_[0];
    sec.alpha_[0] = -1/sec.alpha_[0];

    sec.ck_ = sec.currentTimeStep/(temp1 + sec.psi_[1] + sec.psi_[2]);
  }
  else
  {
    sec.beta_[0] = 1.0 + temp1/sec.psi_[1];
    sec.beta_[1] = -temp1/sec.psi_[1];
    sec.alpha_[0] = 1.0;
    sec.alpha_[1] = -1.0;

    sec.ck_ = sec.currentTimeStep/(temp1 + sec.psi_[1]);
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
// Function      : N_TIA_Gear12::initialize
// Purpose       : Initialize method with initial solution & step-size
// Special Notes : 
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
void N_TIA_Gear12::initialize()
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

  tiaParams.TimeStepLimitedbyBP =  false; 
  
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
      if (sec.currentTime == sec.initialTime)
        currentTimeStep = Xycemin(sec.h0_max_factor_*abs(time_to_stop),sqrt(2.0)/(sec.h0_safety_*dnorm_q));
      else
        currentTimeStep = 0.1* Xycemin(sec.savedTimeStep, abs(time_to_stop));
    } 
    else  // non-time-dependent DAE
    {
      if (sec.currentTime == sec.initialTime)
        currentTimeStep = sec.h0_max_factor_*abs(time_to_stop);
      else
        currentTimeStep = 0.1* Xycemin(sec.savedTimeStep, abs(time_to_stop));
    }
    // choose min of user specified value and our value:
    if (sec.startingTimeStep > 0.0 && (sec.currentTime == sec.initialTime))
      currentTimeStep = Xycemin(sec.startingTimeStep, currentTimeStep);
    // check for maximum step-size:
    double rh = abs(currentTimeStep)*sec.h_max_inv_; 
    if (rh>1.0) currentTimeStep = currentTimeStep/rh;


    // Apply this new stepsize only if it is smaller than the one preceding 
    // the breakpoint, but only do this if this is a non-DCOP breakpoint.
//    if (sec.currentTime != sec.initialTime) // if not DCOP:
//    {
//      sec.currentTimeStep = Xycemin(sec.currentTimeStep, currentTimeStep);
//    }
//    else // else if DCOP:
//    {
      sec.currentTimeStep = currentTimeStep;
//    }
  }

  sec.currentTimeStepRatio = 1.0;
  sec.currentTimeStepSum   = 2.0*sec.currentTimeStep;

  sec.lastTimeStep      = sec.currentTimeStep;
  sec.lastTimeStepRatio = sec.currentTimeStepRatio;
  sec.lastTimeStepSum   = sec.currentTimeStepSum;

  sec.numberSuccessiveFailures = 0;
  sec.stepAttemptStatus        = true;

#ifdef Xyce_VERBOSE_TIME
  if (tiaParams.errorAnalysisOption == 1)
  {
      cout << "ERROROPTION=1:  DeltaT Grow = 2" << "\n" << endl;
      cout << "ERROROPTION=1:  DeltaT Cut = 0.125" << "\n" << endl;
      cout << "ERROROPTION=1:  NL MIN = " << tiaParams.NLmin << "\n" << endl;
      cout << "ERROROPTION=1:  NL MAX = " << tiaParams.NLmax << "\n" << endl;
      cout << "ERROROPTION=1:  DELMAX = " << sec.maxTimeStep << "\n" << endl;
  } 
#endif //Xyce_VERBOSE_TIME
 
//  sec.tolAimFac_ = 0.5;

  sec.nextTime = sec.currentTime + sec.currentTimeStep;

//  if (sec.currentTime == sec.initialTime)
  {
  // x history
  *(ds.xHistory[0]) = *(ds.currSolutionPtr);
  *(ds.xHistory[1]) = *(ds.currSolutionPtr);

  // q history
  *(ds.qHistory[0]) = *(ds.daeQVectorPtr);
//  *(ds.qHistory[1]) = *(ds.daeFVectorPtr);
//  (ds.qHistory[1])->scale(-sec.currentTimeStep);
//  (ds.qHistory[1])->putScalar(0.0);
  *(ds.qHistory[1]) = *(ds.daeQVectorPtr);

  // state history
  *(ds.sHistory[0]) = *(ds.currStatePtr);
  (ds.sHistory[1])->putScalar(0.0);

   // lead current Q compontent history
  *(ds.stoLeadCurrQCompHistory[0]) = *(ds.currStoreLeadCurrQCompPtr);
  *(ds.stoLeadCurrQCompHistory[1]) = *(ds.currStoreLeadCurrQCompPtr);

  }
  // Coefficient initialization 
  sec.numberOfSteps_ = 0;    // number of total time integration steps taken
  sec.currentOrder_ = 1;
//  sec.usedOrder_ = 1;
//  if (sec.currentTime == sec.initialTime)
//  {
  sec.usedOrder_ = 1;
  sec.psi_[0] = sec.currentTimeStep;
  sec.cj_ = 1/sec.psi_[0];
//  }
  sec.nscsco_ = 0;
#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 1)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_Gear12::initialize");
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
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_Gear12::setTwoLevelTimeInfo
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/01/07
//-----------------------------------------------------------------------------
void N_TIA_Gear12::setTwoLevelTimeInfo
  (const N_TIA_TimeIntInfo & tiInfo)
{
  // set initial step-size
  double time_to_stop = sec.stopTime - sec.currentTime;


  // x history
  *(ds.xHistory[0]) = *(ds.currSolutionPtr);
  *(ds.xHistory[1]) = *(ds.currSolutionPtr); 

  // q history
  *(ds.qHistory[0]) = *(ds.daeQVectorPtr);
  *(ds.qHistory[1]) = *(ds.daeQVectorPtr);

  // state history
  *(ds.sHistory[0]) = *(ds.currStatePtr);
  (ds.sHistory[1])->putScalar(0.0); 

  // Coefficient initialization 
  sec.numberOfSteps_ = 0;    // number of total time integration steps taken
  sec.usedOrder_ = 1;
  sec.psi_[0] = sec.currentTimeStep;
  sec.cj_ = 1/sec.psi_[0];
  sec.nscsco_ = 0;

}


//-----------------------------------------------------------------------------
// Function      : N_TIA_Gear12::checkReduceOrder()
// Purpose       : check whether to reduce order independent of local error test
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07 
//-----------------------------------------------------------------------------
void N_TIA_Gear12::checkReduceOrder()
{

}

//-----------------------------------------------------------------------------
// Function      : N_TIA_Gear12::rejectStep()
// Purpose       : code to restore history, choose new order/step-size
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
void N_TIA_Gear12::rejectStep()
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

tiaParams.TimeStepLimitedbyBP = false;

#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";

  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_Gear12::rejectStep");
  }
#endif

  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = (!(tiaParams.constantStepSize) );

  sec.lastAttemptedTimeStep = sec.currentTimeStep;


  double newTimeStep_ = sec.currentTimeStep;
  double rr = 1.0; // step size ratio = new step / old step
  if ((sec.stepAttemptStatus == false) && (adjustStep))
  {
    if (tiaParams.errorAnalysisOption == 1)
    {
       newTimeStep_ = sec.currentTimeStep/8;
    }
    else
    {
      sec.initialPhase_ = false;
      sec.nef_++;
      restoreHistory();
      // restore sec.psi_
  //    for (int i=1;i<=sec.currentOrder_;++i)
  //      sec.psi_[i-1] = sec.psi_[i] - sec.currentTimeStep;

      if (sec.nef_ >= sec.max_LET_fail_)  
      {
        string msg = "N_TIA_Gear12::rejectStep: ";
        msg += "  Maximum number of local error test failures.  ";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
      }

      if ((sec.newtonConvergenceStatus <= 0))
      {
        /// 11/11/05 erkeite:  If the Newton solver fails, don't 
        // rely on the error estimate - it may be full of Nan's.
//        rr = sec.r_min_;
       
        newTimeStep_ = sec.currentTimeStep/8;
//        sec.currentOrder_ = 1;
         sec.currentOrder_ = sec.minOrder_;
  //      if (sec.nef_ > 2) sec.newOrder_ = 1;//consistent with block below.
      }
      else
      {
        // 03/11/04 tscoffe:  Here is the block for choosing order & 
        // step-size when the local error test FAILS (but Newton 
        // succeeded). 
        if (sec.nef_ == 1) // first local error test failure
        {
    //	sec.estOverTol_ 
    sec.Est_ = sec.estOverTol_;
//#ifndef Xyce_USE_Q_NORM
//          rr = sec.r_factor_*pow(sec.r_safety_*(sec.Est_+sec.r_fudge_),-1.0/(sec.newOrder_+1.0));
//          rr = Xycemax(sec.r_min_,Xycemin(sec.r_max_,rr));
//#else
 //         rr = sec.r_factor_*pow(sec.r_safety_*sec.Est_+sec.r_fudge_,-1.0/(sec.newOrder_+1.0)); 

          rr = sec.tolAimFac_/(sec.estOverTol_ + 0.0001);
          rr = pow(rr, 1.0/(sec.currentOrder_+1.0));
          rr = Xycemax(sec.r_min_,Xycemin(sec.r_max_,rr));
//#endif
          newTimeStep_ = rr * sec.currentTimeStep;

//          sec.currentOrder_ = 1; 
//          sec.currentOrder_ = sec.minOrder_;
        }
        else // if (sec.nef_ == 2) // second Dae failure
        {
          rr = sec.r_min_;
          newTimeStep_ = rr * sec.currentTimeStep;
         
//          sec.currentOrder_ = 1;
          sec.currentOrder_ = sec.minOrder_;
        }
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
  }
  else if ((sec.stepAttemptStatus == false) & (!adjustStep))
  {
    string tmp = "  Gear12:rejectStep: Warning: Local error test failed with constant step-size.\n";
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
      tiaParams.TimeStepLimitedbyBP = true;
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
// Function      : N_TIA_Gear12::rejectStepForHabanero
// Purpose       : step rejection, but from an outside program (Habanero API)
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 08/11/09
//-----------------------------------------------------------------------------
void N_TIA_Gear12::rejectStepForHabanero()
{
  restoreHistory();
  sec.setTimeStep(sec.currentTimeStep);
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_Gear12::completeStep()
// Purpose       : code to update history, choose new order/step-size
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07 
//-----------------------------------------------------------------------------
void N_TIA_Gear12::completeStep()
{
  
  tiaParams.TimeStepLimitedbyBP = false;

  sec.numberOfSteps_ ++;
  sec.nef_ = 0;
  sec.lastTime    = sec.currentTime;
  sec.currentTime = sec.nextTime;


#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";

  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_Gear12::completeStep");
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


/*   if (sec.numberOfSteps_ >= 2)
   {
     sec.currentOrder_ = 2;
//     (ds.relErrTolPtr)->putScalar(1e-2);	
   } 
*/ 
    if (tiaParams.errorAnalysisOption == 1)
    {
      
      if (sec.numberOfSteps_ >= 2 &&  sec.maxOrder_ == 2)
      {
        sec.currentOrder_ = 2;   
      } 
     
      rr = 1;

      if (sec.nIterations <= tiaParams.NLmin)
        rr = 2;  
      
      if (sec.nIterations > tiaParams.NLmax)
        rr = 1.0/8;
      
      newTimeStep_ = rr*sec.currentTimeStep;
    }
    else
    {
//#ifndef Xyce_USE_Q_NORM
//      rr = pow(sec.r_safety_*(sec.Est_+sec.r_fudge_),-1.0/(sec.currentOrder_+1.0));
//#else
//      rr = pow(sec.r_safety_*sec.Est_+sec.r_fudge_,-1.0/(sec.currentOrder_+1.0));
//#endif


      rr = sec.tolAimFac_/(sec.estOverTol_ + 0.0001);
      rr = pow(rr, 1.0/(sec.currentOrder_+1.0));

      if (sec.numberOfSteps_ >= 2 && sec.maxOrder_ == 2)
      {
        if (sec.currentOrder_ == 1)
        { 
          sec.currentOrder_ = 2;   
          rr = sec.tolAimFac_/(sec.estOverTol_ + 0.0001);
          rr = pow(rr, 1.0/(sec.currentOrder_+1.0)); 

          if (rr <= 1.05)
          {
//            sec.currentOrder_ = 1; 
            sec.currentOrder_ = sec.minOrder_;
          }
        }
      } 
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

   updateHistory();

   newTimeStep_ = Xycemax(newTimeStep_, sec.minTimeStep);
   newTimeStep_ = Xycemin(newTimeStep_, sec.maxTimeStep);
   
  // 12/01/05 tscoffe:  This is necessary to avoid currentTimeStep == 0 right
  // before a breakpoint.  So I'm checking to see if currentTime is identically
  // equal to stopTime, in which case we are right before a breakpoint and we
  // should not adjust currentStepSize because that would result in
  // currentStepSize == 0.
//  if (sec.currentTime < sec.stopTime)
  if ((sec.stopTime - sec.currentTime) >= sec.minTimeStep) 
  {
    // If the step needs to be adjusted:
    if (adjustStep)
    {
//      newTimeStep_ = Xycemax(newTimeStep_, sec.minTimeStep);
//      newTimeStep_ = Xycemin(newTimeStep_, sec.maxTimeStep);

      double nextTimePt = sec.currentTime + newTimeStep_;

      if (nextTimePt > sec.stopTime)
      {

        sec.savedTimeStep = newTimeStep_;
        
        nextTimePt  = sec.stopTime;
        newTimeStep_ = sec.stopTime - sec.currentTime;
        tiaParams.TimeStepLimitedbyBP = true;
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

//      sec.currentTimeStep = newTimeStep_;
    }
    else // if time step is constant for this step:
    {
      double nextTimePt = sec.currentTime + sec.currentTimeStep;

      if (nextTimePt > sec.stopTime)
      {
        sec.savedTimeStep = sec.currentTimeStep;
        
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

  sec.currentTimeStep = newTimeStep_;
// 11/02/04 tscoffe:  This should be done at the beginning of an integration
//                    stop rather than at the end, so its been moved into
//                    ControlAlgorithm::transientLoop_
  // Update next stop time in StepErrorControl:
  // sec.updateStopTime();
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_Gear12::updateStateDeriv
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/17/05
//-----------------------------------------------------------------------------
void N_TIA_Gear12::updateStateDeriv ()
{

  ds.nextStateDerivPtr->linearCombo(sec.alpha_[0],*ds.nextStatePtr, sec.alpha_[1],*ds.sn0Ptr);


  if (sec.currentOrder_ == 2)
  {
    ds.nextStateDerivPtr->linearCombo(1.0,*ds.nextStateDerivPtr, sec.alpha_[2], *(ds.sHistory[1]));
  }

  ds.nextStateDerivPtr->scale(1.0/sec.currentTimeStep);
  
#ifdef Xyce_DEBUG_TIME
    if (tiaParams.debugLevel > 1)
    {
    cout << "\n next state Ptr: \n" << endl;
    ds.nextStatePtr->printPetraObject();
    cout << endl;

    cout << "\n sn0: \n" << endl;
    ds.sn0Ptr->printPetraObject();
    cout << endl;

    cout << "\n next State Deriv: \n" << endl;
    ds.nextStateDerivPtr->printPetraObject();
    cout << endl;
    
    }
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_OneStep::updateLeadCurrent
// Purpose       : calculates lead currents in store vector with 
//                 the storeLeadCurrQCompVec. 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling, SNL
// Creation Date : 03/22/2013
//-----------------------------------------------------------------------------
void N_TIA_Gear12::updateLeadCurrent ()
{
  ds.nextStoreLeadCurrQCompDerivPtr->linearCombo(sec.alpha_[0],*ds.nextStoreLeadCurrQCompPtr, sec.alpha_[1],*ds.stoQCn0Ptr);

  if (sec.currentOrder_ == 2)
  {
    ds.nextStoreLeadCurrQCompDerivPtr->linearCombo(1.0,*ds.nextStoreLeadCurrQCompDerivPtr, sec.alpha_[2], *(ds.stoLeadCurrQCompHistory[1]));
  }

  ds.nextStoreLeadCurrQCompDerivPtr->scale(1.0/sec.currentTimeStep);

  ds.nextStorePtr->linearCombo(1.0,*ds.nextStorePtr,1.0,*ds.nextStoreLeadCurrQCompDerivPtr);
#ifdef Xyce_DEBUG_TIME
    if (tiaParams.debugLevel > 1)
    {
    cout << "\n next store Ptr: \n" << endl;
    ds.nextStorePtr->printPetraObject();
    cout << endl;
    
    }
#endif

}

//-----------------------------------------------------------------------------
// Function      : N_TIA_Gear12::getInitialQnorm
// Purpose       : Needed by 2-level solves.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/18/07
//-----------------------------------------------------------------------------
void N_TIA_Gear12::getInitialQnorm (N_TIA_TwoLevelError & tle)
{
  tle.q1HistorySum = ds.partialSum_q1();
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_Gear12::setupTwoLevelError
// Purpose       : Needed by 2-level solves.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
void N_TIA_Gear12::setupTwoLevelError(N_TIA_TwoLevelError & tle)
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

