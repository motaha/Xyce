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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_NumericalJacobian.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 04/30/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.67 $
//
// Revision Date  : $Date: 2014/02/24 23:49:15 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------   Standard Includes   ----------
#include <N_UTL_Misc.h>
#include <string>
#include <iostream>

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif


#include <sstream>
// ----------   Xyce Includes   ----------
#include <N_DEV_NumericalJacobian.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceMgr.h>

#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

// ---------- Static Initializations ----------


namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::NumericalJacobian
// Purpose       : block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/30/02
//-----------------------------------------------------------------------------
NumericalJacobian::NumericalJacobian(
  MatrixLoadData & mlData1,
  const SolverState &ss1,
  const ExternData  &ed1,
  const DeviceOptions & do1)
  : mlData(mlData1),
    cols(mlData1.cols),
    vals(mlData1.vals),
    Qvals(mlData1.Qvals),
    val_local(mlData1.val_local),
    Qval_local(mlData1.Qval_local),
    col_local(mlData1.col_local),
    row_local(mlData1.row_local),
    internalFlag(mlData1.internalFlag),
    solState (ss1),
    extData  (ed1),
    devOptions(do1),
    maxCols(10) // guess
{}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::NumericalJacobian
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/30/02
//-----------------------------------------------------------------------------
NumericalJacobian::NumericalJacobian(const NumericalJacobian & right)
  : mlData(right.mlData),
    cols(right.cols),
    vals(right.vals),
    Qvals(right.Qvals),
    val_local(right.val_local),
    Qval_local(right.Qval_local),
    col_local(right.col_local),
    row_local(right.row_local),
    internalFlag(right.internalFlag),
    solState (right.solState),
    extData  (right.extData),
    devOptions(right.devOptions),
    maxCols(right.maxCols)
{

}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::~NumericalJacobian
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/30/00
//-----------------------------------------------------------------------------
NumericalJacobian::~NumericalJacobian()
{

}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::testDAEMatrices
//
// Purpose       : Performs a numerical jacobian test on the dFdx and dQdx
//                 matrices.
//
// Special Notes : This is the main numerical jacoabian test.  There is
//                 another numerical jacobian test, which pre-dates this
//                 one, but is not recommended.  The older one is handled
//                 by the functions preceding this one in this file.
//
//                 The old one attempt to run the code using a purely
//                 numerical jacobian.  It does not work with votlage limiting.
//
//                 The new one (controlled from this function) operates along
//                 side the analytic jacobian, but does not replace it, as
//                 the solvers will still use the analytic jacobian.  This
//                 function will test, on a device-by-device basis, the device
//                 contributions to the analytic jaocian.  That way,
//                 it is possible to test a single device (rather than the whole
//                 problem), and test it at a precise point during the simulation.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/15/06
//-----------------------------------------------------------------------------
bool NumericalJacobian::testDAEMatrices( DeviceInstance & instance, const std::vector<std::string> & nameVec)
{

  // Set up various references to indexing arrays
  const std::vector<int> & devLIDs               = instance.getDevLIDs();
  const std::vector<int> & devStateLIDs          = instance.getStaLIDVec();
  const std::vector< std::vector<int> > & devJacLIDs = instance.getDevJacLIDs();
  const std::vector< std::vector<int> > & jacStamp   = instance.jacobianStamp();


  // Set up references to temporary data structures.
  std::vector< std::vector<double> > & numJacF = mlData.numJac;
  std::vector< std::vector<double> > & saveJacF = mlData.saveJac;
  std::vector< std::vector<double> > & devJacF = mlData.devJac;
  std::vector< std::vector<double> > & diffJacF = mlData.diffJac;
  std::vector< std::vector<double> > & relJacF = mlData.relJac;

  std::vector< std::vector<double> > & numJacQ = mlData.numJacQ;
  std::vector< std::vector<double> > & saveJacQ = mlData.saveJacQ;
  std::vector< std::vector<double> > & devJacQ = mlData.devJacQ;
  std::vector< std::vector<double> > & diffJacQ = mlData.diffJacQ;
  std::vector< std::vector<double> > & relJacQ = mlData.relJacQ;

  std::vector< std::vector<int> > & statusF  = mlData.status;
  std::vector< std::vector<int> > & statusQ  = mlData.statusQ;
  std::vector< std::vector<int> > & stencil = mlData.stencil;

  std::vector<double> & saveF = mlData.saveRHS;
  std::vector<double> & pertF = mlData.pertRHS;
  std::vector<double> & origF = mlData.origRHS;
  std::vector<double> & saveQ = mlData.saveQ;
  std::vector<double> & pertQ = mlData.pertQ;
  std::vector<double> & origQ = mlData.origQ;

  std::vector<double> & saveSoln = mlData.saveSoln;
  std::vector<double> & pertSoln = mlData.pertSoln;
  std::vector<double> & saveCurrSoln = mlData.saveCurrSoln;

  std::vector<double> & saveLastState = mlData.saveLastState;
  std::vector<double> & saveCurrState = mlData.saveCurrState;
  std::vector<double> & saveNextState = mlData.saveNextState;
  std::vector<double> & saveStateDerivs = mlData.saveStateDerivs;

  // set up references to epetra objects.
  N_LAS_Vector & Fvec          = (*extData.daeFVectorPtr);
  N_LAS_Vector & Qvec          = (*extData.daeQVectorPtr);

  N_LAS_Vector & currSol      = (*extData.currSolVectorPtr);
  N_LAS_Vector & nextSol      = (*extData.nextSolVectorPtr);

  N_LAS_Vector & lastSta      = (*extData.lastStaVectorPtr);
  N_LAS_Vector & currSta      = (*extData.currStaVectorPtr);
  N_LAS_Vector & nextSta      = (*extData.nextStaVectorPtr);
  N_LAS_Vector & nextStaDeriv = (*extData.nextStaDerivVectorPtr);

  N_LAS_Matrix & dQdxMat      = (*extData.dQdxMatrixPtr);
  N_LAS_Matrix & dFdxMat      = (*extData.dFdxMatrixPtr);

  int numRows, numCols;
  numCols = devLIDs.size();
  numRows = jacStamp.size();

  int testSize= (numRows>numCols)?numRows:numCols;
  if (testSize > saveF.size())
  {
    mlData.resizeTestJacSolData(testSize);
    mlData.resizeTestJacQData(testSize);
  }

  int numState = devStateLIDs.size();
  if (numState > saveCurrState.size())
  {
    mlData.resizeTestJacStateData(numState);
  }

  int i, j, jCol;

  if(devJacLIDs.empty())
  {
#ifdef Xyce_DEBUG_DEVICE
    Report::UserWarning() << instance.getName() << " does not have jacLIDs available";
#endif

    return true;
  }

  if (instance.getOrigFlag() && numRows > 0 && numRows == jacStamp.size())
  {

    // Zero out all the mlData structures
    saveF.assign(saveF.size(),0.0);
    pertF.assign(pertF.size(),0.0);
    origF.assign(origF.size(),0.0);

    saveQ.assign(saveQ.size(),0.0);
    pertQ.assign(pertQ.size(),0.0);
    origQ.assign(origQ.size(),0.0);

    saveSoln.assign(saveSoln.size(),0.0);
    pertSoln.assign(pertSoln.size(),0.0);
    saveCurrSoln.assign(saveCurrSoln.size(),0.0);

    saveLastState.assign(saveLastState.size(),0.0);
    saveCurrState.assign(saveCurrState.size(),0.0);
    saveNextState.assign(saveNextState.size(),0.0);
    saveStateDerivs.assign(saveStateDerivs.size(),0.0);
    for (i=0;i<numJacF.size();++i)
    {
      numJacF[i].assign(numJacF[i].size(),0.0);
      saveJacF[i].assign(saveJacF[i].size(),0.0);
      devJacF[i].assign(devJacF[i].size(),0.0);
      diffJacF[i].assign(diffJacF[i].size(),0.0);
      statusF[i].assign(statusF[i].size(),-1);
      stencil[i].assign(stencil[i].size(),0);

      numJacQ[i].assign(numJacQ[i].size(),0.0);
      saveJacQ[i].assign(saveJacQ[i].size(),0.0);
      devJacQ[i].assign(devJacQ[i].size(),0.0);
      diffJacQ[i].assign(diffJacQ[i].size(),0.0);
      statusQ[i].assign(statusQ[i].size(),-1);
    }

    // Save Soln, RHS, and State for this device
    bool origFlag = instance.getOrigFlag();
    //for (i=0 ; i<numRows ; ++i)
    int tmpSize= (numRows>numCols)?numRows:numCols;
    double sqrtEta=devOptions.testJac_SqrtEta;
    for (i=0 ; i<tmpSize; ++i)
    {
      saveF[i]      = Fvec[devLIDs[i]];
      saveQ[i]      = Qvec[devLIDs[i]];

      saveSoln[i]     = nextSol[devLIDs[i]];
      saveCurrSoln[i] = currSol[devLIDs[i]];
      pertSoln[i]     = sqrtEta * (1.0 + fabs(saveSoln[i]));
    }

    for (i=0 ; i<numState ; ++i)
    {
      saveLastState[i] = lastSta[devStateLIDs[i]];
      saveCurrState[i] = currSta[devStateLIDs[i]];
      saveNextState[i] = nextSta[devStateLIDs[i]];
      saveStateDerivs[i] = nextStaDeriv[devStateLIDs[i]];
    }

    // Save the original matrix for later:
    for (i=0 ; i<numRows ; ++i)
    {
      jCol = devJacLIDs[i].size();
      for (j=0 ; j<jCol ; ++j)
      {
        double valF = dFdxMat[devLIDs[i]][devJacLIDs[i][j]];
        saveJacF[i][j] = valF;
        double valQ = dQdxMat[devLIDs[i]][devJacLIDs[i][j]];
        saveJacQ[i][j] = valQ;
      }
    }

    // Zeroing needs to be done after all saved values are
    // recorded because there can be multiple references
    // to the same element
    for (i=0 ; i<numRows ; ++i)
    {
      jCol = devJacLIDs[i].size();
      for (j=0 ; j<jCol ; ++j)
      {
        dFdxMat[devLIDs[i]][devJacLIDs[i][j]] = 0;
        dQdxMat[devLIDs[i]][devJacLIDs[i][j]] = 0;
      }
    }

    // Now that the original load has been zeroed out, re-load the
    // analytic contributions, to get the contributions from *just* this
    // device.
    instance.loadDAEdQdx ();
    instance.loadDAEdFdx ();

    for (i=0 ; i<numRows ; ++i)
    {
      devJacF[i].assign(devJacF[i].size(),0.0);
      devJacQ[i].assign(devJacQ[i].size(),0.0);
      stencil[i].assign(stencil[i].size(),0);

      jCol = devJacLIDs[i].size();
      for (j=0 ; j<jCol ; ++j)
      {
        double valF = dFdxMat[devLIDs[i]][devJacLIDs[i][j]];
        devJacF[i][jacStamp[i][j]] = valF;
        double valQ = dQdxMat[devLIDs[i]][devJacLIDs[i][j]];
        devJacQ[i][jacStamp[i][j]] = valQ;
        stencil[i][jacStamp[i][j]] = 1;
      }
    }

    // zero out the RHS, and re-load, so that we have only the
    // elements from one device present.
    for (i=0 ; i<numRows ; ++i)
    {
      Fvec[devLIDs[i]] = 0.0;
      Qvec[devLIDs[i]] = 0.0;
    }
    // re-load for just this instance:
    loadLocalDAEVectors (instance);

    // Save RHS for just this instance:
    for (i=0 ; i<numRows ; ++i)
    {
      origF[i] = Fvec[devLIDs[i]];
      origQ[i] = Qvec[devLIDs[i]];
    }

    for (i=0 ; i<numRows ; ++i)
    {
      jCol = devJacLIDs[i].size();
      for (j=0 ; j<jCol ; ++j)
      {
        statusF[i][jacStamp[i][j]] = -1;
        statusQ[i][jacStamp[i][j]] = -1;
      }
    }

    // This test, of the merged numRows/numCols, should only be needed 1 time.
    // Leaving this off by default for now.  For reasons, see the header
    // for the mergeTest function.
    if (false)
    {
      mergeTest (instance, nameVec);
    }

    // These are the tolerances for determining jacobian agreement
    //      relT: similar to reltol, fractional error that is allowed
    //      absT: similar to abstol, error that is allowed for small derivatives
    double relTol=devOptions.testJac_relTol;
    double absTol=devOptions.testJac_absTol;

    double ndFdx, adFdx, ddFdx, relError_dFdx;
    double ndQdx, adQdx, ddQdx, relError_dQdx;

    bool failedTest = false;
    for (i=0 ; i<numCols ; ++i)
    {
      // Don't bother perturbing gnd.
      if (nameVec[devLIDs[i]] == "gnd") continue;

      for (j=0 ; j<numRows ; ++j)
      {
        nextSol[devLIDs[j]] = saveSoln[j];
        currSol[devLIDs[j]] = saveCurrSoln[j];
        Fvec[devLIDs[j]] = 0.0;
        Qvec[devLIDs[j]] = 0.0;
      }

      for (j=0 ; j<numState ; ++j)
      {
        lastSta[devStateLIDs[j]] = saveLastState[j];
        currSta[devStateLIDs[j]] = saveCurrState[j];
        nextSta[devStateLIDs[j]] = saveNextState[j];
        nextStaDeriv[devStateLIDs[j]] = saveStateDerivs[j];
      }

      // Perturb the solution.
      double dX = pertSoln[i];
      nextSol[devLIDs[i]] += dX;

      // Re-load the F,Q vectors:
      loadLocalDAEVectors (instance);

      for (j=0 ; j<numRows ; ++j)
      {
        pertF[j] = Fvec[devLIDs[j]];
        pertQ[j] = Qvec[devLIDs[j]];
      }

#ifdef Xyce_DEBUG_TESTJAC
      testDebugHead (instance, nameVec, i, dX);
#endif

      for (j=0 ; j<numRows ; ++j)
      {
        // Don't bother with the gnd row.
        if (nameVec[devLIDs[j]] == "gnd") continue;

        // if this derivative is not loaded analytically, don't bother
        // with it.
        if (stencil[j][i]!=1) continue;

        double dF = (pertF[j]-origF[j]);
        double dQ = (pertQ[j]-origQ[j]);
        numJacF[j][i] = dF/dX;
        numJacQ[j][i] = dQ/dX;

        ndFdx = numJacF[j][i];
        adFdx = devJacF[j][i];
        ddFdx = fabs(adFdx-ndFdx);
        relError_dFdx = ddFdx/(relTol*fabs(ndFdx)+absTol);

        diffJacF[j][i] = ddFdx;
        relJacF[j][i] = relError_dFdx;

        ndQdx = numJacQ[j][i];
        adQdx = devJacQ[j][i];
        ddQdx = fabs(adQdx-ndQdx);
        relError_dQdx = ddQdx/(relTol*fabs(ndQdx)+absTol);

        diffJacQ[j][i] = ddQdx;
        relJacQ[j][i] = relError_dQdx;

#ifdef Xyce_DEBUG_TESTJAC
        if (isnan(numJacF[j][i]))
        {
          testDebugOut (instance, nameVec, i, j);
        }
#endif
        // if the device is a Inductor, and IC= has been specified,
        // then skip this term as it is a special case.
        ExtendedString varNameI(nameVec[devLIDs[i]]); varNameI.toUpper();
        ExtendedString varNameJ(nameVec[devLIDs[j]]); varNameJ.toUpper();
        if ( ((solState.dcopFlag) && varNameI[0]=='L' && varNameJ[0]=='L') )
        {
          // For the inductor branch current, the matrix element will be
          // there whether IC= was specified or not.
          if (adFdx == 1 && ndFdx == 0)
          {
            statusF[j][i] = 3;
            statusQ[j][i] = 3;
          }
        }
        // if the device is a capacitor, and it has a branch current,
        // that means that IC= has been used, thus it is a "special case"
        // that should be skipped. The branch current will only be there
        // for IC=.
        else if((!(solState.dcopFlag) && varNameI[0]=='C') && varNameJ[0]=='C')
        {
          statusF[j][i] = 3;
          statusQ[j][i] = 3;
        }

        if ( statusF[j][i] != 3 )
        {
          if (relError_dFdx > 1.0) // failure
          {
            statusF[j][i] = -2;
            failedTest = true;
          }
          else // success
          {
            statusF[j][i] = 1;
          }

          if (relError_dQdx > 1.0) // failure
          {
            statusQ[j][i] = -2;
            failedTest = true;
          }
          else // success
          {
            statusQ[j][i] = 1;
          }
        }
      }

#ifdef Xyce_DEBUG_TESTJAC
      testDebugTail (instance, nameVec);
#endif
    }

#ifdef Xyce_DEBUG_DEVICE
    // Output Jacobians Differences.
    // If debug enabled, always output.
    // If not, only output for failures.
    printJacobian_ (dout(), instance, nameVec, failedTest);
#else
    if (failedTest)
      printJacobian_ (lout(), instance, nameVec, failedTest);
#endif

    // Restore jacobian, RHS for this device
    instance.setOrigFlag(origFlag);
    tmpSize= (numRows>numCols)?numRows:numCols;
    for (i=0 ; i<tmpSize; ++i)
    {
      Fvec[devLIDs[i]] = saveF[i];
      Qvec[devLIDs[i]] = saveQ[i];
      nextSol[devLIDs[i]] = saveSoln[i];
      currSol[devLIDs[i]] = saveCurrSoln[i];
    }
    for (i=0 ; i<numState ; ++i)
    {
      lastSta[devStateLIDs[i]] = saveLastState[i];
      currSta[devStateLIDs[i]] = saveCurrState[i];
      nextSta[devStateLIDs[i]] = saveNextState[i];
      nextStaDeriv[devStateLIDs[i]] = saveStateDerivs[i];
    }


    for (i=0 ; i<numRows ; ++i)
    {
      jCol = devJacLIDs[i].size();
      for (j=0 ; j<jCol ; ++j)
      {
        dFdxMat[devLIDs[i]][devJacLIDs[i][j]] = saveJacF[i][j];
        dQdxMat[devLIDs[i]][devJacLIDs[i][j]] = saveJacQ[i][j];
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::loadLocalDAEVectors
// Purpose       :
// Special Notes : Note this function ignores the B-vector.
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/15/06
//-----------------------------------------------------------------------------
void NumericalJacobian::loadLocalDAEVectors (DeviceInstance & instance)
{
  N_LAS_Vector & currSta      = (*extData.currStaVectorPtr);
  N_LAS_Vector & nextSta      = (*extData.nextStaVectorPtr);
  N_LAS_Vector & nextStaDeriv = (*extData.nextStaDerivVectorPtr);
  N_LAS_Vector & nextSol      = (*extData.nextSolVectorPtr);

  const std::vector<int> & devStateLIDs = instance.getStaLIDVec();
  int numState = devStateLIDs.size();

  instance.updateDependentParameters(nextSol); // this line necessary for expressions
  instance.updatePrimaryState ();

  // Assume backward euler integration, so that the time integrator
  // accessors are not needed.
  //anaIntPtr_->updateDivDiffs(devLIDs, devStateLIDs);
  //anaIntPtr_->updateDerivs(devLIDs, devStateLIDs);
  for (int j=0 ; j<numState ; ++j)
  {
    nextStaDeriv[devStateLIDs[j]] =
      solState.pdt * (nextSta[devStateLIDs[j]]-currSta[devStateLIDs[j]]);
  }

  instance.updateSecondaryState ();
  instance.loadDAEQVector ();
  instance.loadDAEFVector ();

  return;
}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::printJacobian_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/12/06
//-----------------------------------------------------------------------------
void NumericalJacobian::printJacobian_(
  std::ostream &                        os,
  DeviceInstance &                      instance,
  const std::vector<std::string> &      nameVec,
  bool                                  failed)
{
  bool NAflag = false;
  const std::vector<int> & devLIDs               = instance.getDevLIDs();

  // These curly brackets are for scoping only.
  {
    const std::vector< std::vector<double> > & numJac = mlData.numJac;
    const std::vector< std::vector<double> > & anaJac = mlData.devJac;
    const std::vector< std::vector<double> > & diffJac = mlData.diffJac;
    const std::vector< std::vector<double> > & relJac = mlData.relJac;
    const std::vector< std::vector<int> > & stencil = mlData.stencil;
    const std::vector< std::vector<int> > & status = mlData.status;


    os << Xyce::section_divider << std::endl;
    os << "dFdx matrix for " << instance.getName();

    if (failed)
    {
      os << ":  JACOBIAN TEST FAILURE";
    }
    else
    {
      os << ":  JACOBIAN TEST SUCCESS";
    }
    os << " at time = " << solState.currTime;
    os << " at Niter = " << solState.newtonIter;
    os << std::endl;
    os << "       Numerical     Analytic      absDiff       relative Error  Status   Names  (row, col)"<<std::endl;

    int i,j;
    int numCols = devLIDs.size();
    int numRows = (status.size()<numCols)?status.size():numCols;

    for (i = 0; i < numRows; ++i)
    {
      if (nameVec[devLIDs[i]] == "gnd") continue;

      for (j = 0; j < numCols; ++j)
      {
        if (nameVec[devLIDs[j]] == "gnd") continue;

        // if this variable has not been tested for any reason, skip.
        if (status[i][j]==-1) continue;

        // if this derivative is not loaded analytically, skip.
        if (stencil[i][j]!=1) continue;

        // Note: JT=jacobian test is there to make this easy to grep.
        static char tmpChar[128];
        static char prefix[4];

        sprintf(prefix,"%s","FT:");

        if (status[i][j]==-2)
        {
          sprintf(tmpChar,"%s  %12.4e  %12.4e  %12.4e  %12.4e      fail",prefix,
                  numJac[i][j], anaJac[i][j], diffJac[i][j], relJac[i][j]);
        }
        else if (status[i][j]==3)
        {
          NAflag = true;
          sprintf(tmpChar,"%s  %12.4e  %12.4e  %12.4e  %12.4e      NA  ",prefix,
                  numJac[i][j], anaJac[i][j], diffJac[i][j], relJac[i][j]);
        }
        else
        {
          sprintf(tmpChar,"%s  %12.4e  %12.4e  %12.4e  %12.4e          ",prefix,
                  numJac[i][j], anaJac[i][j], diffJac[i][j], relJac[i][j]);
        }

        os << std::string(tmpChar);

        os << "     ("<< nameVec[devLIDs[i]]
            << ", " << nameVec[devLIDs[j]]
            << ") "
            << " row,col=[ " << i << ", " << j << "]" << std::endl;


      }
    }
  }

  const std::vector< std::vector<double> > & numJac = mlData.numJacQ;
  const std::vector< std::vector<double> > & anaJac = mlData.devJacQ;
  const std::vector< std::vector<double> > & diffJac = mlData.diffJacQ;
  const std::vector< std::vector<double> > & relJac = mlData.relJacQ;
  const std::vector< std::vector<int> > & stencil = mlData.stencil;
  const std::vector< std::vector<int> > & status = mlData.statusQ;

  os << "dQdx matrix for " << instance.getName();

  if (failed)
  {
    os << ":  JACOBIAN TEST FAILURE";
  }
  else
  {
    os << ":  JACOBIAN TEST SUCCESS";
  }
  os << " at time = " << solState.currTime;
  os << std::endl;
  os << "       Numerical     Analytic      absDiff       relative Error  Status   Names  (row, col)"<<std::endl;

  int i,j;
  int numCols = devLIDs.size();
  int numRows = (status.size()<numCols)?status.size():numCols;

  for (i = 0; i < numRows; ++i)
  {
    if (nameVec[devLIDs[i]] == "gnd") continue;

    for (j = 0; j < numCols; ++j)
    {
      if (nameVec[devLIDs[j]] == "gnd") continue;

      // if this variable has not been tested for any reason, skip.
      if (status[i][j]==-1) continue;

      // if this derivative is not loaded analytically, skip.
      if (stencil[i][j]!=1) continue;

      // Note: QT=jacobian test is there to make this easy to grep.
      static char tmpChar[128];
      if (status[i][j]==-2)
      {
        sprintf(tmpChar,"QT:  %12.4e  %12.4e  %12.4e  %12.4e      fail",
                numJac[i][j], anaJac[i][j], diffJac[i][j], relJac[i][j]);
      }
      else if (status[i][j]==3)
      {
        NAflag = true;
        sprintf(tmpChar,"QT:  %12.4e  %12.4e  %12.4e  %12.4e      NA  ",
                numJac[i][j], anaJac[i][j], diffJac[i][j], relJac[i][j]);
      }
      else
      {
        sprintf(tmpChar,"QT:  %12.4e  %12.4e  %12.4e  %12.4e          ",
                numJac[i][j], anaJac[i][j], diffJac[i][j], relJac[i][j]);
      }

      os << std::string(tmpChar);

      os << "     ("<< nameVec[devLIDs[i]]
          << ", " << nameVec[devLIDs[j]]
          << ") "
          << " row,col=[ " << i << ", " << j << "]" << std::endl;

    }
  }

  if(NAflag)
    os << " Note:  NA = untestable special case, such as IC=, etc." << std::endl;
  os << Xyce::section_divider << std::endl;

  if (failed)
  {
    if (!(devOptions.testJacWarn))
    {
      Report::UserError() << "Numerical Jacobian test failure" << std::endl
                          << "If you want this failure to be a warning, rather than an error, "
                          << "run with .options device testjacwarn=1 in the netlist.";
    }
    else
    {
      Report::UserWarning() << "Numerical Jacobian test failure";
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::mergeTest
//
// Purpose       : ERK:  This function tests merged rows and cols, to make sure
//                 that they are actually equivalent.  This is adapted from
//                 a version of this diagnostic by Dave Shirley.
//
//                 Another purpose of this test is to avoid testing the same
//                 element (numerical vs. analytic) more than once.  I've
//                 chosen to use a different simpler approach for this, using
//                 a matrix stencil.  The nice thing about stencils is that
//                 they don't require very many if-statements.
//
//                 With the stencil implemented, the only purpose for this test
//                 is essentially to test that the merged rows and cols are
//                 in fact correctly merged.  I think that is a low-probability
//                 bug, so for now I'm not calling this test function, in
//                 the interests of speed.
//
//                 Finally, I'm leaving this test out, because it relies on
//                 "==" style comparisons of floating point numbers, which can
//                 be dangerous and easy to break.
//
//                 However, for historical reasons I'm leaving it in the code.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/14/06
//-----------------------------------------------------------------------------
void NumericalJacobian::mergeTest( DeviceInstance & instance, const std::vector<std::string> & nameVec)
{
  int i,j,k;
  const std::vector<int> & devLIDs               = instance.getDevLIDs();
  const std::vector< std::vector<double> > & devJac = mlData.devJac;
  std::vector< std::vector<int> > & status = mlData.status;
  int numRows, numCols;
  numRows = devLIDs.size();
  numCols = numRows;

  // insure that for any instance this is only checked once.
  if (!(instance.getMergeRowColChecked()))
  {
    for (i=1 ; i<numRows ; ++i)
    {
      for (j=0 ; j<i ; ++j)
      {
        if (devLIDs[i] == devLIDs[j])
        {
          for (k=0 ; k<numCols ; ++k)
          {
            if (devJac[i][k] != devJac[j][k])
            {
              Report::UserWarning() << "In device " + instance.getName() << " different non-zero values in row merge";
            }
            status[i][k] = 2;
          }
          for (k=0 ; k<numRows ; ++k)
          {
            if (devJac[k][i] != devJac[k][j])
            {
              Report::UserWarning() << "In device " + instance.getName() << " different non-zero values in column merge";
            }
            status[k][i] = 2;
          }
        }
      }
    }
    instance.setMergeRowColChecked(true);
  }

}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::testDebugHead
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/14/06
//-----------------------------------------------------------------------------
void NumericalJacobian::testDebugHead( DeviceInstance & instance, const std::vector<std::string> & nameVec, int i, double dX)
{
  const std::vector<int> & devLIDs = instance.getDevLIDs();

  Xyce::dout() << Xyce::section_divider<<std::endl;
  Xyce::dout() << "Perturbing (LID="<<devLIDs[i]<<") " << nameVec[devLIDs[i]] << " by " << dX << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::testDebugOut
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/14/06
//-----------------------------------------------------------------------------
void NumericalJacobian::testDebugOut( DeviceInstance & instance, const std::vector<std::string> & nameVec, int i, int j)
{
  const std::vector<int> & devLIDs = instance.getDevLIDs();
  const std::vector<double> & pertRHS = mlData.pertRHS;
  const std::vector<double> & origRHS = mlData.origRHS;

  const std::vector< std::vector<double> > & numJac = mlData.numJac;
  const std::vector< std::vector<double> > & relJac = mlData.relJac;

  Xyce::dout().width(15); Xyce::dout().precision(7); Xyce::dout().setf(std::ios::scientific);
  Xyce::dout() << "dFdX: ";
  Xyce::dout() << " (" << devLIDs[j] << ", " << devLIDs[i] << ") ";
  Xyce::dout() << numJac[j][i];
  Xyce::dout() << " Forig = " << origRHS[j];
  Xyce::dout() << " Fperturb = " << pertRHS[j];
  double dF = -(pertRHS[j]-origRHS[j]);
  Xyce::dout() << " dF = " << dF;
  Xyce::dout() << " (" << nameVec[devLIDs[j]] << ", " << nameVec[devLIDs[i]] << ") ";
  Xyce::dout() << std::endl;
  Xyce::dout() << "  relative error = " << relJac[j][i] << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::testDebugTail
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/14/06
//-----------------------------------------------------------------------------
void NumericalJacobian::testDebugTail( DeviceInstance & instance, const std::vector<std::string> & nameVec)
{
  Xyce::dout() << Xyce::section_divider<<std::endl;
}

} // namespace Device
} // namespace Xyce
