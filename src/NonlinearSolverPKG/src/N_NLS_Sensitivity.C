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
// Filename       : $RCSfile: N_NLS_Sensitivity.C,v $
//
// Purpose        : Body for the sensitivity class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/30/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.82 $
//
// Revision Date  : $Date $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------   Standard Includes   ----------

#include <N_UTL_Misc.h>
#include <sstream>

// ----------   Xyce Includes   ----------

#include <N_NLS_Sensitivity.h>
#include <N_NLS_Manager.h>

#include <N_LOA_Loader.h>

#include <N_UTL_OptionBlock.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_LAS_Solver.h>
#include <N_LAS_Problem.h>

#include <N_LAS_System.h>
#include <N_LAS_Builder.h>

#include <N_ERH_ErrorMgr.h>
#include <N_ANP_AnalysisInterface.h>

#include <N_TOP_Topology.h>

#include <N_UTL_Expression.h>

#include <N_PDS_ParMap.h>
#ifdef Xyce_PARALLEL_MPI
#include <N_PDS_ParComm.h>
#include <mpi.h>
#else
#include <N_PDS_SerialComm.h>
#endif


// ----------   Static Declarations ----------

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::N_NLS_Sensitivity
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/02
//-----------------------------------------------------------------------------
N_NLS_Sensitivity::N_NLS_Sensitivity (N_NLS_NonLinearSolver & nls,
                                      N_TOP_Topology & topTmp,
                                      N_IO_CmdParse & cp)
    : N_NLS_NonLinearSolver(cp),
      allocateddXVec_ (false),
      debugLevel_(1),
      solutionSize_(0),
      solveDirectFlag_(true),
      solveAdjointFlag_(true),
      outputScaledFlag_(false),
      outputUnscaledFlag_(true),
      maxParamStringSize_(0),
      stdOutputFlag_(true),
      fileOutputFlag_(false),
      dakotaFileOutputFlag_(false),
      numSolves_(0),
      difference(SENS_FWD),
      objFuncGiven_(false),
      objFuncGIDsetup_(false),
      expNumVars_(0),
      expVal_(0.0),
      objFuncString_(""),
      curValue_(0.0),
      objFuncEval_(0.0),
      dOdp_(0.0),
      sqrtEta_(1.0e-8),
      sqrtEtaGiven_(false),
      dOdXVectorPtr_(0),
      lambdaVectorPtr_(0),
      savedRHSVectorPtr_(0),
      savedNewtonVectorPtr_(0),
      origFVectorPtr_(0),
      pertFVectorPtr_(0),
      testFVectorPtr_(0),
      nls_(nls),
      top_(topTmp),
      expPtr_(0),
      numSensParams_(0)
{
  // if the base nonlinear solver class had a copy constructor, I could use
  // that here... maybe I'll set that up later. ERK 11/15/02.
  lasSysPtr_    = nls_.lasSysPtr_;
  anaIntPtr_    = nls_.anaIntPtr_;
  loaderPtr_    = nls_.loaderPtr_;
  rhsVectorPtr_ = nls_.rhsVectorPtr_;

  NewtonVectorPtr_     = nls_.NewtonVectorPtr_;
  lasSolverPtr_        = nls_.lasSolverPtr_;
  jacobianMatrixPtr_   = nls_.jacobianMatrixPtr_;
  nextSolVectorPtrPtr_ = nls_.nextSolVectorPtrPtr_;

  dOdXVectorPtr_        = lasSysPtr_->builder().createVector();
  savedRHSVectorPtr_    = lasSysPtr_->builder().createVector();
  savedNewtonVectorPtr_ = lasSysPtr_->builder().createVector();

  origFVectorPtr_ = lasSysPtr_->builder().createVector();
  pertFVectorPtr_ = lasSysPtr_->builder().createVector();
  testFVectorPtr_ = lasSysPtr_->builder().createVector();

  solutionSize_ = lasSysPtr_->getSolutionSize();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::~N_NLS_Sensitivity
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/02
//-----------------------------------------------------------------------------
N_NLS_Sensitivity::~N_NLS_Sensitivity()
{
  delete dOdXVectorPtr_;
  dOdXVectorPtr_ = 0;

  delete savedRHSVectorPtr_;
  savedRHSVectorPtr_ = 0;

  delete savedNewtonVectorPtr_;
  savedNewtonVectorPtr_ = 0;

  delete origFVectorPtr_;
  origFVectorPtr_ = 0;

  delete pertFVectorPtr_;
  pertFVectorPtr_ = 0;

  delete testFVectorPtr_;
  testFVectorPtr_ = 0;

  if (lambdaVectorPtr_) // only allocated for adjoint solves
  {
    delete lambdaVectorPtr_;
    lambdaVectorPtr_ = 0;
  }

  if (expPtr_)
  {
    delete expPtr_;
    expPtr_ = 0;
  }

  // For all the stuff that is to be deleted in the nonlinear solver
  // base class destructor, just set those pointers to zero because
  // they will have been deleted already.
  //
  // This is the consequence of having the sensitivity class derived
  // off of the nonlinear solver base class, but having it use all
  // the same linear solver objects, etc., as the nonlinear solver used
  // to solve the problem.

  NewtonVectorPtr_     = 0;
  gradVectorPtr_       = 0;
  solWtVectorPtr_      = 0;
  petraOptionBlockPtr_ = 0;
  lasSolverPtr_        = 0;

  // delete contents of dfdp vector:
  for (int iparam=0;iparam<numSensParams_; ++iparam)
  {
    if( dfdpPtrVector_[iparam] != 0)
    {
      delete dfdpPtrVector_[iparam];
      dfdpPtrVector_[iparam] = 0;
    }
    if (dXdpPtrVector_[iparam] != 0)
    {
      delete dXdpPtrVector_[iparam];
      dXdpPtrVector_[iparam] = 0;
    }
  }
  dfdpPtrVector_.clear();
  dXdpPtrVector_.clear();
}



//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::stdOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/18/2014
//-----------------------------------------------------------------------------
void N_NLS_Sensitivity::stdOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities)
{
  // Send the sensitivity information to the screen:
#ifdef Xyce_PARALLEL_MPI
  N_PDS_Comm *pdsCommPtr = pdsMgrPtr_->getPDSComm();
  int myPID = pdsCommPtr->procID();
  if (myPID==0)
#endif
  {
    Xyce::dout() << "\n"<<idString << " Sensitivities of objective function:" 
        << objFuncString_ << std::endl;

    Xyce::dout() << std::setw(maxParamStringSize_)<<"Name";
    Xyce::dout() << "\t"<<std::setw(13)<<"Value";
    Xyce::dout() << "\t"<<std::setw(13)<<"Sensitivity";
    Xyce::dout() << "\t"<<std::setw(13)<<"Normalized"<<std::endl;

    for (int iparam=0; iparam< numSensParams_; ++iparam)
    {
      Xyce::dout() << std::setw(maxParamStringSize_)<<paramNameVec_[iparam];

      Xyce::dout() << "\t" << std::setw(13)<< std::scientific<< std::setprecision(4)
        << paramVals[iparam];

      Xyce::dout() << "\t" << std::setw(13)<< std::scientific<< std::setprecision(4)
        << sensitivities[iparam];

      Xyce::dout() << "\t" << std::setw(13)<< std::scientific<< std::setprecision(4)
        << scaled_sensitivities[iparam] << std::endl;
    }
  }
#ifdef Xyce_PARALLEL_MPI
  pdsCommPtr->barrier();
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::fileOutput
// Purpose       : Dump sensitivity information out to a file.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/18/2014
//-----------------------------------------------------------------------------
void N_NLS_Sensitivity::fileOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities)
{

#ifdef Xyce_PARALLEL_MPI
  N_PDS_Comm *pdsCommPtr = pdsMgrPtr_->getPDSComm();
  int myPID = pdsCommPtr->procID();
  if (myPID==0)
#endif
  {
    std::ostringstream numSolvesOStr;
    numSolvesOStr << numSolves_;
    std::string dodpFileName = netlistFileName_ + numSolvesOStr.str() + "_dodp" + idString +".txt";
    FILE *fp = fopen(dodpFileName.c_str(),"w");
    for (int iparam=0;iparam< numSensParams_; ++iparam)
    {
      fprintf(fp,"\t%16.8e\n", sensitivities[iparam]);
    }
    fclose(fp);
  }
#ifdef Xyce_PARALLEL_MPI
  pdsCommPtr->barrier();
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::dakOutput
// Purpose       : Dump sensitivity information out to a dakota-style file.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/19/2014
//-----------------------------------------------------------------------------
void N_NLS_Sensitivity::dakOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities)
{
#ifdef Xyce_PARALLEL_MPI
  N_PDS_Comm *pdsCommPtr = pdsMgrPtr_->getPDSComm();
  int myPID = pdsCommPtr->procID();
  if (myPID==0)
#endif
  {
    // write a file format that can be used by dakota :
    // Note that currently, this will simply overwrite the same 
    // file every time this function is called.
    std::string dakotaFileName = netlistFileName_ + "_dodp" + idString + "_all.txt";
    FILE *fp2 = fopen(dakotaFileName.c_str(),"w");
    fprintf(fp2,"%16.8e", objFuncEval_ );
    fprintf(fp2,"%s","\n[\n");
    for (int iparam=0;iparam< numSensParams_; ++iparam)
    {
      fprintf(fp2,"\t%16.8e\n", sensitivities[iparam]);
    }
    fprintf(fp2,"%s","]\n");
    fclose(fp2);
  }
#ifdef Xyce_PARALLEL_MPI
  pdsCommPtr->barrier();
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::solve
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/21/02
//-----------------------------------------------------------------------------
int N_NLS_Sensitivity::solve ( 
     std::vector<double> & objectiveVec,
     std::vector<double> & dOdpVec, 
     std::vector<double> & dOdpAdjVec,
     std::vector<double> & scaled_dOdpVec, 
     std::vector<double> & scaled_dOdpAdjVec)
{
  if (!solveDirectFlag_ && !solveAdjointFlag_) return 1;

  dOdpVec_.clear();
  dOdpAdjVec_.clear();

  scaled_dOdpVec_.clear();
  scaled_dOdpAdjVec_.clear();

  // It may now be neccessary to re-load the jacobian and rhs vectors.
  // It is necccessary, for example, if the two-level Newton solver is the
  // solver being used.
  nls_.enableSensitivity ();

  // first get the derivatives of the RHS vector w.r.t. the
  // user-specified optimization parameters.
  calcSensitivities ();

  calcObjFuncDerivs ();

  objectiveVec.clear();
  objectiveVec.push_back(objFuncEval_);

  if (solveDirectFlag_) 
  {
    solveDirect ();

    if (outputUnscaledFlag_)
    {
      dOdpVec = dOdpVec_;
    }

    if (outputScaledFlag_)
    {
      scaled_dOdpVec = scaled_dOdpVec_;
    }
  }

  if (solveAdjointFlag_) 
  {
    solveAdjoint ();
    if (outputUnscaledFlag_)
    {
      dOdpAdjVec = dOdpAdjVec_;
    }

    if (outputScaledFlag_)
    {
      scaled_dOdpAdjVec = scaled_dOdpAdjVec_;
    }
  }

  numSolves_++;

  return 1;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::solveDirect
//
// Purpose       : This function calculates the direct sensitivities for
//                 the user specified parameters.
//
//                 The ultimate goal of this function is to obtain dO/dp,
//                 where O is the objective function, and p is a
//                 user-defined optimization parameter.
//
//                 This is a bit confusing because the DAKOTA folks use
//                 a different naming convention than I have tended to
//                 use.  In Dakota's documentation, dO/dp would be referred
//                 to as df/du.  In Xyce, f is already considered to be the
//                 right-hand-side residual vector.  For clarity, this is
//                 the key between my notation and DAKOTA's:
//                        DAK     Xyce
//                         f       O
//                         y       x
//                         u       p
//                         c       f
//
//                 To obtain dOdp, the device manager is first called, and
//                 told to calculate dfdp, which is the derivative of the
//                 residual vector w.r.t. user-defined params.  It does
//                 this by numerical differentiation.  Then, after that,
//                 this function solves the equation:
//
//                 J dx/dp = df/dp    ->   dx/dp = J^-1 df/dp
//
//                 for each p.  
//
//                 J is the jacobian matrix, so J=df/dx.  (dc/dy)
//
//                 After performing these linear solves, dO/dp is to be
//                 obtained using the chain rule by:
//
//                 dO/dp = - dO/dx * dx/dp  + dO/dp
//
//                 The O in the dO/dp on the left side of this equation
//                 should have a hat over it "^", to indicate that it is
//                 different than the O on the right hand side.
//
//                 Note, this method is best if you have lots of objective
//                 functions, and a small number of parameters.  For adjoint
//                 calculations it is the other way around.
//
//                 11/19/02.
//                 It is assumed (for now) that dO/dp on the right hand side
//                 is zero, i.e., there is no direct
//                 dependence on p by O.  So, for example, if a user
//                 defined parameter to be used is the length of a MOSFET,
//                 the MOSFET length will NOT appear in the analytical
//                 expression for O.
//
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/02
//-----------------------------------------------------------------------------
int N_NLS_Sensitivity::solveDirect ()
{
#ifdef Xyce_DEBUG_NONLINEAR
  if (debugLevel_ > 0)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << "In N_NLS_Sensitivity::solveDirect" << std::endl;
  }
#endif

  int iparam;

  if( !allocateddXVec_ )
  {
    dXdpPtrVector_.resize(numSensParams_,0);
    for (iparam=0;iparam<numSensParams_;++iparam)
      dXdpPtrVector_[iparam] = lasSysPtr_->builder().createVector();
    allocateddXVec_ = true;
  }

  // first save a copy of the rhs vector, in case we want it later.
  savedRHSVectorPtr_->update(1.0, *(rhsVectorPtr_), 0.0);
  savedNewtonVectorPtr_->update(1.0, *(NewtonVectorPtr_), 0.0);

  // Now solve the series of linear systems to get dXdp.
  for (iparam=0; iparam< numSensParams_; ++iparam)
  {
    // copy the current dfdp vector into the f vector data structure.
    rhsVectorPtr_->update(1.0, *(dfdpPtrVector_[iparam]), 0.0);

    lasSolverPtr_->solve();

    // allocate the dxdp vector for this param, and
    // copy the resulting deltax vector into the dxdp data structure.
    (dXdpPtrVector_[iparam])->update(1.0, *(NewtonVectorPtr_), 0.0);

#ifdef Xyce_DEBUG_NONLINEAR
    // do debug output.
    if (debugLevel_ > 0)
    {
      Xyce::dout() << "iparam="<<iparam << "\t" << paramNameVec_[iparam] <<std::endl;
      for (int k = 0; k < solutionSize_; ++k)
      {
        Xyce::dout() << "k = " << std::setw(4) << k
          <<" dXdp = "<< std::setw(11)<< std::scientific
          << std::setprecision(4)<< (*(dXdpPtrVector_[iparam]))[k]
          <<std::endl;
      }

      std::ostringstream filename; 
      filename << netlistFileName_ << "_dxdp";
      filename << std::setw(3) << std::setfill('0') << iparam;
      filename << ".txt";
      dXdpPtrVector_[iparam]->writeToFile(const_cast<char *>(filename.str().c_str()));
    }
#endif

  }// end of for loop

  // Restore the RHS and Newton vectors.
  rhsVectorPtr_->update(1.0, *(savedRHSVectorPtr_),0.0);
  NewtonVectorPtr_->update(1.0, *(savedNewtonVectorPtr_),0.0);

  // Now get the final dOdp's (one for each param).
  for (iparam=0; iparam< numSensParams_; ++iparam)
  {
    double tmp = (-1.0) * dOdXVectorPtr_->dotProduct( (*(dXdpPtrVector_[iparam]))  );
    tmp += dOdp_;

    dOdpVec_.push_back(tmp);
 
    // get scaled value.  dO/dp*(p/100)
    double normalize = paramOrigVals_[iparam]/100.0;
    tmp *= normalize;
    scaled_dOdpVec_.push_back(tmp);
  }

  if (stdOutputFlag_)
  {
    stdOutput(std::string("Direct"), paramOrigVals_, dOdpVec_, scaled_dOdpVec_);
  }

  if (fileOutputFlag_)
  {
    fileOutput(std::string("Direct"), paramOrigVals_, dOdpVec_, scaled_dOdpVec_);
  }

  if (dakotaFileOutputFlag_)
  {
    dakOutput(std::string("Direct"), paramOrigVals_, dOdpVec_, scaled_dOdpVec_);
  }

  return 1;
}



//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::calcObjFuncDerivs
//
// Purpose       : This function assumes an objective function, O, and
//                 uses it to calculate other stuff, including dO'/du,
//                 dO/dy, dO/du.
//
//                 The objective function simply assumes that you are
//                 trying to match a solution variable with a
//                 user-specified number.   (objective=1)
//
//                 Or that you are trying to match your final solution
//                 with a specified solution. (objective=2)
//
// Special Notes : This is mostly for testing purposes.  Usually, the
//                 objective function would be set and managed from
//                 the DAKOTA side.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/15/02
//-----------------------------------------------------------------------------
bool N_NLS_Sensitivity::calcObjFuncDerivs ()
{
  bool bsuccess = true;
  int i;

  dOdXVectorPtr_->putScalar(0.0);

  // first check if the user specified a function via the expressions
  // class:
  if (objFuncGiven_)
  {
    if (!objFuncGIDsetup_)
    {
      // set up the gid's:
      expVarGIDs_.resize( expNumVars_ );
      expVarLocal_.resize( expNumVars_ );

      for (i = 0; i < expNumVars_; ++i)
      {
        std::list<int> svGIDList1, dummyList;
        char type1;
        bool found = top_.getNodeSVarGIDs(NodeID(expVarNames_[i], Xyce::_VNODE), svGIDList1, dummyList, type1);

        bool foundLocal = found;
#ifdef Xyce_PARALLEL_MPI
        N_PDS_Comm *pdsCommPtr = pdsMgrPtr_->getPDSComm();
        {
          // synchronize in parallel
        double found_glob=0.0;
        double found_local = found?1.0:0.0;
        pdsCommPtr->barrier();
        pdsCommPtr->sumAll(&found_local, &found_glob, 1);
        found = (found_glob != 0.0)?true:false;
#ifdef Xyce_DEBUG_NONLINEAR
        Xyce::dout() << "global found = " << found_glob <<std::endl;
#endif
        }
#endif

        bool found2 = false;
        if (!found)// if looking for this as a voltage node failed, try a "device" (i.e. current) node.
        {
          found2 = top_.getNodeSVarGIDs(NodeID(expVarNames_[i], Xyce::_DNODE), svGIDList1, dummyList, type1);
        }

        bool foundLocal2 = found2;
#ifdef Xyce_PARALLEL_MPI
        if (!found)
        {
          // synchronize in parallel
          double found_glob=0.0;
          double found_local = found2?1.0:0.0;
          pdsCommPtr->barrier();
          pdsCommPtr->sumAll(&found_local, &found_glob, 1);
          found2 = (found_glob != 0.0)?true:false;
#ifdef Xyce_DEBUG_NONLINEAR
          Xyce::dout() << "global found2 = " << found_glob <<std::endl;
#endif
        }
#endif

        if (!found && !found2)
        {
          static std::string tmp = "objective function variable not found!\n";
          N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL, tmp);
        }

        int tmpGID = -1;
        if (foundLocal || foundLocal2)
        {
          tmpGID = svGIDList1.front();
        }

#ifdef Xyce_DEBUG_NONLINEAR
        Xyce::dout() << "tmpGID = " << tmpGID <<std::endl;
#endif
        expVarGIDs_[i] = tmpGID;
      }
      objFuncGIDsetup_ = true;
    }

    // obtain the expression variable values:
    expVarVals_.resize (expNumVars_);
    expVarDerivs_.resize (expNumVars_);
    for (i = 0; i < expNumVars_; ++i)
    {
      if (expVarGIDs_[i] != -1)
      {
        expVarVals_[i] =
         (*nextSolVectorPtrPtr_)->getElementByGlobalIndex(expVarGIDs_[i], 0);
      }
      else
      {
        expVarVals_[i] = 0.0;
      }
    }

    //get expression value and partial derivatives
    expPtr_->evaluate( expVal_, expVarDerivs_, expVarVals_ );
    objFuncEval_ = expVal_;
    dOdXVectorPtr_->putScalar(0.0);
    for (i=0;i<expNumVars_;++i)
    {
      int tmpGID = expVarGIDs_[i];
      double tmpDODX = expVarDerivs_[i];

#ifdef Xyce_DEBUG_NONLINEAR
      Xyce::dout() << "i="<<i<<"  gid = " << tmpGID << "  dodx = "<< tmpDODX << std::endl;
#endif

      if (tmpGID != -1)
      {
        dOdXVectorPtr_->setElementByGlobalIndex(tmpGID, tmpDODX, 0);
      }
    }

    // Assuming this is zero:
    dOdp_ = 0.0;
  }
  else // give a warning.
  {
    static std::string tmp = " ***** WARNING: objective function was not specified.\n\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, tmp);
  }

#ifdef Xyce_DEBUG_NONLINEAR
  if (debugLevel_ > 0)
  {
    std::string filename = netlistFileName_ + "_dodx.txt";
    dOdXVectorPtr_->writeToFile(const_cast<char *>(filename.c_str()));
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::solveAdjoint
// Purpose       : Solves first for the vector, lambda, of adjoint variables.
//                 Afterwards, it solves for dO/dp.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/20/02
//-----------------------------------------------------------------------------
int N_NLS_Sensitivity::solveAdjoint ()
{
  if( !allocateddXVec_ )
  {
    dXdpPtrVector_.resize(numSensParams_,0);
    for (int iparam=0;iparam<numSensParams_;++iparam)
      dXdpPtrVector_[iparam] = lasSysPtr_->builder().createVector();
    allocateddXVec_ = true;
  }

  // first save a copy of the rhs vector, in case we want it later.
  savedRHSVectorPtr_->update(1.0, *(rhsVectorPtr_),0.0);
  savedNewtonVectorPtr_->update(1.0, *(NewtonVectorPtr_),0.0);

  lambdaVectorPtr_ = lasSysPtr_->builder().createVector();
  bool useTran = jacobianMatrixPtr_->useTranspose ();
  jacobianMatrixPtr_->setUseTranspose (true);

#ifdef Xyce_DEBUG_NONLINEAR
  if (debugLevel_ > 0)
  {
    std::string matrixFile = netlistFileName_ + "_adjointMatrix.txt";
    jacobianMatrixPtr_->writeToFile(const_cast<char *>(matrixFile.c_str()));
  }
#endif

  // copy the current dOdx vector into the f vector data structure.
  rhsVectorPtr_->update(1.0, *(dOdXVectorPtr_),0.0);

  int status = lasSolverPtr_->solve();
  if (status!=0)
  {
    std::string msg("N_NLS_Sensitivity::solveAdjoint.  Solver failed\n");
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // allocate the dxdp vector for this param, and
  // copy the resulting deltax vector into the dxdp data structure.
  lambdaVectorPtr_->update(1.0, *(NewtonVectorPtr_),0.0);

#ifdef Xyce_DEBUG_NONLINEAR
  if (debugLevel_ > 0)
  {
    std::string filename = netlistFileName_ + "_lambda.txt";
    lambdaVectorPtr_->writeToFile(const_cast<char *>(filename.c_str()));
  }
#endif

  // Now that we have lambda, get the dOdp's by doing dot products of
  // lambda * df/dp.

  // do the final dot products, one for each param.
  for (int iparam=0; iparam< numSensParams_; ++iparam)
  {
    double tmp = -1.0 * lambdaVectorPtr_->dotProduct(*(dfdpPtrVector_[iparam]));
    dOdpAdjVec_.push_back(tmp);

    // get scaled value.  dO/dp*(p/100)
    double normalize = paramOrigVals_[iparam]/100.0;
    tmp *= normalize;
    scaled_dOdpAdjVec_.push_back(tmp);
  }


  // restore the useTranspose flag to the original setting. (probably
  // false)
  jacobianMatrixPtr_->setUseTranspose (useTran);

  // Restore the RHS and Newton vectors.
  rhsVectorPtr_->update(1.0, *(savedRHSVectorPtr_),0.0);
  NewtonVectorPtr_->update(1.0, *(savedNewtonVectorPtr_),0.0);

  // set the sensitivity information to the screen:
  if (stdOutputFlag_)
  {
    stdOutput(std::string("Adjoint"), paramOrigVals_, dOdpAdjVec_, scaled_dOdpAdjVec_);
  }
  if (fileOutputFlag_)
  {
    fileOutput(std::string("Adjoint"), paramOrigVals_, dOdpVec_, scaled_dOdpVec_);
  }

  if (dakotaFileOutputFlag_)
  {
    dakOutput(std::string("Adjoint"), paramOrigVals_, dOdpVec_, scaled_dOdpVec_);
  }

  return 1;
}


//-----------------------------------------------------------------------------
// Function      : Sensitivity::calcSensitivities
//
// Purpose       : This function is to be called after a calculation has
//                 converged.  It calculates vectors of derivatives:  df/dp,
//                 where f is the residual, and p is a varried parameter.
//
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------
bool N_NLS_Sensitivity::calcSensitivities ()
{
  int iparam;
  std::string msg;
  N_LAS_System & lasSys_ = (*lasSysPtr_);

#ifdef Xyce_PARALLEL_MPI
  N_PDS_Comm *pdsCommPtr = pdsMgrPtr_->getPDSComm();
  pdsCommPtr->barrier();
#endif

  paramOrigVals_.clear();
  
  // save a copy of the f vector.
  origFVectorPtr_->update(1.0, *(rhsVectorPtr_), 0.0);

  // Loop over the vector of parameters.  For each parameter, find the
  // device entity (a model or an instance) which corresponds to it, and
  // perform the finite difference calculation.
  std::vector<std::string>::iterator firstParam = paramNameVec_.begin ();
  std::vector<std::string>::iterator lastParam  = paramNameVec_.end ();
  std::vector<std::string>::iterator iterParam;
  for ( iterParam=firstParam, iparam=0;
        iterParam!=lastParam; ++iterParam, ++iparam )
  {

#ifdef Xyce_DEBUG_NONLINEAR
    if (debugLevel_ > 0)
    {
      Xyce::dout() << std::endl << "  Calculating df/dp for: ";
      Xyce::dout() << *iterParam << std::endl;
    }
#endif

    // save a copy of the rhs vector, in case we want it later.
    // (-RHS, as it supports a Newton solve)
    origFVectorPtr_->update(-1.0, *(rhsVectorPtr_), 0.0);

    // now perturb the value of this parameter.
    std::string paramName(*iterParam);
    double paramOrig = 0.0;
    bool found = loaderPtr_->getParamAndReduce(paramName, paramOrig);
    if (!found)
    {
      std::string msg("Sensitivity::calcSensitivities: cannot find parameter ");
      msg += paramName;
      N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::DEV_FATAL, msg);
    }

    double dp = sqrtEta_ * (1.0 + fabs(paramOrig));
    double paramPerturbed = paramOrig;
    paramOrigVals_.push_back(paramOrig);

    if (difference==SENS_FWD)
    {
      paramPerturbed += dp;
    }
    else if (difference==SENS_REV)
    {
      paramPerturbed -= dp;
    }
    else if (difference==SENS_CNT)
    {
      static std::string tmp = "difference=central not supported.\n";
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0, tmp);
    }
    else
    {
      static std::string tmp = "difference not recognized!\n";
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0, tmp);
    }

#ifdef Xyce_DEBUG_NONLINEAR
    if (debugLevel_ > 0)
    {
      Xyce::dout() << std::setw(maxParamStringSize_)<< *iterParam
        << " dp = " << std::setw(11)<< std::scientific<< std::setprecision(4) << dp 
        << " original value = " << std::setw(16)<< std::scientific<< std::setprecision(9) << paramOrig 
        << " modified value = " << std::setw(16)<< std::scientific<< std::setprecision(9) << paramPerturbed 
        <<std::endl;
    }
#endif
    loaderPtr_->setParam (paramName, paramPerturbed);

    // Now that the parameter has been perturbed,
    // calculate the numerical derivative.

    // rhs_(); // this is the same as loadPtr_->loadRHS, except that it keeps statistics
    loaderPtr_->loadRHS();

    // save the pertF vector 
    // (-RHS, as it supports a Newton solve)
    pertFVectorPtr_->update(-1.0, *(rhsVectorPtr_), 0.0);

    // calculate the df/dp vector.  
    double rdp=1/dp;
    dfdpPtrVector_[iparam]->putScalar(0.0);
    dfdpPtrVector_[iparam]->addVec (+1.0, *(pertFVectorPtr_)); //+Fperturb
    dfdpPtrVector_[iparam]->addVec (-1.0, *(origFVectorPtr_)); //-Forig
    dfdpPtrVector_[iparam]->scale(rdp);

#ifdef Xyce_DEBUG_NONLINEAR
    if (debugLevel_ > 0)
    {
      Xyce::dout() << *iterParam << ": ";
      Xyce::dout().width(15); Xyce::dout().precision(7); Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << "deviceSens_dp = " << dp << std::endl;

      for (int k1 = 0; k1 < solutionSize_; ++k1)
      {
        Xyce::dout() << "k = " << std::setw(4) << k1
          <<" fpert = "<< std::setw(11)<< std::scientific<< std::setprecision(4)<< (*(pertFVectorPtr_))[k1]
          <<" forig = "<< std::setw(11)<< std::scientific<< std::setprecision(4)<< (*(origFVectorPtr_))[k1]
          <<" dfdp = "<< std::setw(11)<< std::scientific<< std::setprecision(4)<< (*(dfdpPtrVector_[iparam]))[k1]
          <<std::endl;
      }

      std::ostringstream filename; 
      filename << netlistFileName_ << "_dfdp";
      filename << std::setw(3) << std::setfill('0') << iparam;
      filename << ".txt";
      dfdpPtrVector_[iparam]->writeToFile(const_cast<char *>(filename.str().c_str()));

      filename.str("");
      filename << netlistFileName_ << "_fpert";
      filename << std::setw(3) << std::setfill('0') << iparam;
      filename << ".txt";
      pertFVectorPtr_->writeToFile(const_cast<char *>(filename.str().c_str()));
    }
#endif

    // now reset the parameter and rhs to previous values.
    loaderPtr_->setParam (paramName, paramOrig);
    rhsVectorPtr_->update(-1.0, *(origFVectorPtr_), 0.0);
  }

#ifdef Xyce_PARALLEL_MPI
  pdsCommPtr->barrier();
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::setOptions
//
// Purpose       : This function processes the .SENS line
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/15/02
//-----------------------------------------------------------------------------
bool N_NLS_Sensitivity::setOptions(const N_UTL_OptionBlock& OB)
{
  bool bsuccess = true;
  std::list<N_UTL_Param>::const_iterator iter = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator end   = OB.getParams().end();

  numSensParams_ = 0;
  for ( ; iter != end; ++ iter)
  {
    if (iter->uTag() == "OBJFUNC")
    {
      objFuncString_ = iter->stringValue();
      expPtr_ = new N_UTL_Expression(iter->stringValue());
      objFuncGiven_ = true;
    }
    else if ( std::string( iter->uTag() ,0,5) == "PARAM") // this is a vector
    {
      ExtendedString tag = iter->stringValue();
      tag.toUpper();
      // set up the initial skeleton of the maps:
      ++numSensParams_;
      paramNameVec_.push_back(tag);
      int sz = tag.size();
      if (sz > maxParamStringSize_)
      {
        maxParamStringSize_ = sz;
      }
    }
    else
    {
      Xyce::Report::UserWarning() << iter->uTag() << " is not a recognized sensitivity solver option.\n" << std::endl;
    }
  }

  // parse the expression now, so if there are any errors, they will come
  // up early in the simulation.
  if (objFuncGiven_)
  {
    // setup the names:
    expVarNames_.clear();

    std::vector<std::string> nodes;
    expPtr_->get_names(XEXP_NODE, nodes);
    std::vector<std::string> instances;
    expPtr_->get_names(XEXP_INSTANCE, instances);

    // Make the current (instance) strings all upper case.
    // The topology directory apparently requires this.
    std::vector<std::string>::iterator iter;
    for (iter=instances.begin();iter!=instances.end();++iter)
    {
      ExtendedString tmpString = *iter;
      tmpString.toUpper ();
      *iter  = tmpString;
    }

    expVarNames_.insert(expVarNames_.end(), nodes.begin(), nodes.end());
    expVarNames_.insert(expVarNames_.end(), instances.begin(), instances.end());

    // Order the names in the expression so that it agrees with the order
    // in expVarNames.
    if ( !expVarNames_.empty() )
    {
      expPtr_->order_names( expVarNames_ );
    }

    expNumVars_ = expVarNames_.size();
  }

#ifdef Xyce_DEBUG_NONLINEAR
  if (debugLevel_ > 0)
  {
    std::vector<std::string>::iterator iter;
    std::vector<std::string>::iterator begin = paramNameVec_.begin();
    std::vector<std::string>::iterator end   = paramNameVec_.end  ();

    for (iter=begin;iter!=end;++iter)
    {
      Xyce::dout() << *iter<<std::endl;
    }
  }
#endif

  // now that the number of parameters is known, allocate the dfdp vector:
  dfdpPtrVector_.resize (numSensParams_, 0);
  N_LAS_System & lasSys_ = (*lasSysPtr_);

  for (int iparam=0;iparam<numSensParams_; ++iparam)
  {
    dfdpPtrVector_[iparam] = lasSys_.builder().createVector();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::setSensitivityOptions
// Purpose       : This function processes the .options SENSITIVITY line
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/20/2013
//-----------------------------------------------------------------------------
bool N_NLS_Sensitivity::setSensitivityOptions(const N_UTL_OptionBlock& OB)
{
  bool bsuccess = true;
  std::list<N_UTL_Param>::const_iterator iter = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator end   = OB.getParams().end();

  for ( ; iter != end; ++ iter)
  {
    if (iter->uTag() == "ADJOINT")
    {
      solveAdjointFlag_ = 
        static_cast<bool>(iter->getImmutableValue<bool>());
    }
    else if (iter->uTag() == "DIRECT")
    {
      solveDirectFlag_ = 
        static_cast<bool>(iter->getImmutableValue<bool>());
    }
    else if (iter->uTag() == "OUTPUTSCALED")
    {
      outputScaledFlag_ = 
        static_cast<bool>(iter->getImmutableValue<bool>());
    }
    else if (iter->uTag() == "OUTPUTUNSCALED")
    {
      outputUnscaledFlag_ = 
        static_cast<bool>(iter->getImmutableValue<bool>());
    }
    else if (iter->uTag() == "STDOUTPUT")
    {
      stdOutputFlag_ = 
        static_cast<bool>(iter->getImmutableValue<bool>());
    }
    else if (iter->uTag() == "DAKOTAFILE")
    {
      dakotaFileOutputFlag_ = 
        static_cast<bool>(iter->getImmutableValue<bool>());
    }
    else if (iter->uTag() == "DIAGNOSTICFILE")
    {
      fileOutputFlag_ = 
        static_cast<bool>(iter->getImmutableValue<bool>());
    }
    else if (iter->uTag() == "DIFFERENCE")
    {
      ExtendedString sval=iter->stringValue();
      sval.toUpper();
      if(sval=="FORWARD")
      {
        difference=SENS_FWD;
      }
      else if(sval=="REVERSE")
      {
        difference=SENS_REV;
      }
      else if(sval=="CENTRAL")
      {
        difference=SENS_CNT;
        static std::string tmp = "difference=central not supported.\n";
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0, tmp);
      }
      else
      {
        static std::string tmp = "difference not recognized!\n";
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0, tmp);
      }
    }
    else if (iter->uTag() == "SQRTETA")
    {
      sqrtEta_ = iter->getImmutableValue<double>();
      sqrtEtaGiven_ = true;
    }
#ifdef Xyce_DEBUG_NONLINEAR
    else if (iter->uTag() == "DEBUGLEVEL")
    {
      debugLevel_ = iter->getImmutableValue<int>();
    }
#endif
    else
    {
      Xyce::Report::UserWarning() << iter->uTag() << " is not a recognized sensitivity solver option.\n" << std::endl;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::setTranOptions
//
// Purpose       : Not used yet.  Same as setOptions, but for transient
//                 mode, if and when the sensitivity stuff is ever set
//                 up to work in transient.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/02
//-----------------------------------------------------------------------------
bool N_NLS_Sensitivity::setTranOptions(const N_UTL_OptionBlock& OB)
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::setHBOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool N_NLS_Sensitivity::setHBOptions(const N_UTL_OptionBlock& OB)
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::getMaxNormF
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and MEMS Modeling
// Creation Date : 9/28/2009
//-----------------------------------------------------------------------------
double N_NLS_Sensitivity::getMaxNormF() const
{
  return nls_.getMaxNormF();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::getMaxNormFindex
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and MEMS Modeling
// Creation Date : 9/28/2009
//-----------------------------------------------------------------------------
int N_NLS_Sensitivity::getMaxNormFindex () const
{
  return nls_.getMaxNormFindex ();
}

