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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_ConductanceExtractor.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/03/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.29.2.2 $
//
// Revision Date  : $Date $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------   Standard Includes   ----------

#include <N_UTL_Misc.h>

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

// ----------   Xyce Includes   ----------

#include <N_NLS_ConductanceExtractor.h>
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
#include <N_PDS_Comm.h>

#include <N_IO_CmdParse.h>
#include <N_IO_PkgOptionsMgr.h>

#include <Epetra_CrsMatrix.h>

// ----------   Static Declarations ----------

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConductanceExtractor::N_NLS_ConductanceExtractor
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/03/06
//-----------------------------------------------------------------------------
N_NLS_ConductanceExtractor::N_NLS_ConductanceExtractor (
    N_NLS_NonLinearSolver & nls,
    N_TOP_Topology & topTmp,
    N_IO_CmdParse & cp)
    : commandLine_(cp),
      nls_(nls),
      top_(topTmp),
      solutionSize_(0),
      debugLevel_(0),
      lasSysPtr_(0),
      anaIntPtr_(0),
      loaderPtr_(0),
      rhsVectorPtr_(0),
      dfdvVectorPtr_(0),
      NewtonVectorPtr_(0),
      dxdvVectorPtr_(0),
      matrixDiagonalPtr_(0),
      lasSolverPtr_(0),
      jacobianMatrixPtr_(0),
      nextSolVectorPtrPtr_(0),
      currSolVectorPtrPtr_(0),
      savedRHSVectorPtr_(0),
      savedNewtonVectorPtr_(0),
      gradVectorPtr_(0),
      columnVectorPtr_(0),
      columnMapPtr_(0),
      gidsSetUpFlag_(false)
{
  lasSysPtr_    = nls_.lasSysPtr_;
  anaIntPtr_    = nls_.anaIntPtr_;
  loaderPtr_    = nls_.loaderPtr_;
  rhsVectorPtr_ = nls_.rhsVectorPtr_;
  dfdvVectorPtr_ = nls_.rhsVectorPtr_; // using the same vector for dfdv.

  NewtonVectorPtr_      = nls_.NewtonVectorPtr_;
  dxdvVectorPtr_        = nls_.NewtonVectorPtr_; // using same vector for dxdv.
  lasSolverPtr_         = nls_.lasSolverPtr_;
  jacobianMatrixPtr_    = nls_.jacobianMatrixPtr_;
  nextSolVectorPtrPtr_  = nls_.nextSolVectorPtrPtr_;
  currSolVectorPtrPtr_  = nls_.currSolVectorPtrPtr_;

  // creations
  savedRHSVectorPtr_    = lasSysPtr_->builder().createVector();
  savedNewtonVectorPtr_ = lasSysPtr_->builder().createVector();
  matrixDiagonalPtr_    = lasSysPtr_->builder().createVector();
  solutionSize_         = lasSysPtr_->getSolutionSize();

#ifdef Xyce_PARALLEL_MPI
  // construct column vector used for parallel construction of RHS's
  const Epetra_Map & col_map = jacobianMatrixPtr_->epetraObj().ColMap();
  columnMapPtr_ = new N_PDS_ParMap( &(const_cast<Epetra_Map&>(col_map)),
                                    savedRHSVectorPtr_->pmap()->pdsComm() );
  columnVectorPtr_ = new N_LAS_Vector( *(savedRHSVectorPtr_->pmap()), *columnMapPtr_ );
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConductanceExtractor::~N_NLS_ConductanceExtractor
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/03/06
//-----------------------------------------------------------------------------
N_NLS_ConductanceExtractor::~N_NLS_ConductanceExtractor()
{
  if (savedRHSVectorPtr_)
  {
    delete savedRHSVectorPtr_;
    savedRHSVectorPtr_ = 0;
  }

  if (savedNewtonVectorPtr_)
  {
    delete savedNewtonVectorPtr_;
    savedNewtonVectorPtr_ = 0;
  }

  if (matrixDiagonalPtr_)
  {
    delete matrixDiagonalPtr_;
    matrixDiagonalPtr_ = 0;
  }

  int vecSize = dIdxPtrVector_.size();
  for (int ivec=0; ivec<vecSize ;++ivec)
  {
    if (dIdxPtrVector_[ivec])
    {
      delete dIdxPtrVector_[ivec];
      dIdxPtrVector_[ivec] = 0;
    }
  }

  if( columnVectorPtr_ )
  {
    delete columnVectorPtr_;
    delete columnMapPtr_;
    columnVectorPtr_ = 0;
    columnMapPtr_ = 0;
  }

  // For all the stuff that is to be deleted in the nonlinear solver
  // base class destructor, just set those pointers to zero because
  // they are deleted elsewhere.
  lasSysPtr_            = 0;
  anaIntPtr_            = 0;
  loaderPtr_            = 0;
  rhsVectorPtr_         = 0;
  dfdvVectorPtr_        = 0;
  NewtonVectorPtr_      = 0;
  dxdvVectorPtr_        = 0;
  lasSolverPtr_         = 0;
  jacobianMatrixPtr_    = 0;
  nextSolVectorPtrPtr_  = 0;
  currSolVectorPtrPtr_  = 0;
  solutionSize_         = 0;
  gradVectorPtr_        = 0;
  lasSolverPtr_         = 0;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConductanceExtractor::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool N_NLS_ConductanceExtractor::registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;
  string netListFile = "";
  if (commandLine_.getArgumentValue("netlist") != "")
  {
    netListFile = commandLine_.getArgumentValue("netlist");
  }
  pkgOptMgrPtr_->submitRegistration(
      "CONDUCTANCE", netListFile, new N_NLS_ConductanceExtractor_OptionsReg( this ) );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConductanceExtractor::setOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/03/06
//-----------------------------------------------------------------------------
bool N_NLS_ConductanceExtractor::setOptions(const N_UTL_OptionBlock& OB)
{
  bool bsuccess = true;

  for (list<N_UTL_Param>::const_iterator it_tpL = OB.getParams().begin();
       it_tpL != OB.getParams().end(); ++ it_tpL)
  {
    // no-op for now.
    if (it_tpL->uTag() == "DEBUGLEVEL")
    {
      debugLevel_ = it_tpL->iVal();
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConductanceExtractor::setupIDs_
//
// Purpose       : This function sets up the various GIDs and matrix rows,
//                 that are needed in the extract function.
//
// Special Notes : There are 2 types of IDs needed, which are both associated
//                 with independent voltage sources.  One type is the GID of
//                 the current variable, which is needed for a "get" operation.
//                 The other is the LID of the matrix row for the voltage drop
//                 equation.  In serial, this is the same as the current GID,
//                 but in parallel they will generally be different.  For the
//                 matrix row, the LID is used for a "put" operation, so it
//                 makes more sense for it to be an LID.
//
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 03/03/06
//-----------------------------------------------------------------------------
bool N_NLS_ConductanceExtractor::setupIDs_
  (const map<string,double> & inputMap)
{
  bool bsuccess = true;

  // Setup vectors for each terminal current
  int idSize = inputMap.size();
  currentGIDs_.resize(idSize);
  currentLIDs_.resize(idSize);
  vsrcPosGIDs_.resize(idSize);
  vsrcPosLIDs_.resize(idSize);
  for (int i=0; i<idSize ;++i)
  {
    dIdxPtrVector_.push_back(lasSysPtr_->builder().createVector());
    currentGIDs_[i] = -1;
    currentLIDs_[i] = -1;
    vsrcPosGIDs_[i] = -1;
    vsrcPosLIDs_[i] = -1;
  }

  // Loop over the map, which should contain, in the first
  // arguments, the names of the voltage sources we need.
  map<string,double>::const_iterator iterM = inputMap.begin();
  map<string,double>::const_iterator  endM = inputMap.end  ();
  int i=0;
  for (; iterM != endM; ++i, ++iterM)
  {
    // Note: When passing the sourceName to topology, it must
    // be all CAPS, or topology won't find it.
    ExtendedString src = iterM->first;
    src.toUpper();
    string sourceName = src;
    char type;
    int index;

    // This is to get the IDs for the currents through the
    // voltage sources specified in the map.
    list<int> GIDList, extGIDList;
    top_.getNodeSVarGIDs(NodeID(sourceName,_DNODE), GIDList, extGIDList, type);

    list<int>::iterator iterI;
    if (!(GIDList.empty ()))
    {
      iterI = GIDList.begin();
      currentGIDs_[i] = *iterI;
      currentLIDs_[i] = dIdxPtrVector_[0]->pmap()->globalToLocalIndex(currentGIDs_[i]);

      iterI = extGIDList.begin();
      vsrcPosGIDs_[i] = *iterI;
      vsrcPosLIDs_[i] = dIdxPtrVector_[0]->pmap()->globalToLocalIndex(vsrcPosGIDs_[i]);

      if( vsrcPosLIDs_[i] == -1 )
      {
        string msg = "N_NLS_ConductanceExtractor::setupIDs_";
        msg += " The " + sourceName + " source has the positive node"
         " owned by another processor.  The 2-level solve can't handle that.";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }

      // check that vneg is connected to gnd
      ++iterI;
      int vnegGID = *iterI;
      if (vnegGID != -1)
      {
        string msg = "N_NLS_ConductanceExtractor::setupIDs_";
        msg += " The " + sourceName + " source has the negative node"
         " connected to something other than ground!  The 2-level solve can't handle that.";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }
    }

  }

#ifdef Xyce_DEBUG_CONDUCTANCE
  if (debugLevel_ > 0)
  {
    cout << "current GIDs: " << endl;
    for( int i1=0; i1 < idSize; ++i1 )
    {
      cout << "  currentGIDs_["<<i1<<"] = " << currentGIDs_[i1] << ", currentLIDs_["<<i1<<"] = " << currentLIDs_[i1] << endl;
    }

    cout << "Vsrc pos equation rows: " << endl;
    for( int i1=0; i1 < idSize; ++i1 )
    {
      cout << "  vsrcPosGIDs_["<<i1<<"] = " << vsrcPosGIDs_[i1] << ", vsrcPosLIDs_["<<i1<<"] = " << vsrcPosLIDs_[i1] << endl;
    }
  }
#endif

  gidsSetUpFlag_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConductanceExtractor::setup_dIdX_Vectors_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/05/06
//-----------------------------------------------------------------------------
bool N_NLS_ConductanceExtractor::setup_dIdX_Vectors_ ()
{
  bool bsuccess = true;

  // Set up dIdx's.  These correspond to rows in the Jacobian.
  // There will be one for each Vsrc current.   In general, we are assuming
  // that each "connecting" voltage source is connected to ground at the
  // negative node, and to the rest of the circuit at the positive node,
  // so the node of interest (and thus KCL of interest) is the
  // positive node.  dIdx is a vector containing all the derivatives
  // from the corresponding KCL row, minus the derivative with respect to
  // I.  (I is a solution variable, but needs to be excluded from dIdx,
  // or the various dot products will cancel some terms that they
  // should not).
  int idSize = currentGIDs_.size();

  N_LAS_Vector * currentVec;

  // int iC_row=0;
  for( int iC_row=0; iC_row < idSize; ++iC_row )
  {
#ifdef Xyce_PARALLEL_MPI
    currentVec = columnVectorPtr_;
#else
    currentVec = dIdxPtrVector_[iC_row];
#endif

    currentVec->putScalar(0.0);

    if( currentGIDs_[iC_row] != -1 )
    {
      int iRow = vsrcPosGIDs_[iC_row];
      int rowLength = jacobianMatrixPtr_->getRowLength(iRow);
      int numEntries = rowLength;
      vector<double> coeffs(rowLength, 0.0);
      vector<int> colIndices(rowLength, -1);

      jacobianMatrixPtr_->getRowCopy
        (iRow, rowLength, numEntries, &coeffs[0], &colIndices[0]);

      for (int ic=0;ic<rowLength;++ic)
      {
        // need to exclude entries that are with respect to
        // the the variable, 'I'.
        int gid = colIndices[ic];
        if (gid ==  currentGIDs_[iC_row]) coeffs[ic] = 0.0;
      }

      for (int icol=0;icol<rowLength;++icol)
      {
        double val = coeffs[icol];
        int gid = colIndices[icol];
        if (gid != -1)
          currentVec->setElementByGlobalIndex(gid, val, 0);
      }
    }
    currentVec->fillComplete();

#ifdef Xyce_PARALLEL_MPI
    *(dIdxPtrVector_[iC_row]) = *columnVectorPtr_;
#endif

#ifdef Xyce_DEBUG_CONDUCTANCE
    if (debugLevel_ > 0)
    {
      cout << "\ndIdx[" << iC_row << "]:" << endl;
      dIdxPtrVector_[iC_row]->printPetraObject();
    }
#endif
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConductanceExtractor::extract
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/03/06
//-----------------------------------------------------------------------------
bool N_NLS_ConductanceExtractor::extract (
        const map<string,double> & inputMap,
        vector<double> & outputVector,
        vector< vector<double> > & jacobian )

{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_CONDUCTANCE
  const string dashedline2 = "---------------------";
  if (debugLevel_ > 0)
  {
    cout << dashedline2 << endl;
    cout << "N_NLS_ConductanceExtractor::extract" << endl;
    cout << dashedline2 << endl;
  }
#endif

  if (inputMap.empty() ||
      outputVector.empty() ||
      jacobian.empty() )
  {
    return false;
  }

  // Need the "next" vector.
  N_LAS_Vector * solnVecPtr = *(nextSolVectorPtrPtr_);
  //N_LAS_Vector * solnVecPtr = *(currSolVectorPtrPtr_);

  if (!gidsSetUpFlag_)
  {
    setupIDs_(inputMap);
  }

  // Now that the solve is complete, obtain the currents and conductances
  // that are needed for the upper level nonlinear solve.

  // The currents are owned by the individual sources in the device package:
  // Note: this may need some parallel refactoring.

  // First obtain the currents.
  int idSize = currentGIDs_.size();
  for (int i=0;i<currentGIDs_.size();++i)
  {
    int index1 = currentGIDs_[i];
    if( index1 > -1 )
      outputVector[i] = solnVecPtr->getElementByGlobalIndex(index1);
    else
      outputVector[i] = 0.0;
  }

#ifdef Xyce_PARALLEL_MPI
  //sumAll to get all currents locally
  N_PDS_Comm * comm = dfdvVectorPtr_->pmap()->pdsComm();
  vector<double> tmpVector(idSize,0.0);
  for( int i = 0; i < idSize; ++i )
  {
    tmpVector[i] = outputVector[i];
    outputVector[i] = 0.0;
  }

  comm->sumAll( &(tmpVector[0]), &(outputVector[0]), idSize );
#endif

#ifdef Xyce_DEBUG_CONDUCTANCE
  if (debugLevel_ > 0)
  {
    int itmp;
    for (itmp=0;itmp < outputVector.size();++itmp)
    {
      cout << "currentVector["<< itmp <<"] = " << outputVector[itmp]<<endl;
    }
  }
#endif

  // This function needs to solve one or more linear systems.  The linear
  // systems, (for a direct sensitivity) is:
  // J * dx/dv = df/dv
  //
  // where "v" is an applied voltage, not part of the system of unknowns,
  // which is represented by "x".
  //
  // J is the traditional Jacobian matrix, df/dx.
  //
  // For the direct approach, a different dx/dv is needed for each v, so
  // a different linear system, with a different df/dv will be solved for
  // each applied voltage, v.
  //
  // For an adjoint approach, it could be possible to do a single linear
  // solve, and get the sensitivities to all the v's all at once.
  //
  // Note:  In general, this extraction is done in the context of a
  // multi-level Newton solve, in which the extracted conductances are
  // used one "level" up in the simulation.  The inputs are all voltages,
  // which are applied via voltage sources.  This means that the equation
  // used to complete the system is of the form:
  //
  //   f = V(n) - Vexternal = 0
  //
  //   where V(n) is one of the solution variables in x, and Vexternal
  //   is the voltage applied from above, represented in the above
  //   equations as "v".
  //
  //   Given that this single, simple equation is the entire coupling,
  //   the df/dv derivative, for any Vexternal, is -1.0.

  // Save the old rhs and newton vectors.  (This might not be neccessary).
  // first save a copy of the rhs vector, in case we want it later.
  savedRHSVectorPtr_->putScalar(0.0);
  savedRHSVectorPtr_->addVec(1.0, *(rhsVectorPtr_));

  savedNewtonVectorPtr_->putScalar(0.0);
  savedNewtonVectorPtr_->addVec(1.0, *(NewtonVectorPtr_));

  // Before we try to do any linear solves, check that the Jacobian
  // actually has been loaded.  Sometimes, early in the run, the inner
  // solve will not load the Jacobian, b/c the norm of f is zero, or very
  // small. Typically this will happen if the inner problem is linear,
  // and is entirely driven by the upper level circuit via the connected
  // Vsrc devices.

  jacobianMatrixPtr_->getDiagonal ( *matrixDiagonalPtr_ );
  double diagInfNorm = 0.0;
  matrixDiagonalPtr_->infNorm(&diagInfNorm);

  // if the infinite norm of the diagonal is zero, then the matrix
  // probably hasn't been set up yet.  If that is the case, we should
  // call loadJacobian now.  It should be safe to call this, given
  // that the RHS vector was loaded, no matter what, as part of the
  // nonlinear solve.
  if (diagInfNorm < 1.0e-30)
  {
#ifdef Xyce_DEBUG_CONDUCTANCE
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0,
        "\n\tJacobian for inner problem not loaded.  Forcing a load.\n");
#endif
    nls_.jacobian_();
  }

  // Get dIdx vectors, one for each I-row.
  bool b1 = setup_dIdX_Vectors_(); bsuccess = bsuccess && b1;

  // This loop is over the different applied V's.
  // This loop is also over columns of the small (output) Jacobian.
  map<string,double>::const_iterator iterM = inputMap.begin();
  map<string,double>::const_iterator  endM = inputMap.end  ();
  int iV_col=0;
  for (;iV_col<idSize;++iV_col,++iterM)
  {
    // set up dfdv:
    dfdvVectorPtr_->putScalar(0.0);
    int irow = currentLIDs_[iV_col];
    // note: dfdv = -1.0, but -dfdv needs to go in here.
    if( irow != -1 )
      (*dfdvVectorPtr_)[irow] = 1.0;   // = -dfdv.

    lasSolverPtr_->solve();

    int iC_row=0;
    for (;iC_row<idSize;++iC_row)
    {
      // Get the dot product of dIdx and dxdv, to get dIdv.
      double dIdv = dIdxPtrVector_[iC_row]->dotProduct(*(dxdvVectorPtr_));

#ifdef Xyce_DEBUG_CONDUCTANCE
      if (debugLevel_ > 0)
      {
        string vsrcName = iterM->first;
        printPetraObjects_ (vsrcName);
        cout << "dIdv = " << dIdv << endl;
      }
#endif
      // put dIdV's into the small matrix:
      jacobian[iC_row][iV_col] = dIdv;
    }

  } // cols of output Jacobian (iV_col)

  // Restore the RHS and Newton vectors. (again, this may not be necessary).
  rhsVectorPtr_->putScalar(0.0);
  rhsVectorPtr_->addVec(1.0, *(savedRHSVectorPtr_));

  NewtonVectorPtr_->putScalar(0.0);
  NewtonVectorPtr_->addVec(1.0, *(savedNewtonVectorPtr_));

#ifdef Xyce_VERBOSE_CONDUCTANCE
  printJacobian_ (inputMap,jacobian);
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConductanceExtractor::setupISO_IDs_
//
// Purpose       : This function sets up the various GIDs and matrix rows,
//                 that are needed in the extract function.
//
// Special Notes : This one is specific to iso devices, but uses the
//                 same STL objects as the original setup_IDs_ function.
//
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 03/03/06
//-----------------------------------------------------------------------------
bool N_NLS_ConductanceExtractor::setupISO2_IDs_(const string & isoName)
{
  bool bsuccess = true;

  // Setup vectors for each terminal current
  int idSize = 2; // for now assuming size 2 for iso2 devices.
  currentGIDs_.resize(idSize);
  currentLIDs_.resize(idSize);
  vsrcPosGIDs_.resize(idSize);
  vsrcPosLIDs_.resize(idSize);
  for (int i=0; i<idSize ;++i)
  {
    dIdxPtrVector_.push_back(lasSysPtr_->builder().createVector());
    currentGIDs_[i] = -1;
    currentLIDs_[i] = -1;
    vsrcPosGIDs_[i] = -1;
    vsrcPosLIDs_[i] = -1;
  }

  // Note: When passing the name to topology, it must
  // be all CAPS, or topology won't find it.
  ExtendedString src = isoName;
  src.toUpper();
  string sourceName = src;
  char type;
  int index;

  // This is to get the IDs for the currents through the
  // voltage sources specified in the map.
  list<int> GIDList, extGIDList;
  top_.getNodeSVarGIDs(NodeID(sourceName,_DNODE), GIDList, extGIDList, type);

  if (!(GIDList.empty ()))
  {
    list<int>::iterator gidIter = GIDList.begin();
    list<int>::iterator gidExtIter = extGIDList.begin();

    int i=0;
    for (i=0;i<idSize;++i,++gidIter,++gidExtIter)
    {
      currentGIDs_[i] = *(gidIter);
      currentLIDs_[i] = dIdxPtrVector_[i]->pmap()->globalToLocalIndex(currentGIDs_[i]);
      //cout << "currentGIDs_["<<i<<"]="<< currentGIDs_[i]<<endl;
    }

    // reset i, but keep iterating gidExtIter.
    for (i=0;i<idSize;++i,++gidExtIter)
    {
      vsrcPosGIDs_[i] = *(gidExtIter);
      vsrcPosLIDs_[i] = dIdxPtrVector_[i]->pmap()->globalToLocalIndex(vsrcPosGIDs_[i]);

      //cout << "vsrcPosGIDs_["<<i<<"]="<< vsrcPosGIDs_[i]<<endl;

      if( vsrcPosLIDs_[i] == -1 )
      {
        string msg = "N_NLS_ConductanceExtractor::setupISO_IDs_";
        msg += " The " + sourceName + " source has the positive node"
         " owned by another processor.  The 2-level solve can't handle that.";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }
    }

  }

  // No need to check if the ISO devices are connected to gnd, as
  // they are built that way internally.

#ifdef Xyce_DEBUG_CONDUCTANCE
  if (debugLevel_ > 0)
  {
    cout << "current GIDs: " << endl;
    for( int i1=0; i1 < idSize; ++i1 )
    {
      cout << "  currentGIDs_["<<i1<<"] = " << currentGIDs_[i1] << endl;
    }

    cout << "Vsrc pos equation rows: " << endl;
    for( int i1=0; i1 < idSize; ++i1 )
    {
      cout << "  vsrcPosGIDs_["<<i1<<"] = " << vsrcPosGIDs_[i1] << endl;
    }
  }
#endif

  gidsSetUpFlag_ = true;

  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_NLS_ConductanceExtractor::extract
//
// Purpose       : Slightly different version of the extract function,
//                 intended for ISO-2 devices, rather than Vsrc's.
//
// Special Notes : Unlike the other version of this function, this function
//                 does not need to pull out Vsrc currents.  It only needs
//                 to worry about conductances.  Also, this means that the
//                 function doesn't have to use the solution vector at all.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/03/06
//-----------------------------------------------------------------------------
bool N_NLS_ConductanceExtractor::extract
( const string & isoName, vector< vector<double> > & jacobian )
{
  bool bsuccess = true;

  int idSize = jacobian.size();

#ifdef Xyce_DEBUG_CONDUCTANCE
  const string dashedline2 = "---------------------";
  if (debugLevel_ > 0)
  {
    cout << dashedline2 << endl;
    cout << "N_NLS_ConductanceExtractor::extract - iso2" << endl;
    cout << dashedline2 << endl;
  }
#endif

  if (!gidsSetUpFlag_)
  {
    setupISO2_IDs_(isoName);
  }

  // Set asside the original Newton Vectors:
  savedRHSVectorPtr_->putScalar(0.0);
  savedRHSVectorPtr_->addVec(1.0, *(rhsVectorPtr_));

  savedNewtonVectorPtr_->putScalar(0.0);
  savedNewtonVectorPtr_->addVec(1.0, *(NewtonVectorPtr_));

  // Get dIdx vectors, one for each I-row.
  bool b1 = setup_dIdX_Vectors_(); bsuccess = bsuccess && b1;

  // This loop is over the different applied V's.
  // This loop is also over columns of the small (output) Jacobian.
  int iV_col=0;
  for (;iV_col<idSize;++iV_col)
  {
    // set up dfdv:
    dfdvVectorPtr_->putScalar(0.0);
    int irow = currentLIDs_[iV_col];
    // note: dfdv = -1.0, but -dfdv needs to go in here.
    if( irow != -1 )
      (*dfdvVectorPtr_)[irow] = 1.0;   // = -dfdv.

    lasSolverPtr_->solve();

    int iC_row=0;
    for (;iC_row<idSize;++iC_row)
    {
      // Get the dot product of dIdx and dxdv, to get dIdv.
      double dIdv = dIdxPtrVector_[iC_row]->dotProduct(*(dxdvVectorPtr_));

#ifdef Xyce_DEBUG_CONDUCTANCE
      if (debugLevel_ > 0)
      {
        string tmpName = "unknown";
        printPetraObjects_ (tmpName);
        cout << "dIdv = " << dIdv << endl;
      }
#endif
      // put dIdV's into the small matrix:
      jacobian[iC_row][iV_col] = dIdv;
    }

  } // cols of output Jacobian (iV_col)

  // Restore the RHS and Newton vectors. (again, this may not be necessary).
  rhsVectorPtr_->putScalar(0.0);
  rhsVectorPtr_->addVec(1.0, *(savedRHSVectorPtr_));

  NewtonVectorPtr_->putScalar(0.0);
  NewtonVectorPtr_->addVec(1.0, *(savedNewtonVectorPtr_));

#ifdef Xyce_VERBOSE_CONDUCTANCE
  varMap_["Var 0"] = 0.0;
  varMap_["Var 1"] = 0.0;
  printJacobian_ (varMap_,jacobian);
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConductanceExtractor::printJacobian_
// Purpose       : Prints the small STL jacobian to the screen.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/08/06
//-----------------------------------------------------------------------------
void N_NLS_ConductanceExtractor::printJacobian_
    (const map<string,double> & inputMap,
     vector< vector<double> > & jacobian)
{
  cout.width(15); cout.precision(7); cout.setf(ios::scientific);
  cout << "Output Jacobian/Conductance array: \n";
  cout <<"              ";

  int iE1, iE2;
  int numElectrodes = jacobian.size();
  map<string,double>::const_iterator iterM = inputMap.begin();
  map<string,double>::const_iterator  endM = inputMap.end  ();
  for (iE1 = 0; iE1 < numElectrodes; ++iE1,++iterM)
  {
    cout << "\t"<<iterM->first;
  }
  cout << endl;
  iterM = inputMap.begin();
  for (iE1 = 0; iE1 < numElectrodes; ++iE1,++iterM)
  {
    cout << "\t"<<iterM->first;
    for (iE2 = 0; iE2 < numElectrodes; ++iE2)
    {
      cout << "\t" <<jacobian[iE1][iE2];
    }
    cout << endl;
  }
  cout << endl;

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConductanceExtractor::printPetraObjects_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/08/06
//-----------------------------------------------------------------------------
void N_NLS_ConductanceExtractor::printPetraObjects_ (const string & varName)
{
  cout.width(15); cout.precision(7); cout.setf(ios::scientific);
  string srcName = varName;
  cout << "Info for input voltage: " << srcName << endl;
  cout << "Jacobian:" << endl;
  jacobianMatrixPtr_->printPetraObject();

  // now print out the dxdv vector:
  cout << "dxdv:" << endl;
  dxdvVectorPtr_->printPetraObject();

  cout << "dfdv:" << endl;
  dfdvVectorPtr_->printPetraObject();

  return;
}

