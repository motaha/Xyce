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
// Filename       : $RCSfile: N_NLS_NonLinearSolver.C,v $
//
// Purpose        : Body file which declares an interface common to all
//                  supported nonlinear solver algorithms.  The Manager class
//                  uses this interface to call a concrete algorithm.
//
// Special Notes  : This is the "Strategy" class in the Strategy design
//                  pattern.
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.114.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:48 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

//#include <N_UTL_Misc.h>
#include <N_UTL_OptionBlock.h>
#include <N_NLS_Manager.h>
#include <N_NLS_NonLinearSolver.h>
#include <N_NLS_ConstraintBT.h>
#include <N_NLS_TwoLevelNewton.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_LAS_Solver.h>
#include <N_LAS_Problem.h>
#include <N_LAS_SolverFactory.h>
#include <N_LAS_PrecondFactory.h>

#include <N_LAS_System.h>
#include <N_LAS_Builder.h>

#include <N_ERH_ErrorMgr.h>

#include <N_ANP_AnalysisInterface.h>
#include <N_ANP_AnalysisManager.h>

#include <N_LOA_Loader.h>
#include <N_IO_OutputMgr.h>

#include <N_IO_CmdParse.h>
#include <N_PDS_Manager.h>
#include <N_PDS_ParComm.h>

// Harmonic Balance matrix free stuff
#include <N_NLS_MatrixFreeEpetraOperator.h>

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::N_NLS_NonLinearSolver
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
N_NLS_NonLinearSolver::N_NLS_NonLinearSolver(N_IO_CmdParse &cp)
  : commandLine_(cp),
    nextSolVectorPtrPtr_(0),
    currSolVectorPtrPtr_(0),
    tmpSolVectorPtrPtr_(0),
    rhsVectorPtr_(0),
    jacobianMatrixPtr_(0),
    gradVectorPtr_(0),
    NewtonVectorPtr_(0),
    solWtVectorPtr_(0),
    lasSysPtr_(0),
    lasSolverPtr_(0),
    petraOptionBlockPtr_(0),
    tlnPtr_(0),
    loaderPtr_(0),
    nlpMgrPtr_(0),
    outMgrPtr_(0),
    pdsMgrPtr_(0),
    anaIntPtr_(0),
    netlistFileName_(""),
#ifdef Xyce_DEBUG_NONLINEAR
    outputStepNumber_(0),
    debugTimeFlag_(true),
    contStep_(0),
#endif
    matrixFreeFlag_(false)
{
  N_NLS_NonLinearSolver::resetCountersAndTimers_();

  if (commandLine_.getArgumentValue("netlist") != "")
  {
    netlistFileName_ = commandLine_.getArgumentValue("netlist");
  }
}


//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::~N_NLS_NonLinearSolver
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
N_NLS_NonLinearSolver::~N_NLS_NonLinearSolver()
{
  if (NewtonVectorPtr_)
  {
    delete NewtonVectorPtr_;
    NewtonVectorPtr_ = 0;
  }

  if (gradVectorPtr_)
  {
    delete gradVectorPtr_;
    gradVectorPtr_ = 0;
  }

  if (solWtVectorPtr_)
  {
    delete solWtVectorPtr_;
    solWtVectorPtr_ = 0;
  }

  if (petraOptionBlockPtr_)
  {
    delete petraOptionBlockPtr_;
    petraOptionBlockPtr_ = 0;
  }

  if (lasSolverPtr_)
  {
    delete lasSolverPtr_;
    lasSolverPtr_ = 0;
  }

}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::setPetraOptions
//
// Purpose       : Passes option block corresponding to "LINSOL" onto
//                 nonlinear solver.
//
// Special Notes:  These options are saved to be passed onto the object from
//                 registerLinearSolver() in the initializeAll() function.
//                 This is only called if "NONLIN" is present in the
//                 circuit file.
//
// Return Type   : boolean
//
// - Input Arguments -
//
//    OB         : Option block containing options corresponding to
//                 "LINSOL" in the netlist.
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 11/9/00
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::setPetraOptions(const N_UTL_OptionBlock & OB)
{
  petraOptionBlockPtr_ = new N_UTL_OptionBlock(OB);
  return (petraOptionBlockPtr_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::setDCOPRestartOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 09/17/2007
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::setDCOPRestartOptions(const N_UTL_OptionBlock& OB)
{
  string msg = "DCOP restart options not supported for this solver.  Use nox instead. ";
	N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::setICOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 09/17/2007
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::setICOptions(const N_UTL_OptionBlock& OB)
{
  string msg = ".IC options not supported for this nonlinear solver.  Use nox instead. ";
	N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::setNodeSetOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 09/25/2007
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::setNodeSetOptions(const N_UTL_OptionBlock& OB)
{
  string msg = ".NODESET options not supported for this nonlinear solver.  Use nox instead. ";
	N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::setLocaOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::setLocaOptions (const N_UTL_OptionBlock& OB)
{
  string msg = "N_NLS_NonLinearSolver::setLocaOptions - not implemented for this solver. Use nox instead. ";
	N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return true;
}

//---------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::setTwoLevelLocaOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::setTwoLevelLocaOptions (const N_UTL_OptionBlock& OB)
{
  string msg = "N_NLS_NonLinearSolver::setTwoLevelLocaOptions - not implemented for this solver.  Use nox instead.";
	N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return true;
}

//---------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::setTwoLevelOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//---------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::setTwoLevelOptions (const N_UTL_OptionBlock& OB)
{
  return true;
}

//---------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::setTwoLevelTranOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//---------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::setTwoLevelTranOptions (const N_UTL_OptionBlock& OB)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::registerRHSVector
// Purpose       : 
// Special Notes : 
// Return Type   : boolean
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::registerRHSVector(N_LAS_Vector* tmp_RHSVecPtr)
{
  rhsVectorPtr_ = tmp_RHSVecPtr;
  return (rhsVectorPtr_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::registerLoader
// Purpose       : 
// Special Notes : 
// Return Type   : boolean
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::registerLoader(N_LOA_Loader* tmp_LoaderPtr)
{
  loaderPtr_ = tmp_LoaderPtr;
  return (loaderPtr_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::registerLinearSystem
// Purpose       : 
// Special Notes : 
// Return Type   : boolean
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/09/00
//-----------------------------------------------------------------------------

bool N_NLS_NonLinearSolver::registerLinearSystem(N_LAS_System* tmp_LasSysPtr)
{
  lasSysPtr_ = tmp_LasSysPtr;
  return (lasSysPtr_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::registerPrecondFactory
// Purpose       : 
// Special Notes : 
// Return Type   : boolean
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::registerPrecondFactory(const RefCountPtr<N_LAS_PrecondFactory>& tmp_LasPrecPtr)
{
  lasPrecPtr_ = tmp_LasPrecPtr;
  return (!Teuchos::is_null(lasPrecPtr_));
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::registerParallelMgr
// Purpose       : 
// Special Notes : 
// Return Type   : boolean
// Scope         : public
// Creator       : Eric Keiter, Sandia
// Creation Date : 6/8/2013
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::registerParallelMgr(N_PDS_Manager * pdsMgrPtr)
{
  pdsMgrPtr_ = pdsMgrPtr;
  return (pdsMgrPtr_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::registerAnalysisInterface
// Purpose       : 
// Special Notes : 
// Return Type   : boolean
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/03/00
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::registerAnalysisInterface(N_ANP_AnalysisInterface* tmp_anaIntPtr)
{
  anaIntPtr_ = tmp_anaIntPtr;
  return (anaIntPtr_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::registerOutputMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/23/03
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::registerOutputMgr (N_IO_OutputMgr * outPtr)
{
  outMgrPtr_ = outPtr;
  return (outMgrPtr_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::registerTwoLevelSolver
// Purpose       : This function is called in the event that the two-level
//                 Newton method has been invoked.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/24/02
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::registerTwoLevelSolver
  (N_NLS_TwoLevelNewton * tmp_tlnPtr)
{
  tlnPtr_ = tmp_tlnPtr;
  return (tlnPtr_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::registerParamMgr
// Purpose       : This function is called in the event that the two-level
//                 Newton method has been invoked.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/24/02
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::registerParamMgr
  (N_NLS_ParamMgr * ptr)
{
  nlpMgrPtr_ = ptr;
  return (nlpMgrPtr_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::registerTopology
// Purpose       : This function is called to register a topology object
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, 1437
// Creation Date : 12/3/08
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::registerTopology(N_TOP_Topology * ptr)
{
  topologyRcp_ = rcp( ptr, false);
  return (! Teuchos::is_null(topologyRcp_) );
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::initializeAll
//
// Purpose       : Called after all register and set functions.
//                 Once the various registrations have taken place,
//                 this function sets the remaining pointers.
//
// Special Notes:  This function obtains the solution, temporary solution and
//                 rhs vectors from the LAS system class.
//
// Return Type   : boolean
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/12/00
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::initializeAll()
{
  bool bsuccess = true;


  // Check the registerLinearSystem has been successfully called.
  if (lasSysPtr_ == 0)
    return false;

  // get the temporaries:
  tmpSolVectorPtrPtr_ = lasSysPtr_->getTmpSolVectorPtr();
  bsuccess = bsuccess && (tmpSolVectorPtrPtr_ != 0);

  // get rhs vector:
  rhsVectorPtr_ = lasSysPtr_->getRHSVector();
  bsuccess = bsuccess && (rhsVectorPtr_ != 0);

  // get current solution vectors:
  currSolVectorPtrPtr_ = lasSysPtr_->getCurrSolVectorPtr();
  bsuccess = bsuccess && (currSolVectorPtrPtr_ != 0);

  // get next solution vectors:
  nextSolVectorPtrPtr_ = lasSysPtr_->getNextSolVectorPtr();
  bsuccess = bsuccess && (nextSolVectorPtrPtr_ != 0);

  // get jacobian:
  jacobianMatrixPtr_ = lasSysPtr_->getJacobianMatrix();
  bsuccess = bsuccess && (jacobianMatrixPtr_ != 0);

  N_LAS_Builder & bld_ = lasSysPtr_->builder();

  // create gradient vector:
  //gradVectorPtr_ = lasSysPtr_->builder().createVector();
  gradVectorPtr_ = bld_.createVector();
  bsuccess = bsuccess && (gradVectorPtr_ != 0);

  // create Newton update vector:
  //NewtonVectorPtr_ = lasSysPtr_->builder().createVector();
  NewtonVectorPtr_ = bld_.createVector();
  bsuccess = bsuccess && (NewtonVectorPtr_ != 0);

  // create solution weighting vector:
  //solWtVectorPtr_ = lasSysPtr_->builder().createVector();
  solWtVectorPtr_ = bld_.createVector();
  bsuccess = bsuccess && (solWtVectorPtr_ != 0);

  if( !petraOptionBlockPtr_ )
  {
    petraOptionBlockPtr_ = new N_UTL_OptionBlock();
    petraOptionBlockPtr_->getParams().push_back( N_UTL_Param( "TYPE", "DEFAULT" ) );
  }

  Teuchos::RefCountPtr<N_LAS_Vector> NewtonVectorRCPtr = Teuchos::rcp(NewtonVectorPtr_, false);
  Teuchos::RefCountPtr<N_LAS_Vector> rhsVectorRCPtr = Teuchos::rcp(rhsVectorPtr_, false);

  if (!matrixFreeFlag_)
  {
    Teuchos::RefCountPtr<N_LAS_Matrix> jacobianMatrixRCPtr = Teuchos::rcp(jacobianMatrixPtr_,false);
    // Normal full matrix linear solver options
    lasProblemRCPtr_ = rcp( new N_LAS_Problem( jacobianMatrixRCPtr,
                                    NewtonVectorRCPtr,
                                    rhsVectorRCPtr) );
    // before we make the solver, we need to check and see if this
    // is a parallel binary running on one processor.  If so, we'll
    // need to set extra options so that the correct solver is
    // created
    if( lasSysPtr_->getPDSManager()->getPDSComm()->isSerial() )
    {
      N_UTL_OptionBlock::ParameterList::iterator currentParam = petraOptionBlockPtr_->begin();
      N_UTL_OptionBlock::ParameterList::iterator endParam = petraOptionBlockPtr_->end();
      while( currentParam != endParam )
      {
        if( (currentParam->uTag() == "TYPE") && (currentParam->sVal() == "DEFAULT") )
        {
          currentParam->setVal("KLU");
        }
        currentParam++;
      }
      petraOptionBlockPtr_->getParams().push_back( N_UTL_Param( "TR_PARTITION", 0 ) );

    }
  }
  else
  {
    // Matrix free harmonic balance linear solver options
    // Create MatrixFreeLinearProblem
    Teuchos::RefCountPtr<N_NLS_NonLinearSolver> NonlinearSolverRCPtr = Teuchos::rcp(this, false);
    Teuchos::RefCountPtr<N_NLS_MatrixFreeEpetraOperator>
      matFreeOp = matrixFreeEpetraOperator(
          NonlinearSolverRCPtr,
          NewtonVectorRCPtr,
          rhsVectorRCPtr,
          //lasSysPtr_->builder().getSolutionMap()
          bld_.getSolutionMap()
          );

    Teuchos::RefCountPtr<Epetra_Operator> epetraOperator = Teuchos::rcp_dynamic_cast<Epetra_Operator>(matFreeOp,true);
    // Create N_LAS_Problem
    lasProblemRCPtr_ = rcp( new N_LAS_Problem(
          epetraOperator,
          NewtonVectorRCPtr,
          rhsVectorRCPtr
          )
        );

  }

  lasSolverPtr_ = N_LAS_SolverFactory::create( *petraOptionBlockPtr_,
                                              *lasProblemRCPtr_ , commandLine_);

  // If a preconditioner factory has been provided by the analysis package,
  // use it to generate a preconditioner for the linear solver.
  if (!Teuchos::is_null(lasPrecPtr_)) {
    Teuchos::RefCountPtr<N_LAS_Preconditioner> precond = lasPrecPtr_->create( Teuchos::rcp( lasSysPtr_, false ) );
    lasSolverPtr_->setPreconditioner( precond );
  }

#ifdef Xyce_DEBUG_NONLINEAR
    ostringstream ost;
  ost << "size of solution vector: " << lasSysPtr_->getGlobalSolutionSize();
  ost << endl;
  ost << "size of state vector: " << lasSysPtr_->getGlobalStateSize();
  ost << endl;
  ost << "End of N_NLS_NonLinearSolver::initializeAll\n";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, ost.str());
#endif

  return bsuccess;
}

#ifdef Xyce_DEBUG_NONLINEAR
//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::debugOutput1_
// Purpose       : Write Jacobian to the file matrix.(n).txt, and the residual
//                 to the file rhs.(n).txt
//
// Special Notes : These are objects that will be availabled prior to the
//                 linear solve performed at each Newton step.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/00
//-----------------------------------------------------------------------------
void N_NLS_NonLinearSolver::debugOutput1
   (N_LAS_Matrix & jacobian, N_LAS_Vector & rhs)
{
  int debugLevel = getDebugLevel();
  int newtStep = getNumIterations();
  int screenOutput = getScreenOutputFlag ();
  int contStep = getContinuationStep();
  int paramNumber = getParameterNumber ();

  if (!debugTimeFlag_ || debugLevel < 1) return;

  char filename1[256]; for (int ich = 0; ich < 256; ++ich) filename1[ich] = 0;
  char filename2[256]; for (int ich = 0; ich < 256; ++ich) filename2[ich] = 0;

  if (debugLevel >= 3)
  {
    sprintf(filename1, "matrix_%03d_%03d_%03d_%03d.txt",outputStepNumber_,paramNumber,contStep,newtStep);
  }
  else if (debugLevel == 2)
  {
    sprintf(filename1, "matrix_%03d_%03d.txt",outputStepNumber_,newtStep);
  }
  else
  {
    sprintf(filename1, "matrix_%03d.txt", newtStep);
  }

  jacobian.writeToFile(filename1,false, getMMFormat () );

  if (screenOutput == 1)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0,
                          "\n\t***** Jacobian matrix:\n");
    jacobian.printPetraObject();
  }

  if (debugLevel >= 3)
  {
    sprintf(filename2, "rhs_%03d_%03d_%03d_%03d.txt",outputStepNumber_,paramNumber,contStep,newtStep);
  }
  else if (debugLevel == 2)
  {
    sprintf(filename2, "rhs_%03d_%03d.txt",outputStepNumber_,newtStep);
  }
  else
  {
    sprintf(filename2, "rhs_%03d.txt", newtStep);
  }

  if (screenOutput == 1)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0,
              "\n\t***** RHS vector:\n");

    rhs.printPetraObject();
  }

  rhs.writeToFile(filename2);

#ifdef Xyce_DEBUG_VOLTLIM
  debugOutputJDX_VOLTLIM ();
#endif

  debugOutputDAE ();

}
#endif // Xyce_DEBUG_NONLINEAR

#ifdef Xyce_DEBUG_VOLTLIM
//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::debugOutputJDX_VOLTLIM_
// Purpose       : Write JDX vector to output files.
// Special Notes : This requires a matvec.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/05/05
//-----------------------------------------------------------------------------
void N_NLS_NonLinearSolver::debugOutputJDX_VOLTLIM()
{
  int debugLevel = getDebugLevel();
  int newtStep = getNumIterations();
  int contStep = getContinuationStep();
  int paramNumber = getParameterNumber ();

  char filename1[256]; for (int ich = 0; ich < 256; ++ich) filename1[ich] = 0;
  char filename2[256]; for (int ich = 0; ich < 256; ++ich) filename2[ich] = 0;
  char filename3[256]; for (int ich = 0; ich < 256; ++ich) filename3[ich] = 0;

  if (debugLevel >= 3)
  {
    sprintf(filename1, "jdxVL_%03d_%03d_%03d_%03d.txt",outputStepNumber_,paramNumber,contStep,newtStep);
    sprintf(filename2, "fdxVL_%03d_%03d_%03d_%03d.txt",outputStepNumber_,paramNumber,contStep,newtStep);
    sprintf(filename3, "qdxVL_%03d_%03d_%03d_%03d.txt",outputStepNumber_,paramNumber,contStep,newtStep);
  }
  else if (debugLevel == 2)
  {
    sprintf(filename1, "jdxVL_%03d_%03d.txt",outputStepNumber_,newtStep);
    sprintf(filename2, "fdxVL_%03d_%03d.txt",outputStepNumber_,newtStep);
    sprintf(filename3, "qdxVL_%03d_%03d.txt",outputStepNumber_,newtStep);
  }
  else
  {
    sprintf(filename1, "jdxVL_%03d.txt", newtStep);
    sprintf(filename2, "fdxVL_%03d.txt", newtStep);
    sprintf(filename3, "qdxVL_%03d.txt", newtStep);
  }

  bool Transpose = false;  // if set to true, the matvec does the transpose.

  jdxVLVectorPtr_->putScalar(0.0);
  fdxVLVectorPtr_->putScalar(0.0);
  qdxVLVectorPtr_->putScalar(0.0);

  jacTestMatrixPtr_->matvec( Transpose , *dxVoltlimVectorPtr_, *jdxVLVectorPtr_);
  dFdxTestMatrixPtr_->matvec( Transpose , *dxVoltlimVectorPtr_, *fdxVLVectorPtr_);
  dQdxTestMatrixPtr_->matvec( Transpose , *dxVoltlimVectorPtr_, *qdxVLVectorPtr_);

  jdxVLVectorPtr_->writeToFile(filename1);
  fdxVLVectorPtr_->writeToFile(filename2);
  qdxVLVectorPtr_->writeToFile(filename3);

  if (debugLevel >= 3)
  {
    sprintf(filename1, "jtest_%03d_%03d_%03d_%03d.txt",outputStepNumber_,paramNumber,contStep,newtStep);
    sprintf(filename2, "ftest_%03d_%03d_%03d_%03d.txt",outputStepNumber_,paramNumber,contStep,newtStep);
    sprintf(filename3, "qtest_%03d_%03d_%03d_%03d.txt",outputStepNumber_,paramNumber,contStep,newtStep);
  }
  else if (debugLevel == 2)
  {
    sprintf(filename1, "jtest_%03d_%03d.txt",outputStepNumber_,newtStep);
    sprintf(filename2, "ftest_%03d_%03d.txt",outputStepNumber_,newtStep);
    sprintf(filename3, "qtest_%03d_%03d.txt",outputStepNumber_,newtStep);
  }
  else
  {
    sprintf(filename1, "jtest_%03d.txt", newtStep);
    sprintf(filename2, "ftest_%03d.txt", newtStep);
    sprintf(filename3, "qtest_%03d.txt", newtStep);
  }

  jacTestMatrixPtr_->writeToFile(filename1);
  dFdxTestMatrixPtr_->writeToFile(filename2);
  dQdxTestMatrixPtr_->writeToFile(filename3);

}
#endif

#ifdef Xyce_DEBUG_NONLINEAR
//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::debugOutputDAE
// Purpose       : Write DAE vectors and matrices to output files.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/27/03
//-----------------------------------------------------------------------------
void N_NLS_NonLinearSolver::debugOutputDAE()
{
  int debugLevel = getDebugLevel();
  int newtStep = getNumIterations();
  int contStep = getContinuationStep();
  int paramNumber = getParameterNumber ();

  char filename1[256]; for (int ich = 0; ich < 256; ++ich) filename1[ich] = 0;
  char filename2[256]; for (int ich = 0; ich < 256; ++ich) filename2[ich] = 0;

  char filename4[256]; for (int ich = 0; ich < 256; ++ich) filename4[ich] = 0;
  char filename6[256]; for (int ich = 0; ich < 256; ++ich) filename6[ich] = 0;
  char filename7[256]; for (int ich = 0; ich < 256; ++ich) filename7[ich] = 0;
  char filename8[256]; for (int ich = 0; ich < 256; ++ich) filename8[ich] = 0;
  char filename9[256]; for (int ich = 0; ich < 256; ++ich) filename9[ich] = 0;

  N_LAS_Matrix *dQdx    = lasSysPtr_->getDAEdQdxMatrix ();
  N_LAS_Matrix *dFdx    = lasSysPtr_->getDAEdFdxMatrix ();

  N_LAS_Vector *daeQ    = lasSysPtr_->getDAEQVector();
  N_LAS_Vector *daeF    = lasSysPtr_->getDAEFVector();

  N_LAS_Vector *daeFlim = lasSysPtr_->getdFdxdVpVector ();
  N_LAS_Vector *daeQlim = lasSysPtr_->getdQdxdVpVector ();

  //cout << "In debugOutputDAE" << endl;

  if (debugLevel >= 3)
  {
    sprintf(filename1, "dQdx_%03d_%03d_%03d_%03d.txt"    ,outputStepNumber_,paramNumber,contStep,newtStep);
    sprintf(filename2, "dFdx_%03d_%03d_%03d_%03d.txt"    ,outputStepNumber_,paramNumber,contStep,newtStep);

    sprintf(filename4, "daeQ_%03d_%03d_%03d_%03d.txt"    ,outputStepNumber_,paramNumber,contStep,newtStep);
    sprintf(filename6, "daeF_%03d_%03d_%03d_%03d.txt"    ,outputStepNumber_,paramNumber,contStep,newtStep);

    sprintf(filename8, "daeQlim_%03d_%03d_%03d_%03d.txt"    ,outputStepNumber_,paramNumber,contStep,newtStep);
    sprintf(filename9, "daeFlim_%03d_%03d_%03d_%03d.txt"    ,outputStepNumber_,paramNumber,contStep,newtStep);
  }
  else if (debugLevel >= 2)
  {
    sprintf(filename1, "dQdx_%03d_%03d.txt"    ,outputStepNumber_,newtStep);
    sprintf(filename2, "dFdx_%03d_%03d.txt"    ,outputStepNumber_,newtStep);

    sprintf(filename4, "daeQ_%03d_%03d.txt"    ,outputStepNumber_,newtStep);
    sprintf(filename6, "daeF_%03d_%03d.txt"    ,outputStepNumber_,newtStep);

    sprintf(filename8, "daeQlim_%03d_%03d.txt"    ,outputStepNumber_,newtStep);
    sprintf(filename9, "daeFlim_%03d_%03d.txt"    ,outputStepNumber_,newtStep);
  }
  else
  {
    sprintf(filename1, "dQdx_%03d.txt"    , newtStep);
    sprintf(filename2, "dFdx_%03d.txt"    , newtStep);

    sprintf(filename4, "daeQ_%03d.txt"    , newtStep);
    sprintf(filename6, "daeF_%03d.txt"    , newtStep);

    sprintf(filename8, "daeQlim_%03d.txt"    , newtStep);
    sprintf(filename9, "daeFlim_%03d.txt"    , newtStep);
  }

  // write the matrices:
  dQdx->writeToFile (filename1,false, getMMFormat () );
  dFdx->writeToFile (filename2,false, getMMFormat () );

  // write the vectors:
  daeQ->writeToFile(filename4);
  daeF->writeToFile(filename6);
  daeQlim->writeToFile(filename8);
  daeFlim->writeToFile(filename9);
}
#endif // Xyce_DEBUG_NONLINEAR

#ifdef Xyce_DEBUG_NONLINEAR
//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::debugOutput3_
// Purpose       : Write out the update vector and the new solution.
//
// Special Notes : These are objects that will be availabled *after* the
//                 linear solve performed at each Newton step.  That
//                 differentiates this function from debugOutput1_.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/00
//-----------------------------------------------------------------------------
void N_NLS_NonLinearSolver::debugOutput3
   (N_LAS_Vector & dxVector, N_LAS_Vector & xVector)
{
  int debugLevel = getDebugLevel();
  int nlStep = getNumIterations();
  int contStep = getContinuationStep();
  int paramNumber = getParameterNumber ();

  if (!debugTimeFlag_ || debugLevel < 1) return;

  char filename[256];  for (int ich = 0; ich < 256; ++ich) filename[ich] = 0;

  if (debugLevel >= 3)
  {
    sprintf(filename, "update_%03d_%03d_%03d_%03d.txt",outputStepNumber_,paramNumber,contStep,nlStep);
  }
  else if (debugLevel == 2)
  {
    sprintf(filename, "update_%03d_%03d.txt",outputStepNumber_,nlStep);
  }
  else
  {
    sprintf(filename, "update_%03d.txt", nlStep);
  }
  dxVector.writeToFile(filename);

#if 0
  if (nlParams.getScreenOutputFlag () )
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0,
                         "\n\t***** Update vector:\n\n");
    //searchDirectionPtr_->printPetraObject();
    dxVector.printPetraObject();
  }
#endif


  if (debugLevel >= 3)
  {
    sprintf(filename, "solution_%03d_%03d_%03d_%03d.txt",outputStepNumber_,paramNumber,contStep,nlStep);
  }
  if (debugLevel == 2)
  {
    sprintf(filename,"solution_%03d_%03d.txt",outputStepNumber_,nlStep);
  }
  else
  {
    sprintf(filename, "solution_%03d.txt", nlStep);
  }
  xVector.writeToFile(filename);

#if 0
  if (nlParams.getScreenOutputFlag () )
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0,
                           "\n\t***** Solution vector:\n");
    //x->printPetraObject();
    xVector.printPetraObject();
  }
#endif

}
#endif

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::resetCountersAndTimers_
//
// Purpose       : Resets all the counters and timers in this object.
//
// Scope         : protected
// Creator       : Tamara G. Kolda, SNL, CSMR (8950)
//                 Eric Keiter, SNL, Parallel Computational Sciences. (9233)
// Creation Date : 01/24/02
//-----------------------------------------------------------------------------
void N_NLS_NonLinearSolver::resetCountersAndTimers_()
{
  numJacobianLoads_ = 0;
  numJacobianFactorizations_ = 0;
  numLinearSolves_ = 0;
  numFailedLinearSolves_ = 0;
  numResidualLoads_ = 0;
  totalNumLinearIters_ = 0;
  totalLinearSolveTime_ = 0.0;
  totalResidualLoadTime_ = 0.0;
  totalJacobianLoadTime_ = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::setX0_()
//
// Purpose       : This should be called at the beginning of each nonlinear
//                 iteration. Copies information from nextSolVector (and
//                 related vectors that are important but hidden from the
//                 nonlinear solver) into tmpSolVector.
//
// Return Type   : boolean
// Scope         : protected
// Creator       : Tamara G. Kolda, SNL, CSMR (8950)
//                 Eric Keiter, SNL, Parallel Computational Sciences. (9233)
// Creation Date : 01/24/02
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::setX0_()
{
  anaIntPtr_->equateTmpVectors();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::rhs_()
//
// Purpose       : Calculates the RHS corresponding to the current solution
//                 vector. More specifically, it fills in the content of
//                 RHSVectorPtr_ based on the contents of nextSolVectorPtr_.
//
// Special Notes : The rhsVectorPtr_ is really the NEGATIVE of F(x).
//
// Scope         : private
// Creator       : Tamara G. Kolda, SNL, Compuational Sciences and
//                 Mathematics Research Department
// Creation Date : 06/19/01
//-----------------------------------------------------------------------------

bool N_NLS_NonLinearSolver::rhs_()
{
  loaderPtr_->loadRHS();
  ++numResidualLoads_;
  totalResidualLoadTime_ += loaderPtr_->getResidualTime();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::jacobian_()
//
// Purpose       : Calculates the Jacobian corresponding to the current
//                 solution vector. More specifically, it fills in the
//                 content of jacobianMatrixPtr_ based on the contents of
//                 nextSolVectorPtr_.
//
// Return Type   : boolean
// Scope         : private
// Creator       : Tamara G. Kolda, SNL, Compuational Sciences and
//                 Mathematics Research Department
// Creation Date : 06/19/01
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::jacobian_()
{
  loaderPtr_->loadJacobian();
  ++numJacobianLoads_;
  totalJacobianLoadTime_ += loaderPtr_->getJacobianTime();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::applyJacobian()
//
// Purpose       : Applies the Jacobian corresponding to the current
//                 solution vector.
//
// Return Type   : boolean
// Scope         : private
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::applyJacobian(const N_LAS_Vector& input, N_LAS_Vector& result)
{
  loaderPtr_->applyJacobian(input,result);
  ++numJacobianLoads_;
  totalJacobianLoadTime_ += loaderPtr_->getJacobianTime();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::newton_
//
// Purpose       : Calculates the Newton direction corresponding to the
//                 current RHS and Jacobian matrix.
//
// Return Type   : boolean
//
// Scope         : private
// Creator       : Tamara G. Kolda, SNL, Compuational Sciences and
//                 Mathematics Research Department
// Creation Date : 06/19/01
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::newton_()
{
  int solutionStatus = lasSolverPtr_->solve( false );

  totalLinearSolveTime_ += lasSolverPtr_->solutionTime();
  ++numLinearSolves_;

  if( lasSolverPtr_->isIterative() )
  {
    N_UTL_Param param( "Iterations", 0 );
    lasSolverPtr_->getInfo( param );
    totalNumLinearIters_ += param.iVal();

    if( solutionStatus ) ++numFailedLinearSolves_;
  }
  else
  {
    N_UTL_Param param( "Refactored", 0 );
    lasSolverPtr_->getInfo( param );
    if( param.iVal() ) ++numJacobianFactorizations_;
    if( solutionStatus ) ++numFailedLinearSolves_;
  }

  if( solutionStatus ) return false;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::gradient_()
// Purpose       : Calculates the Gradient direction corresponding to the
//                 current RHS and Jacobian matrix.
//                 Computes gradient using jacobianMatrixPtr_ and
//                 the rhsVectorPtr_. On output, gradVectorPtr_ contains
//                 the gradient of 0.5 * ||F(x)||^2.
//
// Special Notes : The rhsVectorPtr_ is really the NEGATIVE of F(x).
//
// Return Type   : boolean
// Scope         : private
// Creator       : Tamara G. Kolda, SNL, Compuational Sciences and
//                 Mathematics Research Department
// Creation Date : 06/19/01
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::gradient_()
{
  // Compute gradient = Jacobian' * RHS
  bool transpose = true;
  jacobianMatrixPtr_->matvec(transpose, *rhsVectorPtr_, *gradVectorPtr_);

  // We have to scale by -1 because the rhsVectorPtr_ is really the
  // NEGATIVE of F(X). This gives us gradient = Jacobian' * F(X).
  gradVectorPtr_->scale(-1.0);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::getCouplingMode ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Compuational Sciences and
// Creation Date : 12/05/02
//-----------------------------------------------------------------------------
TwoLevelNewtonMode N_NLS_NonLinearSolver::getCouplingMode ()
{
  return FULL_PROBLEM;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::setMatrixFreeFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/29/08
//-----------------------------------------------------------------------------
void N_NLS_NonLinearSolver::setMatrixFreeFlag (bool matrixFreeFlag)
{
  matrixFreeFlag_ = matrixFreeFlag;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::getMatrixFreeFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/29/08
//-----------------------------------------------------------------------------
bool N_NLS_NonLinearSolver::getMatrixFreeFlag ()
{
  return matrixFreeFlag_;
}

#ifdef Xyce_DEBUG_NONLINEAR
//-----------------------------------------------------------------------------
// Function      : N_NLS_NonLinearSolver::setDebugFlags
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric R. Keiter, SNL, Compuational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
void N_NLS_NonLinearSolver::setDebugFlags()
{
  outputStepNumber_ = anaIntPtr_->getAnalysisMgr ()->getStepNumber()+1;

  double currTime = anaIntPtr_->getTime();

  debugTimeFlag_ =
  (currTime       >= getDebugMinTime()  &&
   currTime       <= getDebugMaxTime() ) &&
   (outputStepNumber_ >= getDebugMinTimeStep() &&
    outputStepNumber_ <= getDebugMaxTimeStep());

  bool steadyStateFlag = anaIntPtr_->getAnalysisMgr ()->getSteadyStateFlag ();

  if (tlnPtr_ != 0)
    contStep_ = tlnPtr_->getContStepNumber();
  else
    contStep_ = 0;

  if (steadyStateFlag == true)
  {
    outputStepNumber_ -= 1;
  }
}
#endif

