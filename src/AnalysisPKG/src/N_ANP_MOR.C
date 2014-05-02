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
//-------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: N_ANP_MOR.C,v $
// Purpose       : MOR analysis functions.
// Special Notes :
// Creator       : Ting Mei
// Creation Date :  7/11
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.37.2.1 $
// Revision Date  : $Date: 2014/03/04 23:50:53 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iomanip>


// ----------   Xyce Includes   ----------
#include <N_ANP_AnalysisManager.h>
#include <N_TIA_Assembler.h>
#include <N_LOA_Loader.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_IO_CmdParse.h>
#include <N_TIA_StepErrorControl.h>
#include <N_UTL_Timer.h>
#include <N_UTL_ExpressionData.h>
#include <N_NLS_ReturnCodes.h>

#include <N_UTL_Misc.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockSystemHelpers.h>

#include <N_ANP_MOR.h>
#include <N_ANP_Report.h>
#include <N_LAS_MOROperators.h>

#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>

#ifdef Xyce_PARALLEL_MPI
#include <N_PDS_ParComm.h>
#include <mpi.h>
#else
#include <N_PDS_SerialComm.h>
#endif

// ----------   Other Includes   ----------

#include <Epetra_CrsMatrix.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_MultiVector.h>
#include <Epetra_LinearProblem.h>
#include <Amesos.h>

#ifdef Xyce_BELOS
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresIter.hpp>
#include <BelosDGKSOrthoManager.hpp>
#include <BelosStatusTestMaxIters.hpp>
#include <BelosOutputManager.hpp>
#include <BelosEpetraAdapter.hpp>
#endif

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_LAPACK.hpp>
#include <Teuchos_BLAS.hpp>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : MOR::MOR( AnalysisManager * )
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 5/11
//-----------------------------------------------------------------------------
MOR::MOR( AnalysisManager * anaManagerPtr ) :
  AnalysisBase(anaManagerPtr),
  dcopFlag_(true),
  morEvalSize_(0),
  numPorts_(0),
  stepMult_(0.0),
  fStep_(0.0),
  currentFreq_(0.0),
  s0_(0.0)
{
}

//-----------------------------------------------------------------------------
// Function      : MOR::~MOR()
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 5/11
//-----------------------------------------------------------------------------
MOR::~MOR()
{
}

//-----------------------------------------------------------------------------
// Function      : MOR::setAnalysisParams
// Purpose       :
// Special Notes : These are from the .MOR statement.
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 6/11
//-----------------------------------------------------------------------------
bool MOR::setAnalysisParams(const N_UTL_OptionBlock & paramsBlock)
{
  std::list<N_UTL_Param>::const_iterator it_tp;
  std::list<N_UTL_Param>::const_iterator first = paramsBlock.getParams().begin();
  std::list<N_UTL_Param>::const_iterator last = paramsBlock.getParams().end();
  for (it_tp = first; it_tp != last; ++it_tp)
  {
    if (it_tp->uTag()      == "SIZE")
    {
      tiaParams.ROMsize = it_tp->getImmutableValue<int>();
    }
    else if (it_tp->uTag() == "PORTLIST")
    {
      portList_         = it_tp->getValue<std::vector<std::string> >();
    }
  }

//  tiaParams.debugLevel = 1;
#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    dout() << std::endl
           << section_divider << std::endl
           <<" MOR simulation parameters" << std::endl
           << " size = " << tiaParams.ROMsize << std::endl;
  }
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : MOR::run()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool MOR::run()
{
  bool bsuccess = true;

  bsuccess = bsuccess & init();
  bsuccess = bsuccess & reduceSystem();

  if (tiaParams.morCompOrigTF)
  {
    // Evaluate original system
    bsuccess = bsuccess & evalOrigTransferFunction();
  }

  // Reset the output adapter, just in case the reduced system needs to output
  // its transfer functions.
  outputMgrAdapterRCPtr_->resetOutputMORTF();

  if (tiaParams.morCompRedTF)
  {
    // Evaluate reduced system
    bsuccess = bsuccess & evalRedTransferFunction();
  }

  // if processing the loop failed,
  // then skip finish step
  if( bsuccess )
  {
    bsuccess = bsuccess & finish();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool MOR::init()
{
  bool bsuccess = true;

  // Compute expansion point in transformed space.
  s0_ =  2.0 * M_PI * tiaParams.morExpPoint;

  if (tiaParams.morCompOrigTF || tiaParams.morCompRedTF)
  {
    morEvalSize_ = setupSweepParam_();
  }

  if(dcopFlag_)
  {
    anaManagerRCPtr_->currentMode_ = 0;
  }

  // Get set to do the operating point.
  integrationMethod_ = TIAMethod_NONE;
  wimRCPtr_->createTimeIntegMethod(integrationMethod_);

  stepNumber            = 0;
  doubleDCOPFlag_ = loaderRCPtr_->getDoubleDCOPFlag ();
  doubleDCOPStep_ = tiaParams.firstDCOPStep;

  // set initial guess, if there is one to be set.
  // this setInitialGuess call is to up an initial guess in the
  // devices that have them (usually PDE devices).  This is different than
  // the "intializeProblem" call, which sets IC's.  (initial conditions are
  // different than initial guesses.
  loaderRCPtr_->setInitialGuess (anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr);

  // If available, set initial solution
  inputOPFlag_ =
    outputMgrAdapterRCPtr_->setupInitialConditions( *(anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr),
                        *(anaManagerRCPtr_->getTIADataStore()->flagSolutionPtr));

  // Set a constant history for operating point calculation
  anaManagerRCPtr_->getTIADataStore()->setConstantHistory();
  anaManagerRCPtr_->getTIADataStore()->computeDividedDifferences();
  wimRCPtr_->obtainCorrectorDeriv();

  // solving for DC op
  handlePredictor();
  loaderRCPtr_->updateSources();
  secRCPtr_->newtonConvergenceStatus = nlsMgrRCPtr_->solve();
  anaManagerRCPtr_->getTIADataStore()->stepLinearCombo ();
  gatherStepStatistics_ ();
  secRCPtr_->evaluateStepError ();

  if ( secRCPtr_->newtonConvergenceStatus <= 0)
  {
    std::string msg=
        "Solving for DC operating point failed! Cannot continue MOR analysis.";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

  // Create B matrix stamp
  std::vector<int> tempVec;
  loaderRCPtr_->getBMatrixEntriesforMOR(tempVec, bMatPosEntriesVec_);

  // Create inductor  matrix stamp
//  loaderRCPtr_->getInductorsEntriesforMOR(inductorEntriesVec_);

  // Determine how many global ports their are performing a sumAll()
  int hsize = tempVec.size();
  N_PDS_Comm * pdsComm = (lasSystemRCPtr_->getPDSManager())->getPDSComm();
  pdsComm->sumAll( &hsize, &numPorts_, 1 );

  // Check that the listed ports on the .mor line are the same number as the B matrix entries.
  if ( numPorts_ != (int)portList_.size() )
  {
    std::string msg=
        "Number of specified ports in .MOR line is inconsistent with number of voltage sources.";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

  // Getting the GIDs for the port list in the order specified in the .MOR line
  // This is used to permute the bMatEntriesVec_ in the order specified by the user.
  std::vector<int> gidPosEntries( hsize );
  for (int i=0; i<hsize; ++i)
  {
    std::list<int> svGIDList1, dummyList;
    char type1;
    anaManagerRCPtr_->topoMgrPtr->getNodeSVarGIDs(NodeID(portList_[i],_VNODE), svGIDList1, dummyList, type1);

    // Grab the GID for this port.
    gidPosEntries[i] = svGIDList1.front();
  }

  // Use the base map to get the global IDs
  // Get the parallel maps for the original system
  RCP<N_PDS_Manager> pdsMgrPtr_ = rcp(lasSystemRCPtr_->getPDSManager(), false);
  RCP<N_PDS_ParMap> BaseMap_ = rcp(pdsMgrPtr_->getParallelMap( "SOLUTION" ), false);

  // Find the voltage source corresponding to each port and place the LID in the bMatEntriesVec_.
  bMatEntriesVec_.resize( hsize );
  for (int i=0; i<hsize; ++i)
  {
    int gid = gidPosEntries[i];
    bool found = false;
    for (int j=0; j<hsize; ++j)
    {
      if (gid == BaseMap_->localToGlobalIndex(bMatPosEntriesVec_[j]))
      {
        bMatEntriesVec_[i] = tempVec[j];
        found = true;
        break;
      }
    }
    if (found == false)
    {  std::string msg=
         "Did not find voltage source corresponding to port.";
         N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    }
  }

  if (tiaParams.morCompOrigTF)
  {
    // Resize transfer function matrices
    origH_.shape(numPorts_, numPorts_);
  }
  if (tiaParams.morCompRedTF)
  {
    // Resize transfer function matrices
    redH_.shape(numPorts_, numPorts_);
  }

  // Create C and G matrices from DCOP solution
  anaManagerRCPtr_->getTIADataStore()->daeQVectorPtr->putScalar(0.0);
  anaManagerRCPtr_->getTIADataStore()->daeFVectorPtr->putScalar(0.0);

  anaManagerRCPtr_->getTIADataStore()->dFdxdVpVectorPtr->putScalar(0.0);
  anaManagerRCPtr_->getTIADataStore()->dQdxdVpVectorPtr->putScalar(0.0);
  anaManagerRCPtr_->getTIADataStore()->dQdxMatrixPtr->put(0.0);
  anaManagerRCPtr_->getTIADataStore()->dFdxMatrixPtr->put(0.0);

  loaderRCPtr_->updateState
                ((anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr),
               (anaManagerRCPtr_->getTIADataStore()->currSolutionPtr),
               (anaManagerRCPtr_->getTIADataStore()->lastSolutionPtr),
               (anaManagerRCPtr_->getTIADataStore()->nextStatePtr),
               (anaManagerRCPtr_->getTIADataStore()->currStatePtr),
               (anaManagerRCPtr_->getTIADataStore()->lastStatePtr),
               (anaManagerRCPtr_->getTIADataStore()->nextStorePtr),
               (anaManagerRCPtr_->getTIADataStore()->currStorePtr),
               (anaManagerRCPtr_->getTIADataStore()->lastStorePtr)
               );

  loaderRCPtr_->loadDAEVectors
              ((anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr),
               (anaManagerRCPtr_->getTIADataStore()->currSolutionPtr),
               (anaManagerRCPtr_->getTIADataStore()->lastSolutionPtr),
               (anaManagerRCPtr_->getTIADataStore()->nextStatePtr),
               (anaManagerRCPtr_->getTIADataStore()->currStatePtr),
               (anaManagerRCPtr_->getTIADataStore()->lastStatePtr),
               (anaManagerRCPtr_->getTIADataStore()->nextStateDerivPtr),
               (anaManagerRCPtr_->getTIADataStore()->nextStorePtr),
               (anaManagerRCPtr_->getTIADataStore()->currStorePtr),
               (anaManagerRCPtr_->getTIADataStore()->lastStorePtr),
               (anaManagerRCPtr_->getTIADataStore()->nextStoreLeadCurrQCompPtr),
               (anaManagerRCPtr_->getTIADataStore()->daeQVectorPtr),
               (anaManagerRCPtr_->getTIADataStore()->daeFVectorPtr),
               (anaManagerRCPtr_->getTIADataStore()->dFdxdVpVectorPtr),
               (anaManagerRCPtr_->getTIADataStore()->dQdxdVpVectorPtr) );

  loaderRCPtr_->loadDAEMatrices(anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr,
      anaManagerRCPtr_->getTIADataStore()->nextStatePtr, anaManagerRCPtr_->getTIADataStore()->nextStateDerivPtr,
      anaManagerRCPtr_->getTIADataStore()->nextStorePtr,
      anaManagerRCPtr_->getTIADataStore()->dQdxMatrixPtr,  anaManagerRCPtr_->getTIADataStore()->dFdxMatrixPtr);

  CPtr_ = rcp(anaManagerRCPtr_->getTIADataStore()->dQdxMatrixPtr, false);
  GPtr_ = rcp(anaManagerRCPtr_->getTIADataStore()->dFdxMatrixPtr, false);


  ///  Xyce::dout() << "Branch nodes: " << std::endl;
  ///  for (unsigned int i=0; i < bMatEntriesVec_.size(); ++i)
  ///  {
  ///    Xyce::dout() << "Node " << i << " : " << bMatEntriesVec_[i] << std::endl;
  ///  }

  ///  Xyce::dout() << "Printing GPtr: " << std::endl;
  ///  GPtr_->printPetraObject();
  ///  Xyce::dout() << "Printing CPtr: " << std::endl;
  ///  CPtr_->printPetraObject();


  // Storage for row extraction
  int length=1, numEntries=0;
  std::vector<int> colIndices(length);
  std::vector<double> coeffs(length);

  for (unsigned int i=0; i < bMatEntriesVec_.size(); ++i)
  {
     // Get the number of non-zero entries in this row.
     numEntries = GPtr_->getLocalRowLength( bMatEntriesVec_[i] );
     if ( numEntries != 1 )
     {
       std::string msg=
         "Supposed voltage source row has too many entries! Cannot continue MOR analysis.";
         N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
     }

     // Extract out rows of G based on indices in bMatEntriesVec_.
     GPtr_->getLocalRowCopy(bMatEntriesVec_[i], length, numEntries, &coeffs[0], &colIndices[0]);

     // If the coefficient for this voltage source is positive, make it negative.
     if ( coeffs[0] > 0.0 )
     {
       coeffs[0] *= -1.0;
       GPtr_->putLocalRow(bMatEntriesVec_[i], length, &coeffs[0], &colIndices[0]);
     }
  }


  ///  Xyce::dout() << "Printing GPtr (after scaling): " << std::endl;
  ///  GPtr_->printPetraObject();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::reduceSystem()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 6/4/2012
//-----------------------------------------------------------------------------
bool MOR::reduceSystem()
{
  // At this time, the reduceSystem() method will compute the reduced-order system using
  // multi-port PRIMA.  This method should allow the user multiple algorithms, defined by
  // the tiaParams.morMethod option, which would be best implemented using a factory design.

  bool bsuccess = true;

  // Get the parallel maps for the original system
  RCP<N_PDS_Manager> pdsMgrPtr_;
  pdsMgrPtr_ = rcp(lasSystemRCPtr_->getPDSManager(), false);

  RCP<N_PDS_ParMap> BaseMap_;
  BaseMap_ = rcp(pdsMgrPtr_->getParallelMap( "SOLUTION" ), false);
  //BaseMap_->petraMap()->Print(Xyce::dout());

  // Create matrix for (G + s0 * C)
  if (s0_ != 0.0)
  {
    sCpG_MatrixPtr_ = rcp( new N_LAS_Matrix( &CPtr_->epetraObj(), false ) );
    sCpG_MatrixPtr_->scale( s0_ );
    sCpG_MatrixPtr_->add( *GPtr_ );
  }
  else
  {
    sCpG_MatrixPtr_ = rcp( new N_LAS_Matrix( &GPtr_->epetraObj(), false ) );
  }

  // Create multivector for B and R
  RPtr_ = rcp( new N_LAS_MultiVector( *BaseMap_, numPorts_ ) );
  RPtr_->putScalar( 0.0 );
  BPtr_ = rcp( new N_LAS_MultiVector( *BaseMap_, numPorts_ ) );
  for (unsigned int j=0; j < bMatEntriesVec_.size(); ++j)
  {
    BPtr_->setElementByGlobalIndex( BaseMap_->localToGlobalIndex(bMatEntriesVec_[j]), -1.0, j );
  }

  ///  Xyce::dout() << "Printing out BPtr" << std::endl;
  ///  BPtr_->epetraObj().Print(Xyce::dout());

  ///  Xyce::dout() << "Printing out sCpG" << std::endl;
  ///  (sCpG_MatrixPtr_->epetraObj()).Print(Xyce::dout());

  // Create linear problem for (G + s0 * C)
  origProblem_ = rcp(new Epetra_LinearProblem(&sCpG_MatrixPtr_->epetraObj(), &RPtr_->epetraObj(), &BPtr_->epetraObj() ) );

  // Create solver object for this linear problem, which will be used to generate the projection basis
  Amesos amesosFactory;
  origSolver_ = rcp( amesosFactory.Create( "Klu", *origProblem_ ) );

  // Solve for R = inv(G + s0*C)*B
  // Perform symbolic factorization.
  int linearStatus = origSolver_->SymbolicFactorization();
  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos symbolic factorization exited with error: " << linearStatus << std::endl;
    bsuccess = false;
  }

  // Perform numeric factorization
  linearStatus = origSolver_->NumericFactorization();
  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos numeric factorization exited with error: " << linearStatus << std::endl;
    bsuccess = false;
  }

  // Perform solve for R = inv(G + s0*C)*B
  linearStatus = origSolver_->Solve();
  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos solve exited with error: " << linearStatus << std::endl;
    bsuccess = false;
  }

  ///  Xyce::dout() << "Printing out R" << std::endl;
  ///  (RPtr_->epetraObj()).Print(Xyce::dout());

  // Create an Epetra_Operator object to apply the operator inv(G + s0*C)*C
  RCP<Epetra_Operator> COpPtr_ = rcp( &CPtr_->epetraObj(), false );
  RCP<N_LAS_AmesosGenOp> AOp = rcp( new N_LAS_AmesosGenOp( origSolver_, COpPtr_ ) );

  // Check to see if the requested size of the ROM is valid.
  // An orthogonal basis cannot be generated that is larger than the dimension of the original system
  int kblock = 0;
  if (tiaParams.ROMsize > RPtr_->globalLength())
  {
    kblock = (int)(RPtr_->globalLength() / numPorts_);
    UserWarning(*this) << "Requested reduced-order model dimension is larger than original system dimension, resizing to original system dimension";
  }
  else
  {
    kblock = (int)(tiaParams.ROMsize / numPorts_);
  }
  int k = kblock * numPorts_;

  // Resize the projection matrices
  redG_.shape(k, k);
  redC_.shape(k, k);
  redB_.shape(k, numPorts_);
  redL_.shape(k, numPorts_);

#ifdef Xyce_BELOS
  // ---------------------------------------------------------------------
  // Now use Belos to compute the basis vectors for K_k(inv(G + s0*C)*C, R)
  // ---------------------------------------------------------------------

  // Helpful typedefs for the templates
  typedef double                            ST;
  typedef Epetra_MultiVector                MV;
  typedef Epetra_Operator                   OP;
  typedef Belos::MultiVecTraits<ST,MV>      MVT;

  // Output manager.
  Belos::OutputManager<ST> printer;

  // Status test.
  Belos::StatusTestMaxIters<ST, MV, OP> maxIterTest( kblock );

  // Orthogonalization manager.
  Belos::DGKSOrthoManager<ST, MV, OP> orthoMgr;

  // Linear Problem.
  // Reuse RPtr_.  We need the basis vectors for K(inv(G + s0*C)*C, R)
  N_LAS_MultiVector temp( *BaseMap_, numPorts_ );
  Belos::LinearProblem<ST, MV, OP > problem( AOp,
                                             rcp( &temp.epetraObj(), false ),
                                             rcp( &RPtr_->epetraObj(), false ));
  problem.setProblem();

  // Create parameter list.
  Teuchos::ParameterList params;
  params.set("Num Blocks", kblock);
  params.set("Block Size", numPorts_ );

  // Create Krylov subspace iteration from Belos
  Belos::BlockGmresIter<ST, MV, OP> krylovIter( rcp( &problem, false ),
                                                rcp( &printer, false ),
                                                rcp( &maxIterTest, false ),
                                                rcp( &orthoMgr, false ),
                                                params );

  // Get a matrix to hold the orthonormalization coefficients.
  RCP<Teuchos::SerialDenseMatrix<int,ST> > z
    = rcp( new Teuchos::SerialDenseMatrix<int,ST>( numPorts_, numPorts_ ) );

  // Orthonormalize the Krylov kernel vectors V
  RCP<MV> V = rcp( new Epetra_MultiVector( RPtr_->epetraObj() ) );
  orthoMgr.normalize( *V, z );

  // Set the new state and initialize the solver to use R as the kernel vectors.
  Belos::GmresIterationState<ST,MV> initState;
  initState.V = V;
  initState.z = z;
  initState.curDim = 0;
  krylovIter.initializeGmres(initState);

  // Have the solver iterate until the basis size is computed
  try {
    krylovIter.iterate();
  }
  catch (const Belos::GmresIterationOrthoFailure &e) {
    // This might happen if the basis size is the same as the operator's dimension.
  }
  catch (const std::exception &e) {
    bsuccess = false; // Anything else is a failure.
  }

  // Get the basis vectors back from the iteration object
  Belos::GmresIterationState<ST,MV> newState = krylovIter.getState();

  // Create a view into the first k vectors of newState.V
  std::vector<int> indices( k );
  for (int i=0; i<k; ++i) { indices[i] = i; }
  V = rcp( new Epetra_MultiVector( View, *newState.V, &indices[0], k ) );

  RCP<MV> W = rcp( new Epetra_MultiVector( V->Map(), k ) );
  W = V;
//  Xyce::dout() << "Printing out V" << std::endl;
//  V->Print(Xyce::dout());
//  Xyce::dout() << "Printing out W" << std::endl;
//  W->Print(Xyce::dout());

  // ---------------------------------------------------------------------
  // Sparsify the reduced system, if requested.
  // ---------------------------------------------------------------------
  if (tiaParams.morSparsificationType)
  {

    N_LAS_MultiVector xyceV( &*V, false );
    N_LAS_MultiVector temp2( *BaseMap_, k );

    // G * V
    GPtr_->matvec( false, xyceV, temp2 );
    // V' * G * V
    MVT::MvTransMv( 1.0, xyceV.epetraObj(), temp2.epetraObj(), redG_ );
    //redG_.print( Xyce::dout() );

    // C * V
    CPtr_->matvec( false, xyceV, temp2 );
    // V' * C * V
    MVT::MvTransMv( 1.0, xyceV.epetraObj(), temp2.epetraObj(), redC_ );
    //redC_.print( Xyce::dout() );

    // V' * B
    MVT::MvTransMv( 1.0, xyceV.epetraObj(), BPtr_->epetraObj(), redB_ );
    //redB_.print( Xyce::dout() );

    bsuccess = bsuccess & sparsifyRedSystem_();
  }

  // -----------------------------------------------------------------------------
  // Scale the basis vectors, if requested, before computing the projected system
  // -----------------------------------------------------------------------------

  // Create the scaling vector.
  Teuchos::SerialDenseMatrix<int,ST> Xhatscale( k, 1 );

  int scaleType = tiaParams.morScaleType;

  if ( scaleType == 1 )
  {
    for (int i=0; i<k; ++i)
    {
      Xhatscale( i, 0 ) = tiaParams.morScaleFactor;
    }
  }

//  if ( scaleType == 2 )
//  {

//  }

if ( scaleType == 2 || scaleType == 3  || scaleType ==4)
  {
    Epetra_MultiVector Xmag( V->Map(), 1 );

    if ( scaleType == 2 )
      Xmag.PutScalar( tiaParams.morScaleFactor );

    if ( scaleType == 4)
      Xmag.PutScalar(1.0);

    MVT::MvTransMv( 1.0, *V, Xmag, Xhatscale );

//    Xhatscale.print(Xyce::dout());
    for (int i=0; i<k; ++i)
    {

      if ( scaleType == 2 )
      {
        if ( fabs(Xhatscale( i, 0 )) >  tiaParams.morScaleFactor)
        {
          Xhatscale( i, 0 ) = tiaParams.morScaleFactor / fabs( Xhatscale( i, 0 ) );
        }
        else
        {
          Xhatscale( i, 0 ) =  tiaParams.morScaleFactor1 / fabs( Xhatscale( i, 0 ) );
        }
      }

      if ( scaleType == 3 )
        Xhatscale( i, 0 ) = tiaParams.morScaleFactor * fabs( Xhatscale( i, 0 ) );

      if ( scaleType == 4 )
      {
        if ( fabs(Xhatscale( i, 0 )) > 1.0 )
        {
          Xhatscale( i, 0 ) = tiaParams.morScaleFactor / fabs( Xhatscale( i, 0 ) );
        }
        else
        {
          Xhatscale( i, 0 ) =  tiaParams.morScaleFactor1 / fabs( Xhatscale( i, 0 ) );
        }
      }


    }
  }

//  Xhatscale.print(Xyce::dout());

  // Scale the computed basis vectors before projecting the original system
  Teuchos::SerialDenseMatrix<int,ST> D(k,k);
  if ( scaleType != 0 )
  {
    RCP<MV> Vtemp = rcp( new Epetra_MultiVector( V->Map(), k ) );
//    Teuchos::SerialDenseMatrix<int,ST> D(k,k);
    for (int i=0; i<k; ++i)
    {
      if (Xhatscale( i, 0 ) != 0.0)
        D(i,i) = 1.0/Xhatscale( i, 0 );
      else
        D(i,i) = 1.0;
    }
    Xyce::dout() << " the scaling matrix "  <<  std::endl;

//    D.print(Xyce::dout());

    MVT::MvTimesMatAddMv( 1.0, *V, D, 0.0, *Vtemp );
    V = Vtemp;
  }

//  Xyce::dout() << "Printing out V" << std::endl;
//  V->Print(Xyce::dout());
//  Xyce::dout() << "Printing out W" << std::endl;
//  W->Print(Xyce::dout());

  // ---------------------------------------------------------------------
  // Now use the basis vectors for K_k(inv(G + s0*C)*C, R) to compute
  // the projected system.
  // ---------------------------------------------------------------------
  N_LAS_MultiVector xyceV( &*V, false );
  N_LAS_MultiVector xyceW( &*W, false );
  N_LAS_MultiVector temp2( *BaseMap_, k );


//  redL_ = redB_;

  if (!tiaParams.morSparsificationType)
  {
//    redL_ = redB_;
  // G * V
    GPtr_->matvec( false, xyceV, temp2 );
  // V' * G * V
    MVT::MvTransMv( 1.0, xyceW.epetraObj(), temp2.epetraObj(), redG_ );
    //redG_.print( Xyce::dout() );

  // C * V
    CPtr_->matvec( false, xyceV, temp2 );
  // V' * C * V
    MVT::MvTransMv( 1.0, xyceW.epetraObj(), temp2.epetraObj(), redC_ );
  //redC_.print( Xyce::dout() );

  // V' * B
    MVT::MvTransMv( 1.0, xyceW.epetraObj(), BPtr_->epetraObj(), redB_ );
  //redB_.print( Xyce::dout() );

    MVT::MvTransMv( 1.0, xyceV.epetraObj(), BPtr_->epetraObj(), redL_ );


//  Xyce::dout() << "Printing out redB" << std::endl;
//  redB_.print(Xyce::dout());
//  Xyce::dout() << "Printing out redL" << std::endl;
//  redL_.print(Xyce::dout());
  }
  else
  {
//    RCP<MV> Vtemp = rcp( new Epetra_MultiVector( V->Map(), k ) );

//    redL_.print(Xyce::dout());

    if ( scaleType <= 1 )
    {
      redCPtr_->scale(1.0/tiaParams.morScaleFactor);
      redGPtr_->scale(1.0/tiaParams.morScaleFactor);

      redL_.scale(1.0/tiaParams.morScaleFactor);
//      redL_.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, D, redL_, 0.0 );

//      redL_.print(Xyce::dout());

    }
    else
    {
      std::string msg=
      "MOR::reduceSystem: MOR options sparsificationType=1 can only be used with scaletype=1. Other scale types have not been supported for sparsification!";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    }

  }
//  // ---------------------------------------------------------------------
//   // Sparsify the reduced system, if requested.
//   // ---------------------------------------------------------------------
//   if (tiaParams.morSparsificationType)
//   {

//     N_LAS_MultiVector xyceV( &*V, false );
//     N_LAS_MultiVector temp2( *BaseMap_, k );

//     // G * V
//     GPtr_->matvec( false, xyceV, temp2 );
//     // V' * G * V
//     MVT::MvTransMv( 1.0, xyceV.epetraObj(), temp2.epetraObj(), redG_ );
//     //redG_.print( Xyce::dout() );

//     // C * V
//     CPtr_->matvec( false, xyceV, temp2 );
//     // V' * C * V
//     MVT::MvTransMv( 1.0, xyceV.epetraObj(), temp2.epetraObj(), redC_ );
//     //redC_.print( Xyce::dout() );

//     // V' * B
//     MVT::MvTransMv( 1.0, xyceV.epetraObj(), BPtr_->epetraObj(), redB_ );
//     //redB_.print( Xyce::dout() );

//     bsuccess = bsuccess & sparsifyRedSystem_();
//   }
//
  // ---------------------------------------------------------------------
  // Output the projected system, redG_, redC_, and redB_.
  // ---------------------------------------------------------------------

  if (tiaParams.morSaveRedSys)
  {
    // If sparsification was used, write out the N_LAS_Matrix object instead of the dense object.
    if (tiaParams.morSparsificationType)
      outputMgrAdapterRCPtr_->outputROM ( *redGPtr_, *redCPtr_, redB_, redL_ );
    else
      outputMgrAdapterRCPtr_->outputROM ( redG_, redC_, redB_, redL_ );  // L = B'
  }

#else

  std::string msg=
   "Belos is necessary to compute a reduced-order model!  Please recompile Xyce with Belos enabled!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);

#endif // Xyce_BELOS

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::evalOrigTransferFunction()
// Purpose       : Evaluate original transfer function.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 5/29/2012
//-----------------------------------------------------------------------------
bool MOR::evalOrigTransferFunction()
{
  bool bsuccess = true;

  createOrigLinearSystem_();

  int currentStep = 0;
  int finalStep = morEvalSize_;

  bool stepAttemptStatus;

  while (currentStep < finalStep)
  {
    updateCurrentFreq_(currentStep);

    updateOrigLinearSystemFreq_();

    stepAttemptStatus = solveOrigLinearSystem_();

    currentStep++;

    if (stepAttemptStatus)
    {
      processSuccessfulStep(true);
    }
    else // stepAttemptStatus  (ie do this if the step FAILED)
    {
      processFailedStep();
    }

  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::evalRedTransferFunction()
// Purpose       : Evaluate reduced transfer function.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 5/29/2012
//-----------------------------------------------------------------------------
bool MOR::evalRedTransferFunction()
{
  bool bsuccess = true;

  createRedLinearSystem_();

  int currentStep = 0;
  int finalStep = morEvalSize_;

  bool stepAttemptStatus;

  while (currentStep < finalStep)
  {
    updateCurrentFreq_(currentStep);

    updateRedLinearSystemFreq_();

    stepAttemptStatus = solveRedLinearSystem_();

    currentStep++;

    if (stepAttemptStatus)
    {
      processSuccessfulStep(false);
    }
    else // stepAttemptStatus  (ie do this if the step FAILED)
    {
      processFailedStep();
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::createOrigLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/20/2011
//-----------------------------------------------------------------------------
bool MOR::createOrigLinearSystem_()
{
  bool bsuccess = true;

  RCP<N_PDS_Manager> pdsMgrPtr_;
  pdsMgrPtr_ = rcp(lasSystemRCPtr_->getPDSManager(), false);

  RCP<N_PDS_ParMap> BaseMap_, oBaseMap_;
  RCP<Epetra_CrsGraph> BaseFullGraph_;

  BaseMap_ = rcp(pdsMgrPtr_->getParallelMap( "SOLUTION" ), false);
  oBaseMap_ = rcp(pdsMgrPtr_->getParallelMap( "SOLUTION_OVERLAP_GND"), false);
  BaseFullGraph_ = rcp( pdsMgrPtr_->getMatrixGraph("JACOBIAN"), false );

  int numBlocks = 2;

  std::vector<RCP<N_PDS_ParMap> > blockMaps = createBlockParMaps(numBlocks, *BaseMap_, *oBaseMap_);

  // Create a block vector
  REFBPtr_ = rcp ( new N_LAS_BlockVector ( numBlocks, blockMaps[0], BaseMap_ ) );

  std::vector<std::vector<int> > blockPattern(2);
  blockPattern[0].resize(2);
  blockPattern[0][0] = 0; blockPattern[0][1] = 1;
  blockPattern[1].resize(2);
  blockPattern[1][0] = 0; blockPattern[1][1] = 1;

  int MaxGID = BaseMap_->maxGlobalEntity();
  int offset=1;
  while ( offset <= MaxGID ) offset *= 10;

  RCP<Epetra_CrsGraph> blockGraph = createBlockGraph( offset, blockPattern, *blockMaps[0], *BaseFullGraph_);

  sCpG_REFMatrixPtr_ = rcp ( new N_LAS_BlockMatrix( numBlocks, offset, blockPattern, *blockGraph, *BaseFullGraph_) );

  // Load diagonal blocks of real equivalent form: (G - s0*C)
  sCpG_REFMatrixPtr_->put( 0.0 ); // Zero out whole matrix

  // Add scaled C matrix first if expansion point is not zero.
  if (s0_ != 0.0)
  {
    sCpG_REFMatrixPtr_->block( 0, 0 ).add(*CPtr_);
    sCpG_REFMatrixPtr_->block( 0, 0 ).scale(-s0_);
    sCpG_REFMatrixPtr_->block( 1, 1 ).add(*CPtr_);
    sCpG_REFMatrixPtr_->block( 1, 1 ).scale(-s0_);
  }
  // Add G into diagonal block to obtain (G - s0*C)
  sCpG_REFMatrixPtr_->block( 0, 0 ).add(*GPtr_);
  sCpG_REFMatrixPtr_->block( 1, 1 ).add(*GPtr_);

  // Create solver factory
  Amesos amesosFactory;

  REFXPtr_ = rcp ( new N_LAS_BlockVector ( numBlocks, blockMaps[0], BaseMap_ ) );
  REFXPtr_->putScalar( 0.0 );

  blockProblem_ = rcp(new Epetra_LinearProblem(&sCpG_REFMatrixPtr_->epetraObj(), &REFXPtr_->epetraObj(), &REFBPtr_->epetraObj() ) );

  blockSolver_ = rcp( amesosFactory.Create( "Klu", *blockProblem_ ) );

  int linearStatus = blockSolver_->SymbolicFactorization();

  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos symbolic factorization exited with error: " << linearStatus;
    bsuccess = false;
  }

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : MOR::createRedLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/20/2011
//-----------------------------------------------------------------------------
bool MOR::createRedLinearSystem_()
{
  // Get the reduced system size.
  int k = redG_.numRows();

  // Resize serial dense matrix for real-equivalent form.
  sCpG_redMatrix_.shape(2*k, 2*k);
  sCpG_tmpMatrix_.shape(2*k, 2*k);
  ref_redB_.shape(2*k, numPorts_);

  // First, load the diagonals.
  // Get a view of the primary diagonal block, insert (G - s0*C).
  Teuchos::SerialDenseMatrix<int, double> subMtx( Teuchos::View, sCpG_tmpMatrix_, k, k, 0, 0 );
  subMtx.assign( redC_ );
  subMtx.scale( -s0_ );
  subMtx += redG_;

  // Get a view of the lower k x k diagonal block, insert (G - s0*C).
  Teuchos::SerialDenseMatrix<int, double> subMtx2( Teuchos::View, sCpG_tmpMatrix_, k, k, k, k );
  subMtx2.assign( redC_ );
  subMtx2.scale( -s0_ );
  subMtx2 += redG_;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : MOR::updateOrigLinearSystemFreq_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/20/2011
//-----------------------------------------------------------------------------

bool MOR::updateOrigLinearSystemFreq_()
{
  double omega =  2.0 * M_PI * currentFreq_;

  sCpG_REFMatrixPtr_->block( 0, 1).put( 0.0);
  sCpG_REFMatrixPtr_->block( 0, 1).add(*CPtr_);
  sCpG_REFMatrixPtr_->block( 0, 1).scale(-omega);

  sCpG_REFMatrixPtr_->block(1, 0).put( 0.0);
  sCpG_REFMatrixPtr_->block(1, 0).add(*CPtr_);
  sCpG_REFMatrixPtr_->block(1, 0).scale(omega);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : MOR::updateRedLinearSystemFreq_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/20/2011
//-----------------------------------------------------------------------------

bool MOR::updateRedLinearSystemFreq_()
{
  int k = redG_.numRows();

  double omega =  2.0 * M_PI * currentFreq_;

  // Reset reduced real equivalent form matrix to sCpG_tmpMatrix_
  sCpG_redMatrix_.assign( sCpG_tmpMatrix_ );

  // Insert off diagonals.
  Teuchos::SerialDenseMatrix<int, double> subMtx( Teuchos::View, sCpG_redMatrix_, k, k, 0, k );
  subMtx.assign( redC_ );
  subMtx.scale( -omega );

  Teuchos::SerialDenseMatrix<int, double> subMtx2( Teuchos::View, sCpG_redMatrix_, k, k, k, 0 );
  subMtx2.assign( redC_ );
  subMtx2.scale( omega );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : MOR::solveLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :  Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/2011
//-----------------------------------------------------------------------------

bool MOR::solveOrigLinearSystem_()
{

  bool bsuccess = true;

  // Solve the block problem
  int linearStatus = blockSolver_->NumericFactorization();

  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos numeric factorization exited with error: " << linearStatus;
    bsuccess = false;
  }

  // Loop over number of I/O ports here
  for (unsigned int j=0; j < bMatEntriesVec_.size(); ++j)
  {
    REFBPtr_->putScalar( 0.0 );
    (REFBPtr_->block( 0 ))[bMatEntriesVec_[j]] = -1.0;

    linearStatus = blockSolver_->Solve();
    if (linearStatus != 0)
    {
      Xyce::dout() << "Amesos solve exited with error: " << linearStatus;
      bsuccess = false;
    }

    // Compute transfer function entries for all I/O
    for (unsigned int i=0; i < bMatEntriesVec_.size(); ++i)
    {
       // Populate H for all ports in L
       // L is the same as B and also a set of canonical basis vectors (e_i), so
       // we can pick off the appropriate entries of REFXPtr to place into H.
       origH_(i,j) = std::complex<double>(-(REFXPtr_->block( 0 ))[bMatEntriesVec_[i]],
                                          -(REFXPtr_->block( 1 ))[bMatEntriesVec_[i]]);
    }
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::solveRedLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :  Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/2011
//-----------------------------------------------------------------------------

bool MOR::solveRedLinearSystem_()
{
  bool bsuccess = true;

  int k = redG_.numRows();
  int numPorts = bMatEntriesVec_.size();
  Teuchos::LAPACK<int, double> lapack;

  ref_redB_.putScalar( 0.0 );
  Teuchos::SerialDenseMatrix<int, double> tmpRedB( Teuchos::View, ref_redB_, redB_.numRows(), redB_.numCols(), 0, 0 );
  tmpRedB.assign( redB_ );

  // First factor matrix using LU.
  int info = 0;
  std::vector<int> ipiv( sCpG_redMatrix_.numRows() );
  lapack.GETRF( sCpG_redMatrix_.numRows(), sCpG_redMatrix_.numCols(), sCpG_redMatrix_.values(),
                sCpG_redMatrix_.stride(), &ipiv[0], &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRF: LU factorization failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Next solve for all ports at once using LU factors.
  const char trans = 'N';
  lapack.GETRS( trans, sCpG_redMatrix_.numRows(), ref_redB_.numCols(), sCpG_redMatrix_.values(),
                sCpG_redMatrix_.stride(), &ipiv[0], ref_redB_.values(), ref_redB_.stride(), &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRS: LU solve failed with error: " << info << std::endl;
    bsuccess = false;
  }

  Teuchos::SerialDenseMatrix<int, double> tmpRedImag( numPorts, numPorts ), tmpRedReal( numPorts, numPorts );

  //  Compute real part.
  tmpRedReal.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, redB_, tmpRedB, 0.0 );

  // Compute imaginary part.
  Teuchos::SerialDenseMatrix<int, double> tmpRedB2( Teuchos::View, ref_redB_, redB_.numRows(), redB_.numCols(), k, 0 );
  tmpRedImag.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, redB_, tmpRedB2, 0.0 );

  // Compute transfer function entries for all I/O
  for (unsigned int j=0; j < bMatEntriesVec_.size(); ++j)
  {
    for (unsigned int i=0; i < bMatEntriesVec_.size(); ++i)
    {
      redH_(i,j) = std::complex<double>(tmpRedReal(i,j), tmpRedImag(i,j));
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::sparsifyRedSystem_()
// Purpose       : Use techniques to sparsify the projected system.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 8/20/2012
//-----------------------------------------------------------------------------
bool MOR::sparsifyRedSystem_()
{
  bool bsuccess = true;

  int k = redB_.numRows();

  // Storage vectors
  Teuchos::SerialDenseMatrix<int, double> R(redB_);
  Teuchos::SerialDenseMatrix<int, double> A(redC_);

  // Factor Ghat.
  Teuchos::LAPACK<int, double> lapack;
  Teuchos::SerialDenseMatrix<int, double> tmp_Ghat( redG_ );

  int info = 0;
  std::vector<int> ipiv( tmp_Ghat.numRows() );
  lapack.GETRF( tmp_Ghat.numRows(), tmp_Ghat.numCols(), tmp_Ghat.values(),
                tmp_Ghat.stride(), &ipiv[0], &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRF: LU factorization of Ghat failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Compute reciprocal condition estimate of Ghat.
  const char norm = '1';
  double condGhat = 0.0;
  std::vector<double> condWork( 4*k );
  std::vector<int> condIWork( k );
  lapack.GECON( norm, tmp_Ghat.numRows(), tmp_Ghat.values(), tmp_Ghat.stride(), tmp_Ghat.normOne(), &condGhat, &condWork[0], &condIWork[0], &info );
  // Xyce::dout() << "Condition estimate for Ghat is : " << 1.0/ condGhat << std::endl;
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GECON: Computing condition estimate of Ghat failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Compute A = inv(Ghat)* Chat
  const char trans = 'N';
  lapack.GETRS( trans, tmp_Ghat.numRows(), A.numCols(), tmp_Ghat.values(),
                tmp_Ghat.stride(), &ipiv[0], A.values(), A.stride(), &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRS: LU solve failed with error: " << info << std::endl;
    bsuccess = false;
  }
  //Xyce::dout() << "A = inv(Ghat)*Chat" << std::endl;
  //A.print(Xyce::dout());

  // Compute R = inv(Ghat)* Bhat
  lapack.GETRS( trans, tmp_Ghat.numRows(), R.numCols(), tmp_Ghat.values(),
                tmp_Ghat.stride(), &ipiv[0], R.values(), R.stride(), &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRS: LU solve failed with error: " << info << std::endl;
    bsuccess = false;
  }
  //Xyce::dout() << "R = inv(Ghat)*Bhat" << std::endl;
  //R.print(Xyce::dout());

  // Reduce A to Schur form and compute the eigenvalues.
  // Allocate the work space. This space will be used below for calls to:
  // * GEES  (requires 3*k for real)
  // * TREVC (requires 3*k for real)
  // Furthermore, GEES requires a real array of length k
  //
  int lwork = 8*k;
  std::vector<double> work( lwork );
  std::vector<double> revals( k );
  std::vector<double> ievals( k );
  Teuchos::SerialDenseMatrix<int, double> Q(k, k);

  // Create diagonal tmpB for the eigenproblem (A, tmpB), beta should return as a vector of ones.
  const int ldvl = 1;
  double vl[ ldvl ];
  std::vector<double> beta( k );
  Teuchos::SerialDenseMatrix<int, double> tmpB( k, k );
  for (int i=0; i<k; i++)
    tmpB(i,i) = 1.0;
  lapack.GGEV( 'N', 'V', k, A.values(), A.stride(), tmpB.values(), tmpB.stride(),
               &revals[0], &ievals[0], &beta[0], &vl[0], ldvl, Q.values(), Q.stride(),
               &work[0], lwork, &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GGEV: Computing eigenvectors and values of A = inv(Ghat)*Chat failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Scale the eigenvectors returned from GGEV to have 2-norm = 1
  // They are initially scaled to have inf-norm = 1 from LAPACK
  Teuchos::BLAS<int, double> blas;
  std::vector<int> rowIndex( k );
  int i = 0;

  while( i < k ) {
    if ( ievals[i] != 0.0 )
    {
      rowIndex[i] = 1;
      rowIndex[i+1] = -1;

      // Compute the 2-norm of the complex eigenvector
      double norm_r = blas.NRM2( k, Q[i], 1 );
      double norm_i = blas.NRM2( k, Q[i+1], 1 );
      double norm = lapack.LAPY2( norm_r, norm_i );

      // Scale the complex eigenvector
      for (int j=0; j<k; j++)
      {
        Q(j,i)   /= norm;
        Q(j,i+1) /= norm;
      }
      i += 2;
    }
    else
    {
      rowIndex[i] = 0;

      // Compute the 2-norm of the i-th column
      double norm = blas.NRM2( k, Q[i], 1 );

      // Scale the real eigenvector
      for (int j=0; j<k; j++)
      {
        Q(j,i) /= norm;
      }
      i++;
    }
  }

  // Xyce::dout() << "Eigenvalues of A: " << std::endl;
  // for (int i=0; i<k; ++i)
  //   Xyce::dout() << revals[i] << "\t" << ievals[i] << "\t" << beta[i] << std::endl;

  // Xyce::dout() << "Eigenvectors of A: " << std::endl;
  // Q.print(Xyce::dout());

  // Construct complex eigenvectors
  Teuchos::SerialDenseMatrix<int, std::complex<double> > V(k, k);
  for (int i=0; i<k; i++)
  {
    if (rowIndex[i] == 1)
    {
      for (int j=0; j<k; j++)
      {
        // Insert both conjugate pairs simultaneously.
        V(j,i) = std::complex<double>( Q(j,i), Q(j,i+1) );
        V(j,i+1) = std::complex<double>( Q(j,i), -Q(j,i+1) );
      }
      i++;  // Need to increment an extra value for the conjugate pair.
    }
    else  // The eigenvector is real, copy directly
    {
      for (int j=0; j<k; j++)
      {
        V(j,i) = std::complex<double>( Q(j,i), 0.0 );
      }
    }
  }

  // Factor V
  Teuchos::LAPACK<int, std::complex<double> > lapack_complex;
  lapack_complex.GETRF( V.numRows(), V.numCols(), V.values(), V.stride(), &ipiv[0], &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRF: LU factorization of eigenvectors failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Compute condition estimate of V.
  double condV = 0.0;
  std::vector<std::complex<double> > condCWork( 2*k );
  lapack_complex.GECON( norm, V.numRows(), V.values(), V.stride(), V.normOne(), &condV, &condCWork[0], &condWork[0], &info );
  //Xyce::dout() << "Condition estimate for V is : " << 1.0/condV << std::endl;
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GECON: Computing condition estimate of V failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Compute inv(V)
  std::vector<std::complex<double> > cwork( lwork );
  lapack_complex.GETRI( V.numRows(), V.values(), V.stride(), &ipiv[0], &cwork[0], lwork, &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRI: Computing inverse of V failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Convert R to a complex vector to use the multiply routine for inv(V).
  Teuchos::SerialDenseMatrix<int, std::complex<double> > tmpR( k, numPorts_ );
  for (int j=0; j<numPorts_; j++)
    for (int i=0; i<k; i++)
      tmpR(i,j) = std::complex<double>( R(i,j), 0.0 );

  // Compute inv(V) * R
  Teuchos::SerialDenseMatrix<int, std::complex<double> > tmp_redB(k ,numPorts_);
  tmp_redB.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, V, tmpR, 0.0);

  // Generate compressed real storage of V*R and inv(V)
  redL_ = redB_;  // Store copy of redB_ before destroying it below.
  Teuchos::SerialDenseMatrix<int, double> invV(k, k);
  for (int i=0; i<k; i++)
  {
    if (rowIndex[i] == 1)  // Complex conjugate pair.
    {
      for (int j=0; j<k; j++)
      {
        invV(i, j) = V(i,j).real();
        invV(i+1, j) = V(i,j).imag();
      }
      for (int j=0; j<numPorts_; j++)
      {
        redB_(i, j) = tmp_redB(i, j).real();
        redB_(i+1, j) = tmp_redB(i, j).imag();
      }
      i++;  // Increment by one to skip complex conjugate.
    }
    else  // Eigenvalue is real, so is eigenvector.
    {
      for (int j=0; j<k; j++)
        invV(i, j) = V(i,j).real();
      for (int j=0; j<numPorts_; j++)
        redB_(i, j) = tmp_redB(i, j).real();
    }
  }

  // Factor invV.
  lapack.GETRF( invV.numRows(), invV.numCols(), invV.values(), invV.stride(), &ipiv[0], &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRF: LU factorization of inv(V) failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Compute A = inv(invV)^T*Lhat
  const char trans2 = 'T';
  lapack.GETRS( trans2, invV.numRows(), redL_.numCols(), invV.values(),
                invV.stride(), &ipiv[0], redL_.values(), redL_.stride(), &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRS: LU solve for Lhat failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Xyce::dout() << "Bhat : " << std::endl;
  // redB_.print(Xyce::dout());
  // Xyce::dout() << "Lhat : " << std::endl;
  // redL_.print(Xyce::dout());

  // Create redCPtr_ and redGPtr_
  // Generate Epetra_Map that puts all the values of the reduced system on one processor.
  RCP<N_PDS_Manager> pdsMgrPtr_ = Teuchos::rcp(lasSystemRCPtr_->getPDSManager(), false);
  RCP<N_PDS_ParMap> BaseMap_ = Teuchos::rcp(pdsMgrPtr_->getParallelMap( "SOLUTION" ), false);
  RCP<Epetra_Map> serialRedMap;
  if (BaseMap_->petraMap()->Comm().MyPID() == 0)
    serialRedMap = Teuchos::rcp( new Epetra_Map( k, k, 0, BaseMap_->petraMap()->Comm() ) );
  else
    serialRedMap = Teuchos::rcp( new Epetra_Map( k, 0, 0, BaseMap_->petraMap()->Comm() ) );

  Epetra_CrsMatrix* tmpRedG = new Epetra_CrsMatrix( Copy, *serialRedMap, 1 );
  Epetra_CrsMatrix* tmpRedC = new Epetra_CrsMatrix( Copy, *serialRedMap, 2 );

  // Let processor 0 insert entries into tmpRedG and tmpRedC
  if (BaseMap_->petraMap()->Comm().MyPID() == 0)
  {
    std::vector<int> index(2);
    std::vector<double> val(2);
    for (int i=0; i<k; i++)
    {
      // Insert G, which is a diagonal of ones.
      index[0] = i; val[0] = 1.0;
      tmpRedG->InsertGlobalValues( i, 1, &val[0], &index[0] );

      // Insert C, which is block diagonal, where the blocks represent conjugate eigenvalues.
      if (rowIndex[i] == 1)
      {
        index[0] = i; index[1] = i+1;
        val[0] = revals[i]; val[1] = -ievals[i];
        tmpRedC->InsertGlobalValues( i, 2, &val[0], &index[0] );
      }
      else if (rowIndex[i] == -1)
      {
        index[0] = i-1; index[1] = i;
        val[0] = ievals[i-1]; val[1] = revals[i-1];
        tmpRedC->InsertGlobalValues( i, 2, &val[0], &index[0] );
      }
      else // Must be real.
      {
        index[0] = i; val[0] = revals[i];
        tmpRedC->InsertGlobalValues( i, 1, &val[0], &index[0] );
      }
    }
  }
  // We are done inserting values.
  tmpRedC->FillComplete();
  tmpRedG->FillComplete();

  // Pass tmpRedG and tmpRedC into redGPtr_ and redGPtr_, respectively, which will manage the memory.
  redGPtr_ = rcp( new N_LAS_Matrix( tmpRedG ) );
  redCPtr_ = rcp( new N_LAS_Matrix( tmpRedC ) );

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::processSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool MOR::processSuccessfulStep(bool origSys)
{
  bool bsuccess = true;

  // Output H.
  if (origSys)
  {
    outputMgrAdapterRCPtr_->outputMORTF (origSys, currentFreq_, origH_);
  }
  else
  {
    outputMgrAdapterRCPtr_->outputMORTF (origSys, currentFreq_, redH_);
  }

  if ( !firstDoubleDCOPStep_() )
  {
      stepNumber += 1;
      totalNumberSuccessStepsThisParameter_ += 1;
      totalNumberSuccessfulStepsTaken_ += 1;
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : MOR::processFailedStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool MOR::processFailedStep()
{
  bool bsuccess = true;

  stepNumber += 1;
  morEvalFailures_.push_back(stepNumber);
  totalNumberFailedStepsAttempted_  += 1;
  secRCPtr_->numberSuccessiveFailures += 1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::finish
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool MOR::finish()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_ANALYSIS
  Xyce::dout() << "Calling MOR::finish() outputs!" << std::endl;
#endif

  if (!(morEvalFailures_.empty()))
  {
    bsuccess = false;
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : MOR::handlePredictor
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/24/2013
//-----------------------------------------------------------------------------
bool MOR::handlePredictor()
{
  anaManagerRCPtr_->getTIADataStore()->setErrorWtVector();
  wimRCPtr_->obtainPredictor();
  wimRCPtr_->obtainPredictorDeriv();

  // In case this is the upper level of a 2-level sim, tell the
  // inner solve to do its prediction:
  loaderRCPtr_->startTimeStep ();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : MOR::updateCurrentFreq_(int stepNumber )
// Purpose       :
// Special Notes : Used for MOR analysis classes.
// Scope         : public
// Creator       : Ting Mei, SNL.
// Creation Date :
//-----------------------------------------------------------------------------
bool MOR::updateCurrentFreq_(int stepNumber)
{
  double fStart = tiaParams.morCompFStart;

  if (tiaParams.type=="LIN")
  {
    currentFreq_  = fStart + static_cast<double>(stepNumber)*fStep_;
  }
  else if(tiaParams.type=="DEC" || tiaParams.type=="OCT")
  {
    currentFreq_ = fStart*pow(stepMult_, static_cast<double>(stepNumber) );
  }
  else
  {
    std::string msg=
      "MOR::updateCurrentFreq_: unsupported STEP type";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : MOR::setupSweepParam_
// Purpose       : Processes sweep parameters.
// Special Notes : Used for MOR analysis classes.
// Scope         : public
// Creator       : Ting Mei, SNL.
// Creation Date :
//-----------------------------------------------------------------------------
int MOR::setupSweepParam_()
{
  int fCount = 0;
  double fStart = tiaParams.morCompFStart;
  double fStop = tiaParams.morCompFStop;

#ifdef Xyce_DEBUG_ANALYSIS
  if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
  {
    Xyce::dout() << std::endl << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "MOR::setupSweepParam_" << std::endl;
  }
#endif

  if (tiaParams.morCompType=="LIN")
  {
    if ( tiaParams.morCompNP == 1)
      fStep_ = 0;
    else
      fStep_  = (fStop - fStart)/(tiaParams.morCompNP - 1);

    fCount = tiaParams.morCompNP;
#ifdef Xyce_DEBUG_ANALYSIS
    if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
    {
      Xyce::dout() << "fStep   = " << fStep_  << std::endl;
    }
#endif
  }
  else if(tiaParams.morCompType=="DEC")
  {
    stepMult_ = pow(10,(1/(double)tiaParams.morCompNP));
    fCount   = (int)(floor(fabs(log10(fStart) - log10(fStop)) * tiaParams.morCompNP + 1));
#ifdef Xyce_DEBUG_ANALYSIS
    if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
    {
      Xyce::dout() << "stepMult_ = " << stepMult_  << std::endl;
    }
#endif
  }
  else if(tiaParams.morCompType=="OCT")
  {
    stepMult_ = pow(2,1/(double)(tiaParams.morCompNP));

    // changed to remove dependence on "log2" function, which apparently
    // doesn't exist in the math libraries of FreeBSD or the mingw
    // cross-compilation suite.   Log_2(x)=log_e(x)/log_e(2.0)
    double ln2=log(2.0);
    fCount   = floor(fabs(log(fStart) - log(fStop))/ln2 * tiaParams.morCompNP + 1);
#ifdef Xyce_DEBUG_ANALYSIS
    if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
    {
      Xyce::dout() << "stepMult_ = " << stepMult_  << std::endl;
    }
#endif
  }
  else
  {
    std::string msg=
      "MOR::setupSweepParam: unsupported type";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

  // At this point, pinterval equals the total number of steps
  // for the step loop.
  return fCount;
}

} // namespace Analysis
} // namespace Xyce
