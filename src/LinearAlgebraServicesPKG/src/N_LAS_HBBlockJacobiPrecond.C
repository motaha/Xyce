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
// Filename       : $RCSfile: N_LAS_HBBlockJacobiPrecond.C,v $
//
// Purpose        : Implementation file for the Iterative linear solver
//                  interface.
//
// Special Notes  :
//
// Creator        : Heidi K. Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 11/11/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:45 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

// ---------- Standard Includes ----------

#include <Xyce_config.h>

#include <sstream>

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>

#include <N_LAS_HBBlockJacobiPrecond.h>
#include <N_LAS_HBBlockJacobiEpetraOperator.h>	
#include <N_LAS_HBBuilder.h>
#include <N_LOA_HBLoader.h>
#include <N_MPDE_State.h>

#include <N_LAS_MultiVector.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_Matrix.h>

#include <N_LAS_Problem.h>
#include <N_LAS_Builder.h>
#include <N_LAS_System.h>
#include <N_LOA_Loader.h>

#include <N_UTL_Timer.h>

#include <N_ERH_ErrorMgr.h>
	
#include <N_PDS_SerialComm.h>

#include <Teuchos_RCP.hpp>
#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Amesos.h>

using Teuchos::RCP;
using Teuchos::rcp;

// static class member initializations
// Default preconditioner values
/*
const bool   N_LAS_HBBlockJacobiPrecond::useFactory_default_ = false;
const int    N_LAS_HBBlockJacobiPrecond::overlap_default_    = 0;
const double N_LAS_HBBlockJacobiPrecond::dropTol_default_    = 1.0e-03;
const double N_LAS_HBBlockJacobiPrecond::ilutFill_default_   = 2.0;
const double N_LAS_HBBlockJacobiPrecond::rThresh_default_    = 1.0001;
const double N_LAS_HBBlockJacobiPrecond::aThresh_default_    = 0.0001;
*/

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiPrecond::N_LAS_HBBlockJacobiPrecond
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
N_LAS_HBBlockJacobiPrecond::N_LAS_HBBlockJacobiPrecond()
  : N_LAS_Preconditioner()
{
  setDefaultOptions();
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiPrecond::setDefaultOptions
// Purpose       : resets Ifpack options
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool N_LAS_HBBlockJacobiPrecond::setDefaultOptions()
{
  // Set defaults
/*
  useFactory_ = useFactory_default_;
  dropTol_ = dropTol_default_;
  ilutFill_ = ilutFill_default_;
  rThresh_ = rThresh_default_;
  aThresh_ = aThresh_default_;
  overlap_ = overlap_default_;
*/

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiPrecond::setDefaultOption
// Purpose       : resets Ifpack option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool N_LAS_HBBlockJacobiPrecond::setDefaultOption( const string & option )
{
/*
  if( option == "AZ_athresh" )         aThresh_ = aThresh_default_;
  if( option == "AZ_rthresh" )         rThresh_ = rThresh_default_;
  if( option == "AZ_ilut_fill" )       ilutFill_ = ilutFill_default_;
  if( option == "AZ_drop" )            dropTol_ = dropTol_default_;
  if( option == "AZ_overlap" )         overlap_ = overlap_default_;
  if( option == "use_ifpack_factory" ) useFactory_ = useFactory_default_;
*/

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiPrecond::setOptions
// Purpose       : sets Ifpack options and params from modelblock
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool N_LAS_HBBlockJacobiPrecond::setOptions( const N_UTL_OptionBlock & OB )
{
  // Set the parameters from the list
  list<N_UTL_Param>::const_iterator it_tpL = OB.getParams().begin();
  list<N_UTL_Param>::const_iterator end_tpL = OB.getParams().end();
  for (; it_tpL != end_tpL; ++it_tpL) 
    {
      this->setParam( *it_tpL );
    }
  
  // store for restart of solver_
  if( &OB != options_.get() )
    {
      options_ = Teuchos::rcp( &OB, false );
    }
  
  return STATUS_SUCCESS;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiPrecond::setParam
// Purpose       : sets Ifpack option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool N_LAS_HBBlockJacobiPrecond::setParam( const N_UTL_Param & param )
{
  string tag = param.tag();
  string uTag = param.uTag();

  // Set our copies of these parameters that get passed to the solver in the
  // "iterate" command
/*
  if( tag == "AZ_overlap" )
    overlap_ = param.iVal();
  else if( tag == "AZ_athresh")
    aThresh_ = param.dVal();
  else if( tag == "AZ_rthresh")
    rThresh_ = param.dVal();
  else if( tag == "AZ_drop")
    dropTol_ = param.dVal();
  else if( tag == "AZ_ilut_fill")
    ilutFill_ = param.dVal();
  else if( tag == "use_ifpack_factory")
    useFactory_ = param.iVal();
  else
    return false;
*/

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiPrecond::initGraph
// Purpose       : Initialize the graph using information from the N_LAS_System
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 12/11/08
//-----------------------------------------------------------------------------
bool N_LAS_HBBlockJacobiPrecond::initGraph( const Teuchos::RCP<N_LAS_Problem> & problem )
{
  // Generate the graph of each real equivalent form and then generate 
  // empty linear systems for each frequency.
  RCP<N_LAS_Matrix> appdQdx = rcp( appBuilderPtr_->createMatrix() );
  RCP<N_LAS_Matrix> appdFdx = rcp( appBuilderPtr_->createMatrix() );

  // Relate the conductance and inductance matrices:
  // G(t) = df/dx(x(t))
  // C(t) = dq/dx(x(t))
  // Compute:  (omega*i*j*C_bar + G_bar)^{-1} for i=0,1,...,M,-M,...,-1
  //           using real equivalent form K_1 = [G_bar -i*omega*C_bar; i*omega*C_bar G_bar]
  
  // Generate the new real equivalent graph
  RCP<N_PDS_Comm> pdsCommRCPtr = rcp(new N_PDS_SerialComm);
  int origRows = appdQdx->getLocalNumRows();
  int refRows = 2*origRows;
  epetraMap_ = rcp(new Epetra_Map( refRows, refRows, 0, *(pdsCommRCPtr->petraComm())) );

  // Count up the number of nonzero entries for the 2x2 block matrix.
  std::vector<int> refNNZs(refRows);
  maxRefNNZs_ = 0;
  for ( int i=0; i<origRows; ++i ) {
    refNNZs[i] = appdQdx->getLocalRowLength(i) + appdFdx->getLocalRowLength(i);
    refNNZs[origRows+i] = refNNZs[i];
    if (refNNZs[i] > maxRefNNZs_) maxRefNNZs_ = refNNZs[i];
  }
  epetraGraph_ = rcp(new Epetra_CrsGraph( Copy, *epetraMap_, &refNNZs[0], true ));

  // Put together the indices for each row and insert them into the graph.
  int tmpNNZs=0, tmpNNZs2=0;
  std::vector<double> tmpCoeffs(maxRefNNZs_);
  std::vector<int> refIdxs(maxRefNNZs_), refIdxs2(maxRefNNZs_);

  for ( int i=0; i<origRows; ++i ) {

    // Get the indices for the first block of the matrix (G_bar)
    appdFdx->getRowCopy( i, maxRefNNZs_, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );

    // Get the indices for the third block of the matrix (C_bar)
    appdQdx->getRowCopy( i, maxRefNNZs_, tmpNNZs2, &tmpCoeffs[0], &refIdxs2[0] );
    
    // Insert the indices for the third block into refIdxs, as they are the indices of the second block
    for (int j=0; j<tmpNNZs2; ++j) {
      refIdxs[tmpNNZs+j] = refIdxs2[j]+origRows;
    }
    epetraGraph_->InsertGlobalIndices( i, refNNZs[i], &refIdxs[0] );
  
    // Insert the indices for the first block into refIdxs2, as they are the indices of the fourth block
    for (int j=0; j<tmpNNZs; ++j) {
      refIdxs2[tmpNNZs2+j] = refIdxs[j]+origRows;
    }
    epetraGraph_->InsertGlobalIndices( origRows+i, refNNZs[origRows+i], &refIdxs2[0] );
  }  
  epetraGraph_->FillComplete();

  // Get the Fourier series information and generate the Epetra_LinearSystems.
  RCP<N_LAS_BlockVector> bXt = hbBuilderPtr_->createTimeDomainBlockVector();
  N_ = bXt->blockCount();
  M_ = (int)((N_-1)/2);
 
  // Generate the vectors for the N_ linear problems to be solved on the diagonal.
  epetraMatrix_.resize(N_);
  epetraRHS_.resize(N_);
  epetraSoln_.resize(N_);
  epetraProblem_.resize(N_);
  amesosPtr_.resize(N_); 

  Amesos amesosFactory;  

  for (int i=0; i<N_; ++i) {
    epetraMatrix_[i] = rcp( new Epetra_CrsMatrix( Copy, *epetraGraph_ ) );
    epetraMatrix_[i]->FillComplete();
    epetraRHS_[i] = rcp( new Epetra_MultiVector( *epetraMap_, 1 ) );
    epetraSoln_[i] = rcp( new Epetra_MultiVector( *epetraMap_, 1 ) );
    epetraProblem_[i] = rcp( new Epetra_LinearProblem( &*epetraMatrix_[i], &*epetraRHS_[i], &*epetraSoln_[i] ) ); 
    amesosPtr_[i] = rcp( amesosFactory.Create( "Klu", *epetraProblem_[i] ) );
    amesosPtr_[i]->SymbolicFactorization();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiPrecond::compute
// Purpose       : Initialize the values of the Epetra_LinearSystem's
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 12/11/08
//-----------------------------------------------------------------------------
bool N_LAS_HBBlockJacobiPrecond::initValues( const Teuchos::RCP<N_LAS_Problem> & problem )
{
  // Get the next solution vector in frequency domain and state vector
  N_LAS_BlockVector & bXf = *dynamic_cast<N_LAS_BlockVector*>(lasSysPtr_->getNextSolVector());
  N_LAS_BlockVector & bS = *dynamic_cast<N_LAS_BlockVector*>(lasSysPtr_->getNextStaVector());  
  N_LAS_BlockVector & bdSdt = *dynamic_cast<N_LAS_BlockVector*>(lasSysPtr_->getNextStaDerivVector());  

  // Permute the solution vector into the time domain
  RCP<N_LAS_BlockVector> bXt = hbBuilderPtr_->createTimeDomainBlockVector();
  hbLoaderPtr_->permutedIFT(bXf, &*bXt);

  // Get vectors and matrices for each block 
  RCP<N_LAS_Vector> appVec = rcp( appBuilderPtr_->createVector() );
  RCP<N_LAS_Vector> appStateVec = rcp( appBuilderPtr_->createStateVector() );
  RCP<N_LAS_Vector> appdSdt = rcp( appBuilderPtr_->createStateVector() );

  RCP<N_LAS_Matrix> appdQdx = rcp( appBuilderPtr_->createMatrix() );
  RCP<N_LAS_Matrix> appdFdx = rcp( appBuilderPtr_->createMatrix() );

  // Get matrix to contain sum for dQdx and dFdx
  RCP<N_LAS_Matrix> appdQdxSum = rcp( appBuilderPtr_->createMatrix() );
  RCP<N_LAS_Matrix> appdFdxSum = rcp( appBuilderPtr_->createMatrix() );

  appdQdxSum->put(0.0);
  appdFdxSum->put(0.0);

  // Sum up dQdx and dFdx
  for( int i = 0; i < N_; ++i )
  {
    // Set the fast time.
    state_->fastTime = times_[i];
    devInterfacePtr_->setFastTime( times_[i] );
    appLoaderPtr_->updateSources();

    *appVec = bXt->block(i);
    *appStateVec = bS.block(i);
    *appdSdt = bdSdt.block(i);

    appdQdx->put(0.0);
    appdFdx->put(0.0);

    // Compute dQdx and dFdx for this time point.
    appLoaderPtr_->loadDAEMatrices( &*appVec, &*appStateVec, &*appdSdt, &*appStateVec, &*appdQdx,  &*appdFdx);

    // Add into dQdxSum and dFdxSum
    appdQdxSum->add(*appdQdx);
    appdFdxSum->add(*appdFdx);
  }

  // Average the matrix sum
  appdQdxSum->scale( 1.0/N_ );
  appdFdxSum->scale( 1.0/N_ );  

  // Compute omega
  int timesSize = times_.size();
  double period = times_[timesSize - 1];
  double omega = 2.0 * M_PI/ period;

  // Get the values for each row of appdFdxSum/appdQdxSum and insert them into the matrix.
  int origRows = appdQdx->getLocalNumRows();
  int tmpNNZs=0, tmpNNZs2=0;
  double cplxCoeff = 0.0;
  std::vector<double> tmpCoeffs(maxRefNNZs_);
  std::vector<int> refIdxs(maxRefNNZs_), refIdxs2(maxRefNNZs_);

  for ( int nB=0; nB<N_; ++nB ) {

    // Compute the coefficient on the complex matrix [0,1,...,M,-M,...,-1]
    if (nB <= M_)
      cplxCoeff = nB*omega;
    else 
      cplxCoeff = (-M_+(nB-M_-1))*omega;

    // Insert the entries for all four blocks of the real-equivalent form
    for ( int i=0; i<origRows; ++i ) {
      
      // Load [G_bar -i*omega*C_bar; i*omega*C_bar G_bar]

      // Get the indices for the first block of the matrix (G_bar)
      appdFdxSum->getRowCopy( i, maxRefNNZs_, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );
      
      epetraMatrix_[nB]->ReplaceMyValues( i, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );
      
      // Modify the indices for the fourth block
      for ( int j=0; j<tmpNNZs; ++j ) {
        refIdxs[j] += origRows;
      }
      epetraMatrix_[nB]->ReplaceMyValues( origRows+i, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );

      // Add in complex part of REF matrix if non-zero.
      if (cplxCoeff != 0.0) {

        // Get the indices for the first block of the matrix (G_bar)
        appdQdxSum->getRowCopy( i, maxRefNNZs_, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );
        
        // Insert the second block
        for (int ii=0; ii<tmpNNZs; ++ii) { tmpCoeffs[ii] *= cplxCoeff; }
        epetraMatrix_[nB]->ReplaceMyValues( origRows+i, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );

        // Insert the third block
        for (int ii=0; ii<tmpNNZs; ++ii) { 
          tmpCoeffs[ii] *= -1.0; 
          refIdxs[ii] += origRows;
        }
        epetraMatrix_[nB]->ReplaceMyValues( i, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );
      }
    }
  }    

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiPrecond::compute
// Purpose       : Compute a preconditioner M such that M ~= A^{-1}.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool N_LAS_HBBlockJacobiPrecond::compute()
{
  bool precStatus = true;

  // Compute numeric factorization for each block.
  for ( int i=0; i<N_; ++i ) {
    amesosPtr_[i]->NumericFactorization();
  }

  if ( Teuchos::is_null( epetraPrec_ ) )
    epetraPrec_ = blockJacobiOperator( epetraProblem_, amesosPtr_, hbBuilderPtr_ );
  
  if ( Teuchos::is_null( epetraPrec_ ) )
    return false;

/*
  if (useFactory_) 
  {
    // Build the preconditioner using the values in problem_; numeric factorization.
    IFPACK_CHK_ERR(ifpackPrecond_->Compute());
  }
  else 
  {
    Epetra_CrsMatrix * epetraA = dynamic_cast<Epetra_CrsMatrix*>(problem_->epetraObj().GetMatrix());
    bool transpose = epetraA->UseTranspose();
    
    int factErr = rILUK_->Factor();
    if (factErr < 0) 
      return false;
    
#ifdef Xyce_VERBOSE_LINEAR
    double condest;
    rILUK_->Condest( transpose, condest );
#endif
    
    // Define label for printing out during the solve phase
    ostringstream ost;
    ost << "Ifpack_CrsRiluk Preconditioner: LevelFill = " << ilutFill_ << endl << 
           "                                Overlap = " << overlap_ << endl << 
           "                                Athresh = " << aThresh_ << endl <<
           "                                Rthresh = " << rThresh_ << endl <<
#ifdef Xyce_VERBOSE_LINEAR
           "                                CondEst = " << condest  << endl << 
#endif
           "                                ErrCode = " << factErr  << endl;
    string label = ost.str();
    rILUK_->SetLabel(label.c_str());
    
    rILUK_->SetUseTranspose( transpose );
  }
*/

  return precStatus;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiPrecond::apply
// Purpose       : Calls the actual preconditioner to apply y = M*x.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
int N_LAS_HBBlockJacobiPrecond::apply( N_LAS_MultiVector & x, N_LAS_MultiVector & y )
{
  int precStatus = 0;
  std::cout << " I'm applying the block Jacobi preconditioner!" << std::endl;
  // If there is no preconditioner to apply return a nonzero code
  if( Teuchos::is_null(epetraPrec_) )
    precStatus = -1;
  else
    precStatus = epetraPrec_->Apply( x.epetraObj(), y.epetraObj() );

  return precStatus;
}

