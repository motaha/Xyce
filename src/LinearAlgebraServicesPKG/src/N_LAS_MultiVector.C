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
// Filename       : $RCSfile: N_LAS_MultiVector.C,v $
//
// Purpose        : Implementation file for the Abstract interface to the
//                  multi-vector types (RDP, RSP, CDP or CSP).
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/20/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.101.2.3 $
//
// Revision Date  : $Date: 2013/12/08 17:50:58 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>

#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>
#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>

#include <N_ERH_ErrorMgr.h>

// ---------  Other Includes  -----------

#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

#include <EpetraExt_View_MultiVector.h>
#include <Teuchos_BLAS.hpp>

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::N_LAS_MultiVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
N_LAS_MultiVector::N_LAS_MultiVector(N_PDS_ParMap & map, int numVectors)
:  parallelMap_(&map),
   overlapMap_(&map),
   importer_(0),
   exporter_(0),
   viewTransform_(0),
   isOwned_(true)
{
  if (map.numGlobalEntities() < 0)
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,
    "N_LAS_MultiVector::MultiVector - vector length too short. Vectors must be > 0 in length.");

  // Create a new Petra MultiVector and set the pointer.
  aMultiVector_ = new Epetra_MultiVector( *map.petraBlockMap(), numVectors );

  oMultiVector_ = aMultiVector_;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::N_LAS_MultiVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/06/02
//-----------------------------------------------------------------------------
N_LAS_MultiVector::N_LAS_MultiVector( N_PDS_ParMap & map,
                                      N_PDS_ParMap & ol_map,
                                      int numVectors )
: parallelMap_(&map),
  overlapMap_(&ol_map),
  importer_(0),
  exporter_(0),
  viewTransform_(0),
  isOwned_(true)
{
  if (map.numGlobalEntities() < 0)
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,
    "N_LAS_MultiVector::MultiVector - vector length too short. Vectors must be > 0 in length.");

  // Create a new Petra MultiVector and set the pointer.
  oMultiVector_ = new Epetra_MultiVector(*ol_map.petraBlockMap(), numVectors);

  viewTransform_ = new EpetraExt::MultiVector_View(
    *overlapMap_->petraBlockMap(), *parallelMap_->petraBlockMap());
  aMultiVector_ = &((*viewTransform_)(*oMultiVector_));
#ifdef Xyce_PARALLEL_MPI
  if( map.pdsComm()->numProc() > 1 )
  {
    exporter_ = new Epetra_Export( *overlapMap_->petraBlockMap(),
                                   *parallelMap_->petraBlockMap() );
  }
#endif

  importer_ = new Epetra_Import( *overlapMap_->petraBlockMap(), *parallelMap_->petraBlockMap() );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::N_LAS_MultiVector
// Purpose       : Copy Constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
N_LAS_MultiVector::N_LAS_MultiVector( const N_LAS_MultiVector & right )
: parallelMap_( right.parallelMap_ ),
  overlapMap_( right.overlapMap_ ),
  oMultiVector_( new Epetra_MultiVector( *(right.oMultiVector_) ) ),
  importer_(0),
  exporter_(0),
  viewTransform_(0),
  isOwned_(true)
{
  if (right.aMultiVector_ == right.oMultiVector_)
    aMultiVector_ = oMultiVector_;
  else
  {
    viewTransform_ = new EpetraExt::MultiVector_View( *overlapMap_->petraBlockMap(),
                                                      *parallelMap_->petraBlockMap() );
    aMultiVector_ = &((*viewTransform_)( *oMultiVector_ ));
    if( right.parallelMap_->pdsComm()->numProc() > 1 )
    {
      if( right.exporter_ ) exporter_ = new Epetra_Export( *right.exporter_ );
    }
  }

  if( right.importer_ ) importer_ = new Epetra_Import( *right.importer_ );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::N_LAS_MultiVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
N_LAS_MultiVector::N_LAS_MultiVector( Epetra_MultiVector * overlapMV, Epetra_Map& parMap, bool isOwned )
: parallelMap_(0),
  overlapMap_(0),
  oMultiVector_( overlapMV ),
  importer_(0),
  exporter_(0),
  viewTransform_(0),
  isOwned_(isOwned)
{
  // Make sure there is anything to communicate before creating a transform, importer, or exporter
  if (overlapMV->MyLength() == parMap.NumMyElements())
    aMultiVector_ = oMultiVector_;
  else
  {
    viewTransform_ = new EpetraExt::MultiVector_View( overlapMV->Map(), parMap );
    aMultiVector_ = &((*viewTransform_)(*oMultiVector_));
#ifdef Xyce_PARALLEL_MPI
    if( parMap.Comm().NumProc() > 1 )
      exporter_ = new Epetra_Export( overlapMV->Map(), parMap );
#endif

    importer_ = new Epetra_Import( overlapMV->Map(), parMap );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::N_LAS_MultiVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
N_LAS_MultiVector::N_LAS_MultiVector( Epetra_MultiVector * origMV, bool isOwned )
: parallelMap_(0),
//  aMultiVector_( oMultiVector_ ),
  overlapMap_(0),
  oMultiVector_( origMV ),
  importer_(0),
  exporter_(0),
  viewTransform_(0),
  isOwned_(isOwned)
{
  aMultiVector_ = oMultiVector_;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::operator=
// Purpose       : assignment
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
N_LAS_MultiVector & N_LAS_MultiVector::operator=( const N_LAS_MultiVector & right )
{
  if( (oMultiVector_->Map().NumGlobalElements() == right.oMultiVector_->Map().NumGlobalElements())
     && (oMultiVector_->Map().NumMyElements() == right.oMultiVector_->Map().NumMyElements()) )
  {
    *oMultiVector_ = *right.oMultiVector_;
  }

  if( (aMultiVector_->Map().NumGlobalElements() == right.aMultiVector_->Map().NumGlobalElements())
     && (aMultiVector_->Map().NumMyElements() == right.aMultiVector_->Map().NumMyElements()) )
  {
    *aMultiVector_ = *right.aMultiVector_;
  }
  else
  {
#ifdef Xyce_VERBOSE_LINEAR
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,
	"MultiVector being assigned with different Mapping\n" );
#endif
  }

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::~N_LAS_MultiVector
// Purpose       : Default destructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
N_LAS_MultiVector::~N_LAS_MultiVector()
{
  if( importer_ ) delete importer_;
  if( exporter_ ) delete exporter_;
  if( viewTransform_ ) delete viewTransform_; //destroys of aMultiVector_ as well
  if (isOwned_)
  {
    if( oMultiVector_ ) delete oMultiVector_;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::globalLength
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
int N_LAS_MultiVector::globalLength() const
{
  return aMultiVector_->GlobalLength();
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::localLength
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
int N_LAS_MultiVector::localLength() const
{
  return aMultiVector_->MyLength();
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::numVectors
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
int N_LAS_MultiVector::numVectors() const
{
  return aMultiVector_->NumVectors();
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::dotProduct
// Purpose       : Returns the dot product of "this" vector and another.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
double N_LAS_MultiVector::dotProduct( const N_LAS_MultiVector & y ) const
{
  double result = 0.0;
  int PetraError = aMultiVector_->Dot(*(y.aMultiVector_), &result);

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_MultiVector::dotProduct - ", PetraError );
#endif

  return result;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::scale
// Purpose       : Scales a N_LAS_MultiVector by a constant value.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::scale(const double a)
{
  int PetraError = aMultiVector_->Scale(a);

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_MultiVector::scale - ", PetraError);
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::scale
// Purpose       : Scales a N_LAS_MultiVector by a constant value, but
//                 puts it into this.
//                 this = a*x
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/16/00
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::scale(const double a, const N_LAS_MultiVector &x)
{
  int PetraError = aMultiVector_->Scale(a, *(x.aMultiVector_));

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_MultiVector::scale - ", PetraError);
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::multiply
// Purpose       : Performs element-wise multiplication of two vectors
//                 this = this @ x
//                 where @ represents element-wise multiplication
// Special Notes :
// Scope         : Public
// Creator       : Roger P. Pawlowski, SNL, Parallel Computational Sciences
// Creation Date : 3/24/03
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::multiply(const N_LAS_MultiVector &x)
{
  int PetraError = aMultiVector_->Multiply(1.0, *aMultiVector_,
					   *(x.aMultiVector_), 0.0);

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_MultiVector::scale - ", PetraError);
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::daxpy
// Purpose       : Linear combination of two N_LAS_MultiVectors:
//                 this = y + a*x
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::daxpy(const N_LAS_MultiVector &y, const double a,
                              const N_LAS_MultiVector &x)
{
  int PetraError = aMultiVector_->Update(1.0, *(y.aMultiVector_), a,
                                         *(x.aMultiVector_), 0.0);

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_MultiVector::daxpy - ", PetraError);
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::linearCombo
// Purpose       : Linear combination of two N_LAS_MultiVectors:
//                 this = a*x + b*y
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/08/01
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::linearCombo(const double a, const N_LAS_MultiVector &x,
                                    const double b, const N_LAS_MultiVector &y)
{
  int PetraError = aMultiVector_->Update(a, *(x.aMultiVector_), b,
                                         *(y.aMultiVector_), 0.0);

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_MultiVector::linearCombo - ", PetraError);
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::linearCombo
// Purpose       : Linear combination of three N_LAS_MultiVectors:
//                 this = a*x + b*y + c*z
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/08/01
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::linearCombo(const double a, const N_LAS_MultiVector &x,
                                    const double b, const N_LAS_MultiVector &y,
                                    const double c, const N_LAS_MultiVector &z)
{
  update( a, x, b, y, c, z, 0.0);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::update
// Purpose       :
// Special Notes : ERK. From the epetra documentation:
//
//                 this = s*this + a*A
//
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 02/04/02
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::update( double a, const N_LAS_MultiVector & A,
                                double s )
{
  aMultiVector_->Update( a, *(A.aMultiVector_), s );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::update
// Purpose       :
// Special Notes : ERK.  From the epetra documentation:
//
//                 this = s*this + a*A + b*B
//
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 02/04/02
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::update( double a, const N_LAS_MultiVector & A,
                                double b, const N_LAS_MultiVector & B,
                                double s )
{
  aMultiVector_->Update( a, *(A.aMultiVector_),
                         b, *(B.aMultiVector_),
                         s );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::update
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 02/04/02
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::update( double a, const N_LAS_MultiVector & A,
                                double b, const N_LAS_MultiVector & B,
                                double c, const N_LAS_MultiVector & C,
                                double s )
{
  aMultiVector_->Update( a, *(A.aMultiVector_),
                         b, *(B.aMultiVector_),
                         s );
  aMultiVector_->Update( c, *(C.aMultiVector_),
                         1.0 );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::addVec
// Purpose       : Add multiple of a N_LAS_MultiVector:  this = this + a*y
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::addVec(const double a, const N_LAS_MultiVector &y)
{
  int PetraError = aMultiVector_->Update(a, *(y.aMultiVector_), 1.0);

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_MultiVector::addVec - ", PetraError);
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::gemm
// Purpose       : This routine performs a variety of MM operations.
//                 Variations are due to the fact that A, B and C can be local
//                 replicated or global distributed multi_vectors and that we
//                 may or may not operate with the transpose of A and B.
// Special Notes : See Petra documentation for possible cases.
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::GEMM(const bool transA, const bool transB,
                             const double alpha, const double beta,
                             const N_LAS_MultiVector &A,
                             const N_LAS_MultiVector &B)
{

  char tA, tB;

  transA ? (tA='N') : (tA='T');
  transB ? (tB='N') : (tB='T');

#ifdef Xyce_DEBUG_LINEAR

  const string methodMsg("N_LAS_MultiVector::GEMM - ");


  int PetraError = aMultiVector_->Multiply(tA, tB, alpha,
                         *(A.aMultiVector_), *(B.aMultiVector_), beta);

  cout << "sizes of this:  numvectors = ";
  cout << aMultiVector_->NumVectors ();
  cout << "  mylength = ";
  cout << aMultiVector_->MyLength ();
  cout << endl;


  cout << "sizes of A:  numvectors = ";
  cout << A.aMultiVector_->NumVectors ();

  cout << "  mylength = ";
  cout << A.aMultiVector_->MyLength ();
  cout << endl;

  cout << "sizes of B:  numvectors = ";
  cout << B.aMultiVector_->NumVectors();

  cout << "  mylength = ";
  cout << B.aMultiVector_->MyLength ();
  cout << endl;

  cout << "PetraError = " << PetraError << endl;
  processError(methodMsg, PetraError);

#else

  aMultiVector_->Multiply(tA, tB, alpha,
                          *(A.aMultiVector_), *(B.aMultiVector_), beta);

#endif

}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::lpNorm
// Purpose       : Returns lp norms of each vector in N_LAS_MultiVector
// Special Notes : Only p=1 and p=2 implemented now since this is all Petra
//                 supports.
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
int N_LAS_MultiVector::lpNorm(const int p, double * result) const
{
  int PetraError = -1;
  const string methodMsg("N_LAS_MultiVector::lpNorm - ");

  if (p == 1)
    PetraError = aMultiVector_->Norm1(result);
  else if (p == 2)
    PetraError = aMultiVector_->Norm2(result);
  else
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, methodMsg +
                           "Requested norm is not supported");

#ifdef Xyce_DEBUG_LINEAR
  processError(methodMsg, PetraError);
#endif

  return PetraError;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::infNorm
// Purpose       : Returns infinity norm of each vector in N_LAS_MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/16/01
//-----------------------------------------------------------------------------
int N_LAS_MultiVector::infNorm(double * result) const
{
  int PetraError = aMultiVector_->NormInf(result);

#ifdef Xyce_DEBUG_LINEAR
  const string methodMsg("N_LAS_MultiVector::infNorm - ");
  processError(methodMsg, PetraError);
#endif

  return PetraError;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::infNormIndex
// Purpose       : Returns index of the maximum absolute entry in N_LAS_MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Heidi K. Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 11/21/11
//-----------------------------------------------------------------------------
int N_LAS_MultiVector::infNormIndex(int * index) const
{
  Teuchos::BLAS<int,double> blas;

  int numProcs = aMultiVector_->Comm().NumProc();
  int numVectors = aMultiVector_->NumVectors();
  int myLength = aMultiVector_->MyLength();
  std::vector<int> indexTemp( numVectors, 0 ), indexTempAll( numVectors*numProcs, 0 );
  std::vector<double> doubleTemp( numVectors, 0.0 ), doubleTempAll( numVectors*numProcs, 0.0 );
  double ** pointers = aMultiVector_->Pointers();

  for (int i=0; i < numVectors; i++)
  {
    // Remember that IAMAX returns 1-based indexing, so subtract 1 to get the actual index.
    int jj = blas.IAMAX(myLength, pointers[i], 1) - 1;
    if (jj>-1)
    {
      indexTemp[i] = aMultiVector_->Map().GID(jj);
      doubleTemp[i] = std::abs(pointers[i][jj]);
    }
  }

  // Use the Epetra communicator to gather all the local maximum values and indices
  int result = aMultiVector_->Comm().GatherAll(&indexTemp[0], &indexTempAll[0], numVectors);
  result += aMultiVector_->Comm().GatherAll(&doubleTemp[0], &doubleTempAll[0], numVectors);

  // Compute the global infNorm and index
  for (int i=0; i < numVectors; i++)
  {
    // Remember that IAMAX returns 1-based indexing, so subtract 1 to get the actual index.
    int ii = blas.IAMAX( numProcs, &doubleTempAll[i], numVectors ) - 1;
    index[i] = indexTempAll[ii*numVectors + i];
  }

  return result;
}


//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::wRMSNorm
// Purpose       : Returns weighted root-mean-square of each vector in
//                 N_LAS_MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/12/00
//-----------------------------------------------------------------------------
int N_LAS_MultiVector::wRMSNorm(const N_LAS_MultiVector & weights, double * result) const
{
  int PetraError = aMultiVector_->NormWeighted( *(weights.aMultiVector_), result );

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_MultiVector::wRMSNorm - ", PetraError);
#endif

  return PetraError;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::wMaxNorm
// Purpose       : Returns the weighted
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 03/19/01
//-----------------------------------------------------------------------------
int N_LAS_MultiVector::wMaxNorm(N_LAS_MultiVector & weights,
                                N_LAS_MultiVector & tmpVector,
                                double * result) const
{
  int length  = aMultiVector_->MyLength();
  int numVecs = numVectors();

  for (int i = 0;  i < numVecs; ++i)
    for (int j = 0; j < length; ++j)
      tmpVector[i][j] = fabs((*(this))[i][j]) / weights[i][j];

  int PetraError = tmpVector.aMultiVector_->MaxValue(result);

  return PetraError;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::minValue
// Purpose       : Return the minimum value for the multivector
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/16/00
//-----------------------------------------------------------------------------
int N_LAS_MultiVector::minValue(double * result) const
{
  int PetraError = aMultiVector_->MinValue(result);

  return PetraError;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::maxValue
// Purpose       : Return the maximum value for the multivector
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/16/00
//-----------------------------------------------------------------------------
int N_LAS_MultiVector::maxValue(double * result) const
{
  int PetraError = aMultiVector_->MaxValue(result);

  return PetraError;
}


//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::reduce
// Purpose       : Takes a replicated local multi_vector and performs a global
//                 reduction on each entry.
// Special Notes : No_op with new Petra.
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::reduce()
{
  N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
    "N_LAS_MultiVector::reduce Not available in current PETRA!\n" );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::random
// Purpose       : Generates random numbers drawn from a uniform distribution
//                 on the interval (-1,1) using a multiplicative congruential
//                 generator with modulus 2^31 - 1.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::random()
{
  int PetraError = aMultiVector_->Random();

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_MultiVector::random - ", PetraError);
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::putScalar
// Purpose       : Fills N_LAS_MultiVector with the value "scalar".
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::putScalar(const double scalar)
{
  int PetraError = oMultiVector_->PutScalar(scalar);

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_MultiVector::putScalar - ", PetraError);
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::putScalarIncExt
// Purpose       : Fills N_LAS_MultiVector with the value "scalar".
// Special Notes : includes external indices
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/18/01
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::putScalarIncExt(const double scalar)
{
  if( aMultiVector_ == oMultiVector_ )
  {
    map<int, double>::iterator it_idM = externVectorMap_.begin();
    map<int, double>::iterator end_idM = externVectorMap_.end();
    for (; it_idM != end_idM; ++it_idM)
    {
      it_idM->second = scalar;
    }
  }

  putScalar(scalar);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::addScalar
// Purpose       : Adds to N_LAS_MultiVector with the value "scalar".
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/18/01
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::addScalar(const double scalar)
{
  int length  = aMultiVector_->MyLength();
  int numVecs = numVectors();

  for (int i = 0; i < numVecs; ++i)
    for (int j = 0; j < length; ++j)
      (*aMultiVector_)[i][j] += scalar;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::addScalarIncExt
// Purpose       : Adds to N_LAS_MultiVector with the value "scalar".
// Special Notes : includes external indices
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/18/01
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::addScalarIncExt(const double scalar)
{
  if( aMultiVector_ == oMultiVector_ )
  {
    map<int, double>::iterator it_idM = externVectorMap_.begin();
    map<int, double>::iterator end_idM = externVectorMap_.end ();
    for ( ; it_idM != end_idM; ++it_idM)
    {
      it_idM->second += scalar;
    }
  }

  addScalar(scalar);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::absValue
// Purpose       : Abs value of elements of N_LAS_MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/18/01
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::absValue(const N_LAS_MultiVector & A)
{
  int PetraError = oMultiVector_->Abs(*(A.oMultiVector_));

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_MultiVector::absValue - ", PetraError);
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::absValueIncExt
// Purpose       : Abs value of elements
// Special Notes : includes external indices
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/18/01
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::absValueIncExt(const N_LAS_MultiVector & A)
{
  if( aMultiVector_ == oMultiVector_ )
  {
    map<int, double>::iterator it_idM = externVectorMap_.begin();
    map<int, double>::iterator last = externVectorMap_.end();
    map<int, double>::const_iterator it_idM2 = A.externVectorMap_.begin();

    for ( ; it_idM != last; ++it_idM, ++it_idM2 )
      it_idM->second = fabs(it_idM2->second);
  }

  absValue(A);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::reciprocal
// Purpose       : Reciprocal elements of N_LAS_MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/07/01
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::reciprocal(const N_LAS_MultiVector & A)
{
  int PetraError = oMultiVector_->Reciprocal(*(A.oMultiVector_));

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_MultiVector::reciprocal - ", PetraError);
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::getVectorView
// Purpose       : Const view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
RefCountPtr<const N_LAS_Vector> N_LAS_MultiVector::getVectorView(int index) const
{
  RefCountPtr<const N_LAS_Vector> vec = rcp(new N_LAS_Vector((*oMultiVector_)(index),false), true);
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::getNonConstVectorView
// Purpose       : NonConst view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
RefCountPtr<N_LAS_Vector> N_LAS_MultiVector::getNonConstVectorView(int index)
{
  RefCountPtr<N_LAS_Vector> vec = rcp(new N_LAS_Vector((*oMultiVector_)(index),false), true);
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::operator []
// Purpose       : "[]" operator for N_LAS_MultiVector.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
double *& N_LAS_MultiVector::operator[](int index)
{
  return (*oMultiVector_)[index];
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::operator []
// Purpose       : "[]" operator for N_LAS_MultiVector.
// Special Notes : This version returns a "const double".
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
const double *& N_LAS_MultiVector::operator[](int index) const
{
  return const_cast <const double *&> ((*oMultiVector_)[index]);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::fillComplete
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/29/03
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::fillComplete()
{

#ifdef Xyce_PARALLEL_MPI

  if( exporter_ )
    aMultiVector_->Export( *oMultiVector_, *exporter_, Add );

#endif

}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::vectorImport
// Purpose       : Import using Petra_Import object
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/09/01
//-----------------------------------------------------------------------------
bool N_LAS_MultiVector::vectorImport(N_LAS_MultiVector * vec,
                                     Epetra_Import * importer)
{
  aMultiVector_->Import(*(vec->aMultiVector_), *importer, Insert);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::vectorExport
// Purpose       : Export using Petra_Import object
// Special Notes : Reverse of import function
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/09/01
//-----------------------------------------------------------------------------
bool N_LAS_MultiVector::vectorExport(N_LAS_MultiVector * vec,
                                     Epetra_Import * importer)
{
  aMultiVector_->Export(*(vec->aMultiVector_), *importer, Insert);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::importOverlap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/06/02
//-----------------------------------------------------------------------------
bool N_LAS_MultiVector::importOverlap()
{
  bool flag = false;

  if( importer_ )
    flag = oMultiVector_->Import( *aMultiVector_, *importer_, Insert );

  return flag;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::exportContribution
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/06/02
//-----------------------------------------------------------------------------
bool N_LAS_MultiVector::exportContribution( bool do_sum )
{
  if( exporter_ )
    if( do_sum )
      return aMultiVector_->Export( *oMultiVector_, *exporter_, Add );
    else
      return aMultiVector_->Export( *oMultiVector_, *exporter_, Zero );
  else
    return false;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::writeToFile
// Purpose       : Dumps out the multivector entries to a file.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/19/00
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::writeToFile(char *filename, bool useLIDs )
{
  int numProcs = aMultiVector_->Comm().NumProc();
  int localRank = aMultiVector_->Comm().MyPID();
  int masterRank = 0;

  for( int p = 0; p < numProcs; ++p )
  {
    //A barrier inside the loop so each processor waits its turn.
    aMultiVector_->Comm().Barrier();

    if(p == localRank)
    {
      FILE *file = NULL;

      if(masterRank == localRank)
      {
        //This is the master processor, open a new file.
        file = fopen(filename,"w");

        //Write the RDP_MultiVector dimension n into the file.
        fprintf(file,"%d\n",globalLength());
      }
      else
      {
        //This is not the master proc, open file for appending
        file = fopen(filename,"a");
      }

      //Now loop over the local portion of the RDP_MultiVector.
      int length  = localLength();
      int numVecs = numVectors();

      for (int i = 0; i < numVecs; ++i)
        for (int j = 0; j < length; ++j)
        {
          int loc = aMultiVector_->Map().GID(j);
          if( useLIDs ) loc = j;
          fprintf(file,"%d %d %20.13e\n",i,loc,(*aMultiVector_)[i][j]);
        }
      fclose(file);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::readFromFile
// Purpose       : Reads in the multivector entries from a file.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/29/03
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::readFromFile( char * filename )
{
}

#ifdef Xyce_DEBUG_LINEAR

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::processError
// Purpose       : Concrete implementation which processes Petra (in this case)
//                 error codes taken from the Petra member function returns.
// Special Notes : Petra specific.  NOTE ALSO - this function is currently
//                 within the "Xyce_DEBUG_LINEAR" ifdef and so any calls to
//                 this should also be so bracketed.
// Scope         : Private
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::processError(string methodMsg, int error) const
{

  const string PetraOK("Function returned without warnings or "
                                   "errors.\n");
  const string PetraError("Function returned with an error.\n");

  // Process the error
  switch (error)
  {
  case 0:
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, methodMsg + PetraOK);
    break;

  default:
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, methodMsg +
                           PetraError);
  }

}

#endif

#if 0
//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::getElementByGlobalIndex
// Purpose       : Get element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/6/00
//-----------------------------------------------------------------------------
const double & N_LAS_MultiVector::getElementByGlobalIndex(
  const int & global_index, const int & vec_index) const
{
  if( aMultiVector_ != oMultiVector_ )
    return (*oMultiVector_)[vec_index][overlapMap_->globalToLocalIndex(global_index)];
  else if( !parallelMap_ )
    return (*aMultiVector_)[vec_index][ aMultiVector_->Map().LID(global_index) ];
  else
  {
    int i = parallelMap_->globalToLocalIndex(global_index);

    if (i != -1)
      return ((*aMultiVector_)[vec_index])[i];
    else if ( externVectorMap_.count(global_index) )
      return externVectorMap_[global_index];
    else
    {
      char message[128];
      sprintf(message, "getElementByGlobalIndex: failed to find MultiVector "
              "global index. global_index = %d", global_index);
      string msg(message);

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      return externVectorMap_[-1];
    }
  }
}
#endif

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::getElementByGlobalIndex
// Purpose       : Get element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/6/00
//-----------------------------------------------------------------------------
const double & N_LAS_MultiVector::getElementByGlobalIndex(
  const int & global_index, const int & vec_index) const
{
  if( aMultiVector_ != oMultiVector_ )
    return (*oMultiVector_)[vec_index][overlapMap_->globalToLocalIndex(global_index)];
  else if( !parallelMap_ )
    return (*aMultiVector_)[vec_index][ aMultiVector_->Map().LID(global_index) ];
  else
  {
    int i = parallelMap_->globalToLocalIndex(global_index);

    if (i != -1)
      return ((*aMultiVector_)[vec_index])[i];
    else {
      map<int,double>::const_iterator it = externVectorMap_.find(global_index);
      if (it != externVectorMap_.end())
        return (*it).second;
      else
      {
        char message[128];
        sprintf(message, "getElementByGlobalIndex: failed to find MultiVector "
                "global index. global_index = %d", global_index);
        string msg(message);

        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
        return (*externVectorMap_.find(-1)).second;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::setElementByGlobalIndex
// Purpose       : Set element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/6/00
//-----------------------------------------------------------------------------
bool N_LAS_MultiVector::setElementByGlobalIndex(const int & global_index,
                                                const double & val,
                                                const int & vec_index)
{
  if( aMultiVector_ != oMultiVector_ )
    (*oMultiVector_)[vec_index][overlapMap_->globalToLocalIndex(global_index)] = val;
  else if( !parallelMap_ )
    (*oMultiVector_)[vec_index][ oMultiVector_->Map().LID(global_index) ] = val;
  else
  {
    if (global_index != -1)
    {
      int i = parallelMap_->globalToLocalIndex(global_index);

      if (i != -1)
      {
        ( (*aMultiVector_)[vec_index] )[i] = val;
        return STATUS_SUCCESS;
      }
      else
      {
        string msg = " setElementByGlobalIndex: failed to find MultiVector "
        "global index: ";

        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg, global_index);
        return STATUS_FAILURE;
      }
    }
  }

  return STATUS_SUCCESS;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::setElementByGlobalIndexIncExt
// Purpose       : Set element from vector using global index.
// Special Notes : Allows access to external indices block
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/18/01
//-----------------------------------------------------------------------------
bool N_LAS_MultiVector::setElementByGlobalIndexIncExt(const int & global_index,
                                                      const double & val,
                                                      const int & vec_index)
{
  if( aMultiVector_ != oMultiVector_ )
    (*oMultiVector_)[vec_index][overlapMap_->globalToLocalIndex(global_index)] = val;
  else if( !parallelMap_ )
    (*oMultiVector_)[vec_index][ oMultiVector_->Map().LID(global_index) ] = val;
  else
  {
    if (global_index != -1)
    {
      int i = parallelMap_->globalToLocalIndex(global_index);

      if (i != -1)
      {
        ( (*aMultiVector_)[vec_index] )[i] = val;
        return STATUS_SUCCESS;
      }
      else if(externVectorMap_.count(global_index))
      {
        externVectorMap_[global_index] = val;
        return STATUS_SUCCESS;
      }
      else
      {
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,
           string(" setElementByGlobalIndexIncExt: failed to find MultiVector "
                  "global index "));
        return STATUS_FAILURE;
      }
    }
  }

  return STATUS_SUCCESS;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::sumElementByGlobalIndex
// Purpose       : Set element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/7/00
//-----------------------------------------------------------------------------
bool N_LAS_MultiVector::sumElementByGlobalIndex(const int & global_index,
                                                const double & val,
                                                const int & vec_index)
{
  if( aMultiVector_ != oMultiVector_ )
    (*oMultiVector_)[vec_index][overlapMap_->globalToLocalIndex(global_index)] += val;
  else if( !parallelMap_ )
    (*oMultiVector_)[vec_index][ oMultiVector_->Map().LID(global_index) ] += val;
  else
  {
    if (global_index != -1 )
    {
      int i = parallelMap_->globalToLocalIndex(global_index);

      if (i != -1)
      {
        ( (*aMultiVector_)[vec_index] )[i] += val;
        return STATUS_SUCCESS;
      }
      else
      {
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,
                             string(" sumElementByGlobalIndex: failed to find "
                                    "MultiVector global index "));
        return STATUS_FAILURE;
      }
    }
  }

  return STATUS_SUCCESS;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::sumElementByGlobalIndexIncExt
// Purpose       : Set element from vector using global index.
// Special Notes : Allow access to external indices.
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/18/01
//-----------------------------------------------------------------------------
bool N_LAS_MultiVector::sumElementByGlobalIndexIncExt(const int & global_index,
                                                      const double & val,
                                                      const int & vec_index)
{
  if( aMultiVector_ != oMultiVector_ )
    (*oMultiVector_)[vec_index][overlapMap_->globalToLocalIndex(global_index)] += val;
  else if( !parallelMap_ )
    (*oMultiVector_)[vec_index][ oMultiVector_->Map().LID(global_index) ] += val;
  else
  {
    if (global_index != -1 )
    {
      int i = parallelMap_->globalToLocalIndex(global_index);

      if (i != -1)
      {
        ( (*aMultiVector_)[vec_index] ) [i] += val;
        return STATUS_SUCCESS;
      }
      else if( externVectorMap_.count(global_index) )
      {
        externVectorMap_[global_index] += val;
        return STATUS_SUCCESS;
      }
      else
      {
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,
                             string(" sumElementByGlobalIndexIncExt: failed to"
                                    "find MultiVector global index "));
        return STATUS_FAILURE;
      }
    }
  }

  return STATUS_SUCCESS;
}

//-----------------------------------------------------------------------------
// Function      : epetraVector
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/19/03
//-----------------------------------------------------------------------------
Epetra_Vector * N_LAS_MultiVector::epetraVector( int index ) const
{
  return new Epetra_Vector( View, *aMultiVector_, index );
}

//-----------------------------------------------------------------------------
// Function      : printPetraObject
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/14/00
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::printPetraObject() const
{
  cout << *aMultiVector_;
  cout << *oMultiVector_;

#ifdef Xyce_VERBOSE_LINEAR
  map<int, double>::const_iterator it_idM = externVectorMap_.begin();
  map<int, double>::const_iterator end_idM = externVectorMap_.end();
  if (it_idM != end_idM) cout << "<Extern Vector Map>" << endl;
  for (; it_idM != end_idM; ++it_idM)
  {
    cout << "  " << it_idM->first << "\t" << it_idM->second << endl;
  }
  cout << endl;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::packedByteCount
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/31/01
//-----------------------------------------------------------------------------
int N_LAS_MultiVector::packedByteCount() const
{
  int count = 0;
  count += sizeof(int); //# vectors
  count += sizeof(int); //local size
  count += numVectors() * aMultiVector_->MyLength() * sizeof(double); //local vectors

  count += sizeof(int); //extern map size
  count += externVectorMap_.size() * (sizeof(int)+sizeof(double));

  return count;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::pack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/31/01
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::pack( char * buf, int bsize, int & pos,
			N_PDS_Comm * comm ) const
{
  int nVs = numVectors();
  comm->pack( &nVs, 1, buf, bsize, pos );
  int locLth = aMultiVector_->MyLength();
  comm->pack( &locLth, 1, buf, bsize, pos );

  double * A;
  int stride;
  aMultiVector_->ExtractView( &A, &stride );

  comm->pack( A, locLth * nVs, buf, bsize, pos );

  int size = externVectorMap_.size();
  comm->pack( &size, 1, buf, bsize, pos );

  int loc;
  double val;
  map<int,double>::const_iterator iterXM = externVectorMap_.begin();
  map<int,double>::const_iterator endXM = externVectorMap_.end ();
  for(; iterXM != endXM; ++iterXM )
  {
    loc = iterXM->first;
    val = iterXM->second;
    comm->pack( &loc, 1, buf, bsize, pos );
    comm->pack( &val, 1, buf, bsize, pos );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MultiVector::unpack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/31/01
//-----------------------------------------------------------------------------
void N_LAS_MultiVector::unpack( char * buf, int bsize, int & pos,
			N_PDS_Comm * comm )
{
  int numV, locL;
  comm->unpack( buf, bsize, pos, &numV, 1 );
  comm->unpack( buf, bsize, pos, &locL, 1 );

  if( (numV != numVectors()) || (locL != aMultiVector_->MyLength()) )
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
	"ERROR:: UNPACKING VECTOR FAILED!\n" );

  vector<double> tmpVec( numV * locL );
  comm->unpack( buf, bsize, pos, &(tmpVec[0]), locL * numV );

  int loc = 0;
  for( int i = 0; i < numV; ++i )
    for( int j = 0; j < locL; ++j )
      (*aMultiVector_)[i][j] = tmpVec[loc++];

  int size;
  double val;
  comm->unpack( buf, bsize, pos, &size, 1 );
  externVectorMap_.clear();
  for( int i = 0; i < size; ++i )
  {
    comm->unpack( buf, bsize, pos, &loc, 1 );
    comm->unpack( buf, bsize, pos, &val, 1 );
    externVectorMap_[loc] = val;
  }
}
