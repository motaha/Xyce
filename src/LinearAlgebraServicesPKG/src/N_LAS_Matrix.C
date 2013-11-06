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
// Filename       : $RCSfile: N_LAS_Matrix.C,v $
//
// Purpose        : Implemenation file for the Abstract interface to sparse
//                  matrix type.
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
// Revision Number: $Revision: 1.98.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:45 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <algorithm>

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>

#include <N_PDS_ParMap.h>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>

#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

#include <EpetraExt_View_CrsMatrix.h>
#include <EpetraExt_RowMatrixOut.h>

#include <Epetra_OffsetIndex.h>

#if defined(Xyce_DEBUG_LINEAR) || defined(Xyce_DEBUG_DEVICE)
#include <N_ERH_ErrorMgr.h>
#endif

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::~N_LAS_Matrix
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
N_LAS_Matrix::~N_LAS_Matrix()
{
  if ( isOwned_ ) {
    if( oDCRSMatrix_ != aDCRSMatrix_ )
      if( viewTransform_ ) delete viewTransform_; //destroys aDCRSMatrix_ as well
      else                 delete aDCRSMatrix_;

    if( oDCRSMatrix_ ) delete oDCRSMatrix_;
  }

  if( exporter_ ) delete exporter_;
  if( offsetIndex_ ) delete offsetIndex_;
}


//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::N_LAS_Matrix
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
N_LAS_Matrix::N_LAS_Matrix( N_PDS_ParMap & map, vector<int> & diagArray )
: aDCRSMatrix_(0),
  oDCRSMatrix_(0),
  viewTransform_(0),
  exporter_(0),
  offsetIndex_(0),
  isOwned_(true)
{
  aDCRSMatrix_ = new Epetra_CrsMatrix( Copy, *map.petraMap() , &(diagArray[0]) );
  oDCRSMatrix_ = aDCRSMatrix_;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::N_LAS_Matrix
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
N_LAS_Matrix::N_LAS_Matrix( Epetra_CrsMatrix * origMatrix, bool isOwned )
: aDCRSMatrix_( origMatrix ),
  viewTransform_(0),
  exporter_(0),
  offsetIndex_(0),
  isOwned_(isOwned)
{
  oDCRSMatrix_ = aDCRSMatrix_;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::N_LAS_Matrix
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/21/02
//-----------------------------------------------------------------------------
N_LAS_Matrix::N_LAS_Matrix( Epetra_CrsGraph * overlapGraph,
                            Epetra_CrsGraph * baseGraph )
: aDCRSMatrix_(0),
  oDCRSMatrix_(0),
  viewTransform_(0),
  exporter_(0),
  offsetIndex_(0),
  isOwned_(true)
{
  oDCRSMatrix_ = new Epetra_CrsMatrix( Copy, *overlapGraph );

#ifndef Xyce_PARALLEL_MPI
  viewTransform_ = new EpetraExt::CrsMatrix_View( *overlapGraph, *baseGraph );
  aDCRSMatrix_ = &((*viewTransform_)( *oDCRSMatrix_ ));
#else
  aDCRSMatrix_ = new Epetra_CrsMatrix( Copy, *baseGraph );
  exporter_ = new Epetra_Export( overlapGraph->RowMap(), baseGraph->RowMap() );
  offsetIndex_ = new Epetra_OffsetIndex( *overlapGraph, *baseGraph, *exporter_ );
#endif	
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::fillComplete
// Purpose       : 
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/29/03
//-----------------------------------------------------------------------------
void N_LAS_Matrix::fillComplete()
{
  if( exporter_ )
  {
    aDCRSMatrix_->PutScalar( 0.0 ); //make sure assembled matrix is cleared
    aDCRSMatrix_->Export( *oDCRSMatrix_, *exporter_, Add, offsetIndex_ );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::matvec
// Purpose       : Sparse-matrix vector multiply - multivector version.  This
//                 function forms the product y = Ax where x and y are
//                 multivectors.  If transA is true, multiply by the transpose
//                 of matrix, otherwise just use matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void N_LAS_Matrix::matvec(bool transA, const N_LAS_MultiVector &x,
                          N_LAS_MultiVector &y)
{
  int PetraError = aDCRSMatrix_->Multiply(transA, *(x.aMultiVector_),
					  *(y.aMultiVector_));

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_Matrix::matvec - ", PetraError);
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::filterInfNorm
// Purpose       : Filter out values below inf norm * threshold
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/4/01
//-----------------------------------------------------------------------------
void N_LAS_Matrix::filterInfNorm( const double & threshold )
{
  if( aDCRSMatrix_->Filled() ) return;

  double thres = threshold * aDCRSMatrix_->NormInf();
  int numMyRows = aDCRSMatrix_->NumMyRows();
  int numEntries;
  int * indices;
  double * values;

  for( int i = 0; i < numMyRows; ++i )
  {
    aDCRSMatrix_->ExtractMyRowView( i, numEntries, values, indices );
    for( int j = 0; j < numEntries; ++j )
      if( fabs(values[j]) < thres ) values[j] = 0.0;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::filterRowSum
// Purpose       : Filter out values below row sum * threshold
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/4/01
//-----------------------------------------------------------------------------
void N_LAS_Matrix::filterRowSum( const double & threshold )
{
  if( aDCRSMatrix_->Filled() ) return;

  Epetra_Vector vec( aDCRSMatrix_->RowMap() );
  aDCRSMatrix_->InvRowSums( vec );

  int numMyRows = aDCRSMatrix_->NumMyRows();
  int numEntries;
  int * indices;
  double * values;

  for( int i = 0; i < numMyRows; ++i )
  {
    aDCRSMatrix_->ExtractMyRowView( i, numEntries, values, indices );
    double rowThres = threshold / vec[i];
    for( int j = 0; j < numEntries; ++j )
      if( fabs(values[j]) < rowThres ) values[j] = 0.0;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::put
// Purpose       : Put function for the sparse-matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void N_LAS_Matrix::put( double s )
{
  int PetraError =     aDCRSMatrix_->PutScalar(s);
  PetraError = oDCRSMatrix_->PutScalar(s);

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_Matrix::put - ", PetraError );
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::scale
// Purpose       : Scale function for the sparse-matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void N_LAS_Matrix::scale(double scaleFactor)
{
#ifdef Xyce_PARALLEL_MPI
  aDCRSMatrix_->Scale(scaleFactor);
#endif
  oDCRSMatrix_->Scale(scaleFactor);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::getRowLength
// Purpose       : Returns the number of nonzeroes in the row.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
int N_LAS_Matrix::getRowLength(int row) const
{
  return aDCRSMatrix_->NumAllocatedGlobalEntries(row);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::getRow
// Purpose       : Returns row coefficients and associated column indices.
// Special Notes : Uses Petra's ExtractRowCopy so assumes user has
//               : setup necessary space in arrays.  Could use
//               : ExtractRowView which does not require user to setup space.
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void N_LAS_Matrix::getRow(int row, int &length, double *coeffs, int *colIndices) const
{
  int PetraError = aDCRSMatrix_->ExtractGlobalRowView(row, length, coeffs, colIndices);

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_Matrix::getRow - ", PetraError );
#endif
}


//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::getRowCopy
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/05/06
//-----------------------------------------------------------------------------
void N_LAS_Matrix::getRowCopy
  (int row, int length, int & numEntries, double *coeffs, int *colIndices) const
{
  int PetraError = aDCRSMatrix_->ExtractGlobalRowCopy
       (row, length, numEntries, coeffs, colIndices);

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_Matrix::getRowCopy - ", PetraError );
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::getLocalRowCopy
// Purpose       : 
// Special Notes : 
//               : 
//               : 
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/05/06
//-----------------------------------------------------------------------------
void N_LAS_Matrix::getLocalRowCopy
  (int row, int length, int & numEntries, double *coeffs, int *colIndices) const
{
  int PetraError = aDCRSMatrix_->ExtractMyRowCopy
       (row, length, numEntries, coeffs, colIndices);

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_Matrix::getLocalRowCopy - ", PetraError );
#endif
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::putRow
// Purpose       : Put a row into the sparse matrix.
// Special Notes : Replace already allocated values
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
bool N_LAS_Matrix::putRow(int row, int length, const double *coeffs, const int *colIndices)
{
  int PetraError = oDCRSMatrix_->ReplaceGlobalValues(row, length, coeffs, colIndices);

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_Matrix::putRow - ", PetraError );
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::getLocalNumRows
// Purpose       : Returns the number of nonzeroes in the row.
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/24/06
//-----------------------------------------------------------------------------
int N_LAS_Matrix::getLocalNumRows() const
{
  return aDCRSMatrix_->NumMyRows();
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::getLocalRowLength
// Purpose       : Returns the number of nonzeroes in the row.
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/05/06
//-----------------------------------------------------------------------------
int N_LAS_Matrix::getLocalRowLength(int row) const
{
  return aDCRSMatrix_->NumMyEntries(row);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::putLocalRow
// Function      : N_LAS_Matrix::shirleyPutRow
// Purpose       : Put a row into the sparse matrix.
// Special Notes : Replace already allocated values.  
//                 erkeite: note; unlike putRow, this function uses the
//                 a-matrix.
//
//                 erkeite:  This function was mis-named as putLocalRow.  
//                 The row index being used here needs to be a gid (global), 
//                 not an lid (local).  I have re-named it shirleyPutRow
//                 for now.  (in part because I need a real putLocalRow).
//
//                 If it were truly a local function, it would call the
//                 epetra function: ReplaceMyValues.
//
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/24/06
//-----------------------------------------------------------------------------
//bool N_LAS_Matrix::putLocalRow(int row, int length, double *coeffs, int *colIndices)
bool N_LAS_Matrix::shirleyPutRow(int row, int length, double *coeffs, int *colIndices)
{
  int PetraError = aDCRSMatrix_->ReplaceGlobalValues(row, length, coeffs, colIndices);

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_Matrix::putRow - ", PetraError );
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::putLocalRow
// Purpose       : Put values into a row into the sparse matrix, using local indices.
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/29/12
//-----------------------------------------------------------------------------
bool N_LAS_Matrix::putLocalRow(int row, int length, const double * coeffs,
                                                   const int * colIndices)
{
  double * tmp_c = const_cast<double *>(coeffs);
  int * tmp_i = const_cast<int *>(colIndices);
  int PetraError = aDCRSMatrix_->ReplaceMyValues(row, length, tmp_c, tmp_i);

#if defined(Xyce_DEBUG_LINEAR) || defined(Xyce_DEBUG_DEVICE)
  processError( "N_LAS_Matrix::putLocalRow - ", PetraError );
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::insertRow
// Purpose       : Put a row into the sparse matrix.
// Special Notes : Insert into unallocated locations
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void N_LAS_Matrix::insertRow(int row, int length, double *coeffs, int *colIndices)
{
  int PetraError = oDCRSMatrix_->InsertGlobalValues(row, length, coeffs, colIndices);

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_Matrix::insertRow - ", PetraError );
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::replaceLocalRow
// Purpose       : Replace a row in the sparse matrix.
// Special Notes : replace allocated locations
// Scope         : Public
// Creator       : Todd Coffey, 1414, Heidi Thornquist, 1437
// Creation Date : 01/31/07
//-----------------------------------------------------------------------------
void N_LAS_Matrix::replaceLocalRow(int row, int length, double *coeffs, int *colIndices)
{
  int PetraError = oDCRSMatrix_->ReplaceMyValues(row, length, coeffs, colIndices);

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_Matrix::replaceLocalRow - ", PetraError );
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::getDiagonal
// Purpose       : Return the diagonal entries of the sparse matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void N_LAS_Matrix::getDiagonal( N_LAS_Vector & diagonal ) const
{
  int PetraError = aDCRSMatrix_->ExtractDiagonalCopy( *((*(diagonal.aMultiVector_))(0)) );

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_Matrix::getDiagonal - ", PetraError );
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::replaceDiagonal
// Purpose       : Replace values of diagonal elements
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/05/03
//-----------------------------------------------------------------------------
bool N_LAS_Matrix::replaceDiagonal( const N_LAS_Vector & vec )
{
  Epetra_Vector * eVec = vec.epetraVector();
  int PetraError = aDCRSMatrix_->ReplaceDiagonalValues( *eVec );
  delete eVec;

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_Matrix::replaceDiagonal - ", PetraError );
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::sumIntoDiagonal
// Purpose       : Sums values into diagonal elements
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/05/03
//-----------------------------------------------------------------------------
bool N_LAS_Matrix::sumIntoDiagonal( const N_LAS_Vector & vec )
{
  Epetra_Vector * eVec = vec.epetraVector();
  Epetra_Vector nVec( *eVec );
  aDCRSMatrix_->ExtractDiagonalCopy( nVec );
  nVec.Update( 1.0, *eVec, 1.0 );
  int PetraError = aDCRSMatrix_->ReplaceDiagonalValues( nVec );
  delete eVec;

#ifdef Xyce_DEBUG_LINEAR
  processError( "N_LAS_Matrix::replaceDiagonal - ", PetraError );
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::sumIntoRow
// Purpose       : Sum values into a row into the sparse matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
bool N_LAS_Matrix::sumIntoRow(int row, int length, const double * coeffs,
                                                   const int * colIndices)
{
  double * tmp_c = const_cast<double *>(coeffs);
  int * tmp_i = const_cast<int *>(colIndices);
  int PetraError = oDCRSMatrix_->SumIntoGlobalValues(row, length, tmp_c, tmp_i);

#if defined(Xyce_DEBUG_LINEAR) || defined(Xyce_DEBUG_DEVICE)
  processError( "N_LAS_Matrix::sumIntoRow - ", PetraError );
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::sumIntoLocalRow
// Purpose       : Sum values into a row into the sparse matrix, using local indices.
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/30/10
//-----------------------------------------------------------------------------
bool N_LAS_Matrix::sumIntoLocalRow(int row, int length, const double * coeffs,
                                                   const int * colIndices)
{
  double * tmp_c = const_cast<double *>(coeffs);
  int * tmp_i = const_cast<int *>(colIndices);
  int PetraError = oDCRSMatrix_->SumIntoMyValues(row, length, tmp_c, tmp_i);

#if defined(Xyce_DEBUG_LINEAR) || defined(Xyce_DEBUG_DEVICE)
  processError( "N_LAS_Matrix::sumIntoLocalRow - ", PetraError );
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::returnRawEntryPointer
//
// Purpose       : This function returns a raw double* pointer to a single 
//                 matrix element, specified by the local row,col indices.
//
// Special Notes : This function is much more convenient for developers of the
//                 device package than dealing with local compressed rows and
//                 the offsets required to use the bracket operators. 
//
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/18/2010
//-----------------------------------------------------------------------------
double * N_LAS_Matrix::returnRawEntryPointer (int lidRow, int lidCol)
{
  double * retPtr=0;

  int num_entries;
  int * indices;
  double * values;

  oDCRSMatrix_->ExtractMyRowView( lidRow, num_entries, values, indices );
   
  for( int j = 0; j < num_entries; ++j )
  {
     if (indices[j] == lidCol)
     {
       retPtr = &(values[j]);
       break;
     }
  }

  return retPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::extractLocalRowView
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/18/2010
//-----------------------------------------------------------------------------
int N_LAS_Matrix::extractLocalRowView(int lidRow, int& numEntries, double*& values, int*& indices) const
{
  return oDCRSMatrix_->ExtractMyRowView( lidRow, numEntries, values, indices );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::extractLocalRowView
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/18/2010
//-----------------------------------------------------------------------------
int N_LAS_Matrix::extractLocalRowView(int lidRow, int& numEntries, double*& values) const
{
  return oDCRSMatrix_->ExtractMyRowView( lidRow, numEntries, values);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::sumIntoRowWithTest
// Purpose       : Sum values into a row into the sparse matrix.
// Special Notes : Tests for redundant indices and fixes
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/04/00
//-----------------------------------------------------------------------------
bool N_LAS_Matrix::sumIntoRowWithTest(int row, int length, const double *coeffs,
                                                           const int *colIndices)
{
  int new_length = 0;
  vector<double> new_vals(length);
  vector<int> new_indices(length);

  int i, j;
  bool found;

  for( i = 0; i < length; ++i )
  {
    if( colIndices[i] != -1 )
    {
      found = false;
      for( j = 0; j < new_length; ++j )
        if( new_indices[j] == colIndices[i] )
        {
          new_vals[j] += coeffs[i];
          found = true;
          continue;
        }
      if( !found )
      {
        new_indices[new_length] = colIndices[i];
        new_vals[new_length] = coeffs[i];
        ++new_length;
      }
    }

  }

/*
  cout << "N_LAS_Matrix:sumIntoRowWithTest\n";
  for( int i = 0; i < length; ++i )
    cout << i << " " << colIndices[i] << " " << coeffs[i] << endl;
  cout << "-----------------------\n";
*/

  sumIntoRow( row, new_length, &(new_vals[0]), &(new_indices[0]) );

  return true;
}

#ifdef Xyce_LOCAL_INDEX_MATRIX

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::sumIntoLocalRow
// Purpose       : Sum values into a row into the sparse matrix.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/24/02
//-----------------------------------------------------------------------------
bool N_LAS_Matrix::sumIntoLocalRow( int row,
                                    int length,
                                    const vector<double> & vals,
                                    const vector<int> & indices )
{
  int new_length = 0;
  bool found;
  int i, j;

  vector<double> new_vals(length);
  vector<int> new_indices(length);

  for( i = 0; i < length; ++i )
  {
    if( aDCRSMatrix_->GCID(indices[i]) != -1 )
    {
      found = false;
      for( j = 0; j < new_length; ++j )
        if( new_indices[j] == indices[i] )
        {
          new_vals[j] += vals[i];
          found = true;
          continue;
        }
      if( !found )
      {
        new_indices[new_length] = indices[i];
        new_vals[new_length] = vals[i];
        ++new_length;
      }
    }
  }

/*
  cout << "N_LAS_Matrix:sumIntoLocalRow\n";
  for( int i = 0; i < length; ++i )
    cout << i << " " << indices[i] << " " << vals[i] << endl;
  cout << "-----------------------\n";
  for( int i = 0; i < new_length; ++i )
    cout << i << " " << new_indices[i] << " " << new_vals[i] << endl;
  cout << "-----------------------\n";
*/

  int PetraError = aDCRSMatrix_->SumIntoMyValues( row, new_length, &new_vals[0], &new_indices[0] );

#if defined(Xyce_DEBUG_LINEAR) || defined(Xyce_DEBUG_DEVICE)
  processError( "N_LAS_Matrix::sumIntoLocalRow - ", PetraError );
#endif

  return true;
}

#endif

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::add
// Purpose       : Sums in a matrix contribution
// Special Notes : WARNING: only works if graphs match, no checking
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
void N_LAS_Matrix::add( const N_LAS_Matrix & A )
{
  int NumRows = A.aDCRSMatrix_->NumMyRows();
  int* Indices;
  double* Values;
  int NumIndices;

  for( int i = 0; i < NumRows; ++i )
  {
    A.aDCRSMatrix_->ExtractMyRowView( i, NumIndices, Values, Indices );
    aDCRSMatrix_->SumIntoMyValues( i, NumIndices, Values, Indices );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::linearCombo
// Purpose       : Sums in a matrix contribution
// Special Notes : WARNING: only works if graphs EXACTLY match no checking
//
//                 this = a*A + b*B
//
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 2/13/07
//-----------------------------------------------------------------------------
void N_LAS_Matrix::linearCombo( const double a, const N_LAS_Matrix & A, 
                                const double b, const N_LAS_Matrix & B)
{
  int NumRows = (*aDCRSMatrix_).NumMyRows();
  
  int *aIndices, *bIndices;
  int aNumIndices, bNumIndices;
  double *aValues, *bValues;
  
  for( int i = 0; i < NumRows; ++i ) {
    // Get a view of the i-th row for A and B.
    A.aDCRSMatrix_->ExtractMyRowView( i, aNumIndices, aValues, aIndices );
    B.aDCRSMatrix_->ExtractMyRowView( i, bNumIndices, bValues, bIndices );
 
    // Add in the entries from each matrix.
    for ( int j = 0; j < aNumIndices; ++j )
      (*aDCRSMatrix_)[i][j] = a*aValues[j] + b*bValues[j];      
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::operator[]
// Purpose       : Direct access into matrix rows using local indexing
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/5/02
//-----------------------------------------------------------------------------
double * N_LAS_Matrix::operator[]( int row )
{
  return (*oDCRSMatrix_)[row];
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::operator[] const
// Purpose       : Direct access into matrix rows using local indexing
// Special Notes : const version
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/5/02
//-----------------------------------------------------------------------------
double * const & N_LAS_Matrix::operator[]( int row ) const
{
  return (*oDCRSMatrix_)[row];
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::rowMax
// Purpose       : Returns the maximum absolute value of the entries in the
//                 row.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
double N_LAS_Matrix::rowMax(int row) const
{
#ifdef Xyce_DEBUG_LINEAR
  N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
    "N_LAS_Matrix::rowMax - Not Available with Current Petra\n" );
#endif

  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::writeToFile
// Purpose       : Dumps out the sparse matrix to a file.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/19/00
//-----------------------------------------------------------------------------
void N_LAS_Matrix::writeToFile(char *filename, bool useLIDs, bool mmFormat )
{
  if (!mmFormat)
  {
    int numProcs = aDCRSMatrix_->Comm().NumProc();
    int thisProc = aDCRSMatrix_->Comm().MyPID();
    int masterRank = 0;

    int MaxNumEntries = aDCRSMatrix_->MaxNumEntries();
    vector<int> Indices( MaxNumEntries );
    vector<double> Values( MaxNumEntries );
    int NumEntries;
    int NumMyRows = aDCRSMatrix_->NumMyRows();

    if( !aDCRSMatrix_->Filled() )
    {
      cerr << "N_LAS_Matrix: can't writeToFile unless Filled!" << endl;
      return;
    }

    for( int p = 0; p < numProcs; ++p )
    {
      aDCRSMatrix_->Comm().Barrier();

      if( p == thisProc )
      {
        FILE *file = NULL;

        if( masterRank == thisProc )
        {
          file = fopen( filename, "w" );
          fprintf( file, "%d\n", aDCRSMatrix_->NumGlobalNonzeros() );
        }
        else
          file = fopen( filename, "a" );

        for( int i = 0; i < NumMyRows; ++i )
        {

          if( useLIDs )
          {
            int num_entries;
            int * indices;
            double * values;
            aDCRSMatrix_->ExtractMyRowView( i, num_entries, values, indices );
            for( int j = 0; j < num_entries; ++j )
             fprintf( file, "%d %d %26.18e\n", i, indices[j], values[j] );
          }
          else
          {
            int Row = aDCRSMatrix_->Graph().RowMap().GID(i);

            aDCRSMatrix_->ExtractGlobalRowCopy( Row, MaxNumEntries, NumEntries, &Values[0], &Indices[0] );

            for( int j = 0; j < NumEntries; ++j )
             fprintf( file, "%d %d %26.18e\n", Row, Indices[j], Values[j] );
          }
        }

        fclose( file );
      }
    }
  }
  else
  {
    std::string sandiaReq = "Sandia National Laboratories is a multi-program laboratory managed and operated by Sandia Corporation,\n%";
    sandiaReq += " a wholly owned subsidiary of Lockheed Martin Corporation, for the U.S. Department of Energy's National Nuclear \n%";
    sandiaReq += " Security Administration under contract DE-AC04-94AL85000.\n%\n% Xyce circuit matrix.\n%%";

    EpetraExt::RowMatrixToMatrixMarketFile( filename, *aDCRSMatrix_, sandiaReq.c_str() );
  }
}

#if defined(Xyce_DEBUG_LINEAR) || defined(Xyce_DEBUG_DEVICE)

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::processError
// Purpose       : Concrete implementation which processes Petra (in this case)
//                 error codes taken from the Petra member function returns.
// Special Notes : Petra specific.  NOTE ALSO - this function is currently
//                 within the "Xyce_DEBUG_LINEAR" ifdef and so any calls to
//                 this should also be so bracketed.
// Scope         : Private
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void N_LAS_Matrix::processError(string methodMsg, int error) const
{

  const string PetraError("Function returned with an error.\n");

  // Process the error
  if( error < 0 ) N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, methodMsg + PetraError);

}

#endif

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : output stream pipe operator for N_LAS_Matrix
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/14/00
//-----------------------------------------------------------------------------
void N_LAS_Matrix::printPetraObject() const
{
  cout << *oDCRSMatrix_;
  cout << *aDCRSMatrix_;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::setUseTranspose
// Purpose       : 
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/20/02
//-----------------------------------------------------------------------------
int N_LAS_Matrix::setUseTranspose (bool useTranspose)
{
  return aDCRSMatrix_->SetUseTranspose(useTranspose);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Matrix::useTranspose
// Purpose       : 
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/20/02
//-----------------------------------------------------------------------------
bool N_LAS_Matrix::useTranspose ()
{
  return aDCRSMatrix_->UseTranspose();
}
