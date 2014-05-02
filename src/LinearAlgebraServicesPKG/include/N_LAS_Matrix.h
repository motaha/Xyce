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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_LAS_Matrix.h,v $
//
// Purpose        : Specification file for the Abstract interface to sparse
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
// Revision Number: $Revision: 1.69 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Matrix_h
#define Xyce_N_LAS_Matrix_h

// ---------- Standard Includes ----------

#include <string>
#include <vector>

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>

// ----------  Other Includes   ----------

// ----------  Fwd Declares  -------------

class N_PDS_ParMap;
class N_LAS_MultiVector;
class N_LAS_Vector;

class Epetra_CrsMatrix;
class Epetra_CrsGraph;

class Epetra_Export;
class Epetra_OffsetIndex;

namespace EpetraExt {

class CrsMatrix_View;

}

//-----------------------------------------------------------------------------
// Class         : N_LAS_Matrix
// Purpose       : Abstract interface to sparse matrix type.
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
class N_LAS_Matrix
{

public:

  //Constructors
    N_LAS_Matrix( N_PDS_ParMap & map, std::vector<int> & diagArray);

  N_LAS_Matrix( Epetra_CrsGraph * overlapGraph,
                Epetra_CrsGraph * baseGraph );

  //Constructor from an existing Epetra_CrsMatrix (makes copy of origMatrix)
  N_LAS_Matrix( Epetra_CrsMatrix * origMatrix, bool isOwned = true );

  //Destructor
  virtual ~N_LAS_Matrix();

  // This function needs to be invoked for a transpose solve.
  int setUseTranspose (bool useTranspose);
  bool useTranspose ();

  //Accumulate off processor fill contributions if necessary
  virtual void fillComplete();

  // Sparse-matrix vector multiply - multivector version.  If transA is true,
  // multiply by the transpose of matrix, otherwise just use matrix.
  void matvec(bool transA, const N_LAS_MultiVector & x,
  	N_LAS_MultiVector & y);

  // Add in a matrix contribution
  void add( const N_LAS_Matrix & A );

  void linearCombo ( const double a, const N_LAS_Matrix & A, 
                     const double b, const N_LAS_Matrix & B);

  //Filter to zero out values below threshold times inf norm
  void filterInfNorm(const double & threshold);
  //Filters to zero out values below threshold times the row sum
  void filterRowSum(const double & threshold);

  // Put function for the sparse-matrix.
  virtual void put(double s);

  // Scale the matrix
  void scale(double scaleFactor);

  // Diagonal Operations

  // Get the matrix diagonal (stored in an N_LAS_Vector)
  void getDiagonal(N_LAS_Vector & diagonal) const;
  // Sum a vector's values into the matrix diagonal
  bool sumIntoDiagonal(const N_LAS_Vector & vec);
  bool replaceDiagonal(const N_LAS_Vector & vec);

  // Row Operations

  // Get a row's maximum value
  double rowMax(int row) const;
  // Get a row's length (nonzero)
  int getRowLength(int row) const;
  // Get the non-zero values in a row
  void getRow(int row, int & length, double * coeffs, int * colIndices) const;

  // Get the non-zero values in a row
  void getRowCopy(int row, int length, int & numEntries, double * coeffs, int * colIndices) const;
  void getLocalRowCopy(int row, int length, int & numEntries, double * coeffs, int * colIndices) const;

  int getLocalNumRows() const;
  int getLocalRowLength(int row) const;
  //bool putLocalRow(int row, int length, double * coeffs, int * colIndices);
  bool shirleyPutRow(int row, int length, double * coeffs, int * colIndices);

  // Put a set of values into a row
  bool putRow(int row, int length, const double * coeffs, const int * colIndices);

  // Put a set of values into a row, using local indices
  bool putLocalRow(int row, int length, const double * coeffs, const int * colIndices);

  // Insert a set of values into a row
  void insertRow(int row, int length, double * coeffs, int * colIndices);

  // Replace a set of values into a row
  void replaceLocalRow(int row, int length, double * coeffs, int * colIndices);

  // Sum values into a row into the sparse matrix.
  bool sumIntoRow(int row, int length, const double * coeffs, const int * colIndices);

  // Sum values into a row into the sparse matrix, using local indices
  bool sumIntoLocalRow(int row, int length, const double * coeffs, const int * colIndices);

  // Return a pointer to a single row, col element.
  double * returnRawEntryPointer (int lidRow, int lidCol);

  // get a pointer to the compressed local row.
  int extractLocalRowView(int lidRow, int& numEntries, double*& values, int*& indices) const;

  // get a pointer to the compressed local row.
  int extractLocalRowView(int lidRow, int& numEntries, double*& values) const;

  // Tests for redundant indices first.
  bool sumIntoRowWithTest(int row, int length, const double * coeffs, 
      const int * colIndices);

#ifdef Xyce_LOCAL_INDEX_MATRIX
  // Sum values into a row into the sparse matrix.
  bool sumIntoLocalRow( int row,
                        int length,
                        const std::vector<double> & vals,
                        const std::vector<int> & colIndices);
#endif

  double * operator[]( int row );
  double * const operator[]( int row ) const;

  Epetra_CrsMatrix & epetraObj() { return *aDCRSMatrix_; }

  // Output the matrix to a file
  void writeToFile(char * filename, bool useLIDs = false, bool mmFormat=false );
  // Print the underlying Epetra objects
  virtual void printPetraObject(std::ostream &os) const;

  // Friend in the N_LAS_MultiVector and N_LAS_IterativeSolver classes so their
  // member functions can access our private members.
  friend class N_LAS_MultiVector;
  friend class N_NLS_DampedNewton;

protected:

  // Pointer the Petra multi-vector object.
  Epetra_CrsMatrix * aDCRSMatrix_;

  // Overlapped version of matrix
  Epetra_CrsMatrix * oDCRSMatrix_;

  // Subset View Transform
  EpetraExt::CrsMatrix_View * viewTransform_;

  // Importing Tools
  Epetra_Export * exporter_;
  Epetra_OffsetIndex * offsetIndex_;

private:

  // Default constructor (private)
  N_LAS_Matrix();

  // Copy constructor (private)
  N_LAS_Matrix(const N_LAS_Matrix & right);
  // Assignment operator (private)
  N_LAS_Matrix & operator = (const N_LAS_Matrix & right);

  bool operator == (const N_LAS_Matrix & right) const;
  bool operator != (const N_LAS_Matrix & right) const;

  // isOwned flag
  bool isOwned_;

#if defined(Xyce_DEBUG_LINEAR) || defined(Xyce_DEBUG_DEVICE)
  // Process library error codes.
  void processError(std::string methodMsg, int error) const;
#endif

};

#endif


