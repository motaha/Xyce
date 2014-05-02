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
// Filename       : $RCSfile: N_LAS_MultiVector.h,v $
//
// Purpose        : Specification file for the Abstract interface to the
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
// Revision Number: $Revision: 1.72 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_MultiVector_h
#define Xyce_N_LAS_MultiVector_h

// ---------- Standard Includes ----------

#include <string>
#include <map>
#include <vector>

#include <N_UTL_Xyce.h>
#include <N_PDS_fwd.h>
#include <N_UTL_Packable.h>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

// --------  Forward Declarations --------

class N_PDS_ParMap;
class N_LAS_Vector;

class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Export;
class Epetra_Import;
class Epetra_Map;

namespace EpetraExt {

 class MultiVector_View;

}

//-----------------------------------------------------------------------------
// Class         : N_LAS_MultiVector
// Purpose       : Provides an interface for the multi-vector type
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
class N_LAS_MultiVector : public Packable
{

public:

  // Constructors to map to Petra constructors.
  N_LAS_MultiVector( N_PDS_ParMap & map,
                     int numVectors = 1 );

  N_LAS_MultiVector( N_PDS_ParMap & map,
                     N_PDS_ParMap & ol_map,
                     int numVectors = 1 );

  N_LAS_MultiVector( Epetra_MultiVector * origMV, bool isOwned = true );

  // Constructor takes the oMultiVector and generates the aMultiVector
  N_LAS_MultiVector( Epetra_MultiVector * overlapMV, Epetra_Map& parMap, bool isOwned = true );
  
  // Assignment operator
  N_LAS_MultiVector & operator=(const N_LAS_MultiVector & right);

  // Copy constructor
  N_LAS_MultiVector(const N_LAS_MultiVector & right);

  //Destructor
  virtual ~N_LAS_MultiVector();

  // Returns the dot product of "this" vector and another.
  double dotProduct(const N_LAS_MultiVector & y) const;

  // Scale every entry in the multi-vector by "a"
  void scale(const double a);
  // Scale every entry ([i]) in the multi-vector by "a*x[i]"
  void scale(const double a, const N_LAS_MultiVector & x);
  // Matrix-Matrix multiplication.  this[i] = this[i]*x[i] for each vector
  void multiply(const N_LAS_MultiVector & x);

  // Standard blas DAXPY operation
  void daxpy(const N_LAS_MultiVector & y, const double a,
  	const N_LAS_MultiVector & x);

  // Linear combination with two constants and vectors
  void linearCombo(const double a, const N_LAS_MultiVector & x,
  	const double b, const N_LAS_MultiVector & y);

  // Linear combination with three constants and vectors
  void linearCombo(const double a, const N_LAS_MultiVector & x,
  	const double b, const N_LAS_MultiVector & y, const double c,
  	const N_LAS_MultiVector & z);

  void update(double a, const N_LAS_MultiVector & A, double s);

  void update(double a, const N_LAS_MultiVector & A, double b,
  	const N_LAS_MultiVector & B, double s);

  void update(double a, const N_LAS_MultiVector & A, double b,
  	const N_LAS_MultiVector & B, double c, const N_LAS_MultiVector & C,
  	double s);

  // Compute the l_p norm (e.g., 2-norm is l_2)
  int lpNorm(const int p, double * result) const;

  // Infinity norm
  int infNorm(double * result) const;

  // Infinity norm index (the index pointer should be allocated by the caller)
  int infNormIndex(int * index) const;

  // Weighted root-mean-square norm
  int wRMSNorm(const N_LAS_MultiVector & weights, double * result) const;

  // Weighted max-norm
  int wMaxNorm(N_LAS_MultiVector & weights, N_LAS_MultiVector & tmpVector, double * result) const;

  // maximum value
  int maxValue(double * result) const;

  // minimum value
  int minValue(double * result) const;

  // Add a multiple of a multi-vector (this = a*y).
  void addVec(const double a, const N_LAS_MultiVector & y);

  // Performs a global reduction on each entry of a replicated multi-vector.
  void reduce();

  // Generate random number
  void random();

  // Performs a variety of MM operations.  Variations are due to the fact that
  // A, B and C can be local replicated or global distributed multi-vectors and
  // that we may or may not operate with the transpose of A and B.
  void GEMM(const bool transA, const bool transB, const double alpha,
  	const double beta, const N_LAS_MultiVector & A,
  	const N_LAS_MultiVector & B);

  // Fill vector with constant value.
  virtual void putScalar(const double scalar);

  // Add to vector with constant value.
  virtual void addScalar(const double scalar);

  // Absolute value element-wise for vector
  void absValue(const N_LAS_MultiVector & A);

  // Reciprocal of elements in vector
  void reciprocal(const N_LAS_MultiVector & A);

  // Index operator
  double *& operator[] (int index);
  // Index operator
  const double *& operator[] (int index) const;

  // Vector access function
  RCP<const N_LAS_Vector> getVectorView(int index) const;
  RCP<N_LAS_Vector> getNonConstVectorView(int index);

  // Get the global (across all processors) length of the multi-vector
  int globalLength() const;
  // Get the local (on processor component) length of the multi-vector
  int localLength() const;

  // Get the number of individual vectors in the multi-vector
  int numVectors() const;
  // Get the external vector size
  int externVectorSize() const { return externVectorMap_.size(); }

  // Import/Export capability
  bool vectorImport(N_LAS_MultiVector * vec, Epetra_Import * importer);
  bool vectorExport(N_LAS_MultiVector * vec, Epetra_Import * importer);

  bool importOverlap();
  bool exportContribution( bool do_sum = true );

  // Dump vector entries to file.
  virtual void writeToFile( char * filename, bool useLIDs = false );
  virtual void readFromFile( char * filename );

  // Friend of the N_LAS_Matrix and N_LAS_IterativeSolver classes so their
  // member functions can access our private members.
  friend class N_LAS_Matrix;
  friend class N_NLS_DampedNewton;

  // Get for vector elements by their global index (const version)
  const double & getElementByGlobalIndex(const int & global_index, const int & vec_index = 0) const;

  // Set for vector elements by their global index
  bool setElementByGlobalIndex(const int & global_index, const double & val,
  	const int & vec_index = 0);

  // Sum vector elements by their global index
  bool sumElementByGlobalIndex(const int & global_index, const double & val,
  	const int & vec_index = 0);

  // Clear the external vector map
  void clearExternVectorMap() { externVectorMap_.clear(); }
  // Add an element to the external vector map
  void addElementToExternVectorMap(const int & global_index,
  	const double & value);

  // Print the underlying Petra object.
  virtual void printPetraObject(std::ostream &os) const;

  Packable * instance() const { return new N_LAS_MultiVector(* this); }

  int packedByteCount() const;

  void pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const;
  void unpack(char * buf, int bsize, int & pos, N_PDS_Comm * comm);

  // Get the parallel map associated with this multi-vector
  N_PDS_ParMap * pmap() { return parallelMap_.get(); }
  N_PDS_ParMap * omap() { return overlapMap_.get(); }

  // Get the parallel communicator associated with this multi-vector
  N_PDS_Comm * pdsComm() { return pdsComm_.get(); }

  Epetra_MultiVector & epetraObj() { return *aMultiVector_; }

  Epetra_Vector * epetraVector( int index = 0 ) const;

  //Accumulate off processor fill contributions if necessary
  void fillComplete();

protected:

  // Pointer to the multi-vector's parallel map object
  RCP<N_PDS_ParMap> parallelMap_;

  // Parallel Map for overlapped data
  RCP<N_PDS_ParMap> overlapMap_;

  // Pointer the Petra multi-vector object.
  Epetra_MultiVector * aMultiVector_;

  // Overlapped view of multi-vector
  Epetra_MultiVector * oMultiVector_;

  //Used to distribute dependent data to processors needing it for fills
  Epetra_Import * importer_;

  //Used to accumulate off processor partial contributions to fills
  Epetra_Export * exporter_;

  // Subset View Transform
  EpetraExt::MultiVector_View * viewTransform_;

  // Communicator object, if one is needed.
  // NOTE: This communicator is only created when N_LAS_MultiVector is constructed with an Epetra_MultiVector object.
  RCP<N_PDS_Comm> pdsComm_;

  // isOwned flag
  bool isOwned_;

private:

  // Map containing extern elements from migration
    std::map<int,double> externVectorMap_;

  // Process library error codes.
  void processError(const char *methodMsg, int error) const;
};

//-----------------------------------------------------------------------------
inline void N_LAS_MultiVector::addElementToExternVectorMap(const int &
	global_index, const double & value)
{
  if( !externVectorMap_.count(global_index) )
   externVectorMap_[global_index] = value;
}

#endif

