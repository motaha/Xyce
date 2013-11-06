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
// Filename       : $RCSfile: N_NLS_NOX_Vector.h,v $
//
// Purpose        : Interface to Xyce vectors for NOX.
//
// Special Notes  :
//
// Creator        : Tammy Kolda, NLS, 8950
//
// Creation Date  : 01/31/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.15.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:47 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_Vector_h
#define Xyce_N_NLS_NOX_Vector_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>

// ----------   NOX Includes   ----------

#include "NOX_Abstract_Vector.H"

// ---------- Forward Declarations ----------

class N_LAS_Vector;
class N_LAS_System;

// ---------- Namespace Declarations ----------

// N_NLS namespace is for the Xyce Nonlinear Solver Package
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Class         : N_NLS::NOX::Vector
//
// Purpose       :
//
//      NOX Vector Interface for Xyce vectors.
//
// Creator       : Tammy Kolda, SNL, 8950
//
// Creation Date : 2/1/02
//-----------------------------------------------------------------------------

class Vector : public NOX::Abstract::Vector {

public:

  //---------------------------------------------------------------------------
  // Function      : Vector (constructor)
  //
  // Purpose       : Constructs a NOX-compatiable vector object
  //                 containing the given Xyce-compatible vector.
  //
  // Special Notes : The N_LAS_System pointer is only needed to
  //                 support cloning. If N_LAS_Vector supported the
  //                 copy constructor, this would not be necessary.
  //---------------------------------------------------------------------------
  Vector(N_LAS_Vector& vector, N_LAS_System& lasSys);

  //---------------------------------------------------------------------------
  // Function      : Vector (copy constructor)
  // Purpose       : Constructs a Vector using the source vector.
  //---------------------------------------------------------------------------
  Vector(const Vector& source, NOX::CopyType type = NOX::DeepCopy);

  //---------------------------------------------------------------------------
  // Function      : Destructor
  //
  // Purpose       : Deletes the internal N_LAS_Vector only if that
  //                 vector was created via the copy constructor.
  //---------------------------------------------------------------------------
  ~Vector();

  //---------------------------------------------------------------------------
  // Purpose       : Return the length (i.e., number of entries) of the vector
  //---------------------------------------------------------------------------
  int length() const;

  //---------------------------------------------------------------------------
  // Purpose       : Initialize every entry in the vector to the given value
  //---------------------------------------------------------------------------
  NOX::Abstract::Vector& init(double value);

  //---------------------------------------------------------------------------
  // Purpose       : Compute the element-wise absolute value of source
  //---------------------------------------------------------------------------
  NOX::Abstract::Vector& abs(const Vector& source);
  NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& source);

  //---------------------------------------------------------------------------
  // Purpose       : Copy the source vector
  //---------------------------------------------------------------------------
  NOX::Abstract::Vector& operator=(const Vector& source);
  NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& source);

  //---------------------------------------------------------------------------
  // Purpose       : Compute the element-wise reciprocal of source
  //---------------------------------------------------------------------------
  NOX::Abstract::Vector& reciprocal(const Vector& source);
  NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& source);

  //---------------------------------------------------------------------------
  // Purpose       : Scale this vector by gamma
  //---------------------------------------------------------------------------
  NOX::Abstract::Vector& scale(double gamma);

  //---------------------------------------------------------------------------
  // Purpose       : Scale this vector element-by-element by y
  //---------------------------------------------------------------------------
  NOX::Abstract::Vector& scale(const Vector& y);
  NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& y);

  //---------------------------------------------------------------------------
  // Purpose       : this = alpha * a + gamma * this
  //---------------------------------------------------------------------------
  NOX::Abstract::Vector& update(double alpha, const Vector& a,
				double gamma = 0.0);
  NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a,
				double gamma = 0.0);

  //---------------------------------------------------------------------------
  // Purpose       : this = alpha * a + beta * b + gamma * this
  //---------------------------------------------------------------------------
  NOX::Abstract::Vector& update(double alpha, const Vector& a,
				double beta, const Vector& b,
				double gamma = 0.0);
  NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a,
				double beta, const NOX::Abstract::Vector& b,
				double gamma = 0.0);

  //---------------------------------------------------------------------------
  // Purpose       : this = alpha * a + beta * b + gamma * this
  //---------------------------------------------------------------------------
  NOX::Abstract::Vector& random(bool useSeed=false, int seed=1);

  //---------------------------------------------------------------------------
  // Purpose       : Return a pointer to a new cloned vector of this type.
  //---------------------------------------------------------------------------
  Teuchos::RCP<NOX::Abstract::Vector> 
    clone(NOX::CopyType type = NOX::DeepCopy) const;

  //---------------------------------------------------------------------------
  // Purpose       : Return the norm of this vector
  //---------------------------------------------------------------------------
  double norm(NOX::Abstract::Vector::NormType type = NOX::Abstract::Vector::TwoNorm) const;

  //---------------------------------------------------------------------------
  // Purpose       : Return the norm of weights .* this
  //                 (i.e., elementwise multiplication)
  //---------------------------------------------------------------------------
  double norm(const Vector& weights) const;
  double norm(const NOX::Abstract::Vector& weights) const;

  //---------------------------------------------------------------------------
  // Purpose       : Return this' * y (dot product of this with y)
  //---------------------------------------------------------------------------
  double innerProduct(const Vector& y) const;
  double innerProduct(const NOX::Abstract::Vector& y) const;

  //---------------------------------------------------------------------------
  // Purpose       : Return const pointer to the underlying N_LAS_Vector
  //---------------------------------------------------------------------------
  const N_LAS_Vector& getNativeVectorRef() const {return *vectorPtr_;};

  //---------------------------------------------------------------------------
  // Purpose       : Return pointer to the underlying N_LAS_Vector
  //---------------------------------------------------------------------------
  N_LAS_Vector& getNativeVectorRef() {return *vectorPtr_;};

  N_LAS_Vector* getNativeVectorPtr() {return  vectorPtr_;};

  N_LAS_Vector* getNativeVectorPtr() const {return  vectorPtr_;};

  //---------------------------------------------------------------------------
  // Purpose       : Print the underlying vector
  //---------------------------------------------------------------------------
  void print() const;

  //---------------------------------------------------------------------------
  // Purpose       : Return the pointer to the underlying N_LAS_Vector
  //---------------------------------------------------------------------------
  const N_LAS_Vector& getNativeVectorRef_() const {return *vectorPtr_;};

private:
  // Vector stored by this object
  N_LAS_Vector* vectorPtr_;

  // Only used in copy constructor and only because N_LAS_Vector does
  // not have a copy constructor
  N_LAS_System& lasSys_;

  // True is N_LAS_Vector should be deleted when this object is
  // destructed. This is true for any Vector created by the copy
  // constructor.
  bool doDelete_;

}; // class SharedSystem
} // namespace N_NLS_NOX

#endif // Xyce_N_NLS_NOX_SharedSystem_h

