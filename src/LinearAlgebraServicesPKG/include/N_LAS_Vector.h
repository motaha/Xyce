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
// Filename       : $RCSfile: N_LAS_Vector.h,v $
//
// Purpose        : Specification file for the Abstract interface to the
//                  vector types (RDP, RSP, CDP or CSP).
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Computational Sciences
//
// Creation Date  : 10/13/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.24 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Vector_h
#define Xyce_N_LAS_Vector_h

// ---------- Standard Includes ----------

#include <string>
#include <map>

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>

#include <N_LAS_MultiVector.h>

#ifdef Xyce_INLINED_VECTOR_BRACKET_OPERATOR
 #include <Epetra_MultiVector.h>
#endif

// ----------  Other Includes   ----------

// --------  Forward Declarations --------

class Epetra_Vector;

class N_PDS_ParMap;

//-----------------------------------------------------------------------------
// Class         : N_LAS_Vector
// Purpose       : Provides an abstract interface for the vector type (RDP,
//                 RSP, CDP or CSP).
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Compuational Sciences
// Creation Date : 10/13/00
//-----------------------------------------------------------------------------
class N_LAS_Vector : public N_LAS_MultiVector
  {

public:

  // Constructors to map to Petra constructors.
  N_LAS_Vector(N_PDS_ParMap & map)
  : N_LAS_MultiVector(map, 1)
  {}

  N_LAS_Vector( N_PDS_ParMap & map, N_PDS_ParMap & ol_map )
  : N_LAS_MultiVector( map, ol_map, 1 )
  {}

  //Copy constructor
  N_LAS_Vector( const N_LAS_Vector & right )
  : N_LAS_MultiVector(right)
  {}

  //Copy constructor
  N_LAS_Vector( Epetra_Vector * origV, bool isOwned = true );
  //N_LAS_Vector( Epetra_Vector * origV );

  // Constructor takes the oMultiVector and generates the aMultiVector
  N_LAS_Vector( Epetra_Vector * overlapMV, Epetra_Map& parMap, bool isOwned = true );

  // Destructor
  virtual ~N_LAS_Vector() {}

  // Operation: operator []
#ifndef Xyce_INLINED_VECTOR_BRACKET_OPERATOR
  double & operator[] (int index);
#else
  double & operator[] (int index) { return (*oMultiVector_)[0][index]; }
#endif

  // Operation: operator []
  const double & operator[] (int index) const;

  };

#endif

