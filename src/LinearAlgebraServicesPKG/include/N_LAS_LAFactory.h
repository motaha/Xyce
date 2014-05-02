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
// Filename       : $RCSfile: N_LAS_LAFactory.h,v $
//
// Purpose        : Specification file for the interface for creating linear
//                  algebra objects using the GoF Abstract Factory design
//                  pattern.  It must be used with a compatible linear algebra
//                  package (e.g., Trilinos/Petra to provide the concrete
//                  implementations.
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
// Revision Number: $Revision: 1.20 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_LAFactory_h
#define Xyce_N_LAS_LAFactory_h

// ---------- Standard Includes ----------

#include <vector>

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>

// ----------  Fwd Declares  -------------

class N_PDS_ParMap;

class N_LAS_MultiVector;
class N_LAS_Vector;
class N_LAS_Matrix;

typedef double N_LAS_DataType;

//-----------------------------------------------------------------------------
// Class         : N_LAS_LAFactory
// Purpose       : Provides an interface for creating linear algebra objects
//                 using the GoF Abstract Factory design pattern.
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
class N_LAS_LAFactory
{

 public:

  // Creates a new LAS vector
  static N_LAS_Vector * newVector( N_LAS_DataType type,
                                   N_PDS_ParMap & map );
  static N_LAS_Vector * newVector( N_LAS_DataType type,
                                   N_PDS_ParMap & map,
                                   N_PDS_ParMap & ol_map );

  // Creates a new LAS multi-vector
  static N_LAS_MultiVector * newMultiVector( N_LAS_DataType type,
                                             N_PDS_ParMap & map,
                                             int numVectors );
  static N_LAS_MultiVector * newMultiVector( N_LAS_DataType type,
                                             N_PDS_ParMap & map,
                                             int numVectors,
                                             N_PDS_ParMap & ol_map );

  // Creates a new LAS matrix
  static N_LAS_Matrix * newMatrix( N_LAS_DataType type,
                                   N_PDS_ParMap & map,
  	                           std::vector<int> & diagArray);

 private:

  // The Linear Algebra Services default factory
  N_LAS_LAFactory();

};

#endif

