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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_LAS_LAFactory.C,v $
//
// Purpose        : Implementation file for the interface for creating linear
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
// Revision Number: $Revision: 1.16 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_LAS_LAFactory.h>

#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

//-----------------------------------------------------------------------------
// Function      : N_LAS_LAFactory::newVector
// Purpose       : Concrete implementation of the "newVector" function for
//                 creating a single vector.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 10/16/00
//-----------------------------------------------------------------------------
N_LAS_Vector * N_LAS_LAFactory::newVector( N_LAS_DataType type,
                                           N_PDS_ParMap & map )
{
  N_LAS_Vector * tmp = new N_LAS_Vector( map );
//  tmp->putScalar( type );
  return tmp;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_LAFactory::newVector
// Purpose       : Concrete implementation of the "newVector" function for
//                 creating a single vector.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 10/16/00
//-----------------------------------------------------------------------------
N_LAS_Vector * N_LAS_LAFactory::newVector( N_LAS_DataType type,
                                           N_PDS_ParMap & map,
                                           N_PDS_ParMap & ol_map )
{
  N_LAS_Vector * tmp = new N_LAS_Vector( map, ol_map );
//  tmp->putScalar( type );
  return tmp;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_LAFactory::newMultiVector
// Purpose       : Concrete implementation of the "newMultiVector" function.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
N_LAS_MultiVector * N_LAS_LAFactory::newMultiVector( N_LAS_DataType type,
                                                     N_PDS_ParMap & map,
                                                     int numVectors)
{
  N_LAS_MultiVector * tmp = new N_LAS_MultiVector( map, numVectors );
//  tmp->putScalar( type );
  return tmp;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_LAFactory::newMultiVector
// Purpose       : Concrete implementation of the "newMultiVector" function.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
N_LAS_MultiVector * N_LAS_LAFactory::newMultiVector( N_LAS_DataType type,
                                                     N_PDS_ParMap & map,
                                                     int numVectors,
                                                     N_PDS_ParMap & ol_map )
{
  N_LAS_MultiVector * tmp = new N_LAS_MultiVector( map, ol_map, numVectors );
//  tmp->putScalar( type );
  return tmp;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_LAFactory::newMatrix
// Purpose       : Concrete implementation of the "newMatrix" function.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
N_LAS_Matrix * N_LAS_LAFactory::newMatrix( N_LAS_DataType type,
                                           N_PDS_ParMap & map,
					   std::vector<int> & diagArray)
{
  N_LAS_Matrix * tmp = new N_LAS_Matrix( map, diagArray );
//  tmp->put( type );
  return tmp;
}

