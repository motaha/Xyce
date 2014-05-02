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
// Filename       : $RCSfile: N_LAS_Vector.C,v $
//
// Purpose        : Implementation file for the Abstract interface to the
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
// Revision Number: $Revision $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>
#include <N_LAS_Vector.h>

// ---------  Other Includes  -----------

#include <Epetra_Vector.h>

//-----------------------------------------------------------------------------
// Function      : N_LAS_Vector:::N_LAS_Vector
// Purpose       : constructor 
// Special Notes : 
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
N_LAS_Vector::N_LAS_Vector( Epetra_Vector * origV, bool isOwned )
: N_LAS_MultiVector( dynamic_cast<Epetra_MultiVector *>(origV), isOwned )
{
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Vector:::N_LAS_Vector
// Purpose       : constructor
// Special Notes : 
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
N_LAS_Vector::N_LAS_Vector( Epetra_Vector * overlapV, Epetra_Map& parMap, bool isOwned )
: N_LAS_MultiVector( dynamic_cast<Epetra_MultiVector *>(overlapV), parMap, isOwned )
{
}

#ifndef Xyce_INLINED_VECTOR_BRACKET_OPERATOR

//-----------------------------------------------------------------------------
// Function      : N_LAS_Vector:::operator[]
// Purpose       : "[]" operator for N_LAS_Vector.
// Special Notes : This version returns a "double".
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 10/13/00
//-----------------------------------------------------------------------------
double & N_LAS_Vector::operator[](int index)
{
  return (*oMultiVector_)[0][index];
}

#endif

//-----------------------------------------------------------------------------
// Function      : N_LAS_Vector:::operator[]
// Purpose       : "[]" operator for N_LAS_Vector.
// Special Notes : This version returns a "double".
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 10/13/00
//-----------------------------------------------------------------------------
const double & N_LAS_Vector::operator[](int index) const
{
  return (*oMultiVector_)[0][index];
}

