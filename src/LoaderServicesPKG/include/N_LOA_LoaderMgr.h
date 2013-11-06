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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_LOA_LoaderMgr.h,v $
//
// Purpose        : This file contains definition(s) for the loader
//                  manager class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/03/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.13.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:46 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef  Xyce_LOA_LoaderMgr_H
#define  Xyce_LOA_LoaderMgr_H

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>

// ---------- Enum Definitions ----------

enum Loader_index
{
  _CKTLOADER // 0
};

// ---------- Forward declarations ---------

class N_LOA_Loader;
class N_LOA_CktLoader;
class N_LOA_NonlinearEquationLoader;
class N_LAS_Matrix;
class N_LAS_Vector;

//-----------------------------------------------------------------------------
// Class         : N_LOA_LoaderMgr
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/03/00
//-----------------------------------------------------------------------------
class N_LOA_LoaderMgr
{
public:

  // Default constructor
  N_LOA_LoaderMgr();

  // Destructor
  ~N_LOA_LoaderMgr();

  // Creation functions:

  // Create loader object
  N_LOA_Loader * createLoader(int loaderIndex);

  // Create the circuit loader object
  N_LOA_CktLoader * createCktLoader();

  // Create the DAE  loader object
  N_LOA_NonlinearEquationLoader * createNonlinearEquationLoader();

protected:
private :

  // Pointer to the loader object
  N_LOA_Loader * loaderPtr_;

  // Pointer to the circuit loader object
  N_LOA_CktLoader * cktLoaderPtr_;

  // Pointer to the DAE loader object
  N_LOA_NonlinearEquationLoader * nonlinearEquationLoaderPtr_;

};

#endif



