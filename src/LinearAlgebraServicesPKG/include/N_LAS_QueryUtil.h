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
// Filename       : $RCSfile: N_LAS_QueryUtil.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/09/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.25 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------


#ifndef N_LAS_QueryUtil_h
#define N_LAS_QueryUtil_h 1

// ---------- Standard Includes ----------
#include <string>
#include <list>
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>

// ---------- Forward Declarations ----------

//-----------------------------------------------------------------------------
// Class         : N_LAS_QueryUtil
// Purpose       : Abstract class for Query Utility used by LAS to
//                 determine parameters for linear system instantiation
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/09/00
//-----------------------------------------------------------------------------
class N_LAS_QueryUtil
{
public:

  // Destructor
  virtual ~N_LAS_QueryUtil() {}

  //------ generation of Linear System Data
  virtual void generateRowColData() = 0;

  //------ accessor methods

  virtual int numGlobalRows() const = 0;
  virtual int numLocalRows() const = 0;
  virtual int numExternRows() const = 0;
  virtual int numGlobalExternRows() const = 0;

  virtual int numGlobalStateVars() const = 0;
  virtual int numLocalStateVars() const = 0;
  virtual int numExternStateVars() const = 0;
  virtual int numGlobalExternStateVars() const = 0;

  virtual int numGlobalStoreVars() const = 0;
  virtual int numLocalStoreVars() const = 0;
  virtual int numExternStoreVars() const = 0;
  virtual int numGlobalExternStoreVars() const = 0;

  virtual int numGlobalNZs() const = 0;
  virtual int numLocalNZs() const = 0;

  virtual const std::vector<int> & rowList_GID() const = 0;
  virtual const std::vector< std::pair<int,int> > & rowList_ExternGID() const = 0;

  virtual const std::vector<int> & rowList_StateGID() const = 0;
  virtual const std::vector< std::pair<int,int> > & rowList_ExternStateGID() const = 0;

  virtual const std::vector<int> & rowList_StoreGID() const = 0;
  virtual const std::vector< std::pair<int,int> > & rowList_ExternStoreGID() const = 0;

  virtual const std::vector<int> & rowList_NumNZs() const = 0;

  virtual const std::vector< std::list<int> > & rowList_ColList() const = 0;

  virtual const std::vector<int> & vnodeGIDVec() const = 0;

  virtual const std::vector<int> & vsrcGIDVec() const = 0;

  //virtual const std::vector<int> & noDCPathGIDVec() const = 0;

  //virtual const std::vector<int> & connToOneTermGIDVec() const = 0;
  
  virtual const std::vector<char> & rowList_VarType() const = 0;
};

#endif


