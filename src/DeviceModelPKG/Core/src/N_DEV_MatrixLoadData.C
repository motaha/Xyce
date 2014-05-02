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
// Filename      : $RCSfile: N_DEV_MatrixLoadData.C,v $
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Eric R. Keiter, SNL
//
// Creation Date : 4/30/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.19 $
//
// Revision Date  : $Date: 2014/02/24 23:49:15 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_DeviceInstance.h>

namespace Xyce {
namespace Device {


// ---------  Other Includes  -----------

// ---------  Helper Classes ------------


//-----------------------------------------------------------------------------
// Function      : MatrixLoadData::MatrixLoadData
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/30/02
//-----------------------------------------------------------------------------
MatrixLoadData::MatrixLoadData() :
  isize(0),
  isizeNumJac(0)
{

}

//-----------------------------------------------------------------------------
// Function      : MatrixLoadData::MatrixLoadData
// Purpose       : Copy Constructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/30/02
//-----------------------------------------------------------------------------
MatrixLoadData::MatrixLoadData( MatrixLoadData const & right )
  :
  isize (right.isize),
  isizeNumJac (right.isizeNumJac),
  cols (right.cols),
  vals (right.vals) ,
  Qvals (right.Qvals) ,
  val_local (right.val_local),
  Qval_local (right.Qval_local),
  col_local (right.col_local),
  row_local (right.row_local),
  internalFlag (right.internalFlag),

  numJac(right.numJac),
  saveJac(right.saveJac),
  devJac(right.devJac),
  diffJac(right.diffJac),
  relJac(right.relJac),

  numJacQ(right.numJacQ),
  saveJacQ(right.saveJacQ),
  devJacQ(right.devJacQ),
  diffJacQ(right.diffJacQ),
  relJacQ(right.relJacQ),

  status(right.status),
  stencil(right.stencil),

  statusQ(right.statusQ),

  saveRHS(right.saveRHS),
  pertRHS(right.pertRHS),
  origRHS(right.origRHS),

  saveQ(right.saveQ),
  pertQ(right.pertQ),
  origQ(right.origQ),

  saveSoln(right.saveSoln),
  pertSoln(right.pertSoln),
  saveCurrSoln(right.saveCurrSoln),
  saveLastState(right.saveLastState),
  saveCurrState(right.saveCurrState),
  saveNextState(right.saveNextState),
  saveStateDerivs(right.saveStateDerivs)
{

}

//-----------------------------------------------------------------------------
// Function      : MatrixLoadData::~MatrixLoadData
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/30/02
//-----------------------------------------------------------------------------
MatrixLoadData::~MatrixLoadData()
{

}

//-----------------------------------------------------------------------------
// Function      : MatrixLoadData::initializeAll
// Purpose       : this function initializes the data structures,
//                 if neccessary.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/02/02
//-----------------------------------------------------------------------------
bool MatrixLoadData::initializeAll (int isizeTmp)
{


  if (vals.size() < isizeTmp)
  {
    isize = isizeTmp;
    vals.resize (isize,0.0);
    Qvals.resize (isize,0.0);
    cols.resize (isize,-1);
  }

  if (row_local.size() < isizeTmp)
  {
    isizeNumJac = isizeTmp;
    row_local.resize(isizeNumJac);
    internalFlag.resize(isizeNumJac);

    val_local.resize(isizeNumJac);
    Qval_local.resize(isizeNumJac);
    col_local.resize(isizeNumJac);
    for (int i=0;i<isizeNumJac;++i)
    {
      val_local[i].val.resize(isizeNumJac,0.0);
      Qval_local[i].val.resize(isizeNumJac,0.0);
      col_local[i].col.resize(isizeNumJac,-1);
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : MatrixLoadData::resizeTestJacSolData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 12/14/06
//-----------------------------------------------------------------------------
void MatrixLoadData::resizeTestJacSolData(int size)
{
  numJac.resize(size);
  saveJac.resize(size);
  devJac.resize(size);
  diffJac.resize(size);
  relJac.resize(size);
  status.resize(size);
  stencil.resize(size);

  // Allocate matrix structures
  for (int i=0 ; i<size ; ++i)
  {
    saveJac[i].resize(size,0.0);
    numJac[i].resize(size,0.0);
    devJac[i].resize(size,0.0);
    diffJac[i].resize(size,0.0);
    relJac[i].resize(size,0.0);
    status[i].resize(size,0);
    stencil[i].resize(size,0);
  }

  saveRHS.resize(size,0.0);
  pertRHS.resize(size,0.0);
  origRHS.resize(size,0.0);
  saveSoln.resize(size,0.0);
  pertSoln.resize(size,0.0);
  saveCurrSoln.resize(size,0.0);
}

//-----------------------------------------------------------------------------
// Function      : MatrixLoadData::resizeTestJacQData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 12/14/06
//-----------------------------------------------------------------------------
void MatrixLoadData::resizeTestJacQData(int size)
{
  numJacQ.resize(size);
  saveJacQ.resize(size);
  devJacQ.resize(size);
  diffJacQ.resize(size);
  relJacQ.resize(size);
  statusQ.resize(size);

  // Allocate matrix structures
  for (int i=0 ; i<size ; ++i)
  {
    saveJacQ[i].resize(size,0.0);
    numJacQ[i].resize(size,0.0);
    devJacQ[i].resize(size,0.0);
    diffJacQ[i].resize(size,0.0);
    relJacQ[i].resize(size,0.0);
    statusQ[i].resize(size,0);
  }

  saveQ.resize(size,0.0);
  pertQ.resize(size,0.0);
  origQ.resize(size,0.0);
}

//-----------------------------------------------------------------------------
// Function      : MatrixLoadData::resizeTestJacStateData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 12/14/06
//-----------------------------------------------------------------------------
void MatrixLoadData::resizeTestJacStateData(int size)
{
  saveLastState.resize(size,0.0);
  saveCurrState.resize(size,0.0);
  saveNextState.resize(size,0.0);
  saveStateDerivs.resize(size,0.0);
}

} // namespace Device
} // namespace Xyce
