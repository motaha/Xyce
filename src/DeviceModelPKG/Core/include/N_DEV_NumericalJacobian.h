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
// Filename       : $RCSfile: N_DEV_NumericalJacobian.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 4/30/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.26.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef  N_DEV_NUMERICAL_JACOBIAN_H
#define  N_DEV_NUMERICAL_JACOBIAN_H

// ---------- Standard Includes ----------
#include <list>
#include <vector>
#include <iosfwd>

// ----------   Xyce Includes   ----------
#include <N_UTL_Xyce.h>
#include <N_DEV_fwd.h>

// ----------   Forward Declarations   ----------
class index_pair;

namespace Xyce {
namespace Device {

class colData;
class valData;

//-----------------------------------------------------------------------------
// Class         : NumericalJacobian
// Purpose       :
// Special Notes :
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/30/02
//-----------------------------------------------------------------------------
class NumericalJacobian
{
public:
  NumericalJacobian (
     MatrixLoadData & mlData1,
     const SolverState &ss1,
     const ExternData  &ed1,
     const DeviceOptions & do1);

  NumericalJacobian (const NumericalJacobian & right);

  ~NumericalJacobian ();

  bool testDAEMatrices
  ( DeviceInstance & instance, const std::vector<std::string> & nameVec);

  void loadLocalDAEVectors  (DeviceInstance & instance);

  void printJacobian_
  (  std::ostream &                        os,
     DeviceInstance & instance,
     const std::vector<std::string> & nameVec,
     bool failed);

  void testDebugHead
  ( DeviceInstance & instance,
    const std::vector<std::string> & nameVec,
    int i, double dX);

  void testDebugOut
  ( DeviceInstance & instance,
    const std::vector<std::string> & nameVec,
    int i,int j);

  void testDebugTail
  ( DeviceInstance & instance, const std::vector<std::string> & nameVec);

  void mergeTest ( DeviceInstance & instance, const std::vector<std::string> & nameVec);

public:
  MatrixLoadData & mlData;

  std::vector<int>     & cols;
  std::vector<double>  & vals;
  std::vector<double>  & Qvals;

  std::vector<valData> & val_local;
  std::vector<valData> & Qval_local;
  std::vector<colData> & col_local;
  std::vector<int>     & row_local;
  std::vector<int>     & internalFlag;

  const DeviceOptions & devOptions;
  const SolverState   & solState;
  const ExternData    & extData;

  int maxCols;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::NumericalJacobian N_DEV_NumericalJacobian;

#endif
