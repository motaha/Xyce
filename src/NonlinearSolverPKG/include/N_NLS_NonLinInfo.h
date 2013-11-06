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
// Filename      : $RCSfile: N_NLS_NonLinInfo.h,v $
// Purpose       : 
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 2/11/07
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.4.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:47 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NON_LIN_INFO_H
#define Xyce_N_NLS_NON_LIN_INFO_H

// ---------- Standard Declarations ----------
#include<N_NLS_TwoLevelEnum.h>

// ---------- Forward Declarations ----------

//-----------------------------------------------------------------------------
// Class         : N_NLS_NonLinInfo
// Purpose       : 
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 2/11/07
//-----------------------------------------------------------------------------
class N_NLS_NonLinInfo
{
public:
  N_NLS_NonLinInfo():
  newtonIter(0),
  twoLevelNewtonCouplingMode (FULL_PROBLEM),
  locaFlag(false),
  continuationStep(0),
  firstContinuationParam(false),
  firstSolveComplete(false)
  {};

  virtual ~N_NLS_NonLinInfo() {};

  int newtonIter;
  TwoLevelNewtonMode twoLevelNewtonCouplingMode;

  // LOCA/homotopy related stuff.
  bool locaFlag;
  int continuationStep;
  bool firstContinuationParam;
  bool firstSolveComplete;


};

#endif // Xyce_N_NLS_NON_LIN_INFO_H


