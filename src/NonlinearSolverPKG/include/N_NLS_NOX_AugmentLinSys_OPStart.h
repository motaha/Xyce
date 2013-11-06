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
// Filename       : $RCSfile: N_NLS_NOX_AugmentLinSys_OPStart.h,v $
//
// Purpose        : Concrete class for augmenting the Jacobian for
//                  pseudo-transient solves.
//
// Special Notes  :
//
// Creator        : Dave Shirley, PSSI
//
// Creation Date  : 05/08/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:47 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------


#ifndef Xyce_N_NLS_NOX_AugmentLinSys_OPStart_h
#define Xyce_N_NLS_NOX_AugmentLinSys_OPStart_h

#include <set>

#include <N_UTL_fwd.h>

#include "Teuchos_RefCountPtr.hpp"
#include "N_NLS_NOX_AugmentLinSys.h"
#include "N_PDS_ParMap.h"
#ifdef Xyce_PARALLEL_MPI
#include "N_PDS_ParComm.h"
#endif


class Epetra_MapColoring;

//-----------------------------------------------------------------------------
// Class         : N_NLS::NOX::AugmentLinSysOPStart
// Purpose       :
// Creator       : Roger Pawlowski, SNL, 9233
// Creation Date :
//-----------------------------------------------------------------------------
namespace N_NLS_NOX {

  class AugmentLinSysOPStart : public N_NLS_NOX::AugmentLinSys {

public:

  //! Ctor.
#ifdef Xyce_PARALLEL_MPI
  AugmentLinSysOPStart(Xyce::NodeNamePairMap &, Xyce::NodeNamePairMap &,
      N_PDS_Comm *);
#else
  AugmentLinSysOPStart(Xyce::NodeNamePairMap &, Xyce::NodeNamePairMap &);
#endif

  //! Dtor.
  ~AugmentLinSysOPStart();

  void setProgressVariable(double dummy) {return;}

  void augmentResidual(const N_LAS_Vector * solution, N_LAS_Vector * residual_vector);

  void augmentJacobian(N_LAS_Matrix * jacobian);

 private:

  //! map of specified variables
  Xyce::NodeNamePairMap & op_;
  Xyce::NodeNamePairMap & allNodes_;

  bool skipSet;
  set<int> skipLID;
  set<int> skipGID;

  N_LAS_Vector* residualPtr_;
  const N_LAS_Vector* solutionPtr_;
  int rSize_;
#ifdef Xyce_PARALLEL_MPI
  N_PDS_ParMap * pmap_;
  N_PDS_Comm * pdsCommPtr_;
#endif
};

}

#endif

