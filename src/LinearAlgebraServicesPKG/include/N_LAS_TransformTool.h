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
// Filename       : $RCSfile: N_LAS_TransformTool.h,v $
//
// Purpose        : Constructs Composite Transforms to be applied to
//                  LinearProblems
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 4/2/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.13 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_TransformTool_h
#define Xyce_N_LAS_TransformTool_h

#include <N_UTL_fwd.h>

#include <Teuchos_RCP.hpp>

#include <EpetraExt_Transform_Composite.h>

#include <Epetra_LinearProblem.h>

//-----------------------------------------------------------------------------
// Class         : N_LAS_Transform
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 4/2/03
//-----------------------------------------------------------------------------
struct N_LAS_Transform
: public EpetraExt::Transform_Composite<Epetra_LinearProblem>
{
};

//-----------------------------------------------------------------------------
// Class         : N_LAS_TransformTool
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 4/2/03
//-----------------------------------------------------------------------------
struct N_LAS_TransformTool
{
  // Construction of Composite Transform
  Teuchos::RCP<N_LAS_Transform> operator()( const N_UTL_OptionBlock & options );
};

#endif

