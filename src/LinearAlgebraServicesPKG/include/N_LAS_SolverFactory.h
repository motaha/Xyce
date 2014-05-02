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
// Filename       : $RCSfile: N_LAS_SolverFactory.h,v $
//
// Purpose        : Linear Solver Factory
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/12/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_SolverFactory_h
#define Xyce_N_LAS_SolverFactory_h

// ---------- Standard Includes ----------

#include <string>

#include <N_UTL_fwd.h>
#include <N_IO_fwd.h>

class N_LAS_Solver;
class N_LAS_Problem;

//-----------------------------------------------------------------------------
// Class         : N_LAS_SolverFactory
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 08/12/03
//-----------------------------------------------------------------------------
struct N_LAS_SolverFactory
{
  // Creates a new LAS vector
  static N_LAS_Solver * create( N_UTL_OptionBlock & options,
                                N_LAS_Problem & problem,
                                N_IO_CmdParse & commandLine_);
};

#endif

