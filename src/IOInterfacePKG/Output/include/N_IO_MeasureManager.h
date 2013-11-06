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
// Filename      : $RCSfile: N_IO_MeasureManager.h,v $
//
// Purpose       : This file is a class to manage measure statements in a sim.
//
// Special Notes :
//
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef  Xyce_N_IO_MeasureManager_H
#define Xyce_N_IO_MeasureManager_H

// ---------- Standard Includes ----------

#include <list>
#include <string>
#include <iostream>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>
#include <N_UTL_Xyce.h>
#include <N_IO_fwd.h>
#include <N_UTL_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_LAS_Vector.h>
#include <N_ANP_SweepParam.h>


// ---------- Forward Declarations ----------
class N_IO_MeasureBase;

//-----------------------------------------------------------------------------
// Class         : N_IO_MeasureManager
// Purpose       : This is a manager class for handling measure statements
//                 in a simulation
// Special Notes :
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
class N_IO_MeasureManager
{

public:

  // Default constructor.
  N_IO_MeasureManager( N_IO_OutputMgr & outputMgr );

  // Destructor
  ~N_IO_MeasureManager();

  // Return true if .measure analysis is being performed on any variables.
  bool isMeasureActive() { return (!allMeasuresList_.empty()); }

  // add .measure line from netlist to list of things to measure
  bool addMeasure(const N_UTL_OptionBlock & measureLine);

  // Called during the simulation to update the measure objects held by this class
  // To keep things obvious, use one call for tran and another for DC
  void updateTranMeasures( const double circuitTime, RCP< N_LAS_Vector > solnVecRCP);
  void updateDcMeasures( const vector<N_ANP_SweepParam> & dcParamsVec, RCP< N_LAS_Vector > solnVecRCP);

  void outputResults(std::ostream& outputStream, bool printHeader=true );

  void getMeasureValue (const string &name, double &value, bool &found);

private:
  N_IO_OutputMgr &     outputManager_;

  // package options manager
  RCP<N_IO_PkgOptionsMgr> pkgOptMgrPtr_;

  list<RCP<N_IO_MeasureBase> > allMeasuresList_;
  list<RCP<N_IO_MeasureBase> > activeMeasuresList_;

};

#endif  // Xyce_N_IO_MeasureManager_H
