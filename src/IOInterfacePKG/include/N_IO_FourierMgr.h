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
// Filename      : $RCSfile: N_IO_FourierMgr.h,v $
//
// Purpose       : This file is a class to manage measure statements in a sim.
//
// Special Notes :
//
// Creator       : Heidi Thornquist, SNL, Electrical and Microsystem Modeling
//
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1.2.6 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef  Xyce_N_IO_FourierMgr_H
#define Xyce_N_IO_FourierMgr_H

// ---------- Standard Includes ----------

#include <list>
#include <string>
#include <iostream>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
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

//-----------------------------------------------------------------------------
// Class         : N_IO_FourierMgr
// Purpose       : This is a manager class for handling .four and .fft statements
//                 in a simulation
// Special Notes :
// Creator       : Heidi Thornquist, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
class N_IO_FourierMgr
{

public:

  // Default constructor.
  N_IO_FourierMgr( N_IO_OutputMgr & outputMgr );

  // Destructor
  ~N_IO_FourierMgr();

  // Return true if Fourier analysis is being performed on any variables.
  bool isFourierActive() { return (!freqVector_.empty() && !time_.empty()); }

  // add .four or .fft line from netlist to list of things to perform analysis on.
  bool addFourierAnalysis( const N_UTL_OptionBlock & fourierLine );

  // Called during the simulation to update the fourier objects held by this class
  void updateFourierData( const double circuitTime, RCP< N_LAS_Vector > solnVecRCP );

  void outputResults( std::ostream& outputStream );

private:

  void getLastPeriod_();

  bool interpolateData_();

  void calculateFT_();

  std::ostream& printResult_( std::ostream& os );

  N_IO_OutputMgr &     outputManager_;

  // package options manager
  RCP<N_IO_PkgOptionsMgr> pkgOptMgrPtr_;

  // Store frequencies and solution variables
  int numFreq_, gridSize_;
  bool calculated_;
  std::vector<int> outputVarsPtr_;
  std::vector<double> time_;
  std::vector<double> freqVector_;
  std::vector<double> outputVarsValues_;
  std::vector<std::string> names_;
  list<N_UTL_Param> outputVars_;
  vector< list<N_UTL_Param>::iterator > depSolVarIterVector_;
  vector<int> prdStart_; 
  vector<double> lastPrdStart_;
  vector<double> newTime_, newValues_, mag_, phase_, nmag_, nphase_, freq_, thd_;

};

#endif  // Xyce_N_IO_FourierMgr_H
