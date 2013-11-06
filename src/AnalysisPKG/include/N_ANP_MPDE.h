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
// Filename       : $RCSfile: N_ANP_MPDE.h,v $
//
// Purpose        : MPDE analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Todd Coffey, 1414, Ting Mei 1437
//
// Creation Date  : 07/23/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.6.2.2 $
// Revision Date  : $Date: 2013/10/03 17:23:30 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_MPDE_h
#define Xyce_N_ANP_MPDE_h

// ----------   Xyce Includes   ----------
#include <N_ANP_AnalysisBase.h>

// ---------- Forward Declarations ----------
class N_ANP_AnalysisManager;


//-------------------------------------------------------------------------
// Class         : N_ANP_MPDE
// Purpose       : MPDE analysis class
// Special Notes : 
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class N_ANP_MPDE : public N_ANP_AnalysisBase 
{
public:
  N_ANP_MPDE( N_ANP_AnalysisManager * anaManagerPtr );

  virtual ~N_ANP_MPDE() {};
  
  bool run();
  bool init();
  bool loopProcess();
  bool processSuccessfulDCOP();
  bool processFailedDCOP();
  bool processSuccessfulStep();
  bool processFailedStep();
  bool finish();
  bool resetForStepAnalysis();
  
  bool finalVerboseOutput();

  void printStepHeader();
  void printProgress();

private:
 
  void takeAnIntegrationStep_();

  // a flag to indicate of the simulation is paused
  bool isPaused; 

  double startDCOPtime, endTRANtime; // startTRANtime
  // Timing/loop count info
  RefCountPtr<N_DEV_DeviceInterface> devInterfacePtr_;
  RefCountPtr<N_TOP_Topology> topoMgrPtr_;
  RefCountPtr<N_LOA_NonlinearEquationLoader> nonlinearEquationLoaderPtr_;
  RefCountPtr<N_LAS_Builder> appBuilderPtr_;
};

#endif // Xyce_N_ANP_MPDE_h

