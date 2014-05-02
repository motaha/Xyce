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
// Filename       : $RCSfile: N_DEV_XyceInterface.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/15/05
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.45 $
//
// Revision Date  : $Date: 2014/02/24 23:49:16 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_DEV_XyceInterface.h>
#include <N_CIR_Xyce.h>
#include <N_DEV_DeviceOptions.h>
#include <N_UTL_BreakPoint.h>

#include <N_PDS_SerialParComm.h>
#include <N_PDS_ParComm.h>
#include <N_PDS_MPIComm.h>

#include <N_IO_CmdParse.h>

#include <N_TIA_TimeIntInfo.h>
#include <N_TIA_TwoLevelError.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : XyceInterface::XyceInterface
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/15/05
//-----------------------------------------------------------------------------
XyceInterface::XyceInterface(
  const DeviceOptions & do1,
  const IO::CmdParse &  cp,
  const std::string &   netlist)
  : devOptions_(do1),
    tmpCmdLine_(cp),
    netlistFileName_(netlist),
    XycePtr_(NULL)
{

}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::~XyceInterface
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/15/05
//-----------------------------------------------------------------------------
XyceInterface::~XyceInterface()

{
  if (XycePtr_ != NULL)
  {
    XycePtr_->finishSolvers ();
    XycePtr_->finalize ();
    delete XycePtr_;
    XycePtr_ = NULL;
  }
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::initialize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/01/2006
//-----------------------------------------------------------------------------
bool XyceInterface::initialize(
  N_PDS_Comm *          comm)
{
  int procID = comm->procID();

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0)
  {
    Xyce::dout() << "In XyceInterface::initialize" << std::endl;
  }
#endif

  // First allocate xyce.
#ifdef Xyce_PARALLEL_MPI
  if( comm->isSerial() )
  {
    XycePtr_ = new N_CIR_Xyce( ((dynamic_cast<N_PDS_SerialParComm*>(comm))->parComm()->mpiComm()->comm()) );
  }
  else
  {
    XycePtr_ = new N_CIR_Xyce( ((dynamic_cast<N_PDS_ParComm*>(comm))->mpiComm()->comm()) );
  }
#else
  XycePtr_ = new N_CIR_Xyce();
#endif

  if (procID == 0)
  {
    // Reset the netlist name in the local copy of the command line arguments:
    tmpCmdLine_.setNetlist(netlistFileName_);
  }

#ifdef Xyce_DEBUG__TWOLEVEL
#ifdef Xyce_PARALLEL_MPI
  Xyce::dout() << "\nproc: " << procID;
#endif
  Xyce::dout() << "Command Line Arguments passed to the inner solve:\n";
  for (int i=0;i<tmpCmdLine_.iargs;++i)
  {
#ifdef Xyce_PARALLEL_MPI
    Xyce::dout() << " proc: " << procID;
#endif
    Xyce::dout() << " cargs["<<i<<"] = " << std::string(tmpCmdLine_.cargs[i]) << std::endl;
  }
#endif

  XycePtr_->initialize(tmpCmdLine_.argc(), tmpCmdLine_.argv());

  XycePtr_->startupSolvers();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::simulateStep
// Purpose       :
// Special Notes : The 'input vector' is a vector of voltages corresponding
//                 to the connected nodes.
//
//                 The 'output vector', is a vector of currents corresponding
//                 to the same nodes.
//
//                 The jacobian is mainly conductances - dI/dV.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/01/2006
//-----------------------------------------------------------------------------
bool XyceInterface::simulateStep
      ( const SolverState & solState,
        const std::map<std::string,double> & inputMap,
        std::vector<double> & outputVector,
        std::vector< std::vector<double> > & jacobian,
        N_TIA_TwoLevelError & tlError
      )
{

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0)
  {
    Xyce::dout() << "In XyceInterface::simulateStep" << std::endl;
  }
#endif

  return XycePtr_->simulateStep
      ( solState,
        inputMap,
        outputVector,
        jacobian,
        tlError
      );
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::finalize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/01/2006
//-----------------------------------------------------------------------------
bool XyceInterface::finalize ()
{
#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0)
  {
    Xyce::dout() << "In XyceInterface::finalize" << std::endl;
  }
#endif


  return true;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::run
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/15/05
//-----------------------------------------------------------------------------
bool XyceInterface::run ()
{
#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0)
  {
    Xyce::dout() << "In XyceInterface::run" << std::endl;
  }
#endif

  int i;
  int iargs_loc = 2;
  char **cargs_loc;

  cargs_loc = new char* [iargs_loc];
  for (i=0;i<iargs_loc;++i)
  {
    cargs_loc[i] = new char[128];
  }

  sprintf(cargs_loc[0],"Xyce ");
  sprintf(cargs_loc[1],"%s",netlistFileName_.c_str());

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0)
  {
    Xyce::dout() << "cargs_loc[0] = " << std::string(cargs_loc[0]) << std::endl;
  }
#endif

  XycePtr_->run (iargs_loc, cargs_loc);

  for (i=0;i<iargs_loc;++i)
  {
    delete [] cargs_loc[i];
  }

  delete [] cargs_loc;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::homotopyStepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void XyceInterface::homotopyStepSuccess
      (const std::vector<std::string> & paramNames,
       const std::vector<double> & paramVals)
{
  XycePtr_->homotopyStepSuccess (paramNames, paramVals);
  return;
}


//-----------------------------------------------------------------------------
// Function      : XyceInterface::homotopyStepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/30/06
//-----------------------------------------------------------------------------
void XyceInterface::homotopyStepFailure ()
{
  XycePtr_->homotopyStepFailure ();
  return;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::stepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void XyceInterface::stepSuccess (int analysis)
{
  XycePtr_->stepSuccess (analysis);
  return;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::stepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void XyceInterface::stepFailure (int analysis)
{
  XycePtr_->stepFailure (analysis);
  return;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::getInitialQnorm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool XyceInterface::getInitialQnorm (N_TIA_TwoLevelError & tle)
{
  bool bsuccess = true;
  if (XycePtr_)
  {
    bsuccess = XycePtr_->getInitialQnorm (tle);
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/14/06
//-----------------------------------------------------------------------------
bool XyceInterface::getBreakPoints
    (std::vector<N_UTL_BreakPoint> &breakPointTimes)
{
  bool bsuccess = true;
  if (XycePtr_)
  {
    bsuccess = XycePtr_->getBreakPoints(breakPointTimes);
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::updateStateArrays
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/18/06
//-----------------------------------------------------------------------------
bool XyceInterface::updateStateArrays ()
{
  bool bsuccess = true;
  if (XycePtr_)
  {
    bsuccess = XycePtr_->updateStateArrays ();
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::startTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
bool XyceInterface::startTimeStep ( const N_TIA_TimeIntInfo & tiInfo )
{
  bool bsuccess = true;
  if (XycePtr_)
  {
    bsuccess = XycePtr_->startTimeStep (tiInfo);
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::setInternalParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/18/06
//-----------------------------------------------------------------------------
bool XyceInterface::setInternalParam (std::string & name, double val)
{
  bool bsuccess = true;
  if (XycePtr_)
  {
    bsuccess = XycePtr_->setInternalParam (name, val);
  }
  return bsuccess;
}

} // namespace Device
} // namespace Xyce
