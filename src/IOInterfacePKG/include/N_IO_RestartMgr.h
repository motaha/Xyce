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
// Filename       : $RCSfile: N_IO_RestartMgr.h,v $
//
// Purpose        : Control storage and writeback for restarts
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 7/19/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.21.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_RestartMgr_h
#define Xyce_N_IO_RestartMgr_h

// ---------- Standard Includes ----------
#include <string>
#include <map>
#include <vector>
#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>

// ---------- Forward Declarations ----------
class N_PDS_Manager;

class N_TOP_Topology;

class N_ANP_AnalysisInterface;

class N_IO_CmdParse;

//-----------------------------------------------------------------------------
// Class         : N_IO_RestartMgr
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/19/01
//-----------------------------------------------------------------------------
class N_IO_RestartMgr
{

public:

  // Factory to generate singleton of class
  static N_IO_RestartMgr * factory(N_IO_CmdParse & cp);

  // Destructor
  ~N_IO_RestartMgr() {}

private:

  N_IO_RestartMgr(N_IO_CmdParse & cp);

  // Copy constructor (private)
  N_IO_RestartMgr(const N_IO_RestartMgr & right);

  // Assignment operator
  N_IO_RestartMgr & operator = (const N_IO_RestartMgr & right);

  bool operator == (const N_IO_RestartMgr & right);
  bool operator != (const N_IO_RestartMgr & right);

public:

  // registration functions:

  bool registerTopology(N_TOP_Topology * top) { return (topMgrPtr_ = top); }

  bool registerDeviceInterface(N_DEV_DeviceInterface * dev){return (devIntPtr_ = dev);}

  bool registeranaInt(N_ANP_AnalysisInterface * tia) { return (anaIntPtr_ = tia); }

  bool registerParallelServices(N_PDS_Manager * pds) { return (pdsMgrPtr_ = pds); }

  // Method to register the package options manager
  bool registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr );

  // Restart flag accessor.
  bool isRestart() const { return restartFlag_; }

  bool registerRestartOptions( const N_UTL_OptionBlock & OB );

  bool getRestartIntervals( double & initialInterval,
                            vector< pair<double,double> > & intervalPairs );

  bool registerNodePartitioning( map<string,int> & nodeMap )
  {
    npMap_ = nodeMap;
    return true;
  }

  bool dumpRestartData(const double & time);
  bool restoreRestartData();

  int restartDataSize();

  struct N_IO_RestartMgr_OptionsReg : public N_IO_PkgOptionsReg
  {
    N_IO_RestartMgr_OptionsReg( N_IO_RestartMgr * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->registerRestartOptions( options ); }

    N_IO_RestartMgr * Mgr;
  };

private:

  N_PDS_Manager * pdsMgrPtr_;

  N_TOP_Topology * topMgrPtr_;

  N_DEV_DeviceInterface * devIntPtr_;
  N_ANP_AnalysisInterface * anaIntPtr_;
  // package options manager
  RCP<N_IO_PkgOptionsMgr> pkgOptMgrPtr_;

  bool restartFlag_;
  string restartFileName_;
  string restartJobName_;

  double initialSaveInterval_;
  vector< pair<double,double> > saveIntervalPairs_;

  map<string,int> npMap_;

  bool pack_;

  // command line object
  N_IO_CmdParse & commandLine_;
};

#endif
