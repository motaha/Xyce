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
// Revision Number: $Revision: 1.30.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:42:38 $
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

#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_ANP_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>

class N_PDS_Manager;

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : RestartMgr
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/19/01
//-----------------------------------------------------------------------------
class RestartMgr
{

public:

  // Factory to generate singleton of class
  static RestartMgr * factory(CmdParse & cp);

  // Destructor
  ~RestartMgr() {}

private:

  RestartMgr(CmdParse & cp);

  // Copy constructor (private)
  RestartMgr(const RestartMgr & right);

  // Assignment operator
  RestartMgr & operator = (const RestartMgr & right);

  bool operator == (const RestartMgr & right);
  bool operator != (const RestartMgr & right);

public:

  // registration functions:

  bool registerTopology(Topo::Topology * top) { return (topMgrPtr_ = top); }

  bool registerDeviceInterface(Device::DeviceInterface * dev){return (devIntPtr_ = dev);}

  bool registeranaInt(N_ANP_AnalysisInterface * tia) { return (anaIntPtr_ = tia); }

  bool registerParallelServices(N_PDS_Manager * pds) { return (pdsMgrPtr_ = pds); }

  // Method to register the package options manager
  bool registerPkgOptionsMgr( PkgOptionsMgr *pkgOptPtr );

  // Restart flag accessor.
  bool isRestart() const { return restartFlag_; }

  bool registerRestartOptions( const Util::OptionBlock & OB );

  bool getRestartIntervals( double & initialInterval,
                            std::vector<std::pair<double,double> > & intervalPairs );

  bool registerNodePartitioning(std::map<std::string,int> & nodeMap )
  {
    npMap_ = nodeMap;
    return true;
  }

  bool dumpRestartData(const double & time);
  bool restoreRestartData();

  int restartDataSize();

  struct RestartMgr_OptionsReg : public PkgOptionsReg
  {
    RestartMgr_OptionsReg( RestartMgr * mgr )
      : Mgr(mgr)
    {}

    bool operator()( const Util::OptionBlock & options )
    { return Mgr->registerRestartOptions( options ); }

    RestartMgr * Mgr;
  };

private:

  N_PDS_Manager * pdsMgrPtr_;

  Topo::Topology * topMgrPtr_;

  Device::DeviceInterface * devIntPtr_;
  N_ANP_AnalysisInterface * anaIntPtr_;
  // package options manager
  PkgOptionsMgr *pkgOptMgrPtr_;

  bool restartFlag_;
  std::string restartFileName_;
  std::string restartJobName_;

  double initialSaveInterval_;
  std::vector<std::pair<double,double> > saveIntervalPairs_;

  std::map<std::string,int> npMap_;

  bool pack_;

  // command line object
  CmdParse & commandLine_;
};

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::RestartMgr N_IO_RestartMgr;

#endif
