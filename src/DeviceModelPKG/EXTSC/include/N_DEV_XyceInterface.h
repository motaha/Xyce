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
// Filename       : $RCSfile: N_DEV_XyceInterface.h,v $
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
// Revision Number: $Revision: 1.29 $
//
// Revision Date  : $Date: 2014/02/24 23:49:16 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_XyceInterface_h
#define Xyce_N_DEV_XyceInterface_h

#include <string>
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_ExternCodeInterface.h>
#include <N_UTL_Misc.h>
#include <N_UTL_fwd.h>
#include <N_CIR_fwd.h>
#include <N_IO_CmdParse.h>

class N_TIA_TimeIntInfo;
class N_TIA_TwoLevelError;

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : XyceInterface
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, Parallel Computational Sciences
// Creation Date : 04/15/05
//-----------------------------------------------------------------------------
class XyceInterface : public ExternCodeInterface
{
  public:
    XyceInterface(
      const DeviceOptions & do1,
      const IO::CmdParse & cp,
      const std::string & netlist);

    virtual ~XyceInterface();

  private:
    XyceInterface (const XyceInterface &right);

  public:
    bool initialize(N_PDS_Comm * comm = 0);

    bool simulateStep
      ( const SolverState & solState,
        const std::map<std::string,double> & inputMap,
        std::vector<double> & outputVector,
        std::vector< std::vector<double> > & jacobian,
        N_TIA_TwoLevelError & tlError
      );

    bool finalize ();
    bool run ();

    void homotopyStepSuccess
      (const std::vector<std::string> & paramNames,
       const std::vector<double> & paramVals);

    void homotopyStepFailure ();

    void stepSuccess (int analysis);
    void stepFailure (int analysis);

    bool getBreakPoints (std::vector<N_UTL_BreakPoint> &breakPointTimes);

    bool updateStateArrays ();
    bool startTimeStep ( const N_TIA_TimeIntInfo & tiInfo );
    bool setInternalParam (std::string & name, double val);

    bool getInitialQnorm (N_TIA_TwoLevelError & tle);


  private:

    std::string                 netlistFileName_;
    Circuit::Simulator *        XycePtr_;
    const DeviceOptions &       devOptions_;
    IO::CmdParse                tmpCmdLine_;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::XyceInterface N_DEV_XyceInterface;

#endif

