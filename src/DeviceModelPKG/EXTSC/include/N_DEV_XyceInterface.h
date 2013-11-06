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
// Revision Number: $Revision: 1.23.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:35 $
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
#include <N_IO_CmdParse.h>

// ---------- Forward Declarations ----------
class N_CIR_Xyce;
class N_UTL_BreakPoint;

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
    XyceInterface (DeviceOptions & do1,
        N_IO_CmdParse & cp,
        string & netlist);

    virtual ~XyceInterface();

    bool initialize(
#ifdef Xyce_PARALLEL_MPI
                    N_PDS_Comm * comm = 0
#endif
                   );

    bool simulateStep
      ( const SolverState & solState,
        const map<string,double> & inputMap,
        vector<double> & outputVector,
        vector< vector<double> > & jacobian,
        N_TIA_TwoLevelError & tlError
      );

    bool finalize ();
    bool run ();

    void homotopyStepSuccess
      (const vector<string> & paramNames,
       const vector<double> & paramVals);

    void homotopyStepFailure ();

    void stepSuccess (int analysis);
    void stepFailure (int analysis);

    bool getBreakPoints (vector<N_UTL_BreakPoint> &breakPointTimes);

    bool updateStateArrays ();
    bool startTimeStep ( const N_TIA_TimeIntInfo & tiInfo );
    bool setInternalParam (string & name, double val);

    bool getInitialQnorm (N_TIA_TwoLevelError & tle);

  protected:
  private:
    XyceInterface (const XyceInterface &right);

  public:
  protected:
  private:

    string netlistFileName_;
    N_CIR_Xyce * XycePtr_;
    DeviceOptions & devOptions_;
    N_IO_CmdParse & commandLine_;
    N_IO_CmdParse tmpCmdLine_;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::XyceInterface N_DEV_XyceInterface;

#endif

