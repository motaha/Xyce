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
// Filename       : $RCSfile: N_DEV_ExternCodeInterface.h,v $
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
// Revision Number: $Revision: 1.26 $
//
// Revision Date  : $Date: 2014/02/24 23:49:16 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_ExternCodeInterface_h
#define Xyce_N_DEV_ExternCodeInterface_h

#include <vector>
#include <string>
#include <map>

#include <N_UTL_Xyce.h>
#include <N_DEV_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>

class N_TIA_TimeIntInfo;
class N_TIA_TwoLevelError;

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : N_DEV_ExternCodeInterface
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, Parallel Computational Sciences
// Creation Date : 04/15/05
//-----------------------------------------------------------------------------
class ExternCodeInterface
{
  public:
    ExternCodeInterface ();

    ExternCodeInterface (const ExternCodeInterface &right);
    virtual ~ExternCodeInterface();

    virtual bool initialize(N_PDS_Comm * comm = 0) = 0;

    virtual bool simulateStep
      ( const SolverState & solState,
        const std::map< std::string, double > & inputMap,
        std::vector<double> & outputVector,
        std::vector< std::vector<double> > & jacobian,
        N_TIA_TwoLevelError & tlError
        ) = 0;

    virtual bool finalize () = 0;
    virtual bool run () = 0;

    virtual void homotopyStepSuccess
      (const std::vector<std::string> & paramNames,
       const std::vector<double> & paramVals) = 0;

    virtual void homotopyStepFailure () = 0;

    virtual void stepSuccess (int analysis) = 0;
    virtual void stepFailure (int analysis) = 0;
    virtual bool getBreakPoints (std::vector<N_UTL_BreakPoint> &breakPointTimes) = 0;
    virtual bool updateStateArrays () = 0;
    virtual bool startTimeStep ( const N_TIA_TimeIntInfo & tiInfo ) = 0;
    virtual bool setInternalParam (std::string & name, double val) = 0;

    virtual bool getInitialQnorm (N_TIA_TwoLevelError & tle) = 0;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::ExternCodeInterface N_DEV_ExternCodeInterface;

#endif

