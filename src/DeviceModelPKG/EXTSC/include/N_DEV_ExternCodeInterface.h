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
// Revision Number: $Revision: 1.19.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:35 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_ExternCodeInterface_h
#define Xyce_N_DEV_ExternCodeInterface_h

#include <vector>
#include <string>
#include <map>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_UTL_Xyce.h>

// ---------- Forward Declarations ----------
class N_UTL_BreakPoint;
class N_PDS_Comm;
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

    virtual bool initialize(
#ifdef Xyce_PARALLEL_MPI
                             N_PDS_Comm * comm = 0
#endif
                           ) = 0;

    virtual bool simulateStep
      ( const SolverState & solState,
        const map< string, double > & inputMap,
        vector<double> & outputVector,
        vector< vector<double> > & jacobian,
        N_TIA_TwoLevelError & tlError
        ) = 0;

    virtual bool finalize () = 0;
    virtual bool run () = 0;

    virtual void homotopyStepSuccess
      (const vector<string> & paramNames,
       const vector<double> & paramVals) = 0;

    virtual void homotopyStepFailure () = 0;

    virtual void stepSuccess (int analysis) = 0;
    virtual void stepFailure (int analysis) = 0;
    virtual bool getBreakPoints (vector<N_UTL_BreakPoint> &breakPointTimes) = 0;
    virtual bool updateStateArrays () = 0;
    virtual bool startTimeStep ( const N_TIA_TimeIntInfo & tiInfo ) = 0;
    virtual bool setInternalParam (string & name, double val) = 0;

    virtual bool getInitialQnorm (N_TIA_TwoLevelError & tle) = 0;

  protected:
  private:

  public:
  protected:
  private:

};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::ExternCodeInterface N_DEV_ExternCodeInterface;

#endif

