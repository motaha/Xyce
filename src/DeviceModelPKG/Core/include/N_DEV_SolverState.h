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
// Filename       : $RCSfile: N_DEV_SolverState.h,v $
//
// Purpose        : This is a container class for solver information.
//                  It may occasionally contain stuff that isn't strictly
//                  pertaining to the solver state, but that is its primary
//                  intention.
//
//                  In general, stuff that goes into this class should
//                  be stuff needed by more than one device instance type.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.60.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_SolverState_h
#define Xyce_N_DEV_SolverState_h

#include <map>
#include <string>
#include <vector>

#include <N_UTL_Misc.h>
#include <N_NLS_TwoLevelEnum.h>
#include <N_NLS_NonLinInfo.h>
#include <N_TIA_TimeIntInfo.h>

// Forward declarations
class N_UTL_Expression;

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : SolverState
// Purpose       : Container for current solver data.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/25/02
//-----------------------------------------------------------------------------
class SolverState
{

  public:
    SolverState ();

    void InitializeHomotopyBlockSize(int numBlocks);

  protected:

  private:

  public:
    double pdt;                // alpha/dt
    double currTimeStep;
    double lastTimeStep;
    double currTime;
    double finalTime;
    double startingTimeStep;
    double bpTol;
    double acceptedTime; // for habanero
    double currentOrder, usedOrder; // BNB, integration order for 2-level stamp

    vector<double> timePoints;

    // MPDE stuff
    double currFastTime;
    double finalFastTime;
    double currentHoldTime;  // this is the time the slow sources should use to evalulate
                             // their value during an MPDE initial conditon (typically zero, but
                             // it can be anything).
    bool   mpdeOnFlag;  // This indicates the MPDE phase of the problem (ie not initial condition)
    bool   blockAnalysisFlag; // This indicates an MPDE/HB run.  This is true during both IC and MPDE/HB phase.

    map<string,double> global_params;
    vector<N_UTL_Expression> global_expressions;
  std::vector<std::string>global_exp_names;

    // output flag:
    bool forceFinalOutput;

    bool   doubleDCOPEnabled;   // true if taking 2 DCOP steps for PDE sim.
    int    doubleDCOPStep;      // 0 or 1.  (first or second "double" DCOP).

    int    timeStepNumber;

    // The following "ltra*" data structures are used to track time history for
    // the LTRA device. It requires an independent time history because
    // it can be compacted if the user specifies that option. With no
    // compaction ltraTimeStepNumber will be equal to timeStepNumber+1.
    size_t ltraTimeIndex;
    size_t ltraTimeHistorySize;
    bool ltraDoCompact;
    std::vector<double> ltraTimePoints;

    int    newtonIter;
    int    stepLoopIter;
    int    continuationStepNumber;
    bool   firstContinuationParam;
    bool   firstSolveComplete;

    bool initTranFlag;        // true only on very first(t=0) time step.
    bool beginIntegrationFlag;// true if 1st time step out of breakpoint (incl. t=0)

    bool dcopFlag;         // true if we are in a DCOP calculation
                           //  (sweep, tranop or OP)
    bool inputOPFlag;       // true if starting from a previous OP calculation

    bool transientFlag;    // true if transient analysis(even during tranop)
    bool dcsweepFlag;      // true if DC Sweep or OP calculation.
    bool tranopFlag;       // true if in dcop phase of transient sim.
    bool acopFlag;         // true if in acop phase of ac sim.
    bool PDESystemFlag;    // true if circuit includes a PDE device.

    bool locaEnabledFlag;  // true if LOCA is enabled for DCOP.

    bool initJctFlag;  // true if on the first newton step of the
                         // first dcop solve of the first .STEP iteration.

    bool initFixFlag;  // true if DCOP solve, not first iteration *AND*
                       // any device not converged.  Allows "OFF" to be
                       // applied.

    bool sweepSourceResetFlag;
    bool debugTimeFlag;

    TwoLevelNewtonMode twoLevelNewtonCouplingMode;

    // pde device BC homotopy/ two-level newton
    double pdeAlpha;
    bool PDEcontinuationFlag; // if true, continuation is being used.
    int  maxPDEContinuationSteps;
    int  currPDEContinuationStep; // this may become obsolete...
    int  prevPDEContinuationStep;

    bool chargeHomotopy;
    double chargeAlpha;

    // mosfet homotopy variables:
    bool artParameterFlag;
    vector<double> gainScale;
    double nltermScale;

    bool sizeParameterFlag;
    double sizeScale;
    double previousSizeScale;

    // BJT homotopy variables
    bool bjtArtParameterFlag;

    // 2-level info:
    N_TIA_TimeIntInfo tiInfo;
    N_NLS_NonLinInfo nlInfo;

    // analysis options  (for now bools that say what type of "." line was present in netlist)
    bool ACspecified;
    bool MORspecified;
    bool TRANspecified;
    bool DCspecified;
    bool STEPspecified;
    bool OPspecified;
    bool MPDEspecified;
    bool HBspecified;

    friend ostream& operator<<(ostream& os, const SolverState & ss);

  protected:

  private:

};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::SolverState N_DEV_SolverState;

#endif

