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
// Filename       : $RCSfile: N_IO_MeasureFourier.h,v $
//
// Purpose        : Measure statistics of a simulation variable
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 03/10/2009
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12 $
// Revision Date  : $Date: 2014/02/24 23:49:20 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureFourier_h
#define Xyce_N_IO_MeasureFourier_h

#include <N_IO_MeasureBase.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : Fourier
// Purpose       : Measure statistics of a simulation variable
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 06/05/2013
//-------------------------------------------------------------------------
class Fourier : public Base
{
  public:
    Fourier( const Util::OptionBlock & measureBlock, N_IO_OutputMgr &outputMgr );
    ~Fourier() {};

    void prepareOutputVariables();
    void updateTran( const double circuitTime, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec);
    void updateDC( const std::vector<N_ANP_SweepParam> & dcParamsVec, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec);
    double getMeasureResult();
    std::ostream& printMeasureResult(std::ostream& os);

  private:
    void getLastPeriod_();
    bool interpolateData_();
    void calculateFT_();
    std::string type_;
    int numOutVars_, prdStart_; 
    std::vector<double> outVarValues_, time_, newTime_, newValues_, mag_, phase_, nmag_, nphase_, freq_;
    double period_, lastPrdStart_, thd_;
    bool initialized_, calculated_;

};

} // namespace Measure
} // namespace IO
} // namespace Xyce

typedef Xyce::IO::Measure::Fourier N_IO_MeasureFourier;

#endif
