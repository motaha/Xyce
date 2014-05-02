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
// Filename       : $RCSfile: N_IO_Op.C,v $
//
// Purpose        : Provide tools for accessing output data in parallel or 
//                  serial
//
// Special Notes  :
//
// Creator        : David Baur
//
// Creation Date  : 11/15/2013
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.13.2.1 $
//
// Revision Date  : $Date: 2014/02/25 22:06:24 $
//
// Current Owner  : $Author: rlschie $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_fwd.h>
#include <N_IO_Op.h>
#include <N_DEV_DeviceInterface.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceSensitivities.h>
#include <N_IO_MeasureBase.h>
#include <N_UTL_Algorithm.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : getOpMap
// Purpose       : Return a creator function for a given operator type.
// Special Notes : Creates a static map between types and create functions on 
//                 first call.
// Scope         : namespace Xyce::IO
// Creator       : David Baur, Raytheon
// Creation Date : 11/21/2013
//-----------------------------------------------------------------------------
CreateFunction getOpMap(int op_type)
{
  static OpMap op_map;
  if (op_map.empty()) {
    op_map[Util::INDEX] = OutputMgrCurrentOutputterIndexOp::ReduceOp::create;
    op_map[Util::CONSTANT] = ConstantOp::ReduceOp::create;
    op_map[Util::TEMPERATURE] = OutputMgrTemperatureOp::ReduceOp::create;
    op_map[Util::FREQUENCY] = OutputMgrFrequencyOp::ReduceOp::create;
    op_map[Util::TIME_VAR] = OutputMgrTimeOp::ReduceOp::create;
    op_map[Util::STEP_SWEEP_VAR] = OutputMgrStepSweepOp::ReduceOp::create;
    op_map[Util::DC_SWEEP_VAR] = OutputMgrDCSweepOp::ReduceOp::create;
    op_map[Util::DEVICE_PARAMETER] = DeviceMgrParameterOp::ReduceOp::create;
    op_map[Util::DEVICE_ENTITY_PARAMETER] = DeviceEntityParameterOp::ReduceOp::create;
    op_map[Util::GLOBAL_PARAMETER] = DeviceMgrGlobalParameterOp::ReduceOp::create;
    op_map[Util::OBJECTIVE_FUNCTION] = ObjectiveOp::ReduceOp::create;
    op_map[Util::MEASURE_FUNCTION] = MeasureOp::ReduceOp::create;
    op_map[Util::SOLUTION_VAR] = SolutionOp::ReduceOp::create;
    op_map[Util::SOLUTION_VAR_REAL] = SolutionRealOp::ReduceOp::create;
    op_map[Util::SOLUTION_VAR_IMAG] = SolutionImaginaryOp::ReduceOp::create;
    op_map[Util::SOLUTION_VAR_MAG] = SolutionMagnitudeOp::ReduceOp::create;
    op_map[Util::SOLUTION_VAR_PHASE] = SolutionPhaseOp::ReduceOp::create;
    op_map[Util::SOLUTION_VAR_DB] = SolutionDecibelsOp::ReduceOp::create;
    op_map[Util::VOLTAGE_DIFFERENCE] = VoltageDifferenceOp::ReduceOp::create;
    op_map[Util::VOLTAGE_DIFFERENCE_REAL] = VoltageDifferenceRealOp::ReduceOp::create;
    op_map[Util::VOLTAGE_DIFFERENCE_IMAG] = VoltageDifferenceImaginaryOp::ReduceOp::create;
    op_map[Util::VOLTAGE_DIFFERENCE_MAG] = VoltageDifferenceMagnitudeOp::ReduceOp::create;
    op_map[Util::VOLTAGE_DIFFERENCE_PHASE] = VoltageDifferencePhaseOp::ReduceOp::create;
    op_map[Util::VOLTAGE_DIFFERENCE_DB] = VoltageDifferenceDecibelsOp::ReduceOp::create;
    op_map[Util::STATE_VAR] = StateOp::ReduceOp::create;
    op_map[Util::STORE_VAR] = StoreOp::ReduceOp::create;
  }
  return op_map[op_type];
}

//-----------------------------------------------------------------------------
// Function      : ReduceSum::reduce
// Purpose       : perform an MPI sum reduction operation across all processors.
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
ReduceSum::reduce(Parallel::Machine comm, complex result)
{
  Parallel::AllReduce(comm, MPI_SUM, &result, 1);

  return result;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrCurrentOutputterIndexOp::get
// Purpose       : get the "index" from the currently active outputter
// Special Notes : the index is basically the line number of output, starting
//                 at zero for the first line
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
OutputMgrCurrentOutputterIndexOp::get(const OutputMgrCurrentOutputterIndexOp &op)
{
  return op.outputMgr_.getCurrentOutputterIndex();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrTimeOp::get
// Purpose       : get the current simulation time being output
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
OutputMgrTimeOp::get(const OutputMgrTimeOp &op)
{
  return op.outputMgr_.getTime();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrFrequencyOp::get
// Purpose       : get the current frequency being output
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
OutputMgrFrequencyOp::get(const OutputMgrFrequencyOp &op)
{
  return op.outputMgr_.getFrequency();
}


//-----------------------------------------------------------------------------
// Function      : OutputMgrTemperatureOp::get
// Purpose       : get the current simulation temperature
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
OutputMgrTemperatureOp::get(const OutputMgrTemperatureOp &op)
{
  return op.outputMgr_.getTemperature();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrStepSweepOp::get
// Purpose       : get the current value of the step parameter being swept.
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
OutputMgrStepSweepOp::get(const OutputMgrStepSweepOp &op)
{
  return op.outputMgr_.getStepSweep(op.index_);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrDCSweepOp::get
// Purpose       : get the current value of the DC voltage being swept.
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
OutputMgrDCSweepOp::get(const OutputMgrDCSweepOp &op)
{
  return op.outputMgr_.getDCSweep(op.index_);
}

//-----------------------------------------------------------------------------
// Function      : ObjectiveOp::get
// Purpose       : get the current value of the objective function
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
ObjectiveOp::get(const ObjectiveOp &op)
{
  complex result(0.0, 0.0);

  if (op.objective_.var1.empty() && op.objective_.var2.empty())
  {
    // saving single values(potentially all values) as there isn't
    // any external data to tell us what needs to be saved
    result = complex(op.objective_.save(op.realSolutionVector_, op.stateVector_, op.storeVector_), 0.0);
  }
  else
  {
    // use user supplied external data to save only the simulation
    // results we really need.
    double v1 = 0.0;
    double v2 = 0.0;
    // Util::Param param;
    // if (!objective.var1.empty())
    // {
    //   param.setTag(objective.var1);
    //   v1 = getPrintValue(outputMgr_, param, realSolutionVector_, stateVector_, storeVector_);
    // }
    // if (!objective.var2.empty())
    // {
    //   param.setTag(objective.var2);
    //   v2 = getPrintValue(outputMgr_, param, realSolutionVector_);
    // }
    result = complex(op.objective_.save(v1, v2, op.realSolutionVector_, op.stateVector_, op.storeVector_), 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgrGlobalParameterOp::get
// Purpose       : get the current value of a global param from the device
//                 package
// Special Notes : It is inappropriate for the Device package to be in charge
//                 of global params, but that's where they are right now.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
DeviceMgrGlobalParameterOp::get(const DeviceMgrGlobalParameterOp &op)
{
  return op.deviceInterface_.getGlobalPar(op.deviceParameterName_);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgrParameterOp::get
// Purpose       : get the current value of a device parameter from the device
//                 package
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
DeviceMgrParameterOp::get(const DeviceMgrParameterOp &op)
{
  return op.deviceInterface_.getParamNoReduce(op.deviceParameterName_);
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntityParameterOp::get
// Purpose       : get the current value of a device parameter from a device
//                 entity
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
DeviceEntityParameterOp::get(const DeviceEntityParameterOp &op)
{
  double result;

  const_cast<Device::DeviceEntity &>(op.deviceEntity_).getParam(op.deviceParameterName_, result);

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionOp::get
// Purpose       : get the current value of a solution vector element
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionOp::get(const SolutionOp &op)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1) {
    result = complex((*op.realSolutionVector_)[op.index_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionRealOp::get
// Purpose       : get a solution variable in preparation for computing real
//                 part
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionRealOp::get(const SolutionRealOp &op)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1) {
    result = complex((*op.realSolutionVector_)[op.index_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionRealOp::eval
// Purpose       : take the real part of a complex number
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionRealOp::eval(complex result)
{
  return result.real();
}


//-----------------------------------------------------------------------------
// Function      : SolutionImaginaryOp::get
// Purpose       : get a solution variable in preparation for computing 
//                 imaginary part
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionImaginaryOp::get(const SolutionImaginaryOp &op)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1) {
    result = complex((*op.realSolutionVector_)[op.index_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionImaginaryOp::eval
// Purpose       : take the imaginary part of a complex number
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionImaginaryOp::eval(complex result)
{
  return result.imag();
}


//-----------------------------------------------------------------------------
// Function      : SolutionMagnitudeOp::get
// Purpose       : get the magnitude of a solution vector element
// Special Notes : Actually just gets the solution variable.  Does not
//                 take the magnitude.  eval does that.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionMagnitudeOp::get(const SolutionMagnitudeOp &op)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex((*op.realSolutionVector_)[op.index_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionMagnitudeOp::eval
// Purpose       : take the magnitude of a solution vector element
// Special Notes : Actually just takes the magnitude of a given complex 
//                 value, does NOT access the solution vector itself. 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionMagnitudeOp::eval(complex result)
{
  return std::abs(result);
}


//-----------------------------------------------------------------------------
// Function      : SolutionPhaseOp::get
// Purpose       : get the phase of a solution vector element
// Special Notes : Actually just gets the solution variable.  Does not
//                 take the phase.  eval does that.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionPhaseOp::get(const SolutionPhaseOp &op)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex((*op.realSolutionVector_)[op.index_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionPhaseOp::eval
// Purpose       : compute the phase of a solution vector element
// Special Notes : Actually just computes the phase of a given complex 
//                 value, does NOT access the solution vector itself. 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionPhaseOp::eval(complex result)
{
  return std::arg(result);
}


//-----------------------------------------------------------------------------
// Function      : SolutionDecibelsOp::get
// Purpose       : get the magnitude (in dB) of a solution vector element
// Special Notes : Actually just gets the solution variable.  Does not
//                 find the magnitude.  eval does that.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionDecibelsOp::get(const SolutionDecibelsOp &op)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex((*op.realSolutionVector_)[op.index_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionDecibelsOp::eval
// Purpose       : compute the magnitude (in dB) of a solution vector element
// Special Notes : Actually just computes the magnitude of a given complex 
//                 value, does NOT access the solution vector itself. 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionDecibelsOp::eval(complex result)
{
  return 20.0*std::log10(abs(result));
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like V(A,B)
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceOp::get(const VoltageDifferenceOp &op)
{
  complex result(0.0, 0.0);

  if (op.index1_ != -1)
  {
    result = complex((*op.realSolutionVector_)[op.index1_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index1_]);
  }
  if (op.index2_ != -1)
  {
    result -= complex((*op.realSolutionVector_)[op.index2_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index2_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceRealOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like V(A,B)
// Special Notes : Computes the difference, but does not take the real part.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceRealOp::get(const VoltageDifferenceRealOp &op)
{
  complex result(0.0, 0.0);

  if (op.index1_ != -1)
  {
    result = complex((*op.realSolutionVector_)[op.index1_], op.realSolutionVector_ == 0 ? 0.0 : (*op.realSolutionVector_)[op.index1_]);
  }

  if (op.index2_ != -1)
  {
    result -= complex((*op.realSolutionVector_)[op.index2_], op.realSolutionVector_ == 0 ? 0.0 : (*op.realSolutionVector_)[op.index2_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceRealOp::eval
// Purpose       : take the real part of a voltage difference
// Special Notes : must "get" the difference first
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceRealOp::eval(complex result)
{
  return result.real();
}


//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceImagOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like V(A,B)
// Special Notes : Computes the difference, but does not take the imaginary 
//                 part.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceImaginaryOp::get(const VoltageDifferenceImaginaryOp &op)
{
  complex result(0.0, 0.0);

  if (op.index1_ != -1)
  {
    result = complex((*op.realSolutionVector_)[op.index1_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index1_]);
  }

  if (op.index2_ != -1)
  {
    result -= complex((*op.realSolutionVector_)[op.index2_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index2_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceImagOp::get
// Purpose       : get the imaginary part of a voltage difference
// Special Notes : Must "get" the difference first.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceImaginaryOp::eval(complex result)
{
  return result.imag();
}


//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceMagnitudeOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like V(A,B)
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceMagnitudeOp::get(const VoltageDifferenceMagnitudeOp &op)
{
  complex result(0.0, 0.0);

  if (op.index1_ != -1)
  {
    result = complex((*op.realSolutionVector_)[op.index1_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index1_]);
  }
  if (op.index2_ != -1)
  {
    result -= complex((*op.realSolutionVector_)[op.index2_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index2_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceMagnitudeOp::eval
// Purpose       : Compute the magnitude of a voltage difference
// Special Notes : Must "get" difference first
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceMagnitudeOp::eval(complex result)
{
  return std::abs(result);
}


//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceMagnitudeOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like V(A,B)
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferencePhaseOp::get(const VoltageDifferencePhaseOp &op)
{
  complex result(0.0, 0.0);

  if (op.index1_ != -1)
  {
    result = complex((*op.realSolutionVector_)[op.index1_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index1_]);
  }
  if (op.index2_ != -1)
  {
    result -= complex((*op.realSolutionVector_)[op.index2_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index2_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferencePhaseOp::eval
// Purpose       : Compute the phase of a voltage difference
// Special Notes : Must "get" difference first
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferencePhaseOp::eval(complex result)
{
  return std::arg(result);
}


//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceMagnitudeOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like V(A,B)
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceDecibelsOp::get(const VoltageDifferenceDecibelsOp &op)
{
  complex result(0.0, 0.0);

  if (op.index1_ != -1)
  {
    result = complex((*op.realSolutionVector_)[op.index1_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index1_]);
  }
  if (op.index2_ != -1)
  {
    result -= complex((*op.realSolutionVector_)[op.index2_], op.imaginarySolutionVector_ == 0 ? 0.0 : (*op.imaginarySolutionVector_)[op.index2_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceMagnitudeOp::eval
// Purpose       : Compute the magnitude of a voltage difference (in dB)
// Special Notes : Must "get" difference first
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceDecibelsOp::eval(complex result)
{
  return 20.0*std::log10(std::abs(result));
}

//-----------------------------------------------------------------------------
// Function      : StateOp::get
// Purpose       : Get a value out of the state vector.
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
StateOp::get(const StateOp &op)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op.stateVector_ == 0 ? 0.0 : (*op.stateVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StoreOp::get
// Purpose       : Get a value out of the Store vector.
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
StoreOp::get(const StoreOp &op)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op.storeVector_ == 0 ? 0.0 : (*op.storeVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : MeasureOp::get
// Purpose       : Get a .measure result
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
MeasureOp::get(const MeasureOp &op)
{
  complex result(const_cast<Measure::Base &>(op.measure_).getMeasureResult(), 0.0);

  return result;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionOp::ExpressionOp
// Purpose       : Constructor for expression Op
// Special Notes : Takes pre-constructed expression as second argument
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
ExpressionOp::ExpressionOp(const std::string &name, Util::Expression &expression, const OutputMgr &output_manager)
  : Base(name),
    expressionData_(expression, const_cast<OutputMgr &>(output_manager)),
    outputMgr_(output_manager)
{
  init();
}

//-----------------------------------------------------------------------------
// Function      : ExpressionOp::ExpressionOp
// Purpose       : Constructor for expression Op
// Special Notes : Takes string as second argument.  expressionData_ 
//                 constructor will process the string into an expression
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
ExpressionOp::ExpressionOp(const std::string &name, const std::string &expression, const OutputMgr &output_manager)
  : Base(name),
    expressionData_(expression, const_cast<OutputMgr &>(output_manager)),
    outputMgr_(output_manager)
{
  init();
}

//-----------------------------------------------------------------------------
// Function      : ExpressionOp::init
// Purpose       : initialize an expression
// Special Notes : runs the ExpressionData::setup method to resolve
//                 symbols
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void ExpressionOp::init()
{
  int numUnresolvedSymbols = expressionData_.setup();

#ifdef Xyce_PARALLEL_MPI
  outputMgr_.getCommPtr()->barrier();
  int minNumUnresolvedSymbols=0;
  outputMgr_.getCommPtr()->minAll(&numUnresolvedSymbols, &minNumUnresolvedSymbols, 1);
  numUnresolvedSymbols = minNumUnresolvedSymbols;
#endif

  if (numUnresolvedSymbols > 0 )
  {
    Report::UserFatal0() << "Can't resolve all symbols in expression variable "
                         << expressionData_.getExpression();
  }
  // set numUnresolvedStringsChecked flag to true so we only do this extra work once.
  expressionData_.setUnresolvedStringsChecked(true);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionOp::get
// Purpose       : evaluate an expression
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
ExpressionOp::get(const ExpressionOp &op)
{
  complex result(op.expressionData_.evaluate(op.realSolutionVector_, op.stateVector_, op.storeVector_, op.imaginarySolutionVector_), 0.0);

  return result;
}



//-----------------------------------------------------------------------------
// Function      : makeOp
// Purpose       : given a parameter list iterator, construct an Op for
//                 the item so defined
// Special Notes : 
// Scope         : global (Xyce::IO namespace, but not in any header)
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
Util::Operator *makeOp(
  const OutputMgr &                     output_manager,
  Util::ParameterList::iterator &       param_it)
{
  Util::Operator *new_op = 0;

  std::string rank;
  {
    std::ostringstream oss;
    oss << output_manager.getProcID();
    rank = oss.str();
  }

  // var_type should be I, V, N, B, H as in I(name), V(name), or V(name, name)
  std::string var_type = (*param_it).tag();

  if (var_type == "TEMP")
  {
    new_op  = new OutputMgrTemperatureOp(var_type, output_manager);
  }

  else if (var_type == "TIME")
  {
    new_op = new OutputMgrTimeOp(var_type, output_manager);
  }

  else if (var_type == "FREQUENCY")
  {
    new_op = new OutputMgrFrequencyOp(var_type, output_manager);
  }

  else if (var_type == "INDEX")
  {
    new_op = new OutputMgrCurrentOutputterIndexOp(var_type, output_manager);
  }

  else if (var_type == "GLOBAL_PARAMETER")
  {
    new_op = new DeviceMgrGlobalParameterOp((*param_it).stringValue(), output_manager.getDeviceInterface(), (*param_it).stringValue());
    new_op->addArg((*param_it).stringValue());
  }

  // Expression
  else if (Util::hasExpressionTag(*param_it)) 
  {
    if ((*param_it).getType() == Util::EXPR)
    {
      new_op = new ExpressionOp(var_type, (*param_it).getValue<Util::Expression>(), output_manager);
    }
    else if( ((*param_it).getType() == Util::DBLE) || ((*param_it).getType() == Util::INT) )
    {
      new_op = new ConstantOp(var_type, (*param_it).getImmutableValue<double>());
    }
    else
    {
      new_op = new ExpressionOp(var_type, (*param_it).tag(), output_manager);
    }
  }

  // Step Sweep
  if (!new_op)
  {
    for (size_t i = 0; i < output_manager.getStepParamVec().size(); ++i)
    {
      if (var_type == output_manager.getStepParamVec()[i].name)
      {
        new_op = new OutputMgrStepSweepOp(var_type, output_manager, i);
        break;
      }
    }
  }

  // DC Sweep
  if (!new_op)
  {
    for (size_t i = 0; i < output_manager.getDCParamVec().size(); ++i)
    {
      if (var_type == output_manager.getDCParamVec()[i].name)
      {
        new_op = new OutputMgrStepSweepOp(var_type, output_manager, i);
        break;
      }
    }
  }

  // Objective Function
  if (!new_op)
  {
    if (output_manager.getObjectiveMap().find(var_type) != output_manager.getObjectiveMap().end())
    {
      Objective &objective = output_manager.getObjective(var_type);
      new_op = new ObjectiveOp(var_type, objective);
    }
  }

  // Measurement
  if (!new_op)
  {
    const Measure::Base *measure = output_manager.getMeasureManager().find(var_type);
    if (measure) 
    {
      new_op = new MeasureOp(var_type, *measure);
    }
  }


  // Voltage references, current references and "node" references take up more 
  // than one element of the parameter list.  The first element of the list has
  // a tag indicating the type (V, VI, VR, etc.)  with a value
  // indicating the number of "args".  It is followed in the list by
  // that many names of nodes (in the qTAG of the param).
  // Here, we turn that construction back into a name as the user would
  // have specified it.

  // Get arguments and build pretty name
  std::string name = var_type;
  std::vector<std::string> args;
  if (var_type[0] == 'V' ||var_type[0] == 'I' ||var_type[0] == 'N') 
  {
    std::ostringstream oss;
    oss << var_type << "(";
    int arg_count = (*param_it).getImmutableValue<int>();
    for (int i = 0; i < arg_count; ++i) 
    {
      ++param_it;
      if (i != 0)
        oss << ",";
      oss << (*param_it).tag();
      args.push_back((*param_it).tag());
    }
    oss << ")";
    name = oss.str();
  }

  // Node (e.g. "N(A)" on the .print line, usually internal nodes of devices,
  // but could also be named elements of state or store vector)
  // If the name of a solution var, we support N, NI, NR, NM, etc., otherwise
  // just support "N"
  if (!new_op && var_type[0] == 'N' && args.size() == 1) 
  {
    NodeNamePairMap::const_iterator it = output_manager.getAllNodes().find(args[0]);
    if (it != output_manager.getAllNodes().end()) 
    {
      int index = (*it).second.first;
      if (var_type == "N" )
      {
        new_op = new SolutionOp(name, index);
      }
      else if (var_type == "NR" )
      {
        new_op = new SolutionOp(name, index);
      }
      else if (var_type == "NI" )
      {
        new_op = new SolutionImaginaryOp(name, index);
      }
      else if (var_type == "NM" )
      {
        new_op = new SolutionMagnitudeOp(name, index);
      }
      else if (var_type == "NP" )
      {
        new_op = new SolutionPhaseOp(name, index);
      }
      else if (var_type == "NDB" )
      {
        new_op = new SolutionDecibelsOp(name, index);
      }
    }
    else 
    {
      it = output_manager.getStateNodes().find(args[0]);
      if (it != output_manager.getStateNodes().end()) 
      {
        new_op = new StateOp(name, (*it).second.first);
      }
      else 
      {
        it = output_manager.getStoreNodes().find(args[0]);
        if (it != output_manager.getStoreNodes().end())
        {
          new_op = new StoreOp(name, (*it).second.first);
        }
        else 
        {
          new_op = new Util::UndefinedOp(var_type);
        }
      }
    }

    if (new_op)
      new_op->addArg(args[0]);
  }

  // Voltage
  if (!new_op && var_type[0] == 'V' && args.size() > 0) 
  {

    // Solution variable
    if (args.size() == 1) 
    {
      int index = -1;
      NodeNamePairMap::const_iterator it = output_manager.getAllNodes().find(args[0]);
      if (it != output_manager.getAllNodes().end())
        index = (*it).second.first;

      if (var_type == "V" )
      {
        new_op = new SolutionOp(name, index);
      }
      else if (var_type == "VR" )
      { 
       new_op = new SolutionRealOp(name, index);
      }
      else if (var_type == "VI" )
      {
        new_op = new SolutionImaginaryOp(name, index);
      }
      else if (var_type == "VM" )
      {
        new_op = new SolutionMagnitudeOp(name, index);
      }
      else if (var_type == "VP" )
      {
        new_op = new SolutionPhaseOp(name, index);
      }
      else if (var_type == "VDB" )
      {
        new_op = new SolutionDecibelsOp(name, index);
      }

      if (new_op)
        new_op->addArg(args[0]);
    }

    // Volatage Difference
    if (args.size() == 2) 
    {
      int index1 = -1;
      int index2 = -1;

      NodeNamePairMap::const_iterator it = output_manager.getAllNodes().find(args[0]);
      if (it != output_manager.getAllNodes().end())
        index1 = (*it).second.first;

      it = output_manager.getAllNodes().find(args[1]);
      if (it != output_manager.getAllNodes().end())
        index2 = (*it).second.first;

      if (var_type == "V" )
      {
        new_op = new VoltageDifferenceOp(name, index1, index2);
      }
      else if (var_type == "VR" )
      {
        new_op = new VoltageDifferenceRealOp(name, index1, index2);
      }
      else if (var_type == "VI" )
      {
        new_op = new VoltageDifferenceImaginaryOp(name, index1, index2);
      }
      else if (var_type == "VM" )
      {
        new_op = new VoltageDifferenceMagnitudeOp(name, index1, index2);
      }
      else if (var_type == "VP" )
      {
        new_op = new VoltageDifferencePhaseOp(name, index1, index2);
      }
      else if (var_type == "VDB" )
      {
        new_op = new VoltageDifferenceDecibelsOp(name, index1, index2);
      }

      if (new_op)
        new_op->addArgs(args.begin(), args.end());
    }
  }

  // branch currents or Lead currents 
  if (!new_op && var_type[0] == 'I') 
  {
    // Node name could be circuit_context:DeviceTypeDeviceName while internally it should be DeviceType:circuit_context:DeviceName.
    std::string modifiedName;
    std::string::size_type lastColonInName = args[0].find_last_of(":");
    if ((lastColonInName != std::string::npos) && (lastColonInName + 1 < args[0].length()))
    {
      std::string::iterator deviceName = args[0].begin() + lastColonInName+1;
      std::string::iterator namePrefixEnd = args[0].begin() + lastColonInName;
      modifiedName.append(deviceName, deviceName + 1);
      modifiedName.append(":");
      modifiedName.append(args[0].begin(), namePrefixEnd + 1);
      modifiedName.append(deviceName + 1, args[0].end());
    }
    else
    {
      modifiedName = args[0];
    }

    // could be a device lead current "DEV_I" or a branch current.
    // so we don't have to duplicate solution vars(branch currents) in the
    // store vector, look for each type.
    std::string store_name = modifiedName + ":DEV_" + var_type;  // if it is in the state/store vec.

    // this if block allows for spaces in YPDE names as in I1(YPDE NAME)
    // we try to find devices based on store_name in the following blocks of code,
    // so do this modification now.
    std::string::size_type space = store_name.find_first_of(" ");
    if (space != std::string::npos)
    {
      if (space == 4 && store_name.substr(0, 4) == "YPDE")
      {
        store_name.replace(4, 1, "%");
        store_name.insert(1, "%");
      }
    }

    // Search store
    NodeNamePairMap::const_iterator it = output_manager.getStoreNodes().find(store_name);
    if (it != output_manager.getStoreNodes().end())
    {
      new_op = new StoreOp(name, (*it).second.first);
    }

    // Search solution
    if (!new_op)
    {
      std::string solution_name = modifiedName + "_BRANCH";         // if it is in the solution vec.
      it = output_manager.getAllNodes().find(solution_name);
      if (it != output_manager.getAllNodes().end())
      {
        if (var_type == "I" )
        {
          int index = (*it).second.first;
          new_op = new SolutionOp(name, index);
        }
      }
    }

    // Special case for PDE and old device lead currents
    if (!new_op)
    {
      // this is confusing.  While the solution vector and state/store vector's
      // use maps with modified device names(specifically where the device
      // type is always first as in "D:subcircuitname:devciename" as apposed to
      // "subcircuitname:Ddevicename", the device manager does not use the modified
      // device name to find a device.  So, when we set up the device name below
      // use the "nodeName" rather than the "modifiedDeviceName".
      store_name = args[0] + ":DEV_" + var_type;
      // have to repeat this check for spaces as in I(YPDE NAME)
      std::string::size_type space = store_name.find_first_of(" ");
      if (space != std::string::npos)
      {
        if (space == 4 && store_name.substr(0, 4) == "YPDE")
        {
          store_name.replace(4, 1, "%");
          store_name.insert(1, "%");
        }
      }

      if (output_manager.getDeviceInterface().findParam(store_name)) 
      {
        new_op = new DeviceMgrParameterOp(var_type, output_manager.getDeviceInterface(), store_name);
      }

      std::string ppde("Y%PDE%" + store_name);
      if (!new_op && output_manager.getDeviceInterface().findParam(ppde)) 
      {
        {
          new_op = new DeviceMgrParameterOp(var_type, output_manager.getDeviceInterface(), ppde);
        }
      }
    }
    if (new_op)
      new_op->addArg(args[0]);
  }

  // Device Sensitivity
  if (!new_op)
  {
    const Device::DeviceEntity *device_entity = output_manager.getDeviceInterface().getDeviceEntity(var_type);
    if (device_entity) 
    {
      std::string param_name = Util::paramNameFromFullParamName(var_type);
      new_op = new DeviceEntityParameterOp(var_type, *device_entity, param_name);
    }
  }

  // Global Parameter
  if (!new_op)
  {
    const double *result = output_manager.getDeviceInterface().findGlobalPar(var_type);
    if (result) 
    {
      // Refactor: [DGB] Should use the address.
      new_op = new DeviceMgrGlobalParameterOp(var_type, output_manager.getDeviceInterface(), var_type);
    }
  }

  // Last chance, try as a device parameter
  if (!new_op)
  {
    if (output_manager.getDeviceInterface().findParam(var_type))
    {
      new_op = new DeviceMgrParameterOp(var_type, output_manager.getDeviceInterface(), var_type);
    }
  }

  if (!new_op) 
  {
    new_op = new Util::UndefinedOp(var_type);
    new_op->addArgs(args.begin(), args.end());
  }
  
  return new_op;
}

//-----------------------------------------------------------------------------
// Function      : makeOps
// Purpose       : given an output manager, a parameter list 
//                (defined by begin and end iterators), and a back_insert
//                iterator for an OpList, create all the the Ops with makeOp,
//                synchronize across processors, validate them,  and if
//                all is well, put them onto the end of the OpList using the
//                back_inserter.
//                 
// Special Notes : 
// Scope         : global, Xyce::IO namespace
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void makeOps(const OutputMgr &output_manager, const NetlistLocation &netlist_location, Util::ParameterList::iterator begin, Util::ParameterList::iterator end, std::back_insert_iterator<Util::OpList> inserter)
{
  Util::OpList ops;

  for (Util::ParameterList::iterator it = begin; it != end; ++it)
    ops.push_back(makeOp(output_manager, it));

  // Resolve parallel here
  std::vector<int> op_type;
  for (Util::OpList::const_iterator it = ops.begin(); it != ops.end(); ++it)
    op_type.push_back((*it)->opType());

  Parallel::Machine comm = output_manager.getCommPtr()->comm();
  Parallel::AllReduce(comm, MPI_MAX, op_type);

  // Validate ops and report errors
  std::vector<int>::const_iterator op_type_it = op_type.begin();
  for (std::vector<Util::Operator *>::iterator it = ops.begin(); it != ops.end(); ++it, ++op_type_it) 
  {
    if ((*it)->isType<Util::UndefinedOp>()) 
    {
      std::string name = (*it)->getName();
      const std::vector<std::string> &arg_list = (*it)->getArgs();
      if (!arg_list.empty()) 
      {
        name += "(";
        for (std::vector<std::string>::const_iterator it = arg_list.begin(); it != arg_list.end(); ++it) 
        {
          if (it != arg_list.begin())
            name += ",";
          name += (*it);
        }
        name += ")";
      }
      if ((*op_type_it) == Util::UNDEFINED)
        Report::UserError().at(netlist_location) << "Function or variable " << name << " is not defined";
      else 
      {
        CreateFunction f = getOpMap(*op_type_it);
        delete *it;
        (*it) = f(name);
      }
    }
    else if ((*it)->opType() != (*op_type_it))
      Report::UserError().at(netlist_location) << "Differing types for " << (*it)->getName() << " discovered across processors";
  }

  std::copy(ops.begin(), ops.end(), inserter);
}

//-----------------------------------------------------------------------------
// Function      : getValue
// Purpose       : Given a communicator, an op, and solution/state/store
//                 vectors, evaluate the op and return its value
// Special Notes : This effectively replaces the old "getPrintValue" method,
//                 and is used throughout the Outputter class.
// Scope         : global, Xyce::IO namespace
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex getValue(
  Parallel::Machine             comm,
  const Util::Operator &        op,
  const N_LAS_Vector *          real_solution_vector,
  const N_LAS_Vector *          imaginary_solution_vector,
  const N_LAS_Vector *          state_vector,
  const N_LAS_Vector *          store_vector)
{

  op.setSolutionVector(real_solution_vector);
  op.setSolutionImagVector(imaginary_solution_vector);
  op.setStateVector(state_vector);
  op.setStoreVector(store_vector);

  return op(comm);
}

} // namespace IO
} // namespace Xyce
