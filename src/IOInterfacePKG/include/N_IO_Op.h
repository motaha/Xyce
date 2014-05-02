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
// Filename       : $RCSfile: N_IO_Op.h,v $
//
// Purpose        : Generate global id structures and proc maps
//                  and distribute nodes to processors
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.9.2.2 $
//
// Revision Date  : $Date: 2014/02/26 20:42:38 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_Op_h
#define Xyce_N_IO_Op_h

#include <iterator>

#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_IO_Measure_fwd.h>
#include <N_UTL_fwd.h>

#include <N_IO_OutputMgr.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>
#include <N_UTL_Op.h>
#include <N_UTL_ExpressionData.h>

namespace Xyce {
namespace IO {

struct ReduceSum
{
  static complex reduce(Parallel::Machine comm, complex result);
};

class ConstantOp : public Util::Op<ConstantOp, Util::ReduceNone, Util::EvalNoop>
{
public:
  ConstantOp(const std::string &name, complex value)
    : Base(name),
      value_(value)
  {}

  virtual ~ConstantOp()
  {}

  static complex get(const ConstantOp &op)
  {
    return op.value_;
  }

  const complex       value_;
};


class OutputMgrCurrentOutputterIndexOp : public Util::Op<OutputMgrCurrentOutputterIndexOp, Util::ReduceNone, Util::EvalNoop>
{
public:
  OutputMgrCurrentOutputterIndexOp(const std::string &name, const OutputMgr &output_manager)
    : Base(name),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrCurrentOutputterIndexOp()
  {}

  static complex get(const OutputMgrCurrentOutputterIndexOp &op);

  const OutputMgr &   outputMgr_;
};

class OutputMgrTimeOp : public Util::Op<OutputMgrTimeOp, Util::ReduceNone, Util::EvalNoop>
{
public:
  OutputMgrTimeOp(const std::string &name, const OutputMgr &output_manager)
    : Base(name),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrTimeOp()
  {}

  static complex get(const OutputMgrTimeOp &op);
    
  const OutputMgr &   outputMgr_;
};

class OutputMgrFrequencyOp : public Util::Op<OutputMgrFrequencyOp, Util::ReduceNone, Util::EvalNoop>
{
public:
  OutputMgrFrequencyOp(const std::string &name, const OutputMgr &output_manager)
    : Base(name),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrFrequencyOp()
  {}

  static complex get(const OutputMgrFrequencyOp &op);

  const OutputMgr &   outputMgr_;
};


class OutputMgrTemperatureOp : public Util::Op<OutputMgrTemperatureOp, Util::ReduceNone, Util::EvalNoop>
{
public:
  OutputMgrTemperatureOp(const std::string &name, const OutputMgr &output_manager)
    : Base(name),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrTemperatureOp()
  {}

  static complex get(const OutputMgrTemperatureOp &op);
    
  const OutputMgr &   outputMgr_;
};

class OutputMgrStepSweepOp : public Util::Op<OutputMgrStepSweepOp, Util::ReduceNone, Util::EvalNoop>
{
public:
  OutputMgrStepSweepOp(const std::string &name, const OutputMgr &output_manager, int index)
    : Base(name),
      outputMgr_(output_manager),
      index_(index)
  {}
    
  virtual ~OutputMgrStepSweepOp()
  {}

  static complex get(const OutputMgrStepSweepOp &op);
    
  const int           index_;
  const OutputMgr &   outputMgr_;
};

class OutputMgrDCSweepOp : public Util::Op<OutputMgrDCSweepOp, Util::ReduceNone, Util::EvalNoop>
{
public:
  OutputMgrDCSweepOp(const std::string &name, const OutputMgr &output_manager, int index)
    : Base(name),
      outputMgr_(output_manager),
      index_(index)
  {}
    
  virtual ~OutputMgrDCSweepOp()
  {}
    
  static complex get(const OutputMgrDCSweepOp &op);
    
  const int           index_;
  const OutputMgr &   outputMgr_;
};

class ObjectiveOp : public Util::Op<ObjectiveOp, ReduceSum, Util::EvalNoop>
{
public:
  ObjectiveOp(const std::string &name, Objective &objective)
    : Base(name),
      objective_(objective)
  {}
    
  virtual ~ObjectiveOp()
  {}

  static complex get(const ObjectiveOp &op);
    
  Objective &         objective_;
};

class DeviceMgrGlobalParameterOp : public Util::Op<DeviceMgrGlobalParameterOp, Util::ReduceNone, Util::EvalNoop>
{
public:
  DeviceMgrGlobalParameterOp(const std::string &name, const Device::DeviceInterface &device_interface, const std::string &device_parameter_name)
    : Base(name),
      deviceInterface_(device_interface),
      deviceParameterName_(device_parameter_name)
  {}

  virtual ~DeviceMgrGlobalParameterOp()
  {}

  static complex get(const DeviceMgrGlobalParameterOp &op);

  const std::string                   deviceParameterName_;
  const Device::DeviceInterface &     deviceInterface_;
};

class DeviceMgrParameterOp : public Util::Op<DeviceMgrParameterOp, Util::ReduceNone, Util::EvalNoop>
{
public:
  DeviceMgrParameterOp(const std::string &name, const Device::DeviceInterface &device_interface, const std::string &device_parameter_name)
    : Base(name),
      deviceInterface_(device_interface),
      deviceParameterName_(device_parameter_name)
  {}


  virtual ~DeviceMgrParameterOp()
  {}

  static complex get(const DeviceMgrParameterOp &op);

  const std::string                   deviceParameterName_;
  const Device::DeviceInterface &     deviceInterface_;
};

class DeviceEntityParameterOp : public Util::Op<DeviceEntityParameterOp, ReduceSum, Util::EvalNoop>
{
public:
  DeviceEntityParameterOp(const std::string &name, const Device::DeviceEntity &device_entity, const std::string &device_parameter_name)
    : Base(name),
      deviceEntity_(device_entity),
      deviceParameterName_(device_parameter_name)
  {}
    

  virtual ~DeviceEntityParameterOp()
  {}

  static complex get(const DeviceEntityParameterOp &op);

  const std::string                   deviceParameterName_;
  const Device::DeviceEntity &        deviceEntity_;
};

class SolutionOp : public Util::Op<SolutionOp, ReduceSum, Util::EvalNoop>
{
public:
  SolutionOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}
    

  virtual ~SolutionOp()
  {}

  static complex get(const SolutionOp &op);

  const int           index_;
};

class SolutionRealOp : public Util::Op<SolutionRealOp, ReduceSum, Util::EvalNoop>
{
public:
  SolutionRealOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}
    

  virtual ~SolutionRealOp()
  {}

  static complex get(const SolutionRealOp &op);
  static complex eval(complex result);

  const int           index_;
};

class SolutionImaginaryOp : public Util::Op<SolutionImaginaryOp, ReduceSum>
{
public:
  SolutionImaginaryOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SolutionImaginaryOp()
  {}

  static complex get(const SolutionImaginaryOp &op);
  static complex eval(complex result);

  const int           index_;
};

class SolutionMagnitudeOp : public Util::Op<SolutionMagnitudeOp, ReduceSum>
{
public:
  SolutionMagnitudeOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SolutionMagnitudeOp()
  {}

  static complex get(const SolutionMagnitudeOp &op);
  static complex eval(complex result);

  const int           index_;
};

class SolutionPhaseOp : public Util::Op<SolutionPhaseOp, ReduceSum>
{
public:
  SolutionPhaseOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SolutionPhaseOp()
  {}

  static complex get(const SolutionPhaseOp &op);
  static complex eval(complex result);

  const int           index_;
};

class SolutionDecibelsOp : public Util::Op<SolutionDecibelsOp, ReduceSum>
{
public:
  SolutionDecibelsOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SolutionDecibelsOp()
  {}

  static complex get(const SolutionDecibelsOp &op);
  static complex eval(complex result);

  const int           index_;
};

class VoltageDifferenceOp : public Util::Op<VoltageDifferenceOp, ReduceSum, Util::EvalNoop>
{
public:
  VoltageDifferenceOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferenceOp()
  {}

  static complex get(const VoltageDifferenceOp &op);

  const int           index1_;
  const int           index2_;
};

class VoltageDifferenceRealOp : public Util::Op<VoltageDifferenceRealOp, ReduceSum, Util::EvalNoop>
{
public:
  VoltageDifferenceRealOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferenceRealOp()
  {}

  static complex get(const VoltageDifferenceRealOp &op);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

class VoltageDifferenceImaginaryOp : public Util::Op<VoltageDifferenceImaginaryOp, ReduceSum>
{
public:
  VoltageDifferenceImaginaryOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferenceImaginaryOp()
  {}

  static complex get(const VoltageDifferenceImaginaryOp &op);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

class VoltageDifferenceMagnitudeOp : public Util::Op<VoltageDifferenceMagnitudeOp, ReduceSum>
{
public:
  VoltageDifferenceMagnitudeOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferenceMagnitudeOp()
  {}

  static complex get(const VoltageDifferenceMagnitudeOp &op);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

class VoltageDifferencePhaseOp : public Util::Op<VoltageDifferencePhaseOp, ReduceSum>
{
public:
  VoltageDifferencePhaseOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferencePhaseOp()
  {}

  static complex get(const VoltageDifferencePhaseOp &op);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

class VoltageDifferenceDecibelsOp : public Util::Op<VoltageDifferenceDecibelsOp, ReduceSum>
{
public:
  VoltageDifferenceDecibelsOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferenceDecibelsOp()
  {}

  static complex get(const VoltageDifferenceDecibelsOp &op);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

class StateOp : public Util::Op<StateOp, ReduceSum, Util::EvalNoop>
{
public:
  StateOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~StateOp()
  {}

  static complex get(const StateOp &op);

  const int           index_;
};

class StoreOp : public Util::Op<StoreOp, ReduceSum, Util::EvalNoop>
{
public:
  StoreOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~StoreOp()
  {}

  static complex get(const StoreOp &op);

  const int           index_;
};

class MeasureOp : public Util::Op<MeasureOp, Util::ReduceNone, Util::EvalNoop>
{
public:
  MeasureOp(const std::string &name, const Measure::Base &measure)
    : Base(name),
      measure_(measure)
  {}

  virtual ~MeasureOp()
  {}

  static complex get(const MeasureOp &op);

  const Measure::Base &       measure_;
};

class ExpressionOp : public Util::Op<ExpressionOp, Util::ReduceNone, Util::EvalNoop>
{
public:
  ExpressionOp(const std::string &name, Util::Expression &expression, const OutputMgr &output_manager);
  ExpressionOp(const std::string &name, const std::string &expression, const OutputMgr &output_manager);

  virtual ~ExpressionOp()
  {}

  void init();

  static complex get(const ExpressionOp &op);

  mutable Util::ExpressionData        expressionData_;
  const OutputMgr &                   outputMgr_;
};

} // namespace IO

namespace Util {

template<>
struct OpType<IO::ConstantOp> 
{
  enum {type = CONSTANT};
};

template<>
struct OpType<IO::OutputMgrCurrentOutputterIndexOp> 
{
  enum {type = INDEX};
};

template<>
struct OpType<IO::OutputMgrTimeOp> 
{
  enum {type = TIME_VAR};
};

template<>
struct OpType<IO::OutputMgrFrequencyOp> 
{
  enum {type = FREQUENCY};
};

template<>
struct OpType<IO::OutputMgrTemperatureOp> 
{
  enum {type = TEMPERATURE};
};

template<>
struct OpType<IO::OutputMgrStepSweepOp> 
{
  enum {type = STEP_SWEEP_VAR};
};

template<>
struct OpType<IO::OutputMgrDCSweepOp> 
{
  enum {type = DC_SWEEP_VAR};
};

template<>
struct OpType<IO::ObjectiveOp> 
{
  enum {type = OBJECTIVE_FUNCTION};
};

template<>
struct OpType<IO::DeviceMgrGlobalParameterOp> 
{
  enum {type = GLOBAL_PARAMETER};
};

template<>
struct OpType<IO::DeviceMgrParameterOp> 
{
  enum {type = DEVICE_PARAMETER};
};

template<>
struct OpType<IO::DeviceEntityParameterOp> 
{
  enum {type = DEVICE_ENTITY_PARAMETER};
};

template<>
struct OpType<IO::SolutionOp> 
{
  enum {type = SOLUTION_VAR};
};

template<>
struct OpType<IO::SolutionRealOp> 
{
  enum {type = SOLUTION_VAR_REAL};
};

template<>
struct OpType<IO::SolutionImaginaryOp> 
{
  enum {type = SOLUTION_VAR_IMAG};
};

template<>
struct OpType<IO::SolutionMagnitudeOp> 
{
  enum {type = SOLUTION_VAR_MAG};
};

template<>
struct OpType<IO::SolutionPhaseOp> 
{
  enum {type = SOLUTION_VAR_PHASE};
};

template<>
struct OpType<IO::SolutionDecibelsOp> 
{
  enum {type = SOLUTION_VAR_DB};
};

template<>
struct OpType<IO::VoltageDifferenceOp> 
{
  enum {type = VOLTAGE_DIFFERENCE};
};

template<>
struct OpType<IO::VoltageDifferenceRealOp> 
{
  enum {type = VOLTAGE_DIFFERENCE_REAL};
};

template<>
struct OpType<IO::VoltageDifferenceImaginaryOp> 
{
  enum {type = VOLTAGE_DIFFERENCE_IMAG};
};

template<>
struct OpType<IO::VoltageDifferenceMagnitudeOp> 
{
  enum {type = VOLTAGE_DIFFERENCE_MAG};
};

template<>
struct OpType<IO::VoltageDifferencePhaseOp> 
{
  enum {type = VOLTAGE_DIFFERENCE_PHASE};
};

template<>
struct OpType<IO::VoltageDifferenceDecibelsOp> 
{
  enum {type = VOLTAGE_DIFFERENCE_DB};
};

template<>
struct OpType<IO::StateOp> 
{
  enum {type = STATE_VAR};
};

template<>
struct OpType<IO::StoreOp> 
{
  enum {type = STORE_VAR};
};

template<>
struct OpType<IO::MeasureOp> 
{
  enum {type = MEASURE_FUNCTION};
};

template<>
struct OpType<IO::ExpressionOp> 
{
  enum {type = EXPRESSION};
};

} // namespace Util

namespace IO {

typedef Util::Operator *(*CreateFunction)(const std::string &name);
typedef std::map<int, CreateFunction> OpMap;

complex getValue(
   Parallel::Machine             comm,
   const Util::Operator &        op,
   const N_LAS_Vector *          real_solution_vector,
   const N_LAS_Vector *          imaginary_solution_vector,
   const N_LAS_Vector *          state_vector,
   const N_LAS_Vector *          store_vector);

void makeOps(const OutputMgr &output_manager, const NetlistLocation &netlist_location, ParameterList::iterator begin, ParameterList::iterator end, std::back_insert_iterator<Util::OpList> inserter);

inline void makeOps(const OutputMgr &output_manager, ParameterList::iterator begin, ParameterList::iterator end, std::back_insert_iterator<Util::OpList> inserter)
{
  makeOps(output_manager, NetlistLocation(), begin, end, inserter);
}


} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_Op_h
