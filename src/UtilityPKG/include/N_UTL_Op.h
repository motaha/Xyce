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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_Op.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 5/22/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_Op_h
#define Xyce_N_UTL_Op_h

#include <typeinfo>
#include <string>

#include <N_UTL_fwd.h>
#include <N_PDS_fwd.h>
#include <N_LAS_Vector.h>

namespace Xyce {
namespace Util {

enum {
  UNDEFINED,
  INDEX,
  CONSTANT,
  TEMPERATURE,
  FREQUENCY,
  TIME_VAR,
  STEP_SWEEP_VAR,
  DC_SWEEP_VAR,
  EXPRESSION,
  DEVICE_PARAMETER,
  DEVICE_ENTITY_PARAMETER,
  GLOBAL_PARAMETER,
  OBJECTIVE_FUNCTION,
  MEASURE_FUNCTION,
  SOLUTION_VAR,
  SOLUTION_VAR_REAL,
  SOLUTION_VAR_IMAG,
  SOLUTION_VAR_MAG,
  SOLUTION_VAR_PHASE,
  SOLUTION_VAR_DB,
  VOLTAGE_DIFFERENCE,
  VOLTAGE_DIFFERENCE_REAL,
  VOLTAGE_DIFFERENCE_IMAG,
  VOLTAGE_DIFFERENCE_MAG,
  VOLTAGE_DIFFERENCE_PHASE,
  VOLTAGE_DIFFERENCE_DB,
  STATE_VAR,
  STORE_VAR
};

template <class T>
struct OpType;

template<class T, class R, class E>
class Op;

class Operator
{
  public:
    Operator(const std::string &name)
      : name_(name)
    {}

    virtual ~Operator()
    {}

    template <class U>
    bool isType() const 
    {
      return type() == typeid(U);
    }

    const std::string &getName() const 
    {
      return name_;
    }

    complex operator()(Parallel::Machine comm) const 
    {
      return evaluate(comm);
    }

    virtual const std::type_info &type() const = 0;

    virtual complex evaluate(Parallel::Machine comm) const = 0;

    virtual int opType() const = 0;

    void addArg(const std::string &arg) 
    {
      argList_.push_back(arg);
    }

    template <class II>
    void addArgs(II begin, II end) 
    {
      argList_.assign(begin, end);
    }

    const std::vector<std::string> &getArgs() 
    {
      return argList_;
    }
    
    void setSolutionVector(const N_LAS_Vector *solution_vector) const 
    {
      realSolutionVector_ = solution_vector;
    }

    void setSolutionImagVector(const N_LAS_Vector *solution_vector) const 
    {
      imaginarySolutionVector_ = solution_vector;
    }

    void setStateVector(const N_LAS_Vector *solution_vector) const 
    {
      stateVector_ = solution_vector;
    }

    void setStoreVector(const N_LAS_Vector *solution_vector) const 
    {
      storeVector_ = solution_vector;
    }

  protected:
    const std::string                   name_;
    std::vector<std::string>            argList_;
    mutable const N_LAS_Vector *        realSolutionVector_;
    mutable const N_LAS_Vector *        imaginarySolutionVector_;
    mutable const N_LAS_Vector *        stateVector_;
    mutable const N_LAS_Vector *        storeVector_;
};

template <class T, class R, class E = T>
class ReduceOp_ : public Operator
{
  public:
    static Util::Operator *create(const std::string &name) 
    {
      return new ReduceOp_(name);
    }

    ReduceOp_(const std::string &name)
      : Operator(name)
    {}

    virtual ~ReduceOp_()
    {}

    virtual const std::type_info &type() const 
    {
      return typeid(T);
    }

    virtual int opType() const 
    {
      return OpType<T>::type;
    }

  protected:
    complex get_() const 
    {
      return complex(0.0, 0.0);
    }

    complex reduce_(Parallel::Machine comm, complex x) const 
    {
      return R::reduce(comm, x);
    }

    complex eval_(complex x) const 
    {
      return E::eval(x);
    }

    virtual complex evaluate(Parallel::Machine comm) const 
    {
      return eval_(reduce_(comm, get_()));
    }
};

template<class T, class R, class E>
class Op : public Operator
{
  public:
    typedef Op<T, R, E> Base;
    typedef ReduceOp_<T, R, E> ReduceOp;

    Op(const std::string &name)
      : Operator(name)
    {}

    virtual ~Op()
    {}

    virtual const std::type_info &type() const 
    {
      return typeid(T);
    }

    virtual int opType() const 
    {
      return OpType<T>::type;
    }

  protected:
    complex get_() const 
    {
      return T::get(static_cast<const T &>(*this));
    }

    complex reduce_(Parallel::Machine comm, complex x) const 
    {
      return R::reduce(comm, x);
    }

    complex eval_(complex x) const 
    {
      return E::eval(x);
    }

    virtual complex evaluate(Parallel::Machine comm) const 
    {
      return eval_(reduce_(comm, get_()));
    }
};

struct ReduceNone
{
    static complex reduce(Parallel::Machine comm, complex result) 
    {
      return result;
    }
};

struct EvalNoop 
{
    static complex eval(complex result) 
    {
      return result;
    }
};

class UndefinedOp : public Util::Op<UndefinedOp, ReduceNone, EvalNoop>
{
  public:
    UndefinedOp(const std::string &name = "Undefined")
      : Base(name)
    {}

    virtual ~UndefinedOp()
    {}

    static complex get(const UndefinedOp &op) 
    {
      return complex(0.0, 0.0);
    }
};
template<>
struct OpType<UndefinedOp> 
{
    enum {type = UNDEFINED};
};

typedef std::vector<Operator *> OpList;

} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_Op_h
