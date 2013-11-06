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
// Filename       : $RCSfile: N_UTL_Param.h,v $
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
// Revision Number: $Revision: 1.55.2.4 $
//
// Revision Date  : $Date: 2013/10/03 17:23:52 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef  N_UTL_PARAM_H
#define  N_UTL_PARAM_H

// ---------- Forward Declarations ----------

class N_UTL_ParamData;
class N_UTL_Expression;

// ---------- Standard Includes ----------
#include <string>
#include <vector>
#include <iosfwd>
#ifdef HAVE_STRCASECMP
#include <string.h>
#endif

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>
#include <N_UTL_Packable.h>
#include <N_UTL_ExpressionData.h>

enum { STR, DBLE, INT, LNG, EXPR, BOOL, STR_VEC, INT_VEC, DBLE_VEC, DBLE_VEC_IND, COMPOSITE };


// this enum is used to define what the context of the N_UTL_Param object
// is within the simulation for evaluation of the parameter.   For example
// "TEMPERATURE" indicates that this N_UTL_Param object referes to the
// current temperature defined by the simulator, and no additional
// information is needed.  If the context is "VOLTAGE" then to evaluate the
// N_UTL_Param object one will also need and index into the solution
// vector.  Thus, this context tells one how to use the extIndex1 and
// extIndex2 vars to assign a value to the N_UTL_Param object.  (two
// indicies are needed for voltage differences.  otherwise just about
// everything else needs only one extra index.
//
// One wrinkle in this is that the list of N_UTL_Param objects stores
// I( name ), V(name) and N(name) as two sequential N_UTL_Param ojects one for the
// I, V or N and the following one for the "name".
// Likewise V( name1, name2 ) uses three N_UTL_Param objects, one for V
// and two for name1 and name2.
// These following N_UTL_Param objects which hold names are needed once to
// find indicies but not after that.  Thus there is an enum NODE_OR_DEVICE_NAME
// to categorize those N_UTL_Param objects that can be dropped after context
// resolution.

typedef enum {
  UNDEFINED,
  INDEX,
  CONSTANT,
  TEMPERATURE,
  FREQUENCY,
  TIME_VAR,
  STEP_SWEEP_VAR,
  DC_SWEEP_VAR,
  EXPRESSION,
  VOLTAGE_DIFFERENCE,
  NODE_OR_DEVICE_NAME,
  SOLUTION_VAR,
  STATE_VAR,
  STORE_VAR,
  DEVICE_PARAMETER,
  GLOBAL_PARAMETER,
  OBJECTIVE_FUNCTION,
  MEASURE_FUNCTION,
  SOLUTION_VAR_REAL,               // these are for solution variables under AC analysis where
  SOLUTION_VAR_IMAG,               // solution is complex and could be printed as real, imaginary,
  SOLUTION_VAR_MAG,                // magnitude, phase and Db = (20*log10(magnitude))
  SOLUTION_VAR_PHASE,
  SOLUTION_VAR_DB,
  VOLTAGE_DIFFERENCE_REAL,
  VOLTAGE_DIFFERENCE_IMAG
} SimulatorVariableContext;


// ----------   Other Includes   ----------

class N_UTL_ParamData
{
  public:

    N_UTL_ParamData()
      : tag_(string()), type_(-1), simVarContext(UNDEFINED), extIndex1(-1), extIndex2(-1) {}

    N_UTL_ParamData(const string & t, const string & v)
      : tag_(t), type_(STR), simVarContext(UNDEFINED), extIndex1(-1), extIndex2(-1)
    {
      vals.stringVal_ = new string(v);
    }
    N_UTL_ParamData(const string & t, const double & v)
      : tag_(t), type_(DBLE), simVarContext(UNDEFINED), extIndex1(-1), extIndex2(-1)
    {
      vals.dbleVal_ = v;
    }
    N_UTL_ParamData(const string & t, const int & v)
      : tag_(t), type_(INT), simVarContext(UNDEFINED), extIndex1(-1), extIndex2(-1)
    {
      vals.intVal_ = v;
    }
    N_UTL_ParamData(const string & t, const long & v)
      : tag_(t), type_(LNG), simVarContext(UNDEFINED), extIndex1(-1), extIndex2(-1)
    {
      vals.lngVal_ = v;
    }
    N_UTL_ParamData(const string & t, const bool & v)
      : tag_(t), type_(BOOL), simVarContext(UNDEFINED), extIndex1(-1), extIndex2(-1)
    {
      vals.boolVal_ = v;
    }

    N_UTL_ParamData(const string & t, const vector<string> & v);
    N_UTL_ParamData(const string & t, const vector<double> & v);
    N_UTL_ParamData(const string & t, const vector<int> & v);
    N_UTL_ParamData(const string & t, const N_UTL_Expression & v);

    ~N_UTL_ParamData();

    void freeAllocatedStorage();


    string sReturnVal;
    string tag_;
    int type_;
    union
    {
      N_UTL_Expression *expression;
      string * stringVal_;
      double dbleVal_;
      int    intVal_;
      long   lngVal_;
      bool   boolVal_;
      vector <string> * stringVec_;
      vector <double> * dbleVec_;
      vector <int> * intVec_;
      N_UTL_Expression * exprVal_;
    } vals;

    // The simVarContext is used to determine how to evaluate
    // the value of this parameter during the simulation.  If the
    // context requires additional data for the evaluation (like
    // an offset or index into a vector), that can be stored in
    // extIndex1 or extIndex2.  For the more complex case of an expression,
    // the additional data can be stored in the expressionDataRCP object.
    // This keeps the categorization and evaluation of the N_UTL_Param object
    // compact and easy to re-evaluate every time step.
    SimulatorVariableContext simVarContext;
    int extIndex1;
    int extIndex2;
    RCP<N_UTL_ExpressionData> expressionDataRCP;
    string qualifiedParameterOrFunctionName;
    double fixedValue;
};

//-----------------------------------------------------------------------------
// Class         : N_UTL_Param
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
class N_UTL_Param : public Packable
{
public:

  // Default constructor.
  N_UTL_Param();

  // Constructors
  N_UTL_Param(const string & t, const string & v);
  N_UTL_Param(const string & t, const double & v);
  N_UTL_Param(const string & t, const int & v);
  N_UTL_Param(const string & t, const long & v);
  N_UTL_Param(const string & t, const bool & v);
  N_UTL_Param(const string & t, const char * v);
  N_UTL_Param(const string & t, const vector<string> & v);
  N_UTL_Param(const string & t, const vector<double> & v);
  N_UTL_Param(const string & t, const vector<int> & v);
  N_UTL_Param(const string & t, const N_UTL_Expression & v);

  // Copy constructor.
  N_UTL_Param(N_UTL_Param const & rhsParam);

  // Assignment operator.
  N_UTL_Param & operator = (N_UTL_Param const & rhsParam);

  // deepCompare -- compare TAG and Value if TAGS are the same
  bool deepCompare(N_UTL_Param const & rhsParam) const;

  // Equality operator.
  bool operator == (N_UTL_Param const & rhsParam) const;

  // None-equality operator.
  bool operator != (N_UTL_Param const & rhsParam) const;

  // Less-than operator.
  bool operator < (N_UTL_Param const & rhsParam) const;

  // Greater-than operator.
  bool operator > (N_UTL_Param const & rhsParam) const;

  // Destructor
  virtual ~N_UTL_Param();

  // Methods to set reset value.
  N_UTL_Param & set(const string & tag, const N_UTL_Param &);
  N_UTL_Param & set(const string & tag, const string & val);
  N_UTL_Param & set(const string & tag, const double & val);
  N_UTL_Param & set(const string & tag, const int & val);
  N_UTL_Param & set(const string & tag, const long & val);
  N_UTL_Param & set(const string & tag, const bool & val);
  N_UTL_Param & set(const string & tag, const char * val);
  N_UTL_Param & set(const string & tag, const vector<string> & val);
  N_UTL_Param & set(const string & tag, const vector<double> & val);
  N_UTL_Param & set(const string & tag, const vector<int> & val);
  N_UTL_Param & set(const string & tag, const N_UTL_Expression & val);

  void setTag(const string & tag);

  void setVal(const N_UTL_Param &);
  void setVal(const string & val);
  void setVal(const double & val);
  void setVal(const int & val);
  void setVal(const long & val);
  void setVal(const bool & val);
  void setVal(const char * val);
  void setVal(const vector<string> & val);
  void setVal(const vector<double> & val);
  void setVal(const vector<int> & val);
  void setVal(const N_UTL_Expression & val);

  // Methods to get tag.
  const string & tag() const;
  const string uTag() const;
  const string lTag() const;

  // Methods to get the "val" in the desired format.
  const string & sVal() const;
  const double & dVal() const;
  const int & iVal() const;
  const long & lVal() const;
  bool bVal() const;
  vector<string> * sVecPtr();
  const vector<string> & sVecVal() const;
  vector<double> * dVecPtr();
  const vector<double> & dVecVal() const;
  vector<int> * iVecPtr();
  const vector<int> & iVecVal() const;
  N_UTL_Expression * ePtr();
  N_UTL_Expression * ePtr() const;
  const N_UTL_Expression & eVal() const;

  // Special accessor methods.
  const string usVal() const;
  const string lsVal() const;

  const int & getType() const;

  // Method for checking if the parameter is expression valued.
  bool hasExpressionValue() const;

  // Method for checking if the parameter is expression valued.
  bool hasExpressionTag () const;

  // Method to check if the parameter value is enclosed in double quotes.
  bool isQuoted();

  // Method to check whether a parameter has a legal real or integer numeric
  // value.
  bool isNumeric() const;
  bool isInteger() const;
  bool isBool() const;

  // Methods for working with time dependency of parameters.
  void setTimeDependent( bool const& timeDependent );
  bool isTimeDependent() const;

  // Method for accessing context data
  SimulatorVariableContext getSimContext() const {return data_->simVarContext;};
  int getExtraIndex1() const {return data_->extIndex1;};
  int getExtraIndex2() const {return data_->extIndex2;};
  RCP<N_UTL_ExpressionData> getExpressionDataPointer() const {return data_->expressionDataRCP;};
  string getQualifiedParameterOrFunctionName() const {return data_->qualifiedParameterOrFunctionName; };
  double getFixedValue() const {return data_->fixedValue; };

  // method for setting context and extra data.  Only need three combinations
  // as we figure out the context and extra data at the same time and then set it.
  void setSimContextAndData( SimulatorVariableContext aSimContext);
  void setSimContextAndData( SimulatorVariableContext aSimContext, int extraIndex1);
  void setSimContextAndData( SimulatorVariableContext aSimContext, int extraIndex1, int extraIndex2 );
  void setSimContextAndData( SimulatorVariableContext aSimContext, RCP<N_UTL_ExpressionData> expDataRcp);
  void setSimContextAndData( SimulatorVariableContext aSimContext, string parameterName );
  void setSimContextAndData( SimulatorVariableContext aSimContext, double aValue );

  // Method for outputting for debugging
  virtual void print();

  // Packing Functionality
  virtual Packable * instance() const;

  // Counts bytes needed to pack block.
  virtual int packedByteCount() const;

  // Packs OptionBlock into char buffer using MPI_PACK.
  virtual void pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const;

  // Unpacks OptionBlock from char buffer using MPI_UNPACK.
  virtual void unpack(char * pB, int bsize, int & pos, N_PDS_Comm * comm);

  friend ostream & operator << (ostream & os, const N_UTL_Param & p);

private:

  // Pointer to parameter data.
  N_UTL_ParamData * data_;
};


//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::operator==
// Purpose       : "==" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
inline bool N_UTL_Param::operator==( N_UTL_Param const& rhsParam ) const
{
#ifndef HAVE_STRCASECMP
  return ( uTag() == rhsParam.uTag() );
#else
  return (
  strcasecmp(
  ( data_->tag_.c_str() ),
  ( rhsParam.data_->tag_.c_str() ) )==0
      );
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::operator<
// Purpose       : "<" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
inline bool N_UTL_Param::operator<( N_UTL_Param const& rhsParam ) const
{
#ifndef HAVE_STRCASECMP
  return ( uTag() < rhsParam.uTag() );
#else
  return (
    strcasecmp(
    ( data_->tag_.c_str() ),
    ( rhsParam.data_->tag_.c_str() ) )<0
        );
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::operator>
// Purpose       : ">" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
inline bool N_UTL_Param::operator>( N_UTL_Param const& rhsParam ) const
{
#ifndef HAVE_STRCASECMP
  return ( uTag() > rhsParam.uTag() );
#else
  return (
    strcasecmp(
    ( data_->tag_.c_str() ),
    ( rhsParam.data_->tag_.c_str() ) )>0
        );
#endif
}

//-----------------------------------------------------------------------------
// Function      : upperTag
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
inline const string N_UTL_Param::uTag() const
{
  return ExtendedString( data_->tag_ ).toUpper();
}

//-----------------------------------------------------------------------------
// Function      : lowerTag
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
inline const string N_UTL_Param::lTag() const
{
  return ExtendedString( data_->tag_ ).toLower();
}

#endif
