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

//-------------------------------------------------------------------------
// Filename      : $RCSfile: N_UTL_Param.C,v $
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Robert Hoekstra, SNL
//
// Creation Date : 5/10/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.82.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:52 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <stdio.h>
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_UTL_Param.h>
#include <N_UTL_Expression.h>
#include <N_PDS_Comm.h>
#include <N_ERH_ErrorMgr.h>

// ---------  Other Includes  -----------

// ---------  Helper Functions  -----------


//-----------------------------------------------------------------------------
// Function      : N_UTL_ParamData::N_UTL_ParamData
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, Electrical Modeling
// Creation Date : 12/23/10
//-----------------------------------------------------------------------------
N_UTL_ParamData::N_UTL_ParamData(const string & t, const vector<string> & v)
  : tag_(t), type_(STR_VEC), simVarContext(UNDEFINED), extIndex1(-1), extIndex2(-1)
{
  vals.stringVec_ = new vector<string>(v);
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_ParamData::N_UTL_ParamData
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, Electrical Modeling
// Creation Date : 12/23/10
//-----------------------------------------------------------------------------
N_UTL_ParamData::N_UTL_ParamData(const string & t, const vector<double> & v)
  : tag_(t), type_(DBLE_VEC), simVarContext(UNDEFINED), extIndex1(-1), extIndex2(-1)
{
  vals.dbleVec_ = new vector<double>(v);
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_ParamData::N_UTL_ParamData
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Christy Warrender, Cognitive Modeling
// Creation Date : 3/16/11
//-----------------------------------------------------------------------------
N_UTL_ParamData::N_UTL_ParamData(const string & t, const vector<int> & v)
  : tag_(t), type_(INT_VEC), simVarContext(UNDEFINED), extIndex1(-1), extIndex2(-1)
{
  vals.intVec_ = new vector<int>(v);
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_ParamData::N_UTL_ParamData
// Purpose       :
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
N_UTL_ParamData::N_UTL_ParamData(const string & t, const N_UTL_Expression & v)
      : tag_(t), type_(EXPR), simVarContext(UNDEFINED), extIndex1(-1), extIndex2(-1)
{
  vals.exprVal_ = new N_UTL_Expression(v);
}


//-----------------------------------------------------------------------------
// Function      : N_UTL_ParamData::~N_UTL_ParamData
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, Electrical Modeling
// Creation Date : 12/23/10
//-----------------------------------------------------------------------------
void N_UTL_ParamData::freeAllocatedStorage()
{
  // delete any allocated space held by this object
  if (type_ == STR)
    delete vals.stringVal_;
  if (type_ == STR_VEC )
    delete vals.stringVec_;
  if (type_ == DBLE_VEC )
    delete vals.dbleVec_;
  if (type_ == INT_VEC )
    delete vals.intVec_;
  if (type_ == EXPR)
    delete vals.exprVal_;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_ParamData::~N_UTL_ParamData
// Purpose       :
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
N_UTL_ParamData::~N_UTL_ParamData()
{
  freeAllocatedStorage();
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::N_UTL_Param
// Purpose       : Constructor with no inititalization
// Special Notes :
// Scope         :
// Creator       : Lon Waters, SNL
// Creation Date : 5/16/01
//-----------------------------------------------------------------------------
N_UTL_Param::N_UTL_Param()
  : data_( new N_UTL_ParamData )
{
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::N_UTL_Param
// Purpose       : Constructor initialized to string
// Special Notes :
// Scope         :
// Creator       : Lon Waters, SNL
// Creation Date : 5/16/01
//-----------------------------------------------------------------------------
N_UTL_Param::N_UTL_Param(const string & t, const string & v)
  : data_( new N_UTL_ParamData( t, v ) )
{
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::N_UTL_Param
// Purpose       : Constructor initialized to double
// Special Notes :
// Scope         :
// Creator       : Lon Waters, SNL
// Creation Date : 5/16/01
//-----------------------------------------------------------------------------
N_UTL_Param::N_UTL_Param(const string & t, const double & v)
  : data_( new N_UTL_ParamData( t, v ) )
{
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::N_UTL_Param
// Purpose       : Constructor initialized to int
// Special Notes :
// Scope         :
// Creator       : Lon Waters, SNL
// Creation Date : 5/16/01
//-----------------------------------------------------------------------------
N_UTL_Param::N_UTL_Param(const string & t, const int & v)
  : data_( new N_UTL_ParamData( t, v ) )
{
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::N_UTL_Param
// Purpose       : Constructor initialized to long
// Special Notes :
// Scope         :
// Creator       : Derek Barnes, SNL
// Creation Date : 8/01/01
//-----------------------------------------------------------------------------
N_UTL_Param::N_UTL_Param(const string & t, const long & v)
  : data_( new N_UTL_ParamData( t, v ) )
{
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::N_UTL_Param
// Purpose       : Constructor initialized to bool
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/14/06
//-----------------------------------------------------------------------------
N_UTL_Param::N_UTL_Param(const string & t, const bool & v)
  : data_( new N_UTL_ParamData( t, v ) )
{
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::N_UTL_Param
// Purpose       : Constructor initialized to char *
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/14/06
//-----------------------------------------------------------------------------
N_UTL_Param::N_UTL_Param(const string & t, const char * v)
  : data_( new N_UTL_ParamData( t, string(v) ) )
{
}


//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::N_UTL_Param
// Purpose       : Constructor initialized to vector<string>
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, Electrical Modeling, 1445
// Creation Date : 12/23/10
//-----------------------------------------------------------------------------
N_UTL_Param::N_UTL_Param(const string & t, const vector<string> & v)
  : data_( new N_UTL_ParamData( t, v ) )
{
}


//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::N_UTL_Param
// Purpose       : Constructor initialized to vector<double>
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, Electrical Modeling, 1445
// Creation Date : 12/23/10
//-----------------------------------------------------------------------------
N_UTL_Param::N_UTL_Param(const string & t, const vector<double> & v)
  : data_( new N_UTL_ParamData( t, v ) )
{
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::N_UTL_Param
// Purpose       : Constructor initialized to vector<int>
// Special Notes :
// Scope         :
// Creator       : Christy Warrender, Cognitive Modeling, 1462
// Creation Date : 3/17/11
//-----------------------------------------------------------------------------
N_UTL_Param::N_UTL_Param(const string & t, const vector<int> & v)
  : data_( new N_UTL_ParamData( t, v ) )
{
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::N_UTL_Param
// Purpose       : Constructor initialized to expression
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/08/04
//-----------------------------------------------------------------------------
N_UTL_Param::N_UTL_Param(const string & t, const N_UTL_Expression & v)
  : data_( new N_UTL_ParamData( t, v ) )
{
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::N_UTL_Param
// Purpose       : Copy Constructor
// Special Notes :
// Scope         :
// Creator       : Lon Waters, SNL
// Creation Date : 5/17/01
//-----------------------------------------------------------------------------
N_UTL_Param::N_UTL_Param(N_UTL_Param const& rhsParam)
{
  data_ = new N_UTL_ParamData;

  data_->tag_ = rhsParam.data_->tag_;
  data_->type_ = rhsParam.data_->type_;
  data_->simVarContext = rhsParam.data_->simVarContext;
  data_->extIndex1 = rhsParam.data_->extIndex1;
  data_->extIndex2 = rhsParam.data_->extIndex2;
  data_->expressionDataRCP = rhsParam.data_->expressionDataRCP;
  data_->qualifiedParameterOrFunctionName = rhsParam.data_->qualifiedParameterOrFunctionName;
  data_->fixedValue = rhsParam.data_->fixedValue;

  if (data_->type_ == STR)
    data_->vals.stringVal_ = new string(*(rhsParam.data_->vals.stringVal_));
  else if (data_->type_ == DBLE)
    data_->vals.dbleVal_ = rhsParam.data_->vals.dbleVal_;
  else if (data_->type_ == INT)
    data_->vals.intVal_ = rhsParam.data_->vals.intVal_;
  else if (data_->type_ == LNG)
    data_->vals.lngVal_ = rhsParam.data_->vals.lngVal_;
  else if (data_->type_ == BOOL)
    data_->vals.boolVal_ = rhsParam.data_->vals.boolVal_;
  else if (data_->type_ == EXPR)
    data_->vals.exprVal_ = new N_UTL_Expression(*(rhsParam.data_->vals.exprVal_));
  else if (data_->type_ == STR_VEC)
    data_->vals.stringVec_ = new vector<string>(*(rhsParam.data_->vals.stringVec_));
  else if (data_->type_ == DBLE_VEC)
    data_->vals.dbleVec_ = new vector<double>(*(rhsParam.data_->vals.dbleVec_));
  else if (data_->type_ == INT_VEC)
    data_->vals.intVec_ = new vector<int>(*(rhsParam.data_->vals.intVec_));

}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::operator=
// Purpose       : assignment operator
// Special Notes :
// Scope         :
// Creator       : Lon Waters, SNL
// Creation Date : 05/18/01
//-----------------------------------------------------------------------------
N_UTL_Param & N_UTL_Param::operator=(N_UTL_Param const& rhsParam)
{
  // if this currently contains a string or an expression, then
  // we're going to create a new one.  So we need to delete
  // the old one first.
  data_->freeAllocatedStorage();
  /*
  if( data_->type_ == STR )
  {
    delete data_->vals.stringVal_;
  }
  else if( data_->type_ == EXPR )
  {
    delete data_->vals.exprVal_;
  }
  */

  data_->tag_ = rhsParam.data_->tag_;
  data_->type_ = rhsParam.data_->type_;
  data_->simVarContext = rhsParam.data_->simVarContext;
  data_->extIndex1 = rhsParam.data_->extIndex1;
  data_->extIndex2 = rhsParam.data_->extIndex2;
  data_->expressionDataRCP = rhsParam.data_->expressionDataRCP;
  data_->qualifiedParameterOrFunctionName = rhsParam.data_->qualifiedParameterOrFunctionName;
  data_->fixedValue = rhsParam.data_->fixedValue;

  if (data_->type_ == STR)
    data_->vals.stringVal_ = new string(*(rhsParam.data_->vals.stringVal_));
  else if (data_->type_ == DBLE)
    data_->vals.dbleVal_ = rhsParam.data_->vals.dbleVal_;
  else if (data_->type_ == INT)
    data_->vals.intVal_ = rhsParam.data_->vals.intVal_;
  else if (data_->type_ == LNG)
    data_->vals.lngVal_ = rhsParam.data_->vals.lngVal_;
  else if (data_->type_ == BOOL)
    data_->vals.boolVal_ = rhsParam.data_->vals.boolVal_;
  else if (data_->type_ == EXPR)
    data_->vals.exprVal_ = new N_UTL_Expression(*(rhsParam.data_->vals.exprVal_));
  else if (data_->type_ == STR_VEC)
    data_->vals.stringVec_ = new vector<string>(*(rhsParam.data_->vals.stringVec_));
  else if (data_->type_ == DBLE_VEC)
    data_->vals.dbleVec_ = new vector<double>(*(rhsParam.data_->vals.dbleVec_));
  else if (data_->type_ == INT_VEC)
    data_->vals.intVec_ = new vector<int>(*(rhsParam.data_->vals.intVec_));

   return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::~N_UTL_Param
// Purpose       : Destructor
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
N_UTL_Param::~N_UTL_Param()
{
  if ( data_ != NULL) delete data_;
  data_ = NULL;
}

//-----------------------------------------------------------------------------
// Function      : upperString
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
const string N_UTL_Param::usVal() const
{
  if (data_->type_ == STR)
    return ExtendedString( *data_->vals.stringVal_ ).toUpper();
  return string("");
}

//-----------------------------------------------------------------------------
// Function      : getLowerString
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
const string N_UTL_Param::lsVal() const
{
  if (data_->type_ == STR)
    return ExtendedString( *data_->vals.stringVal_ ).toLower();
  return string("");
}

//-----------------------------------------------------------------------------
// Function      : get type
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Lon Waters, SNL
// Creation Date : 5/16/01
//-----------------------------------------------------------------------------
const int & N_UTL_Param::getType() const
{
  return data_->type_;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::operator!=
// Purpose       : "!=" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
bool N_UTL_Param::operator!=( N_UTL_Param const& rhsParam ) const
{
  return !(this->operator==(rhsParam));
}


//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::deepCompare
// Purpose       : Comapre two N_UTL_Parms deeper than the equity operator
// Special Notes : The equity operator just checks that the TAGS are equal
//                 this compares the value if that makes sense (i.e. INT to INT
//                 or STRING to STRING, but not INT to STRING)
// Scope         : public
// Creator       : Rich Schiek, Electrical Systems Modeling
// Creation Date : 2/06/2012
//-----------------------------------------------------------------------------
// deepCompare -- compare TAG and Value if TAGS are the same
bool N_UTL_Param::deepCompare(N_UTL_Param const & rhsParam) const
{
  bool match=false;


#ifndef HAVE_STRCASECMP
  bool tagMatch = ( uTag() == rhsParam.uTag() );
#else
  bool tagMatch = ( strcasecmp(( data_->tag_.c_str() ), ( rhsParam.data_->tag_.c_str() ) )==0 );
#endif

  if( tagMatch )
  {
    // tag name matches, so check type
    if( getType() == rhsParam.getType() )
    {
      // type matches so it's sensible to check value
      if( sVal() == rhsParam.sVal())
      {
        match=true;
      }
    }
  }

  return match;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::set
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/17/06
//-----------------------------------------------------------------------------
N_UTL_Param & N_UTL_Param::set( const string & tag, const N_UTL_Param & p )
{
  data_->tag_ = tag;
  setVal(p);
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::set
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
N_UTL_Param & N_UTL_Param::set( const string & tag, const string & val )
{
  data_->tag_ = tag;
  setVal(val);
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::set
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
N_UTL_Param & N_UTL_Param::set( const string & tag, const double & val )
{
  data_->tag_ = tag;
  setVal(val);
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::set
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
N_UTL_Param & N_UTL_Param::set( const string & tag, const int & val )
{
  data_->tag_ = tag;
  setVal(val);
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::set
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
N_UTL_Param & N_UTL_Param::set( const string & tag, const long & val )
{
  data_->tag_ = tag;
  setVal(val);
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::set
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/14/06
//-----------------------------------------------------------------------------
N_UTL_Param & N_UTL_Param::set( const string & tag, const bool & val )
{
  data_->tag_ = tag;
  setVal(val);
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::set
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/14/06
//-----------------------------------------------------------------------------
N_UTL_Param & N_UTL_Param::set( const string & tag, const char * val )
{
  data_->tag_ = tag;
  setVal(string(val));
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::set
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, Electrical Modeling
// Creation Date : 12/23/10
//-----------------------------------------------------------------------------
N_UTL_Param & N_UTL_Param::set(const string & tag, const vector<string> & val)
{
  data_->tag_ = tag;
  setVal(val);
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::set
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, Electrical Modeling
// Creation Date : 12/23/10
//-----------------------------------------------------------------------------
N_UTL_Param & N_UTL_Param::set(const string & tag, const vector<double> & val)
{
  data_->tag_ = tag;
  setVal(val);
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::set
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Christy Warrender, Cognitive Modeling
// Creation Date : 3/16/11
//-----------------------------------------------------------------------------
N_UTL_Param & N_UTL_Param::set(const string & tag, const vector<int> & val)
{
  data_->tag_ = tag;
  setVal(val);
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::set
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/08/04
//-----------------------------------------------------------------------------
N_UTL_Param & N_UTL_Param::set( const string & tag, const N_UTL_Expression & val )
{
  data_->tag_ = tag;
  setVal(val);
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::setTag
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
void N_UTL_Param::setTag( const string & tag )
{
  data_->tag_ = tag;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::setVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/17/06
//-----------------------------------------------------------------------------
void N_UTL_Param::setVal( const N_UTL_Param & p )
{
  int typ = p.getType();

  if (typ == STR)
    setVal(p.sVal());
  else if (typ == DBLE)
    setVal(p.dVal());
  else if (typ == INT)
    setVal(p.iVal());
  else if (typ == LNG)
    setVal(p.lVal());
  else if (typ == BOOL)
    setVal(p.bVal());
  else if (typ == EXPR)
    setVal(p.eVal());
  else
  {
    string msg = "N_UTL_Param::setVal: unsupported type in setval(param)";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::setVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
void N_UTL_Param::setVal( const string & val )
{
  if (data_->type_ == EXPR)
    delete data_->vals.exprVal_;
  if (data_->type_ == STR)
    *(data_->vals.stringVal_) = val;
  else
    data_->vals.stringVal_ = new string (val);
  data_->type_ = STR;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::setVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
void N_UTL_Param::setVal( const double & val )
{
  if (data_->type_ == EXPR)
    delete data_->vals.exprVal_;
  if (data_->type_ == STR)
    delete data_->vals.stringVal_;
  data_->type_ = DBLE;
  data_->vals.dbleVal_ = val;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::setVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
void N_UTL_Param::setVal( const int & val )
{
  if (data_->type_ == EXPR)
    delete data_->vals.exprVal_;
  if (data_->type_ == STR)
    delete data_->vals.stringVal_;
  data_->type_ = INT;
  data_->vals.intVal_ = val;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::setVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
void N_UTL_Param::setVal( const long & val )
{
  if (data_->type_ == EXPR)
    delete data_->vals.exprVal_;
  if (data_->type_ == STR)
    delete data_->vals.stringVal_;
  data_->type_ = LNG;
  data_->vals.lngVal_ = val;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::setVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/14/06
//-----------------------------------------------------------------------------
void N_UTL_Param::setVal( const bool & val )
{
  if (data_->type_ == EXPR)
    delete data_->vals.exprVal_;
  if (data_->type_ == STR)
    delete data_->vals.stringVal_;
  data_->type_ = BOOL;
  data_->vals.boolVal_ = val;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::setVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/14/06
//-----------------------------------------------------------------------------
void N_UTL_Param::setVal( const char * val )
{
  string str(val);
  setVal(string(val));
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::setVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, Electrical Modeling
// Creation Date : 12/23/10
//-----------------------------------------------------------------------------
void N_UTL_Param::setVal(const vector<string> & val)
{
  // delete any allocated space held by this object
  data_->freeAllocatedStorage();
  data_->vals.stringVec_ = new vector<string>(val);
  data_->type_ = STR_VEC;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::setVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, Electrical Modeling
// Creation Date : 12/23/10
//-----------------------------------------------------------------------------
void N_UTL_Param::setVal(const vector<double> & val)
{
  // delete any allocated space held by this object
  data_->freeAllocatedStorage();
  data_->vals.dbleVec_ = new vector<double>(val);
  data_->type_ = DBLE_VEC;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::setVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Christy Warrender, Cognitive Modeling
// Creation Date : 3/16/11
//-----------------------------------------------------------------------------
void N_UTL_Param::setVal(const vector<int> & val)
{
  // delete any allocated space held by this object
  data_->freeAllocatedStorage();
  data_->vals.intVec_ = new vector<int>(val);
  data_->type_ = INT_VEC;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::setVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
void N_UTL_Param::setVal( const N_UTL_Expression & val )
{
  // delete any allocated space held by this object
  data_->freeAllocatedStorage();
  data_->vals.exprVal_ = new N_UTL_Expression (val);
  data_->type_ = EXPR;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::tag
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
const string & N_UTL_Param::tag() const
{
  return data_->tag_;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::sVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
const string & N_UTL_Param::sVal() const
{
  char buf[32];
  string rVal;

  if (data_->type_ == STR)
  {
    return *(data_->vals.stringVal_);
  }
  else if (data_->type_ == INT)
  {
    sprintf(buf, "%d", data_->vals.intVal_);
    rVal = buf;
  }
  else if (data_->type_ == LNG)
  {
    sprintf(buf, "%ld", data_->vals.lngVal_);
    rVal = buf;
  }
  else if (data_->type_ == DBLE)
  {
    sprintf(buf, "%g", data_->vals.dbleVal_);
    rVal = buf;
  }
  else if (data_->type_ == BOOL)
  {
    if (data_->vals.boolVal_)
      rVal = "TRUE";
    else
      rVal = "FALSE";
  }
  else if (data_->type_ == STR_VEC)
  {
    rVal = "STR_VEC";
  }
  else if (data_->type_ == DBLE_VEC)
  {
    rVal = "DBLE_VEC";
  }
  else if (data_->type_ == INT_VEC)
  {
    rVal = "INT_VEC";
  }
  else if (data_->type_ == EXPR)
  {
    rVal = data_->vals.exprVal_->get_expression();
  }
  else
  {
    rVal = "";
  }

//DNS: this used to return a static string, which did not work for multithreaded runs
  data_->sReturnVal = rVal;
  return data_->sReturnVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::dVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
const double & N_UTL_Param::dVal() const
{
  double val;

  if (data_->type_ != DBLE)
  {
    if (data_->type_ == STR)
    {
      //ExtendedString tmp( *data_->vals.stringVal_ );
      string & tmp = ( *data_->vals.stringVal_ );
      //if (tmp.isValue())
      if (N_UTL::isValue(tmp))
      {
        //val = tmp.Value();
        val = N_UTL::Value(tmp);
        delete data_->vals.stringVal_;
      }
      else
      {
        string msg("N_UTL_Param::dVal: attempt to assign value for: ");
        msg += data_->tag_;
        msg += " from string: ";
        msg += tmp;
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg );
      }
    }
    else if (data_->type_ == INT)
    {
      val = data_->vals.intVal_;
    }
    else if (data_->type_ == LNG)
    {
      val = data_->vals.lngVal_;
    }
    else if (data_->type_ == BOOL)
    {
      string msg("N_UTL_Param::dVal: attempt to assign value for: ");
      msg += data_->tag_;
      msg += " from bool";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg );
    }
    else if (data_->type_ == EXPR)
    {
      // Only if this param expression is truely constant can it be converted to a double
      // else it is a fatal error, in the parser most likely
      if (data_->vals.exprVal_->num_vars() == 0)
      {
        (void) data_->vals.exprVal_->evaluateFunction (val);
        delete data_->vals.exprVal_;
      }
      else
      {
        string msg("N_UTL_Param::dVal: attempt to evaluate expression: ");
        msg += data_->vals.exprVal_->get_expression();
        msg += ", which contains unknowns";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg );
      }
    }
    else
    {
      val = 0;
    }
    data_->vals.dbleVal_ = val;
    data_->type_ = DBLE;
  }
  return data_->vals.dbleVal_;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::iVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
const int & N_UTL_Param::iVal() const
{
  int val;
  double dVal;

  if (data_->type_ != INT)
  {
    if (data_->type_ == STR)
    {
      //ExtendedString tmp( *data_->vals.stringVal_ );
      string & tmp = ( *data_->vals.stringVal_ );
      //if (tmp.isInt())
      if (N_UTL::isInt(tmp))
      {
        //val = tmp.Ival();
        val = N_UTL::Ival(tmp);
        delete data_->vals.stringVal_;
      }
      //else if (tmp.isValue())
      else if (N_UTL::isValue(tmp))
      {
        //val = static_cast<int>(tmp.Value());
        val = static_cast<int>(N_UTL::Value(tmp));
        delete data_->vals.stringVal_;
      }
      else
      {
        string msg("N_UTL_Param::iVal: attempt to assign value for: ");
        msg += data_->tag_;
        msg += " from string: ";
        msg += tmp;
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg );
      }
    }
    else if (data_->type_ == DBLE)
    {
      val = static_cast<int> (data_->vals.dbleVal_);
    }
    else if (data_->type_ == LNG)
    {
      val = static_cast<int> (data_->vals.lngVal_);
    }
    else if (data_->type_ == BOOL)
    {
      string msg("N_UTL_Param::iVal: attempt to assign value for: ");
      msg += data_->tag_;
      msg += " from bool";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg );
    }
    else if (data_->type_ == EXPR)
    {
      if (data_->vals.exprVal_->num_vars() == 0)
      {
        (void) data_->vals.exprVal_->evaluateFunction (dVal);
        val = static_cast<int> (dVal);
        delete data_->vals.exprVal_;
      }
      else
      {
        string msg("N_UTL_Param::iVal: attempt to evaluate expression: ");
        msg += data_->vals.exprVal_->get_expression();
        msg += ", which contains unknowns";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg );
      }
    }
    else
    {
      val = 0;
    }
    data_->vals.intVal_ = val;
    data_->type_ = INT;
  }
  return data_->vals.intVal_;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::lVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
const long & N_UTL_Param::lVal() const
{
  long val;
  double dVal;

  if (data_->type_ != LNG)
  {
    if (data_->type_ == STR)
    {
      //ExtendedString tmp( *data_->vals.stringVal_ );
      string & tmp = ( *data_->vals.stringVal_ );
      //if (tmp.isInt())
      if (N_UTL::isInt(tmp))
      {
        //val = tmp.Ival();
        val = N_UTL::Ival(tmp);
        delete data_->vals.stringVal_;
      }
      //else if (tmp.isValue())
      else if (N_UTL::isValue(tmp))
      {
        //val = static_cast<int>(tmp.Value());
        val = static_cast<int>(N_UTL::Value(tmp));
        delete data_->vals.stringVal_;
      }
      else
      {
        string msg("N_UTL_Param::iVal: attempt to assign value for: ");
        msg += data_->tag_;
        msg += " from string: ";
        msg += tmp;
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg );
      }
    }
    else if (data_->type_ == DBLE)
    {
      val = static_cast<long> (data_->vals.dbleVal_);
    }
    else if (data_->type_ == INT)
    {
      val = static_cast<long> (data_->vals.intVal_);
    }
    else if (data_->type_ == BOOL)
    {
      string msg("N_UTL_Param::lVal: attempt to assign value for: ");
      msg += data_->tag_;
      msg += " from bool";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg );
    }
    else if (data_->type_ == EXPR)
    {
      if (data_->vals.exprVal_->num_vars() == 0)
      {
        (void) data_->vals.exprVal_->evaluateFunction (dVal);
        val = static_cast<long> (dVal);
        delete data_->vals.exprVal_;
      }
      else
      {
        string msg("N_UTL_Param::lVal: attempt to evaluate expression: ");
        msg += data_->vals.exprVal_->get_expression();
        msg += ", which contains unknowns";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg );
      }
    }
    else
    {
      val = 0;
    }
    data_->vals.lngVal_ = val;
    data_->type_ = LNG;
  }
  return data_->vals.lngVal_;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::bVal
// Purpose       : Return booleen value of param
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/17/04
//-----------------------------------------------------------------------------
bool N_UTL_Param::bVal() const
{
  bool rVal;
  if (data_->type_ == DBLE)
    rVal = (data_->vals.dbleVal_ != 0);
  else if (data_->type_ == INT)
  {
    rVal = (data_->vals.intVal_ != 0);
  }
  else if (data_->type_ == LNG)
  {
    rVal = (data_->vals.lngVal_ != 0);
  }
  else if (data_->type_ == BOOL)
  {
    rVal = data_->vals.boolVal_;
  }
  else if (data_->type_ == STR)
  {
    //ExtendedString tmp( *data_->vals.stringVal_ );
    string & tmp = ( *data_->vals.stringVal_ );
    //if (tmp.isBool())
    if (N_UTL::isBool(tmp))
    {
      //rVal = tmp.Bval();
      rVal = N_UTL::Bval(tmp);
    }
    else
    {
      string msg("N_UTL_Param::bVal: attempt to assign value for: ");
      msg += data_->tag_;
      msg += " from string: ";
      msg += tmp;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg );
    }
  }
  else if (data_->type_ == EXPR)
  {
    if (data_->vals.exprVal_->num_vars() == 0)
    {
      double dVal;
      (void) data_->vals.exprVal_->evaluateFunction (dVal);
      rVal = (dVal != 0);
    }
    else
    {
      string msg("N_UTL_Param::lVal: attempt to evaluate expression: ");
      msg += data_->vals.exprVal_->get_expression();
      msg += ", which contains unknowns";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg );
    }
  }
  return rVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::sVecPtr
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, Electrical Modeling
// Creation Date : 12/23/10
//-----------------------------------------------------------------------------
vector<string> * N_UTL_Param::sVecPtr()
{
  if (data_->type_ != STR_VEC)
  {
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, string ("N_UTL_Param::sVecPtr: attempt to return vector<string> pointer for non-vector<string> param"));
  }
  return data_->vals.stringVec_;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::sVecVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, Electrical Modeling
// Creation Date : 12/23/10
//-----------------------------------------------------------------------------
const vector<string> & N_UTL_Param::sVecVal() const
{
  if (data_->type_ != STR_VEC)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, string ("N_UTL_Param::sVecVal: attempt to return vector<string> for non vector<string> param"));
  }
  return *(data_->vals.stringVec_);
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::dVecPtr
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, Electrical Modeling
// Creation Date : 12/23/10
//-----------------------------------------------------------------------------
vector<double> * N_UTL_Param::dVecPtr()
{
  if (data_->type_ != DBLE_VEC)
  {
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, string ("N_UTL_Param::dVecPtr: attempt to return vector<double> pointer for non-vector<double> param"));
  }
  return data_->vals.dbleVec_;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::dVecVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, Electrical Modeling
// Creation Date : 12/23/10
//-----------------------------------------------------------------------------
const vector<double> & N_UTL_Param::dVecVal() const
{
  if (data_->type_ != DBLE_VEC)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, string ("N_UTL_Param::dVecVal: attempt to return vector<double> for non vector<double> param"));
  }
  return *(data_->vals.dbleVec_);
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::iVecPtr
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Christy Warrender, Cognitive Modeling
// Creation Date : 3/16/11
//-----------------------------------------------------------------------------
vector<int> * N_UTL_Param::iVecPtr()
{
  if (data_->type_ != INT_VEC)
  {
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, string ("N_UTL_Param::iVecPtr: attempt to return vector<int> pointer for non-vector<int> param"));
  }
  return data_->vals.intVec_;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::iVecVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Christy Warrender, Cognitive Modeling
// Creation Date : 3/16/11
//-----------------------------------------------------------------------------
const vector<int> & N_UTL_Param::iVecVal() const
{
  if (data_->type_ != INT_VEC)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, string ("N_UTL_Param::iVecVal: attempt to return vector<int> for non vector<int> param"));
  }
  return *(data_->vals.intVec_);
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::eVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/17/04
//-----------------------------------------------------------------------------
N_UTL_Expression * N_UTL_Param::ePtr()
{
  if (data_->type_ != EXPR)
  {
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, string ("N_UTL_Param::ePtr: attempt to return expression pointer for non-expression param"));
  }
  return data_->vals.exprVal_;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::eVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/17/04
//-----------------------------------------------------------------------------
N_UTL_Expression * N_UTL_Param::ePtr() const
{
  if (data_->type_ != EXPR)
  {
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, string ("N_UTL_Param::ePtr: attempt to return expression pointer for non-expression param"));
  }
  return data_->vals.exprVal_;
}


//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::eVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/17/04
//-----------------------------------------------------------------------------
const N_UTL_Expression & N_UTL_Param::eVal() const
{
  if (data_->type_ != EXPR)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, string ("N_UTL_Param::eVal: attempt to return expression for non-expression param"));
  }
  return *(data_->vals.exprVal_);
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::hasExpressionValue
// Purpose       : Determine if the N_UTL_Param value is an expression.
// Special Notes : This checks the value, not the tag.
// Scope         :
// Creator       : Lon Waters, SNL
// Creation Date : 10/08/01
//-----------------------------------------------------------------------------
bool N_UTL_Param::hasExpressionValue() const
{
  if ( data_->type_ == EXPR ) {
    return true;
  }
  if ( data_->type_ != STR )
  {
    return false;
  }
  else if ( (*data_->vals.stringVal_)[0] != '{' ||
            (*data_->vals.stringVal_)[data_->vals.stringVal_->size()-1] != '}' )
  {
    return false;
  }
  else
  {
    // The value is a string surrounded by braces: {...}
    return true;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::hasExpressionTag
// Purpose       : Determine if the N_UTL_Param value is an expression.
// Special Notes : This checks the tag.
// Scope         :
// Creator       : Eric Keiter, SNL, Computational Sciences
// Creation Date : 08/15/04
//-----------------------------------------------------------------------------
bool N_UTL_Param::hasExpressionTag () const
{
  if ( data_->tag_[0] != '{' ||
       data_->tag_[data_->tag_.size()-1] != '}' )
  {
    return false;
  }
  else
  {
    // The value is a string surrounded by braces: {...}
    return true;
  }
}

//----------------------------------------------------------------------------
// Function       : N_UTL_Param::isQuoted
// Purpose        : Return true if the parameter value is enclosed in double
//                  quotes.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/20/2002
//----------------------------------------------------------------------------
bool N_UTL_Param::isQuoted()
{
  if ( data_->type_ != STR )
  {
    return false;
  }
  else if ( (*data_->vals.stringVal_)[0] != '"' ||
            (*data_->vals.stringVal_)[data_->vals.stringVal_->size()-1] != '"' )
  {
    return false;
  }
  else
  {
    // The value is a string surrounded by double quotes: "..."
    return true;
  }
}

//----------------------------------------------------------------------------
// Function       : N_UTL_Param::isNumeric
// Purpose        : Checks the value of the parameter to see if it is a
//                  legal real or integer numeric value.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 12/06/2002
//----------------------------------------------------------------------------
bool N_UTL_Param::isNumeric() const
{
  // Only do the checking if the parameter is string valued, return true
  // if it is already real or integer valued.
  if ( data_->type_ == DBLE || data_->type_ == INT || data_->type_ == LNG )
    return true;
  if ( data_->type_ == EXPR || data_->type_ == BOOL)
    return false;
  else if ( data_->type_ == STR)
    return N_UTL::isValue( (*data_->vals.stringVal_) );
    //return ExtendedString(*data_->vals.stringVal_).isValue();
  else
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,string("N_UTL_Param::isNumeric: unknown type"));
  }
  return true;
}

//----------------------------------------------------------------------------
// Function       : N_UTL_Param::isInteger
// Purpose        : Checks the value of the parameter to see if it is a
//                  legal integer numeric value.
// Special Notes  :
// Scope          :
// Creator        : Dave Shirley
// Creation Date  : 03/10/2006
//----------------------------------------------------------------------------
bool N_UTL_Param::isInteger() const
{
  // Only do the checking if the parameter is string valued, return true
  // if it is already real or integer valued.
  if ( data_->type_ == INT || data_->type_ == LNG )
    return true;
  else if ( data_->type_ == EXPR || data_->type_ == BOOL)
    return false;
  else if ( data_->type_ == DBLE)
    return true;
  else if ( data_->type_ == STR)
    return N_UTL::isInt(*data_->vals.stringVal_);
    //return ExtendedString(*data_->vals.stringVal_).isInt();
  else
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, string("N_UTL_Param::isInteger: unknown type"));
  }
  return true;
}

//----------------------------------------------------------------------------
// Function       : N_UTL_Param::isBool
// Purpose        : Checks the value of the parameter to see if it is a
//                  legal bool value.
// Special Notes  :
// Scope          :
// Creator        : Dave Shirley
// Creation Date  : 03/10/2006
//----------------------------------------------------------------------------
bool N_UTL_Param::isBool() const
{
  // Only do the checking if the parameter is string valued, return true
  // if it is already real or integer valued.
  if ( data_->type_ == DBLE || data_->type_ == INT || data_->type_ == LNG || data_->type_ == BOOL)
    return true;
  if ( data_->type_ == EXPR)
    return false;
  else if ( data_->type_ == STR)
    return N_UTL::isBool(*data_->vals.stringVal_);
    //return ExtendedString(*data_->vals.stringVal_).isBool();
  else
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,string("N_UTL_Param::isBool: unknown type"));
  }
  return true;
}

//----------------------------------------------------------------------------
// Function       : N_UTL_Param::setTimeDependent
// Purpose        : Set the flag indicating this is a time dependent parameter
//                  to the given value.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/19/2002
//----------------------------------------------------------------------------
void N_UTL_Param::setTimeDependent( bool const& timeDependent )
{
  string exp;

  if (data_->type_ == EXPR && timeDependent)
    return;
  if (data_->type_ != EXPR && !timeDependent)
    return;
  if (data_->type_ != STR)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL,string("N_UTL_Param::setTimeDependent: attempt to convert non-string to expression"));
  }
  if (!timeDependent)
    return;
  exp = *(data_->vals.stringVal_);
  delete data_->vals.stringVal_;
  data_->type_ = EXPR;
  data_->vals.exprVal_ = new N_UTL_Expression(exp);
  return;
}

//----------------------------------------------------------------------------
// Function       : N_UTL_Param::isTimeDependent
// Purpose        : Return the value of the flag indicating whether the
//                  parameter is time dependent.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/19/2002
//----------------------------------------------------------------------------
bool N_UTL_Param::isTimeDependent() const
{
  if (data_->type_ == EXPR)
    return true;
  return false;
}


//----------------------------------------------------------------------------
// Function       : N_UTL_Param::setSimContextAndData
// Purpose        : method for setting context and extra data.  Only need
//                  three combinations as we figure out the context and
//                  extra data at the same time and then set it.
// Special Notes  :
// Scope          :
// Creator        : Rich Schiek, Electrical Systems Modeling, SNL
// Creation Date  : 11/06/2012
//----------------------------------------------------------------------------
void N_UTL_Param::setSimContextAndData( SimulatorVariableContext aSimContext)
{
  data_->simVarContext = aSimContext;
};

//----------------------------------------------------------------------------
// Function       : N_UTL_Param::setSimContextAndData
// Purpose        : method for setting context and extra data.  Only need
//                  three combinations as we figure out the context and
//                  extra data at the same time and then set it.
// Special Notes  :
// Scope          :
// Creator        : Rich Schiek, Electrical Systems Modeling, SNL
// Creation Date  : 10/04/2012
//----------------------------------------------------------------------------
void N_UTL_Param::setSimContextAndData( SimulatorVariableContext aSimContext, int extraIndex1)
{
  data_->simVarContext = aSimContext;
  data_->extIndex1=extraIndex1;
};

void N_UTL_Param::setSimContextAndData( SimulatorVariableContext aSimContext, int extraIndex1, int extraIndex2 )
{
  data_->simVarContext = aSimContext;
  data_->extIndex1=extraIndex1;
  data_->extIndex2=extraIndex2;
};

void N_UTL_Param::setSimContextAndData( SimulatorVariableContext aSimContext, RCP<N_UTL_ExpressionData> expDataRcp)
{
  data_->simVarContext = aSimContext;
  data_->expressionDataRCP = expDataRcp;
}

void N_UTL_Param::setSimContextAndData( SimulatorVariableContext aSimContext, string aName)
{
  data_->simVarContext = aSimContext;
  data_->qualifiedParameterOrFunctionName = aName;
}

void N_UTL_Param::setSimContextAndData( SimulatorVariableContext aSimContext, double aValue )
{
  data_->simVarContext = aSimContext;
  data_->fixedValue = aValue;
}


//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::instance()
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 6/28/01
//-----------------------------------------------------------------------------
Packable * N_UTL_Param::instance() const
{
  return new N_UTL_Param();
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::print
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley
// Creation Date : 12/06/05
//-----------------------------------------------------------------------------
void N_UTL_Param::print()
{
  cout << "Parameter: " << uTag() << "  Type: ";
  if (data_->type_ == STR)
  {
    cout << "String";
    cout << "Value: '" <<  *(data_->vals.stringVal_) << "'";
  }
  else if (data_->type_ == INT)
  {
    cout << "Int";
    cout << "Value: " <<  data_->vals.intVal_;
  }
  else if (data_->type_ == LNG)
  {
    cout << "Long";
    cout << "Value: " << data_->vals.lngVal_;
  }
  else if (data_->type_ == BOOL)
  {
    cout << "Bool";
    cout << "Value: " << data_->vals.boolVal_;
  }
  else if (data_->type_ == DBLE)
  {
    cout << "Double";
    cout << "Value: " << data_->vals.dbleVal_;
  }
  else if (data_->type_ == EXPR)
  {
    cout << "Expression";
    cout << "Value: " << data_->vals.exprVal_->get_expression();
  }
  else
  {
    cout << "Unsupported type";
  }
  cout << " context ";

  switch (getSimContext())
  {
    case UNDEFINED:
      cout << "UNDEFINED";
      break;
    case CONSTANT:
      cout << "CONSTANT";
      break;
    case TEMPERATURE:
      cout << "TEMPERATURE";
      break;
    case TIME_VAR:
      cout << "TIME_VAR";
      break;
    case STEP_SWEEP_VAR:
      cout << "STEP_SWEEP_VAR";
      break;
    case DC_SWEEP_VAR:
      cout << "DC_SWEEP_VAR";
      break;
    case EXPRESSION:
      cout << "EXPRESSION";
      break;
    case VOLTAGE_DIFFERENCE:
      cout << "VOLTAGE_DIFFERENCE";
      break;
    case NODE_OR_DEVICE_NAME:
      cout << "NODE_OR_DEVICE_NAME";
      break;
    case SOLUTION_VAR:
      cout << "SOLUTION_VAR";
      break;
    case STATE_VAR:
      cout << "STATE_VAR";
      break;
    case STORE_VAR:
      cout << "STORE_VAR";
      break;
    case DEVICE_PARAMETER:
      cout << "DEVICE_PARAMETER";
      break;
    case OBJECTIVE_FUNCTION:
      cout << "OBJECTIVE_FUNCTION";
      break;
    case MEASURE_FUNCTION:
      cout << "MEASURE_FUNCTION";
      break;
  }
  cout << endl;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::packedByteCount
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
int N_UTL_Param::packedByteCount() const
{

  int byteCount = 0;

  //tag info
  byteCount += data_->tag_.length() + sizeof(int);

  //type
  byteCount += sizeof(int);

  //value info
  switch( data_->type_ )
  {
    case -1:   break;
    case STR:  byteCount += data_->vals.stringVal_->length() + sizeof(int);
               break;
    case DBLE: byteCount += sizeof(double);
               break;
    case BOOL:
    case INT:  byteCount += sizeof(int);
               break;
    case LNG:  byteCount += sizeof(long);
               break;
    case EXPR: byteCount += data_->vals.exprVal_->get_expression().length() + sizeof(int);
               break;
    case STR_VEC: byteCount += sizeof(int);
                  for (int i=0; i<(int)data_->vals.stringVec_->size(); i++)
                  {  byteCount += (*data_->vals.stringVec_)[i].length() + sizeof(int); }
                  break;
  }

  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::pack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
void N_UTL_Param::pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const
{

  int length;
  string tmp;

  //pack tag
  length = data_->tag_.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( data_->tag_.c_str(), length, buf, bsize, pos );

  //pack type
  comm->pack( &(data_->type_), 1, buf, bsize, pos );

  //pack value
  switch( data_->type_ )
  {
    case -1:   break;
    case STR:  length = data_->vals.stringVal_->length();
               comm->pack( &length, 1, buf, bsize, pos );
               comm->pack( data_->vals.stringVal_->c_str(), length, buf, bsize, pos );
               break;

    case DBLE: comm->pack( &(data_->vals.dbleVal_), 1, buf, bsize, pos );
               break;

    case INT:  comm->pack( &(data_->vals.intVal_), 1, buf, bsize, pos );
               break;

    case BOOL: int i;
               if (data_->vals.boolVal_)
                 i = 1;
               else
                 i = 0;
               comm->pack( &i, 1, buf, bsize, pos );
               break;

    case LNG:  comm->pack( &(data_->vals.lngVal_), 1, buf, bsize, pos );
               break;

    case EXPR: tmp = data_->vals.exprVal_->get_expression();
               length = tmp.length();
               comm->pack( &length, 1, buf, bsize, pos );
               comm->pack( tmp.c_str(), length, buf, bsize, pos );
               break;

    case STR_VEC: length = (int)data_->vals.stringVec_->size();
                  comm->pack( &length, 1, buf, bsize, pos );
                  for (int i=0; i<(int)data_->vals.stringVec_->size(); i++)
                  {
                    length = (*data_->vals.stringVec_)[i].length();
                    comm->pack( &length, 1, buf, bsize, pos );
                    comm->pack( (*data_->vals.stringVec_)[i].c_str(), length, buf, bsize, pos );
                  }
                  break;

    default:   string msg = "N_UTL_Param::pack: unknown type";
               N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::unpack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
void N_UTL_Param::unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm )
{

  int length;

  //unpack tag
  comm->unpack( pB, bsize, pos, &length, 1 );

  data_->tag_ = string( (pB+pos), length );
  pos += length;

  //unpack type
  comm->unpack( pB, bsize, pos, &(data_->type_), 1 );

  switch( data_->type_ )
  {
    case -1:   break;
    case STR:  comm->unpack( pB, bsize, pos, &length, 1 );
               data_->vals.stringVal_ = new string( (pB+pos), length );
               pos += length;
               break;

    case DBLE: comm->unpack( pB, bsize, pos, &(data_->vals.dbleVal_), 1 );
               break;

    case INT:  comm->unpack( pB, bsize, pos, &(data_->vals.intVal_), 1 );
               break;

    case BOOL: int i;
               comm->unpack( pB, bsize, pos, &i, 1 );
               if (i == 0)
                 data_->vals.boolVal_ = false;
               else
                 data_->vals.boolVal_ = true;
               break;

    case LNG:  comm->unpack( pB, bsize, pos, &(data_->vals.lngVal_), 1 );
               break;

    case EXPR: comm->unpack( pB, bsize, pos, &length, 1 );
               data_->vals.exprVal_ = new N_UTL_Expression (string( (pB+pos), length ));
               pos += length;
               break;

    case STR_VEC: comm->unpack( pB, bsize, pos, &length, 1 );
                  data_->vals.stringVec_ = new vector<string>( length );
                  for (int i=0; i<(int)data_->vals.stringVec_->size(); i++)
                  {
                    comm->unpack( pB, bsize, pos, &length, 1 );
                    (*data_->vals.stringVec_)[i] = string( (pB+pos), length );
                    pos += length;
                  }
                  break;

    default:   string msg = "N_UTL_Param::unpack: unknown type";
               N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg );
  }
}


//-----------------------------------------------------------------------------
// Function      : N_UTL_Param::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
ostream & operator << (ostream & os, const N_UTL_Param & p)
{
  os << p.tag() << "\t";

  switch (p.getType())
  {
    case STR:
      os << " STR\t";
      os << p.sVal();
      break;
    case DBLE:
      os << "DBLE\t";
      os << p.dVal();
      break;
    case INT:
      os << " INT\t";
      os << p.iVal();
      break;
    case LNG:
      os << " LNG\t";
      os << p.lVal();
      break;
    case BOOL:
      os << "BOOL\t";
      os << p.bVal();
      break;
    case EXPR:
      os << "EXPR\t";
      os << p.sVal();
      break;

    case STR_VEC:
      os << "STR_VEC\t";
      {
        // extra scope from brackets keeps numElements local to this part of the case statement
        // otherwise one would get compiler warnings that its initialization could be bypassed
        // by later parts of the case statement
        int numElements = p.sVecVal().size();
        for(int i=0; i<numElements; i++)
        {
          os << p.sVecVal()[i] << " ";
        }
      }
      break;
    case INT_VEC:
      os << "INT_VEC\t";
      break;
    case DBLE_VEC:
      os << "DBLE_VEC\t";
      break;
    case DBLE_VEC_IND:
      os << "DBLE_VEC_IND\t";
      break;
    case COMPOSITE:
      os << "COMPOSITE\t";
      break;
  }

  os << endl;

  return os;
}
