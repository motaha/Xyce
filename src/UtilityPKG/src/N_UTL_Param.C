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
// Revision Number: $Revision: 1.99 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
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

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Function      : upperString
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
std::string Param::usVal() const
{
  std::string s;

  if (data_->enumType() == STR) {
    s = getValue<std::string>();
    toUpper(s);
  }

  return s;
}

//-----------------------------------------------------------------------------
// Function      : getLowerString
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
std::string Param::lsVal() const
{
  std::string s;

  if (data_->enumType() == STR) {
    s = getValue<std::string>();
    toLower(s);
  }

  return s;
}

//-----------------------------------------------------------------------------
// Function      : Param::deepCompare
// Purpose       : Comapre two N_UTL_Parms deeper than the equity operator
// Special Notes : The equity operator just checks that the TAGS are equal
//                 this compares the value if that makes sense (i.e. INT to INT
//                 or STRING to STRING, but not INT to STRING)
// Scope         : public
// Creator       : Rich Schiek, Electrical Systems Modeling
// Creation Date : 2/06/2012
//-----------------------------------------------------------------------------
// deepCompare -- compare TAG and Value if TAGS are the same
bool Param::deepCompare(Param const & rhsParam) const
{
  return equal_nocase(tag_, rhsParam.tag_)
    && getType() == rhsParam.getType()
    && stringValue() == rhsParam.stringValue();
}

//-----------------------------------------------------------------------------
// Function      : Param::sVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
std::string Param::stringValue() const
{
  std::ostringstream oss;
  
  if (data_->enumType() == STR)
  {
    oss << getValue<std::string>();
  }
  else if (data_->enumType() == INT)
  {
    oss << getValue<int>();
  }
  else if (data_->enumType() == LNG)
  {
    oss << getValue<long>();
  }
  else if (data_->enumType() == DBLE)
  {
    oss << getValue<double>();
  }
  else if (data_->enumType() == BOOL)
  {
    oss << (getValue<bool>() ? "TRUE" : "FALSE");
  }
  else if (data_->enumType() == STR_VEC)
  {
    oss << "STR_VEC";
  }
  else if (data_->enumType() == DBLE_VEC)
  {
    oss << "DBLE_VEC";
  }
  else if (data_->enumType() == INT_VEC)
  {
    oss << "INT_VEC";
  }
  else if (data_->enumType() == EXPR)
  {
    oss << getValue<Expression>().get_expression();
  }

  return oss.str();
}

//-----------------------------------------------------------------------------
// Function      : Param::getImmutableValue<double>
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
template<>
double Param::getImmutableValue<double>() const
{
  double val;

  if (data_->enumType() != DBLE)
  {
    if (data_->enumType() == STR)
    {
      const std::string & tmp = getValue<std::string>();
      if (isValue(tmp))
      {
        val = Value(tmp);
      }
      else
      {
        Report::DevelFatal() << "Cannot convert " << tmp << " to double for expression " <<  tag_;
      }
    }
    else if (data_->enumType() == INT)
    {
      val = getValue<int>();
    }
    else if (data_->enumType() == LNG)
    {
      val = getValue<long>();
    }
    else if (data_->enumType() == BOOL)
    {
      Report::DevelFatal() << "Cannot convert boolean to double for expression " <<  tag_;
    }
    else if (data_->enumType() == EXPR)
    {
      // Only if this param expression is truely constant can it be converted to a double
      // else it is a fatal error, in the parser most likely
      Expression &expression = const_cast<Expression &>(getValue<Expression>());
      if (expression.num_vars() == 0)
      {
        expression.evaluateFunction(val);
      }
      else
      {
        Report::UserFatal() << "Attempt to evaluate expression " << expression.get_expression() << ", which contains unknowns";
      }
    }
    else
    {
      val = 0;
    }

    const_cast<Param &>(*this).setVal(val);
  }
  
  return getValue<double>();
}

//-----------------------------------------------------------------------------
// Function      : Param::iVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
template<>
int Param::getImmutableValue<int>() const
{
  int val;
  double dVal;

  if (data_->enumType() != INT)
  {
    if (data_->enumType() == STR)
    {
      const std::string & tmp = getValue<std::string>();
      if (isInt(tmp))
      {
        val = Ival(tmp);
      }
      else if (isValue(tmp))
      {
        val = static_cast<int>(Value(tmp));
      }
      else
      {
        Report::DevelFatal() << "Cannot convert " << tmp << " to integer for expression " <<  tag_;
      }
    }
    else if (data_->enumType() == DBLE)
    {
      val = static_cast<int> (getValue<double>());
    }
    else if (data_->enumType() == LNG)
    {
      val = static_cast<int> (getValue<long>());
    }
    else if (data_->enumType() == BOOL)
    {
      Report::DevelFatal() << "Cannot convert boolean to integer for expression " <<  tag_;
    }
    else if (data_->enumType() == EXPR)
    {
      Expression &expression = const_cast<Expression &>(getValue<Expression>());

      if (expression.num_vars() == 0)
      {
        expression.evaluateFunction(dVal);
        val = dVal;
      }
      else
      {
        Report::UserFatal() << "Attempt to evaluate expression " << expression.get_expression() << ", which contains unknowns";
      }
    }
    else
    {
      val = 0;
    }
    const_cast<Param &>(*this).setVal(val);
  }
  return getValue<int>();
}

//-----------------------------------------------------------------------------
// Function      : Param::lVal
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
template<>
long Param::getImmutableValue<long>() const
{
  long val;
  double dVal;

  if (data_->enumType() != LNG)
  {
    if (data_->enumType() == STR)
    {
      const std::string & tmp = getValue<std::string>();
      if (isInt(tmp))
      {
        val = Ival(tmp);
      }
      else if (isValue(tmp))
      {
        val = static_cast<long>(Value(tmp));
      }
      else
      {
        Report::DevelFatal() << "Cannot convert " << tmp << " to long integer for expression " <<  tag_;
      }
    }
    else if (data_->enumType() == DBLE)
    {
      val = static_cast<long> (getValue<double>());
    }
    else if (data_->enumType() == INT)
    {
      val = static_cast<long> (getValue<int>());
    }
    else if (data_->enumType() == BOOL)
    {
      Report::DevelFatal() << "Cannot convert boolean to long integer for expression " <<  tag_;
    }
    else if (data_->enumType() == EXPR)
    {
      Expression &expression = const_cast<Expression &>(getValue<Expression>());

      if (expression.num_vars() == 0)
      {
        expression.evaluateFunction (dVal);
        val = dVal;
      }
      else
      {
        Report::UserFatal() << "Attempt to evaluate expression " << expression.get_expression() << ", which contains unknowns";
      }
    }
    else
    {
      val = 0;
    }
    const_cast<Param &>(*this).setVal(val);
  }
  return getValue<long>();
}

//-----------------------------------------------------------------------------
// Function      : Param::bVal
// Purpose       : Return booleen value of param
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/17/04
//-----------------------------------------------------------------------------
template<>
bool Param::getImmutableValue<bool>() const
{
  bool rVal;
  if (data_->enumType() == DBLE)
  {
    rVal = (getValue<double>() != 0.0);
  }
  else if (data_->enumType() == INT)
  {
    rVal = (getValue<int>() != 0);
  }
  else if (data_->enumType() == LNG)
  {
    rVal = (getValue<long>() != 0);
  }
  else if (data_->enumType() == BOOL)
  {
    rVal = getValue<bool>();
  }
  else if (data_->enumType() == STR)
  {
    const std::string & tmp = getValue<std::string>();
    if (Util::isBool(tmp))
    {
      rVal = Bval(tmp);
    }
    else
    {
      Report::DevelFatal() << "Cannot convert " << tmp << " to boolean for expression " <<  tag_;
    }
  }
  else if (data_->enumType() == EXPR)
  {
    Expression &expression = const_cast<Expression &>(getValue<Expression>());

    if (expression.num_vars() == 0)
    {
      double dVal;
      expression.evaluateFunction (dVal);
      rVal = (dVal != 0);
    }
    else
    {
      Report::UserError() << "Attempt to evaluate expression " << expression.get_expression() << ", which contains unknowns";
    }
  }
  return rVal;
}

//-----------------------------------------------------------------------------
// Function      : Param::hasExpressionValue
// Purpose       : Determine if the Param value is an expression.
// Special Notes : This checks the value, not the tag.
// Scope         :
// Creator       : Lon Waters, SNL
// Creation Date : 10/08/01
//-----------------------------------------------------------------------------
bool Param::hasExpressionValue() const
{
  return data_->enumType() == EXPR
    || (data_->enumType() == STR && hasExpressionTag(getValue<std::string>()));
}

//----------------------------------------------------------------------------
// Function       : Param::isQuoted
// Purpose        : Return true if the parameter value is enclosed in double
//                  quotes.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/20/2002
//----------------------------------------------------------------------------
bool Param::isQuoted()
{
  if ( data_->enumType() == STR )
  {
    const std::string &s = getValue<std::string>();
    if (s[0] == '"' && s[s.size() - 1] == '"' )
    {
      return true;
    }
  }

  return false;
}

//----------------------------------------------------------------------------
// Function       : Param::isNumeric
// Purpose        : Checks the value of the parameter to see if it is a
//                  legal real or integer numeric value.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 12/06/2002
//----------------------------------------------------------------------------
bool Param::isNumeric() const
{
  // Only do the checking if the parameter is string valued, return true
  // if it is already real or integer valued.
  if ( data_->enumType() == DBLE || data_->enumType() == INT || data_->enumType() == LNG )
    return true;
  if ( data_->enumType() == EXPR || data_->enumType() == BOOL)
    return false;
  else if ( data_->enumType() == STR)
    return isValue(getValue<std::string>());
  else
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,std::string("Param::isNumeric: unknown type"));
  }
  return true;
}

//----------------------------------------------------------------------------
// Function       : Param::isInteger
// Purpose        : Checks the value of the parameter to see if it is a
//                  legal integer numeric value.
// Special Notes  :
// Scope          :
// Creator        : Dave Shirley
// Creation Date  : 03/10/2006
//----------------------------------------------------------------------------
bool Param::isInteger() const
{
  // Only do the checking if the parameter is string valued, return true
  // if it is already real or integer valued.
  if ( data_->enumType() == INT || data_->enumType() == LNG )
    return true;
  else if ( data_->enumType() == EXPR || data_->enumType() == BOOL)
    return false;
  else if ( data_->enumType() == DBLE)
    return true;
  else if ( data_->enumType() == STR)
    return isInt(getValue<std::string>());
  else
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, std::string("Param::isInteger: unknown type"));
  }
  return true;
}

//----------------------------------------------------------------------------
// Function       : Param::isBool
// Purpose        : Checks the value of the parameter to see if it is a
//                  legal bool value.
// Special Notes  :
// Scope          :
// Creator        : Dave Shirley
// Creation Date  : 03/10/2006
//----------------------------------------------------------------------------
bool Param::isBool() const
{
  // Only do the checking if the parameter is string valued, return true
  // if it is already real or integer valued.
  if ( data_->enumType() == DBLE || data_->enumType() == INT || data_->enumType() == LNG || data_->enumType() == BOOL)
    return true;
  if ( data_->enumType() == EXPR)
    return false;
  else if ( data_->enumType() == STR)
    return Util::isBool(getValue<std::string>());
  else
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,std::string("Param::isBool: unknown type"));
  }
  return true;
}

//----------------------------------------------------------------------------
// Function       : Param::setTimeDependent
// Purpose        : Set the flag indicating this is a time dependent parameter
//                  to the given value.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/19/2002
//----------------------------------------------------------------------------
void Param::setTimeDependent( bool timeDependent )
{
  std::string exp;

  if (data_->enumType() == EXPR && timeDependent)
    return;
  if (data_->enumType() != EXPR && !timeDependent)
    return;
  if (data_->enumType() != STR)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL,std::string("Param::setTimeDependent: attempt to convert non-string to expression"));
  }
  if (!timeDependent)
    return;
  exp = getValue<std::string>();
  setVal(Expression(exp));
}

//----------------------------------------------------------------------------
// Function       : Param::isTimeDependent
// Purpose        : Return the value of the flag indicating whether the
//                  parameter is time dependent.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/19/2002
//----------------------------------------------------------------------------
bool Param::isTimeDependent() const
{
  return data_->enumType() == EXPR;
}


//-----------------------------------------------------------------------------
// Function      : Param::instance()
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 6/28/01
//-----------------------------------------------------------------------------
Packable * Param::instance() const
{
  return new Param();
}

//-----------------------------------------------------------------------------
// Function      : Param::packedByteCount
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
int Param::packedByteCount() const
{

  int byteCount = 0;

  //tag info
  byteCount += tag_.length() + sizeof(int);

  //type
  byteCount += sizeof(int);

  //value info
  switch( data_->enumType() )
  {
    case -1:
      break;
    case STR:
      byteCount += getValue<std::string>().length() + sizeof(int);
      break;
    case DBLE:
      byteCount += sizeof(double);
      break;
    case BOOL:
    case INT:
      byteCount += sizeof(int);
      break;
    case LNG:
      byteCount += sizeof(long);
      break;
    case EXPR:
      byteCount += getValue<Expression>().get_expression().length() + sizeof(int);
      break;
    case STR_VEC:
      byteCount += sizeof(int);
      for (size_t i = 0; i < getValue<std::vector<std::string> >().size(); i++) {
        byteCount += getValue<std::vector<std::string> >()[i].length() + sizeof(int);
      }
      break;
    case DBLE_VEC:
      byteCount += ( sizeof(int) + getValue<std::vector<double> >().size() * sizeof(double) );
      break;
  }

  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : Param::pack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
void Param::pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const
{

  int length;
  std::string tmp;

  //pack tag
  length = tag_.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( tag_.c_str(), length, buf, bsize, pos );

  //pack type
  int enum_type = data_->enumType();

  comm->pack( &enum_type, 1, buf, bsize, pos );

  //pack value
  switch( data_->enumType() )
  {
    case -1:
      break;

    case STR:
      length = getValue<std::string>().length();
      comm->pack( &length, 1, buf, bsize, pos );
      comm->pack( getValue<std::string>().c_str(), length, buf, bsize, pos );
      break;

    case DBLE:
      comm->pack( &(getValue<double>()), 1, buf, bsize, pos );
      break;

    case INT:
      comm->pack( &(getValue<int>()), 1, buf, bsize, pos );
      break;

    case BOOL:
    {
      int i;
      if (getValue<bool>())
        i = 1;
      else
        i = 0;
      comm->pack( &i, 1, buf, bsize, pos );
    }
      break;

    case LNG:
      comm->pack( &(getValue<long>()), 1, buf, bsize, pos );
      break;

    case EXPR:
      tmp = getValue<Expression>().get_expression();
      length = tmp.length();
      comm->pack( &length, 1, buf, bsize, pos );
      comm->pack( tmp.c_str(), length, buf, bsize, pos );
      break;

    case STR_VEC:
    {
      const std::vector<std::string> &string_vector = getValue<std::vector<std::string> >();
      length = (int) string_vector.size();
      comm->pack( &length, 1, buf, bsize, pos );
      for (int i=0; i < (int) string_vector.size(); i++)
      {
        length = string_vector[i].length();
        comm->pack( &length, 1, buf, bsize, pos );
        comm->pack( string_vector[i].c_str(), length, buf, bsize, pos );
      }
    }

      break;

    case DBLE_VEC:
    {
      const std::vector<double> &double_vector = getValue<std::vector<double> >();
      length = (int) double_vector.size();
      comm->pack( &length, 1, buf, bsize, pos );
      comm->pack( &double_vector[0], length, buf, bsize, pos );
    }

      break;

    default:   Report::DevelFatal() << "Param::pack: unknown type " << data_->enumType();
  }
}

//-----------------------------------------------------------------------------
// Function      : Param::unpack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
void Param::unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm )
{

  int length;
  int vector_size;

  //unpack tag
  comm->unpack( pB, bsize, pos, &length, 1 );

  tag_ = std::string( (pB+pos), length );
  pos += length;

  //unpack type
  int enum_type = -1;
  comm->unpack( pB, bsize, pos, &enum_type, 1 );

  switch (enum_type)
  {
    case -1:
      break;

    case STR:
      comm->unpack( pB, bsize, pos, &length, 1 );
      setVal(std::string( (pB+pos), length ));
      pos += length;
      break;

    case DBLE: 
    {
      double d;
      
      comm->unpack( pB, bsize, pos, &d, 1 );
      setVal(d);
    }
    break;

    case INT: 
    {
      int i;
      comm->unpack( pB, bsize, pos, &i, 1 );
      setVal(i);
    }
    break;

    case BOOL:
      int i;
      comm->unpack( pB, bsize, pos, &i, 1 );
      if (i == 0)
        setVal(false);
      else
        setVal(true);
      break;

    case LNG:
    {
      long l;
      comm->unpack( pB, bsize, pos, &l, 1 );
      setVal(l);
    }
    break;

    case EXPR:
      comm->unpack( pB, bsize, pos, &length, 1 );
      setVal(Expression(std::string( (pB+pos), length )));
      pos += length;
      break;

    case STR_VEC:
    {

      comm->unpack( pB, bsize, pos, &vector_size, 1 );
      setVal(std::vector<std::string>());
      std::vector<std::string> &x = getValue<std::vector<std::string> >();
      x.reserve(vector_size);

      for (int i=0; i< vector_size; i++)
      {
        comm->unpack( pB, bsize, pos, &length, 1 );
        x.push_back(std::string( (pB+pos), length ));
        pos += length;
      }
    }

    break;

    case DBLE_VEC:
    {
      comm->unpack( pB, bsize, pos, &vector_size, 1 );
      setVal(std::vector<double>());
      std::vector<double> &x = getValue<std::vector<double> >();
      x.resize(vector_size, 0.0);
      comm->unpack( pB, bsize, pos, &(x[0]), vector_size );
    }

    break;

    default:
      std::string msg = "Param::unpack: unknown type";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg );
  }
}


//-----------------------------------------------------------------------------
// Function      : Param::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
std::ostream &operator<<(std::ostream & os, const Param & p)
{
  os << p.tag() << "\t";

  switch (p.getType())
  {
    case STR:
      os << " STR\t" << p.stringValue();
      break;
    case DBLE:
      os << "DBLE\t"<< p.getValue<double>();
      break;
    case INT:
      os << " INT\t" << p.getValue<int>();
      break;
    case LNG:
      os << " LNG\t" << p.getValue<long>();
      break;
    case BOOL:
      os << "BOOL\t" << p.getValue<bool>();
      break;
    case EXPR:
      os << "EXPR\t" << p.stringValue();
      break;

    case STR_VEC: 
    {
      os << "STR_VEC\t";
      int numElements = p.getValue<std::vector<std::string> >().size();
      for(int i=0; i<numElements; i++)
      {
        os << p.getValue<std::vector<std::string> >()[i] << " ";
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

  os << std::endl;

  return os;
}

} // namespace Util
} // namespace Xyce
