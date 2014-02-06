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
// Filename       : $RCSfile: N_UTL_Expression.C,v $
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Eric R. Keiter, SNL
//
// Creation Date : 04/17/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.111.4.4 $
//
// Revision Date  : $Date: 2013/12/03 23:30:12 $
//
// Current Owner  : $Author: rlschie $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


#include <N_UTL_Misc.h>
#include <N_DEV_Const.h>

// ---------- Standard Includes ----------

#include <iterator>
#include <iostream>
#include <iomanip>

#include <sstream> 
// ----------   Xyce Includes   ----------
#include <N_UTL_Expression.h>
#include <N_UTL_ExpressionInternals.h>

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::N_UTL_Expression
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------

N_UTL_Expression::N_UTL_Expression( const string & exp )
{
#ifdef _OMP
  #pragma omp critical
#endif
  {
    expPtr_ = Teuchos::rcp(new N_UTL_ExpressionInternals(exp));
  }
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::N_UTL_Expression
// Purpose       : Copy Constructor
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
N_UTL_Expression::N_UTL_Expression( const N_UTL_Expression & right)
{
#ifdef _OMP
  #pragma omp critical
#endif
  {
    // important distinction here.  If we're calling the copy 
    // constructor, we don't need to copy the underlying ExpressionInternals 
    // object.  We would if people wanted multiple versions of an expression 
    // for factoring or parsing.  
    // Deep Copy Not Needed
    // expPtr_ = Teuchos::rcp(new N_UTL_ExpressionInternals( *(right.expPtr_)) );
    // Shallow Copy
    expPtr_ = right.expPtr_;
  }
  return;
}


//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::~N_UTL_Expression
// Purpose       : Destructor
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
N_UTL_Expression::~N_UTL_Expression ()
{
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::set
// Purpose       : Set the value of the expression to a string
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool N_UTL_Expression::set ( const string & exp )
{
  bool retVal = false;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->set (exp);
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::get_names
// Purpose       : Returns the names of input quantities by type
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
void N_UTL_Expression::get_names(int const & type, vector<string> & names ) const
{
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    expPtr_->get_names(type,names);
  }
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::get_type
// Purpose       : Finds the type of an input quantity name
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int N_UTL_Expression::get_type ( const string & var )
{
  int retVal=0;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->get_type (var);
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::get_lead_designator
// Purpose       : Returns the lead designator (" ", "B", "S", "E", "D", "")
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, SNL
// Creation Date : 01/15/13
//-----------------------------------------------------------------------------
char N_UTL_Expression::get_lead_designator ( const string & var )
{
  int retVal=0;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->get_lead_designator (var);
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::make_constant
// Purpose       : Convert a 'string' placeholder into a constant
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool N_UTL_Expression::make_constant (const string & var, const double & val)
{
  bool retVal=false;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->make_constant (var,val);
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::make_var
// Purpose       : Convert a 'string' placeholder into a variable
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool N_UTL_Expression::make_var (string const & var)
{
  bool retVal=false;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->make_var(var);
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::differentiate
// Purpose       : Form the analytic derivative trees for all variables
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int N_UTL_Expression::differentiate ()
{
  int retVal=0;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->differentiate ();
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::set_var
// Purpose       : Sets the value of an input quantity
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool N_UTL_Expression::set_var ( const string & var,
                                 const double & val)
{
  bool retVal=false;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->set_var (var, val);
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::set_vars
// Purpose       : Sets the values of all input quantities
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool N_UTL_Expression::set_vars ( const vector<double> & vals )
{
  bool retVal=false;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->set_vars ( vals );
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::get_expression
// Purpose       : Returns a string of the expression
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
string N_UTL_Expression::get_expression ()
{
  string retVal;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->get_expression ();
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::get_derivative
// Purpose       : Returns a string of a derivative
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
string N_UTL_Expression::get_derivative ( string const & var )
{
  string retVal;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->get_derivative ( var );
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::get_num
// Purpose       : Returns the number of input quantities of a requested type
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int N_UTL_Expression::get_num(int const & type)
{
  int retVal=0;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->get_num(type);
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::evaluate
// Purpose       : Evaluate expression and derivatives using provided input values
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int N_UTL_Expression::evaluate ( double & exp_r,
                                 vector<double> & deriv_r,
                                 vector<double> & vals )
{
  int retVal=0;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->evaluate ( exp_r, deriv_r, vals );
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::evaluateFunction
// Purpose       : Evaluate expression using provided input values.  
// Special Notes : This is for cases in which the user does not need 
//                 the derivatives.
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date : 04/14/08
//-----------------------------------------------------------------------------
int N_UTL_Expression::evaluateFunction ( double & exp_r, vector<double> & vals )
{
  int retVal=0;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->evaluateFunction ( exp_r, vals );
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::evaluate
// Purpose       : Evaluate expression and derivatives using stored input values
// Special Notes :
// Scope         : private
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int N_UTL_Expression::evaluate ( double & exp_r,
                                 vector<double> & deriv_r)
{
  int retVal=0;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->evaluate ( exp_r, deriv_r);
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::evaluateFunction
// Purpose       : Evaluate expression using stored input values
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int N_UTL_Expression::evaluateFunction ( double & exp_r )
{
  int retVal=0;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->evaluateFunction ( exp_r );
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::set_sim_time
// Purpose       : Set 'time' special variable in expression
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool N_UTL_Expression::set_sim_time(double const & time)
{
  bool retVal=false;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->set_sim_time(time);
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::set_temp
// Purpose       : Set 'temp' special variable in expression
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool N_UTL_Expression::set_temp(double const & tempIn)
{
  bool retVal=false;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->set_temp(tempIn);
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::set_accepted_time
// Purpose       : Set accepted time for converged soltion
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
void N_UTL_Expression::set_accepted_time()
{
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    expPtr_->set_accepted_time();
  }
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::get_break_time
// Purpose       : Returns next breakpoint time
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
double N_UTL_Expression::get_break_time()
{
  double retVal=0.0;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->get_break_time();
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::get_break_time_i
// Purpose       : Returns next breakpoint time
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
double N_UTL_Expression::get_break_time_i()
{
  double retVal=0.0;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->get_break_time_i();
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::get_input
// Purpose       : Return expression input string
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
const string & N_UTL_Expression::get_input ()
{
  return expPtr_->get_input ();
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::order_names
// Purpose       : Put input quantity names in a particular order (used for
//                 replace_func which requires identical ordering for expression
//                 and user defined function
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int N_UTL_Expression::order_names(vector<string> const & new_names)
{
  int retVal=0;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->order_names(new_names);
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::replace_func
// Purpose       : Replace user defined function with its definition in expression
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int N_UTL_Expression::replace_func (string const & func_name,
                                   N_UTL_Expression & func_def,
                                   int numArgs)
{
  int retVal=0;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->replace_func (func_name, *(func_def.expPtr_), numArgs);
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::replace_var
// Purpose       : Replace a variable usage with a parsed sub-expression
// Special Notes : This is used for subcircuit parameters that cannot be
//                 fully resolved to a constant because they have global
//                 parameter usage.
// Scope         :
// Creator       : Thomas Russo, SNL
// Creation Date : 08/10/2010
//-----------------------------------------------------------------------------
int N_UTL_Expression::replace_var (string const & var_name,
                                   N_UTL_Expression & subexpr )
{
  int retVal=0;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->replace_var (var_name, *(subexpr.expPtr_));
  }
  return retVal;
}


//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::replace_name
// Purpose       : Change the name of an input quantity
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool N_UTL_Expression::replace_name ( const string & old_name,
                                      const string & new_name)
{
  bool retVal=false;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->replace_name ( old_name, new_name);
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::getNumDdt
// Purpose       : Return the number of ddt() calls in the expression
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int N_UTL_Expression::getNumDdt ()
{
  int retVal=0;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_-> getNumDdt ();
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::getDdtVals
// Purpose       : Return the most recent arguments of ddt() in the expression
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
void N_UTL_Expression::getDdtVals ( vector<double> & vals )
{
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    expPtr_-> getDdtVals ( vals );
  }
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::setDdtDerivs
// Purpose       : Set the evaluated value of the ddt functions
// Special Notes : This is normally done with derivative values from the
//                 time integration package
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
void N_UTL_Expression::setDdtDerivs ( vector<double> & vals )
{
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    expPtr_->setDdtDerivs ( vals );
  }
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::num_vars
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int N_UTL_Expression::num_vars ()
{
  int retVal=0;
#ifdef _OMP
  #pragma omp critical
#endif 
  {
    retVal = expPtr_->num_vars();
  }
  return retVal;
}


//-----------------------------------------------------------------------------
// Function      : Expression::isTimeDependent
// Purpose       : Return true if expression is either explicitly or implicitly
//                 time dependent
// Special Notes : The ExpressionInternals::isTimeDependent method only returns
//                 true if the expression is implicitly time dependent.
//
//                 It is impossible at this time for indirect time dependence
//                 through global_params to be detected through this method.
//
// Scope         :
// Creator       : Richard Schiek, SNL
// Creation Date : 10/07/2013
//-----------------------------------------------------------------------------
bool N_UTL_Expression::isTimeDependent() const
{
  bool implicitTimeDep = expPtr_->isTimeDepedent();
  bool explicitTimeDep = false;
  std::vector<std::string> specials;
  expPtr_->get_names(XEXP_SPECIAL, specials);
  if (!specials.empty())
  {
    explicitTimeDep=(std::find(specials.begin(), specials.end(), "TIME") != specials.end());
  }
  return (implicitTimeDep || explicitTimeDep);
}

#ifdef Xyce_DEBUG_EXPRESSION
//-----------------------------------------------------------------------------
// Function      : N_UTL_Expression::dumpParseTree
// Purpose       : Dump out the parse tree for an expression
// Special Notes : 
// Scope         :
// Creator       : Tom Russo, SNL
// Creation Date : 9/9/10
//-----------------------------------------------------------------------------
void N_UTL_Expression::dumpParseTree()
{
  expPtr_->dumpParseTree();
}
#endif
