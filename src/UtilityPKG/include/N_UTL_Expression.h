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
// Filename       : $RCSfile: N_UTL_Expression.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 04/17/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $$
//
// Revision Date  : $$
//
// Current Owner  : $$
//-----------------------------------------------------------------------------

#ifndef N_UTL_Expression_H
#define N_UTL_Expression_H

// ---------- Standard Includes ----------
#include <vector>
#include <string>
#include <list>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

#include <N_UTL_Misc.h>
#include <N_UTL_Interface_Enum_Types.h>
#include <iosfwd>


// ---------- Forward Declarations ----------
class N_UTL_ExpressionInternals;

//-----------------------------------------------------------------------------
// Class         : N_UTL_Expression
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
class N_UTL_Expression
{

public:

  N_UTL_Expression (string const & exp = string());
  N_UTL_Expression (const N_UTL_Expression &);
  //  N_UTL_Expression& operator=(const N_UTL_Expression& right) ;
  ~N_UTL_Expression (void);

  bool set (string const & exp);
  void get_names (int const & type, vector < string > & names) const;
  int get_type (string const & var);
  char get_lead_designator (string const & var);
  bool make_constant (string const & var, double const & val);
  bool make_var (string const & var);

  int differentiate();

  bool set_var (const string &, const double &);
  bool set_vars (const vector < double > &);

  string get_expression (void);
  string get_derivative(string const & var);
  int get_num(int const & type);

  int evaluate (double &result, vector < double > &derivs, vector < double > &vals);
  int evaluateFunction (double &result, vector < double > &vals);

  int evaluate (double &result, vector < double > &derivs);
  int evaluateFunction (double &result);

  bool set_sim_time (double const & time);
  bool set_temp (double const & temp);
  void set_accepted_time ();
  double get_break_time (void);
  double get_break_time_i (void);
  const string & get_input (void);
  int order_names (vector < string > const & new_names);
  int replace_func (string const & func_name, N_UTL_Expression & func_def, int numArgs);
  bool replace_name (const string & old_name, const string & new_name);
  int replace_var (string const & var_name, N_UTL_Expression & subexpr);
  int getNumDdt();
  void getDdtVals (vector<double> &);
  void setDdtDerivs (vector<double> &);
  int num_vars ();
#ifdef Xyce_DEBUG_EXPRESSION
  void dumpParseTree();
#endif

private:

  RefCountPtr<N_UTL_ExpressionInternals> expPtr_;

};
#endif // N_UTL_EXPRESSION_H

