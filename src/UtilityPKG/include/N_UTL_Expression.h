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

namespace Xyce {
namespace Util {

class ExpressionInternals;

//-----------------------------------------------------------------------------
// Class         : Expression
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
class Expression
{

public:

  Expression (std::string const & exp = std::string());
  Expression (const Expression &);
  //  Expression& operator=(const Expression& right) ;
  ~Expression (void);

  bool parsed() const;
    
  bool set (std::string const & exp);
  void get_names (int const & type, std::vector< std::string > & names) const;
  int get_type (std::string const & var);
  char get_lead_designator (std::string const & var);
  bool make_constant (std::string const & var, double const & val);
  bool make_var (std::string const & var);

  int differentiate();

  bool set_var (const std::string &, const double &);
  bool set_vars (const std::vector< double > &);

  std::string get_expression (void) const;
  std::string get_derivative(std::string const & var);
  int get_num(int const & type);

  int evaluate (double &result, std::vector< double > &derivs, std::vector< double > &vals);
  int evaluateFunction (double &result, std::vector< double > &vals);

  int evaluate (double &result, std::vector< double > &derivs);
  int evaluateFunction (double &result);

  bool set_sim_time (double const & time);
  bool set_temp (double const & temp);
  void set_accepted_time ();
  double get_break_time (void);
  double get_break_time_i (void);
  const std::string & get_input (void);
  int order_names (std::vector< std::string > const & new_names);
  int replace_func (std::string const & func_name, Expression & func_def, int numArgs);
  bool replace_name (const std::string & old_name, const std::string & new_name);
  int replace_var (std::string const & var_name, Expression & subexpr);
  int getNumDdt();
  void getDdtVals (std::vector<double> &);
  void setDdtDerivs (std::vector<double> &);
  int num_vars() const;
  bool isTimeDependent() const;
  void dumpParseTree();

private:

  RefCountPtr<ExpressionInternals> expPtr_;

};

} // namespace Util
} // namespace Xyce

typedef Xyce::Util::Expression N_UTL_Expression;

#endif // N_UTL_EXPRESSION_H

