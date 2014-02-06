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
// Filename       : $RCSfile: N_UTL_ExpressionInternals.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Shirley, PSSI
//
// Creation Date  : 06/07/01
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

#ifndef N_UTL_ExpressionInternals_H
#define N_UTL_ExpressionInternals_H

// ---------- Standard Includes ----------
#include <vector>
#include <string>
#include <list>

#include <N_UTL_Misc.h>
#include <N_UTL_Interface_Enum_Types.h>
#include <iosfwd>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Other Includes   ----------

// ---------- Structure definitions ----------

class ExpressionNode
{
  public:
   ExpressionNode () :
    valueIndex(0), type(0), constant(0), 
    funcname(""), eval_value(0), eval_num(0),
    funcnum(0)
   {
    operands.clear();
    state.clear();
    fptr.unary = NULL; 
    fptr.binary = NULL;
   };

   ~ExpressionNode () {};

  public:
    int type;                                   // One of EXPR_*
    vector <ExpressionNode *> operands;         // list of pointers to operands
    double constant;                            // If EXPR_CONSTANT, or if a 
                                                // table entry then
                                                // use this as a cache for the 
                                                // entry (which is
                                                // required to be constant)
    vector <double> state;                      // state data for operation 
                                                // (e.g. sdt or ddt)
    int valueIndex;                             // If EXPR_VAR, EXPR_ARGUMENT, 
                                                // or (EXPR_FUNCTION and 
                                                // EXPR_F_USER), use as
                                                // index into vars.  If 
                                                // EXPR_FUNCTION and
                                                // funcnum == EXPR_F_TABLE then
                                                // this is used as a cache to 
                                                // record the section that the
                                                // table last evaluated in
    string funcname;                            // If EXPR_FUNCTION, name of 
                                                // function,
                                                // If EXPR_VAR, name of var
    double eval_value;                          // Evaluated value at current 
                                                // eval_num
    int eval_num;                               // Evaluation number
    int funcnum;                                // ... one of EXPR_F_*
    union                                       // ... and pointer to the function
    {
      double (*unary)(double);
      double (*binary)(double, double);
    } fptr;
} ;

class ExpressionElement
{
  public:
    ExpressionElement () :
      token(0),
      type(0),
      name(""),
      number(0),
      node(NULL)
    {};

    ~ExpressionElement () {};

  public:
    int token;                 // see enum TOK_LIST
    int type;                  // see enum TYP_NODES 
    string name;
    double number;
    ExpressionNode *node;
}; 

class table_state 
{
  public:
    ExpressionElement table;   // pointer to table INPparseNode
    double **tab1;             // address of pointer to start of independent 
                               // table values
    double **tab2;             // address of pointer to start of dependent 
                               // table values
    int last_segment;          // segment used for last accepted point */
    int neighbor;              // distance to neighbor segment (= 2 for 
                               // derivatives, 1 for expressions)
    table_state *next;
};

// ---------- Enumerations ----------

enum EXPR_OPS
{
  EXPR_PLACEHOLDER,      //  0       Used during parsing.  This could be a 
                         //          device or node name
                         //          as in V() or I() or it could just be a 
                         //          string encountered in an expression which 
                         //          will later be turned into a
                         //          variable or constant.
  EXPR_OR,               //  1       Binary functions are 1 - 15
  EXPR_XOR,              //  2
  EXPR_AND,              //  3
  EXPR_EQUAL,            //  4
  EXPR_NOTEQ,            //  5
  EXPR_GREAT,            //  6
  EXPR_GREATEQ,          //  7
  EXPR_LESS,             //  8
  EXPR_LESSEQ,           //  9
  EXPR_PLUS,             // 10
  EXPR_MINUS,            // 11
  EXPR_TIMES,            // 12
  EXPR_DIVIDE,           // 13
  EXPR_REMAINDER,        // 14
  EXPR_POWER,            // 15
  EXPR_FUNCTION,         // 16       Function, type according to EXPR_F_FUNCS
  EXPR_CONSTANT,         // 17       Constant
  EXPR_VAR               // 18       Variable, type according to EXPR_T_TYPES 
};

enum EXPR_F_FUNCS                 // These differentiate between different types
                                  // of EXPR_FUNCTIONs
{
  EXPR_F_ABS,             //  0
  EXPR_F_ACOS,            //  1
  EXPR_F_ACOSH,           //  2
  EXPR_F_ASIN,            //  3
  EXPR_F_ASINH,           //  4
  EXPR_F_ATAN,            //  5
  EXPR_F_ATANH,           //  6
  EXPR_F_COS,             //  7
  EXPR_F_COSH,            //  8
  EXPR_F_DDT,             //  9
  EXPR_F_DDX,             // 10
  EXPR_F_EXP,             // 11
  EXPR_F_IF,              // 12
  EXPR_F_LN,              // 13
  EXPR_F_LOG,             // 14
  EXPR_F_NOT,             // 15
  EXPR_F_RAND,            // 16
  EXPR_F_SDT,             // 17
  EXPR_F_SGN,             // 18
  EXPR_F_SIN,             // 19
  EXPR_F_SINH,            // 20
  EXPR_F_SQRT,            // 21
  EXPR_F_TABLE,           // 22
  EXPR_F_F_TABLE,         // 23
  EXPR_F_R_TABLE,         // 24
  EXPR_F_TAN,             // 25
  EXPR_F_TANH,            // 26
  EXPR_F_UMINUS,          // 27
  EXPR_F_URAMP,           // 28
  EXPR_F_USER,            // 29
  EXPR_F_AGAUSS,          // 30
  EXPR_F_INT,             // 31           
  EXPR_F_SCHEDULE         // 32           
};

enum EXPR_T_TYPES                 // These differentiate between different types
                                  // of EXPR_VARs
{
  EXPR_T_NODE=10,         // 10
  EXPR_T_STRING,          // 11
  EXPR_T_INSTANCE,        // 12
  EXPR_T_SPECIAL,         // 13
  EXPR_T_VARIABLE,        // 14
  EXPR_T_FUNCTION         // 15
};

enum TOK_LIST
{
  TOK_END,             //  0
  TOK_OR,              //  1
  TOK_XOR,             //  2
  TOK_AND,             //  3
  TOK_EQUAL,           //  4
  TOK_NOTEQ,           //  5
  TOK_GREAT,           //  6
  TOK_GREATEQ,         //  7
  TOK_LESS,            //  8
  TOK_LESSEQ,          //  9
  TOK_PLUS,            // 10
  TOK_MINUS,           // 11
  TOK_TIMES,           // 12
  TOK_DIVIDE,          // 13
  TOK_REMAINDER,       // 14
  TOK_POWER,           // 15
  TOK_NOT,             // 16
  TOK_UMINUS,          // 17
  TOK_LPAREN,          // 18
  TOK_RPAREN,          // 19
  TOK_VALUE,           // 20
  TOK_COMMA,           // 21
  TOK_LBRACE,          // 22
  TOK_RBRACE,          // 23
  TOK_SPACE            // 24
};

enum TYP_NODES        // These are for expression element types
{
  TYP_NUM,             //  0
  TYP_STRING,          //  1
  TYP_PNODE            //  2
};

enum EXPR_WARNINGS   // These indicate improper operations which are likely 
                     // inconsequential
{
  EXPR_WARNING = 100,    // 100
  EXPR_OR_WARNING,       // 101
  EXPR_XOR_WARNING,      // 102
  EXPR_AND_WARNING,      // 103
  EXPR_SUM_WARNING,      // 104
  EXPR_DIFF_WARNING,     // 105
  EXPR_TIMES_WARNING,    // 106
  EXPR_DIVIDE_WARNING,   // 107
  EXPR_POWER_WARNING,    // 108
  EXPR_NOT_WARNING       // 109
};

enum EXPR_ERRORS   // These indicate possible loss of precision
{
  EXPR_ERROR = 200,      // 200
  EXPR_REMAINDER_ERROR,  // 201
  EXPR_POWER_ERROR,      // 202
  EXPR_COSH_ERROR,       // 203
  EXPR_EXP_ERROR,        // 204
  EXPR_SINH_ERROR,       // 205
  EXPR_SQRT_ERROR       // 206
};

enum EXPR_FATALS   // These indicate ill-defined operations, and wrong answers
{
  EXPR_FATAL = 300,      // 300
  EXPR_DIVIDE_FATAL,     // 301
  EXPR_REMAINDER_FATAL,  // 302
  EXPR_LOG_FATAL,        // 303
  EXPR_TABLE_FATAL,      // 304
  EXPR_NODERIV_FATAL,    // 305
  EXPR_ARGUMENT_FATAL,   // 306
  EXPR_FUNCTION_FATAL    // 307
};

// Some miscellaneous definitions:

#define EXPR_HUGE (1.e+50)
#define N_INT_STATE     5

double EXPRor(double arg1, double arg2);
double EXPRxor(double arg1, double arg2);
double EXPRand(double arg1, double arg2);
double EXPRequal(double arg1, double arg2);
double EXPRnoteq(double arg1, double arg2);
double EXPRgreat(double arg1, double arg2);
double EXPRgreateq(double arg1, double arg2);
double EXPRless(double arg1, double arg2);
double EXPRlesseq(double arg1, double arg2);
double EXPRplus(double arg1, double arg2);
double EXPRminus(double arg1, double arg2);
double EXPRtimes(double arg1, double arg2);
double EXPRdivide(double arg1, double arg2);
double EXPRremainder(double arg1, double arg2);
double EXPRpower(double arg1, double arg2);

double EXPRabs (double arg);
double EXPRacos (double arg);
double EXPRacosh (double arg);
double EXPRasin (double arg);
double EXPRasinh (double arg);
double EXPRatan (double arg);
double EXPRatanh (double arg);
double EXPRcos (double arg);
double EXPRcosh (double arg);
double EXPRexp (double arg);
double EXPRln (double arg);
double EXPRlog (double arg);
double EXPRnot (double arg);
double EXPRsgn (double arg);
double EXPRsin (double arg);
double EXPRsinh (double arg);
double EXPRsqrt (double arg);
double EXPRtan (double arg);
double EXPRtanh (double arg);
double EXPRuminus (double arg);
double EXPRuramp (double arg);
double EXPRint (double arg);

//-----------------------------------------------------------------------------
// Class         : N_UTL_ExpressionInternals
// Purpose       :
// Special Notes :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/20/08
//-----------------------------------------------------------------------------
class N_UTL_ExpressionInternals
{

public:

  N_UTL_ExpressionInternals (string const & exp = string());
  N_UTL_ExpressionInternals (const N_UTL_ExpressionInternals & right);
  //  N_UTL_ExpressionInternals& operator=(const N_UTL_ExpressionInternals& right) ;
  ~N_UTL_ExpressionInternals (void);

  bool set (string const & exp);
  void get_names (int const & type, vector < string > & names);
  int get_type (string const & var);
  char get_lead_designator ( const string & var );
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
  double get_break_time (void) {if(time_index == -2) return 0; else return(get_break_time_i());};
  double get_break_time_i (void);
  const string & get_input (void);
  int order_names (vector < string > const & new_names);
  int replace_func (string const & func_name, 
                    N_UTL_ExpressionInternals & func_def, int numArgs);
  int replace_var (string const & var_name, 
                   N_UTL_ExpressionInternals & subexpr);
  bool replace_name (const string & old_name, const string & new_name);
  int getNumDdt();
  void getDdtVals (vector<double> &);
  void setDdtDerivs (vector<double> &);
  inline ExpressionNode *get_tree() {return tree_;};
  inline int num_vars () {return num_N_+num_I_+num_lead_+num_string_+num_special_+num_var_+num_func_;};
  bool isTimeDepedent() const {return timeDependent_;};
#ifdef Xyce_DEBUG_EXPRESSION
  void dumpParseTree();
#endif

private:

  string Input_;
  bool differentiated_;
  bool ddxProcessed_;

  int num_N_;
  int num_I_;
  int num_lead_;
  int num_string_;
  int num_special_;
  int num_var_;
  int num_func_;
  double sim_time_;
  int time_index;
  bool timeDependent_;
  bool breakpointed_;

  int numVars_;
  vector <int> varTypes_;                     // array of types of variables
  vector <string> varValues_;                 // array of values of variables
  string leadDesignator_;                     // lead designator for current variables
  vector <double> var_vals_;                  // Values of variables, nodes, instances

  ExpressionNode *tree_;                      // The real stuff
  vector <ExpressionNode *> derivs_;          // The derivative parse trees

  vector <ExpressionNode *> free_list_;       // List of ExpressionNodes to be freed
                                              // in destructor
  vector <ExpressionNode *> done_list_;
  vector <ExpressionElement *> ee_list_;      // List of ExpressionElements to be freed
                                              // at end of constructor
  vector <ExpressionNode *> breaks_;
  //vector <bool> refine_break_;
  ExpressionNode *PThead_;
  int Rmode_;
  int ind_replace_;
  string Rstring_;
  double Rcval_;
  int curr_magic_;
  int curr_num_;
  bool values_changed_;

  // functions
  int find_num_ (const string &);

  void set_nums_ ();
  void create_vars_ ();
  void copy_elements_ (list <ExpressionElement *> &to, list <ExpressionElement *> *from);
  void copy_element_ (ExpressionElement *to, ExpressionElement *from);
  ExpressionNode *makepnode_ (ExpressionElement *elem);
  ExpressionNode *mkfnode_ (const string & fname, int num_args, vector <ExpressionNode *> args);
  ExpressionNode *mkfnode_ (const string & fname, int num_args, ExpressionNode *n);
  ExpressionNode *mksnode_ (const string & name);
  ExpressionNode *mkcon_ (double value);
  ExpressionNode *mkb_ (int type, ExpressionNode *left, ExpressionNode *right);
  ExpressionNode *mkf_(int type, ExpressionNode *arg);
  ExpressionNode *newExpressionNode_ ();
  void deleteExpressionNode_ (ExpressionNode *p);
  ExpressionElement *newExpressionElement_ ();
  void RpTree_ (ExpressionNode * pt, ostringstream & s);
  string varStr_ (int i);
  ExpressionNode * PTcheck_(ExpressionNode *p);
  ExpressionNode * diffDDX_(ExpressionNode *p);
  ExpressionNode * PTdiffDDX_(ExpressionNode *p);
  ExpressionNode * com_expr_ (ExpressionNode *c, ExpressionNode *p);
  ExpressionNode * Differentiate_ (ExpressionNode *arg, int varnum);

  void Rconvert_ (ExpressionNode & node);
  void RcountDDT_ (ExpressionNode & node);
  void RgetDDT_ (ExpressionNode & node, vector<double> & vals);
  void RsetDDT_ (ExpressionNode & node, vector<double> & vals);
  void convert_to_constant_ (int i, double c_value);
  void convert_to_variable_ (int i);

  // method to support copy constructor:
  ExpressionNode * copy_exprNode_ (ExpressionNode *n);

  // Methods to support breakpoints:
  void simplify_ (ExpressionNode & node);
  void get_breaks_ (ExpressionNode & node);
  bool dependent_ (ExpressionNode & node, int ind);
  bool dependent_other_ (ExpressionNode & node, int ind);
  bool arithmatic_ (ExpressionNode & node);

  // Methods to suupport order_names:
  void Rmap_ (ExpressionNode & node, int mode, vector<int> &nmap);
  int EXPRaddDummyString_ (string & dummy);

  // Methods to support replace_func:
  void addNode_ (ExpressionNode *n, int ind, ExpressionNode *f, 
                 N_UTL_ExpressionInternals & func_expr,
                 int na_func,
                 vector<ExpressionNode *> operands);
  void Nreplace_ (ExpressionNode *n, ExpressionNode *f, 
                  N_UTL_ExpressionInternals & func_expr,
                  int na_func,
                  vector<ExpressionNode *> operands);
  int Freplace_ (ExpressionNode *n, string const & func_name,
                 N_UTL_ExpressionInternals & func_expr,
                 int na_func);
  void RemoveFentry_ (ExpressionNode *n, int new_ind, int old_ind);

  // method to support replace_var
  int Vreplace_ (ExpressionNode *n, string const & varName,
                 N_UTL_ExpressionInternals & subexpr);

  // Methods to support evaluate:
  void EXPReval_ (ExpressionNode & node, double & res, vector<double> &vals);
  void clear_eval_num_ (ExpressionNode *n);

  // Miscellaneous utility routines
  void compactLine_(string & inputLine, string & compactedLine);
  void tokenize_(string & inputLine, list<ExpressionElement *> & tokenList);
  void convertPolyToExpr_(list<ExpressionElement *> & tokenList);
  void standardizeTable_(list<ExpressionElement *> & tokenList);
#ifdef Xyce_DEBUG_EXPRESSION
  // Debugging methods

  void dumpParseTree_(ExpressionNode *tree, int indentLevel=0);
  void indentWithDashes_(int level);
#endif

  // Catch attempts to use operator= at compile time, by 
  // defining this null-op as a private member.
  // This doesn't work because parameters are often put inside of
  // STL objects, which rely on operator=.
  //N_UTL_ExpressionInternals& operator=(const N_UTL_ExpressionInternals& right) ;

};
#endif // N_UTL_EXPRESSION_H
