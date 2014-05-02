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
// Filename       : $RCSfile: N_UTL_ExpressionData.C,v $
//
// Purpose        : Handle data for one expression object
//
// Special Notes  :
//
// Creator        : Richard Schiek
//
// Creation Date  : 8/24/2009
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.28 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>
#include <N_UTL_Misc.h>

// ---------- Standard Includes ----------
#include <iostream>
#include <fstream>

#include <sstream>

// ---------- Xyce Includes -------------
#include <N_LAS_Vector.h>
#include <N_UTL_ExpressionData.h>
#include <N_UTL_Param.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_UTL_Expression.h>

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Function      : ExpressionData::ExpressionData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------
ExpressionData::ExpressionData (
    const std::string & expression,
    IO::OutputMgr &    output_manager
) :
  expPtr_(0),
  outputManager_(output_manager),
  expression_(expression),
  setupCalled(false),
  val(0.0),
  numUnresolvedStrings(0),
  numUnresolvedStringsChecked(false)
{}


//-----------------------------------------------------------------------------
// Function      : ExpressionData::ExpressionData
// Purpose       : constructor -- used to take an existing expression pointer
//                 not a string that can be converted to an expression
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------
ExpressionData::ExpressionData (
  const Expression &  expression,
  IO::OutputMgr &    output_manager)
  : expPtr_(new Expression( expression )),
    outputManager_(output_manager),
    expression_(expPtr_->get_expression()),
    setupCalled(false),
    val(0.0),
    numUnresolvedStrings(0),
    numUnresolvedStringsChecked(false)
{}


//-----------------------------------------------------------------------------
// Function      : ExpressionData::~ExpressionData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------
ExpressionData::~ExpressionData ()
{
  delete expPtr_;

  for (Util::OpList::iterator it = expressionVars_.begin(); it != expressionVars_.end(); ++it)
    delete *it;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionData::evaluate
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------
double ExpressionData::evaluate (const N_LAS_Vector * solnVecPtr, 
                                 const N_LAS_Vector * stateVecPtr, 
                                 const N_LAS_Vector * stoVecPtr,  
                                 const N_LAS_Vector * solnVecImagPtr)
{
  int i;

  if (!setupCalled)
  {
    numUnresolvedStrings = setup();
  }

  if( solnVecPtr != NULL )
  {
    // loop over expressionVars_ to get all the values.
    varVals.clear();
    for (OpList::const_iterator it = expressionVars_.begin(); it != expressionVars_.end(); ++it)
      varVals.push_back(IO::getValue(outputManager_.getCommPtr()->comm(), *(*it), solnVecPtr, solnVecImagPtr, stateVecPtr, stoVecPtr ).real());

    // STD and DDT are implicitly time dependent.  Check the underlying expression 
    // for time dependence, and if it is get the current time with getPrgetImmutableValue<int>()
    // and set it in the underlying expression 
    if( expPtr_->isTimeDependent() )
    {
      expPtr_->set_sim_time(outputManager_.getTime());
    }

    // now get expression value and partial derivatives.
    // Note that for these purposes (.PRINT output) we only
    // really need the value, not the derivs.
    val = 0.0;
    if (expPtr_)
    {
      expPtr_->evaluateFunction ( val, varVals );
    }
  }

  return val;
}

//-----------------------------------------------------------------------------
// Namespace     : Unnamed
// Purpose       : file-local scoped methods and data
// Special Notes : just the declaration, definition at end of file
// Creator       : Tom Russo, SNL
// Creation Date : 11/27/2013
//-----------------------------------------------------------------------------
namespace {
void convertNodalComputation(std::string &nodalComputation, 
                             std::list<Param> &paramList);
} // namespace (unnamed)

//-----------------------------------------------------------------------------
// Function      : ExpressionData::setup
//
// Purpose       : Manages all the basic setup for this class.
//
// Special Notes : Originally, all this stuff was in the constructor,
//                 but it needs to happen after topology is
//                 completely done setting itself up, and the
//                 constructor call was too early for that.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/10/04
//-----------------------------------------------------------------------------
int ExpressionData::setup()
{
#ifdef Xyce_DEBUG_EXPRESSION
  Xyce::dout() << "In ExpressionData::setup" << std::endl;
#endif

  setupCalled = true;

  // we return the number of unresolved strings on exit.  
  // we can't just exit with an error if there are unresolved strings 
  // because in parallel such strings could be resolved on another 
  // processor.  Thus the number of unresolved strings is saved in 
  // numUnresolvedStrings and the class owner can check that at least 
  // one processor in parallel has 0 for numUnresolvedStrings.
  int returnVal = 0;

  // allocate expression pointer if we need to
  if( expPtr_ == NULL )
  {
    expPtr_ = new Expression(expression_);
  }

  // clear out the STL containers:
  varNames.clear();
  varVals.clear();

  // query the expression object for all of its dependent vars.
  expPtr_->get_names(XEXP_ALL, varNames);

  // this varNames vec is a list of string representations of all of
  // the vars in the expression.  use this list to make a list of
  // Param objects which we can then call the
  // IO::OutputMgr::setParamContextType() and then use
  // getPrintValue_() to get the results.  This reuses code that is
  // already in place to handle solution, state, store and lead
  // currents.

  std::list<Param> param_list;
  
  std::vector<std::string>::iterator iterNode = varNames.begin();
  std::vector<std::string>::iterator iterEnd = varNames.end();
  while(  iterNode != iterEnd )
  {
    // based on the type of variable, create the needed Param
    // objects for setParmContextType_ to work.
    int varType = expPtr_->get_type( *iterNode );
    switch (varType)
    {
      case XEXP_NODAL_COMPUTATION:
        // deconstruct the string and turn it into params, push back into
        // param_list
        convertNodalComputation(*iterNode,param_list);
        break;
      case XEXP_NODE:
        param_list.push_back( Param( "V" , 1 ) );
        param_list.push_back( Param( *iterNode , 0.0 ) );
        break;
      case XEXP_INSTANCE:
      case XEXP_LEAD:
        {
          char leadDesignator = expPtr_->get_lead_designator( *iterNode);
          std::string currentName("I");
          if( (leadDesignator != 0) && (leadDesignator!=' ') )
          {
            currentName = currentName + leadDesignator;
          }
          param_list.push_back( Param( currentName , 1 ) );
          param_list.push_back( Param( *iterNode , 0.0 ) );
        }
        break;
      case XEXP_SPECIAL:
      case XEXP_STRING:
        param_list.push_back( Param( *iterNode , 0.0 ) );
        break;
      case XEXP_VARIABLE:
        // this case is a global param that must be resolved at each use because
        // it can change during a simulation
        param_list.push_back( Param( "GLOBAL_PARAMETER" , *iterNode ) );
        break;
      default:
        {
          std::string errorMsg("Can't find context for expression variable ");
          errorMsg += *iterNode;
          errorMsg += " in full expression ";
          errorMsg += expression_;
          errorMsg += " Please check to ensure this parameter is correct and set in your netlist.";
          N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_WARNING_0,errorMsg);
          // increment the number of unresolved strings reported in the
          // return value.  
          ++returnVal;
          // Just issue the warning for now.  Trying to resolve the 
          // string in param_list will cause problems in parallel 
          // as not all expressions on all procs will have the same 
          // unresolved vars.  The class owner will have to check that
          // at least one of the processors with this expression has
          // all the strings resolved. 
          // param_list.push_back( Param( *iterNode , 0.0 ) );
        }
    }
    iterNode++;
  }

  IO::makeOps(outputManager_, param_list.begin(), param_list.end(), std::back_inserter(expressionVars_));
  
  varVals.reserve( expressionVars_.size() );

  return returnVal;
}
//-----------------------------------------------------------------------------
// Namespace     : Unnamed
// Purpose       : file-local scoped methods and data
// Special Notes : 
// Creator       : Tom Russo, SNL
// Creation Date : 11/27/2013
//-----------------------------------------------------------------------------
namespace {

//-----------------------------------------------------------------------------
// Function      : convertNodalComputation
// Purpose       : given a nodal expression string (e.g. "VM(A,B)"),
//                 construct the set of Params that makeOps would expect for
//                 it
// Special Notes : 
//
// Scope         : file-local
// Creator       : Tom Russo
// Creation Date : 11/27/2013
//-----------------------------------------------------------------------------
void convertNodalComputation(std::string &nodalComputation, 
                             std::list<Param> &paramList)
{
  std::list<Param> tempParamList;
  
  std::size_t firstParen = nodalComputation.find_first_of("(");
  std::size_t lastParen = nodalComputation.find_first_of("(");
  // the length of the name of the param is actually equal to the position
  // of the first paren
  std::string compName=nodalComputation.substr(0,firstParen);
  std::string args=nodalComputation.substr(firstParen+1,nodalComputation.length()-compName.length()-2);

#ifdef Xyce_DEBUG_EXPRESSION
  Xyce::dout() << "Processing nodalComputation : " << nodalComputation 
               << Util::push<< std::endl;
  Xyce::dout() << "name of computation: " << compName << std::endl;
  Xyce::dout() << "args: " << args << Util::push << std::endl;
#endif // Xyce_DEBUG_EXPRESSION

  std::size_t firstComma=args.find_first_of(",");
  while (firstComma != std::string::npos)
  {
    std::string arg = args.substr(0,firstComma);
    std::size_t argsLength = args.length();
    args = args.substr(firstComma+1,argsLength-arg.length()-1);
    firstComma = args.find_first_of(",");
    tempParamList.push_back(Param(arg,0.0));
#ifdef Xyce_DEBUG_EXPRESSION
    Xyce::dout() << "arg " << arg << std::endl;
#endif
  }

  tempParamList.push_back(Param(args,0.0));

#ifdef Xyce_DEBUG_EXPRESSION
  Xyce::dout() << "Remaining arg " << args << std::endl;
  Xyce::dout() << "There were " << tempParamList.size() << " args." 
               << Util::pop << std::endl;
  Xyce::dout() << Util::pop << std::endl;
#endif

  paramList.push_back(Param(compName,static_cast<double>(tempParamList.size())));
  std::copy (tempParamList.begin(), tempParamList.end(), std::back_inserter(paramList)); 

} // namespace (unnammed)

}
} // namespace Util
} // namespace Xyce
