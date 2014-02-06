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
// Revision Number: $Revision: 1.13.2.6 $
//
// Revision Date  : $Date: 2013/12/03 23:30:12 $
//
// Current Owner  : $Author: rlschie $
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
#include<N_UTL_Expression.h>


//-----------------------------------------------------------------------------
// Function      : N_UTL_ExpressionData::N_UTL_ExpressionData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------
N_UTL_ExpressionData::N_UTL_ExpressionData (
    const string & expression1,
    N_IO_OutputMgr * outputMgrPtr
) :
  expression (expression1),
  expPtr (NULL),
  outputMgrPtr_(outputMgrPtr),
  setupCalled(false),
  val(0.0),
  numVars(0),
  numSpecialVars(0),
  numUnresolvedStrings(0),
  numUnresolvedStringsChecked(false),
  procID_(0)
{
}


//-----------------------------------------------------------------------------
// Function      : N_UTL_ExpressionData::N_UTL_ExpressionData
// Purpose       : constructor -- used to take an existing expression pointer
//                 not a string that can be converted to an expression
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------
N_UTL_ExpressionData::N_UTL_ExpressionData (
    const N_UTL_Expression * expressionPointer,
    N_IO_OutputMgr * outputMgrPtr
) :
  expPtr (NULL),
  outputMgrPtr_(outputMgrPtr),
  setupCalled(false),
  val(0.0),
  numVars(0),
  numSpecialVars(0),
  numUnresolvedStrings(0),
  numUnresolvedStringsChecked(false),
  procID_(0)
{
  // copy the expression object because the pointer passed in
  // is owned by another object.
  expPtr = new N_UTL_Expression( *expressionPointer );
  // get a string version of the expression to fill in
  // this data.
  expression = expPtr->get_expression();
}


//-----------------------------------------------------------------------------
// Function      : N_UTL_ExpressionData::~N_UTL_ExpressionData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------
N_UTL_ExpressionData::~N_UTL_ExpressionData ()
{
  if (expPtr)
  {
    delete expPtr;
  }

  varNames.clear();
  specialVarNames.clear();
  currentParam.clear();
  varGIDs.clear();
  varVals.clear();
  expressionVars_.clear();
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_ExpressionData::evaluate
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------
double N_UTL_ExpressionData::evaluate (const N_LAS_Vector * solnVecPtr, const N_LAS_Vector *stateVecPtr, const N_LAS_Vector * stoVecPtr)
{
  int i;

  if (!setupCalled)
  {
    numUnresolvedStrings = setup ();
    setupCalled = true;
  }

#ifdef Xyce_DEBUG_IO
  cout << "N_UTL_ExpressionData::evaluate" << endl;
#endif

  if( solnVecPtr != NULL )
  {
    // loop over expressionVars_ to get all the values.
    list<N_UTL_Param>::iterator iterCur = expressionVars_.begin();
    list<N_UTL_Param>::iterator iterEnd = expressionVars_.end();
    int index=0;
    while( iterCur != iterEnd )
    {
      varVals[index] = outputMgrPtr_->getPrintValue( iterCur, solnVecPtr, stateVecPtr, stoVecPtr );
      index++;
      iterCur++;
    }
    
    // STD and DDT are implicitly time dependent.  Check the underlying expression 
    // for time dependence, and if it is get the current time with getPrgetImmutableValue<int>()
    // and set it in the underlying expression 
    if( expPtr->isTimeDependent() )
    {
      expPtr->set_sim_time(outputMgrPtr_->getCircuitTime());
    }

    // now get expression value and partial derivatives.
    // Note that for these purposes (.PRINT output) we only
    // really need the value, not the derivs.
    val = 0.0;
    if (expPtr)
    {
      expPtr->evaluateFunction ( val, varVals );
    }
  }

  return val;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_ExpressionData::setup
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
int N_UTL_ExpressionData::setup ()
{
  // we return the number of unresolved strings on exit.  
  // we can't just exit with an error if there are unresolved strings 
  // because in parallel such strings could be resolved on another 
  // processor.  Thus the number of unresolved strings is saved in 
  // numUnresolvedStrings and the class owner can check that at least 
  // one processor in parallel has 0 for numUnresolvedStrings.
  int returnVal = 0;
  // setup the string pN.
  ostringstream ost;
  ost << procID_;
  string pN (ost.str());

#ifdef Xyce_DEBUG_IO
  cout << endl << "N_UTL_ExpressionData::setup" << endl;
#endif

  // allocate expression pointer if we need to
  if( expPtr == NULL )
  {
    expPtr = new N_UTL_Expression(expression);
  }

  // clear out the STL containers:
  varNames.clear();
  currentParam.clear();
  varGIDs.clear();
  varVals.clear();
  expressionVars_.clear();

  // query the expression object for all of it's dependent vars.
  expPtr->get_names(XEXP_ALL, varNames);

  // this varNames vec is a list of string representations of all of the vars in the expression.
  // use this list to make a list of N_UTL_Param objects which we can then
  // call the N_IO_OutputMgr::setParamContextType_() and then use getPrintValue_()
  // to get the results.  This reuses code that is already in place to handle solution,
  // state, store and lead currents.

  vector<string>::iterator iterNode = varNames.begin();
  vector<string>::iterator iterEnd = varNames.end();
  while( iterNode != iterEnd )
  {
    // based on the type of variable, create the needed N_UTL_Param
    // objects for setParmContextType_ to work.
    int varType = expPtr->get_type( *iterNode );
#ifdef Xyce_DEBUG_IO
    std::cout << "N_UTL_ExpressionData::setup " << *iterNode << endl;
#endif
    switch (varType)
    {
      case XEXP_NODE:
        expressionVars_.push_back( N_UTL_Param( "V" , 1 ) );
        expressionVars_.push_back( N_UTL_Param( *iterNode , 0.0 ) );
#ifdef Xyce_DEBUG_IO
        std::cout << " Type was XEXP_NODE " << endl;
#endif
        break;
      case XEXP_INSTANCE:
      case XEXP_LEAD:
        {
          char leadDesignator = expPtr->get_lead_designator( *iterNode);
          string currentName("I");
          if( (leadDesignator != 0) && (leadDesignator!=' ') )
          {
            currentName = currentName + leadDesignator;
          }
          expressionVars_.push_back( N_UTL_Param( currentName , 1 ) );
          expressionVars_.push_back( N_UTL_Param( *iterNode , 0.0 ) );
        }
#ifdef Xyce_DEBUG_IO
        std::cout << " Type was XEXP_LEAD " << endl;
#endif
        break;
      case XEXP_SPECIAL:
      case XEXP_STRING:
        expressionVars_.push_back( N_UTL_Param( *iterNode , 0.0 ) );
#ifdef Xyce_DEBUG_IO
        std::cout << " Type was XEXP_STRING " << endl;
#endif
        break;
      case XEXP_VARIABLE:
        // this case is a global param that must be resolved at each use because
        // it can change during a simulation
        expressionVars_.push_back( N_UTL_Param( "GLOBAL_PARAMETER" , *iterNode ) );
        break;
      default:
        {
          string errorMsg("Can't find context for expression variable ");
          errorMsg += *iterNode;
          errorMsg += " in full expression ";
          errorMsg += expression;
          errorMsg += " Please check to ensure this parameter is correct and set in your netlist.";
          N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_WARNING_0,errorMsg);
          // increment the number of unresolved strings reported in the
          // return value.  
          returnVal++;
          // Just issue the warning for now.  Trying to resolve the 
          // string in expressionVars_ will cause problems in parallel 
          // as not all expressions on all procs will have the same 
          // unresolved vars.  The class owner will have to check that
          // at least one of the processors with this expression has
          // all the strings resolved. 
          // expressionVars_.push_back( N_UTL_Param( *iterNode , 0.0 ) );
        }
    }
    iterNode++;
  }

  list<N_UTL_Param>::iterator iterParam = expressionVars_.begin();
  list<N_UTL_Param>::iterator paramEnd = expressionVars_.end();
  while( iterParam != paramEnd )
  {
    // check if *iterParam has its context set.  if not then call setParamContextType_
    SimulatorVariableContext theParamsContext = iterParam->getSimContext();
    if( theParamsContext == UNDEFINED )
    {
      bool contextFound = outputMgrPtr_->setParamContextType_(iterParam);
      if( !contextFound )
      {
        string msg("Can't find context for expression variable ");
        msg += iterParam->tag() + " " + iterParam->tag() + " " + iterParam->tag();
        // msg is not quite complete at this stage.  it could be that it has two tags
        // need better error report in this case.
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0,msg);
      }
    }
    // in seting the context of the expressionVars_, some of them will no
    // longer be needed. for example V( a ) is stored in two parameters one for "V"
    // and the second for "a".  After setParamContextType_ is called, the first
    // param will have all the data needed to evalueate V(a) and the second one will
    // be labeled as NODE_OR_DEVICE_NAME.  IF it is, then we can safely erase it
    // from the params
    if( iterParam->getSimContext() == NODE_OR_DEVICE_NAME )
    {
      // this erases the element pointed to iterParam and then points
      // iterParam at the next value.  Thus we don't need a ++iterParam.
      iterParam = expressionVars_.erase( iterParam );
    }
    else
    {
      ++iterParam;
    }
  }

  // now resize the varVals vec
  varVals.resize( expressionVars_.size() );
  return returnVal;
}

