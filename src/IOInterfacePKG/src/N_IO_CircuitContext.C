//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Filename       : $RCSfile: N_IO_CircuitContext.C,v $
//
// Purpose        : Define the circuit "context".
//
// Special Notes  :
//
// Creator        : Lon Waters
//
// Creation Date  : 01/21/2003
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revsion$
//
// Revision Date   : $Date: 2013/10/03 17:23:43 $
//
// Current Owner  : $Author: tvrusso $
//----------------------------------------------------------------------------

#include <Xyce_config.h>
// ---------- Standard Includes ---------

#include <N_UTL_Misc.h>

#include <iostream>

#ifdef HAVE_ALGORITHM
#include <algorithm>
#else
#ifdef HAVE_ALGO_H
#include <algo.h>
#else
#error Must have either <algorithm> or <algo.h>!
#endif
#endif

#include <sstream>


// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceBlock.h>

#include <N_ERH_ErrorMgr.h>

#include <N_IO_CircuitContext.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_ParameterBlock.h>
#include <N_IO_DeviceBlock.h>

#include <N_UTL_Expression.h>
#include <N_UTL_Misc.h>

#include <N_PDS_Comm.h>
#include <N_IO_CircuitMetadata.h>

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::N_IO_CircuitContext
// Purpose        : Constructor
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 01/21/2003
//----------------------------------------------------------------------------
N_IO_CircuitContext::N_IO_CircuitContext( N_IO_CircuitMetadata & md,
                        list<N_IO_CircuitContext*> & cL,
                        N_IO_CircuitContext *& ccPtr)
  : name_(""),
    deviceCount_(0),
    resolved_(false),
    subcircuitPrefix_(""),
    parentContextPtr_(NULL),
    resolvedParams_(md),
    resolvedGlobalParams_(md),
    metadata_(md),
    contextList_(cL),
    currentContextPtr_(ccPtr)
{
  if (currentContextPtr_ == NULL)
  {
    currentContextPtr_ = this;
  }
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::~N_IO_CircuitContext
// Purpose        : Destructor
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 01/21/2003
//----------------------------------------------------------------------------
N_IO_CircuitContext::~N_IO_CircuitContext()
{
  // Delete each subcircuit context in this context.
  map< string, N_IO_CircuitContext * >::iterator itsc = circuitContextTable_.begin();
  map< string, N_IO_CircuitContext * >::const_iterator itsc_end = circuitContextTable_.end();

  for ( ; itsc != itsc_end; ++itsc )
  {
    delete itsc->second;
  }

  circuitContextTable_.clear();


  // We delete all our own stored model pointers.  Because we also
  // delete all our subcircuit contexts, *their* destructors take care of
  // deleting their model pointers.
  map<string, N_IO_ParameterBlock*>::iterator modelIter =
    models_.begin();
  map<string, N_IO_ParameterBlock*>::iterator modelIterEnd =
    models_.end();
  for (; modelIter != modelIterEnd; ++modelIter)
    delete modelIter->second;
  models_.clear();
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CircuitContext::operator=
// Purpose       : assignment operator
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/18/2006
//-----------------------------------------------------------------------------
N_IO_CircuitContext & N_IO_CircuitContext ::operator=
   (const N_IO_CircuitContext & right)
{
    currentContextPtr_ = right.currentContextPtr_;
    parentContextPtr_ = right.parentContextPtr_;
    contextList_ = right.contextList_;
    name_ = right.name_;

    deviceCount_ = right.deviceCount_;
    instanceList_ = right.instanceList_;
    instanceErrorInfo_ = right.instanceErrorInfo_;

    nodeList_ = right.nodeList_;
    subcircuitParameters_ = right.subcircuitParameters_;
    circuitContextTable_ = right.circuitContextTable_;
    models_ = right.models_;

    unresolvedParams_ = right.unresolvedParams_;
    unresolvedGlobalParams_ = right.unresolvedGlobalParams_;
    unresolvedFunctions_ = right.unresolvedFunctions_;
    mutualInductances_ = right.mutualInductances_;

    sharedInductorTable_ = right.sharedInductorTable_;
    allCoupledInductors_ = right.allCoupledInductors_;
    allIndexedMIs_ = right.allIndexedMIs_;
    kLines_ = right.kLines_;

    subcircuitPrefix_ = right.subcircuitPrefix_;
    nodeMap_ = right.nodeMap_;
    resolved_ = right.resolved_;
    resolvedParams_ = right.resolvedParams_;
    resolvedGlobalParams_ = right.resolvedGlobalParams_;
    resolvedFunctions_ = right.resolvedFunctions_;
    metadata_ = right.metadata_;

  return *this;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::beginSubcircuitContext
// Purpose        : Add a subcircuit context to the current context.
// Special Notes  : This routine sets the current context to a newly
//                  created N_IO_CircuitContext object. This context
//                  remains the current context until it is explicitly
//                  terminated with a call to endSubcircuitContext.
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 01/24/2003
//----------------------------------------------------------------------------
void N_IO_CircuitContext::beginSubcircuitContext(
    string const& netlistFileName,
    vector<N_IO_SpiceSeparatedFieldTool::StringToken> & subcircuitLine)

{
#ifdef Xyce_DEBUG_IO
  cout << "N_IO_CircuitContext::beginSubcircuitContext" << endl;
  cout << "*******************************************" << endl;
  cout << " subcircuit file name: " << netlistFileName << endl;
  cout << "  size of line vector: " << subcircuitLine.size() << endl;
#endif

  // preprocessing hacketry
  vector<N_IO_SpiceSeparatedFieldTool::StringToken>::iterator iterLine =
   subcircuitLine.begin() + 2;
  while( iterLine != subcircuitLine.end() )
  {
    if ( iterLine->string_ == "(" || iterLine->string_ == ")" )
    {
      subcircuitLine.erase( iterLine );
    }
    else
    {
      iterLine++;
    }
  }

  // Create a new circuit context for the subcircuit.
  N_IO_CircuitContext* subcircuitContextPtr =
    new N_IO_CircuitContext(metadata_, contextList_, currentContextPtr_);

  // Set the parent context, save the current context and reset it to the
  // newly created context.
  subcircuitContextPtr->parentContextPtr_ = currentContextPtr_;
  contextList_.push_front(currentContextPtr_);
  currentContextPtr_ = subcircuitContextPtr;

  // Extract the subcircuit data from subcircuitLine.
  int numFields = subcircuitLine.size();

  if (numFields < 3)
  {
    string msg("Not enough data on .subckt line");
    if (numFields > 1)
    {
      msg += " for subcircuit name: " + subcircuitLine[1].string_ + "\n";
    }
    else
    {
      msg += ", no subcircuit name given\n";
    }
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
        netlistFileName, subcircuitLine[0].lineNumber_);
  }

#ifdef Xyce_DEBUG_IO
  cout << "  Subcircuit name: " << subcircuitLine[1].string_ << endl;
  cout << "  Number of fields on subckt line = " << numFields << endl;
#endif

  // Extract and set the subcircuit name.
  ExtendedString field ( subcircuitLine[1].string_ );
  field.toUpper();
  subcircuitContextPtr->setName(field);

  // Extract the subcircuit external nodes.
  int i;
  for (i = 2; i < numFields; ++i)
  {
    ExtendedString fieldES (  subcircuitLine[i].string_ );
    fieldES.toUpper();
    // Exit the loop if there are netlist parameters
    if (fieldES == "PARAMS:")
    {
       break;
    }
    else if (i < numFields-1 && subcircuitLine[i+1].string_ == "=")
    {
      --i;
      break;
    }
    subcircuitContextPtr->nodeList_.push_back(fieldES);
  }

  // Extract the subcircuit parameters.
  N_UTL_Param parameter("","");

  ++i; // Advance to start of parameters.
  ExtendedString fieldES( "" );
  while (i+2 < numFields)
  {
    fieldES =  subcircuitLine[i].string_;
    fieldES.toUpper();
    if (!fieldES.possibleParam())
    {
      string msg("Parameter name ");
      msg += parameter.uTag();
      msg += " contains illegal character(s)";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
          netlistFileName, subcircuitLine[i].lineNumber_);
    }
    parameter.setTag( fieldES );

    if ( (parameter.uTag() == "TEMP") || (parameter.uTag() == "VT") ||
         (parameter.uTag() == "GMIN") || (parameter.uTag() == "TIME") )
    {
      string msg("Parameter name ");
      msg += parameter.uTag();
      msg += " not allowed in subcircuit parameter list\n";
      msg += " for subcircuit ";
      msg += getName() + "\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
          netlistFileName, subcircuitLine[i].lineNumber_);
    }

    if ( subcircuitLine[i+1].string_ != "=" )
    {
      string msg("Equal sign required between parameter and value in\n");
      msg += "PARAM list for subcircuit ";
      msg += getName() + "\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
          netlistFileName, subcircuitLine[0].lineNumber_);
    }

    i+=2; // Advance past "=" sign

    fieldES =  subcircuitLine[i].string_;
    fieldES.toUpper();
    parameter.setVal(fieldES);
    subcircuitContextPtr->subcircuitParameters_.push_back(parameter);
#ifdef Xyce_DEBUG_IO
    cout << " Adding parameter " << parameter.uTag() << " with value " << fieldES << endl;
#endif
    ++i;
  }

  // check for truncated params: list
  if ( i < numFields )
  {
    string msg("Parameter list error in subcircuit ");
    msg += getName() + "\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
    netlistFileName, subcircuitLine[0].lineNumber_);
  }

#ifdef Xyce_DEBUG_IO
  cout << "End of N_IO_CircuitContext::beginSubcircuitContext" << endl;
  cout << "*******************************************" << endl;
#endif
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::endSubcircuitContext
// Purpose        : End the current context, push it onto the previous context's
//                  list of contexts, and reset the current context pointer to
//                  the previous context.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 01/24/2003
//----------------------------------------------------------------------------
void N_IO_CircuitContext::endSubcircuitContext()
{
#ifdef Xyce_DEBUG_IO
  cout << "N_IO_CircuitContext::endSubcircuitContext" << endl;
#endif
  // Get the previous context from the stack.
  if ( ! contextList_.empty() )
  {
    // Add the current context to the previous context's list of contexts.
    ( contextList_.front() )->circuitContextTable_[ currentContextPtr_->name_ ] =
     currentContextPtr_;

    // Reset the current context to the previous context.
    currentContextPtr_ = contextList_.front();

    // remove from stack
    contextList_.pop_front();
  }

  else
  {
    string msg("Error encountered while traversing subcircuits.\n");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg );
  }
#ifdef Xyce_DEBUG_IO
  cout << "End of N_IO_CircuitContext::endSubcircuitContext" << endl;
#endif
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::addModel
// Purpose        : Add a model to the current context.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/03/2003
//----------------------------------------------------------------------------
void N_IO_CircuitContext::addModel(N_IO_ParameterBlock * modelPtr)
{
  N_IO_ParameterBlock* tmpModel;
  if (findModel(modelPtr->getName(), tmpModel))
  {
    string msg("Reading model named ");
    msg += modelPtr->getName();
    msg += " in the ";
    if (getCurrentContextName() == "")
      msg += "main circuit";
    else
    {
      msg += "subcircuit ";
      msg += getCurrentContextName();
    }
    msg += " and found one or more models previously defined in this scope.";
    msg += "\n\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg );
  }

  currentContextPtr_->models_[modelPtr->getName()] = modelPtr;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::addParams
// Purpose        : Add a set of .PARAM parameters to the current context.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/03/2003
//----------------------------------------------------------------------------
void N_IO_CircuitContext::addParams(N_IO_OptionBlock const& param)
{
  int numberOfParams = param.getNumberOfParameters();

  N_UTL_Param parameter;
  for (int i = 0; i < numberOfParams; ++i)
  {
    parameter = param.getParameter(i);
    resolveQuote(parameter);
    currentContextPtr_->unresolvedParams_.push_back(parameter);
  }
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::addGlobalParams
// Purpose        : Add a set of .GLOBAL_PARAM parameters to the current context.
// Special Notes  :
// Scope          :
// Creator        : Dave Shirley, PSSI
// Creation Date  : 11/17/2005
//----------------------------------------------------------------------------
void N_IO_CircuitContext::addGlobalParams(N_IO_OptionBlock const& param)
{
  int numberOfParams = param.getNumberOfParameters();

  N_UTL_Param parameter;
  for (int i = 0; i < numberOfParams; ++i)
  {
    parameter = param.getParameter(i);
    resolveQuote(parameter);
    currentContextPtr_->unresolvedGlobalParams_.push_back(parameter);
  }
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::addGlobalNode
// Purpose        : Add a global node from a .GLOBAL line
// Special Notes  :
// Scope          :
// Creator        : Dave Shirley, PSSI
// Creation Date  : 10/02/2009
//----------------------------------------------------------------------------
void N_IO_CircuitContext::addGlobalNode( string & gnode)
{
  currentContextPtr_->globalNodes_.insert(gnode);
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::resolveQuote
// Purpose        : Resolve quoted parameters as soon as they are encountered.
// Special Notes  : This avoids every processor in a parallel run trying to
//                  access the same file.  Also exit handling does not work on
//                  error for parallel runs otherwise.
// Scope          :
// Creator        : Dave Shirley, PSSI
// Creation Date  :
//----------------------------------------------------------------------------
void N_IO_CircuitContext::resolveQuote(N_UTL_Param & parameter)
{
  if (parameter.isQuoted())
  {
    // The parameter is time dependent with its time history defined by the set
    // of time-value pairs in the file given by the value of the parameter.
    // Open and read these values, using them to build a "Table" expression
    // for the value of the parameter.
    ifstream paramDataIn;
    string parameterFile(parameter.sVal().substr(1,parameter.sVal().size()-2));
    paramDataIn.open(parameterFile.c_str(), ios::in);
    if ( !paramDataIn.is_open() )
    {
      string msg("Could not find file " + parameterFile + "\n");
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg );
    }

    string table("table(time");
    string time;
    string value;
    while ( paramDataIn >> time )
    {
      if ( paramDataIn >> value )
      {
        table += "," + time + "," + value;
      }
      else
      {
        string msg("Reached end of file in ");
        msg += parameterFile + " while expecting another value\n";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg );
      }
    }

    paramDataIn.close();

    if( table.size() <= 10 ) // the length of "table(time" from above
    {
      string msg("Failed to successfully read contents of " + parameterFile + "\n");
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg );
    }

    table += ")";

    parameter.setVal( N_UTL_Expression(table) );
    return;
  }
}
//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::addFunction
// Purpose        : Add a .FUNC function to the current context.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/03/2003
//----------------------------------------------------------------------------
void N_IO_CircuitContext::addFunction(N_IO_FunctionBlock const& function)
{
  currentContextPtr_->unresolvedFunctions_.push_back(function);
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::resolve
// Purpose        : Resolve parameters in the current context.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/04/2003
//----------------------------------------------------------------------------
bool N_IO_CircuitContext::resolve( vector<N_DEV_Param> const& subcircuitInstanceParams )
{
  vector<N_UTL_Param> retryParams;
  vector<N_IO_FunctionBlock> retryFunctions;

  if (currentContextPtr_->subcircuitParameters_.empty() &&
      subcircuitInstanceParams.empty() &&
      currentContextPtr_->resolved_)
  {
    // If there are no subcircuit parameters and this context has already
    // been resolved, then no work to do.
    return true;
  }

#ifdef Xyce_DEBUG_IO
  cout << " N_IO_CircuitContext::resolve have something to do..." << endl;
#endif

  currentContextPtr_->resolved_ = true;

  // Clear resolved params and funcs.
  currentContextPtr_->resolvedParams_.clearParameters();
  currentContextPtr_->resolvedGlobalParams_.clearParameters();
  currentContextPtr_->resolvedFunctions_.clear();


  vector<N_UTL_Param> asYetUnresolvedSubcircuitParameters=currentContextPtr_->subcircuitParameters_;
  vector<N_UTL_Param> asYetUnresolvedParameters=currentContextPtr_->unresolvedParams_;
  vector<N_UTL_Param> asYetUnresolvedGlobalParameters=currentContextPtr_->unresolvedGlobalParams_;
  vector<N_IO_FunctionBlock> asYetUnresolvedFunctions=currentContextPtr_->unresolvedFunctions_;

  bool resolvedSomethingThisLoop=true;
  bool somethingLeftToDo= (!(asYetUnresolvedSubcircuitParameters.empty()) ||
                           !(asYetUnresolvedParameters.empty()) ||
                           !(asYetUnresolvedGlobalParameters.empty()) ||
                           !(asYetUnresolvedFunctions.empty())  ||
                           !(subcircuitInstanceParams.empty()) );

  while (resolvedSomethingThisLoop && somethingLeftToDo)
  {
    resolvedSomethingThisLoop=false;
    somethingLeftToDo=false;
    retryParams.clear();
    retryFunctions.clear();

    // Add subcircuitParameters_ to the set of resolved parameters.
    // currentContextPtr_->resolvedParams_.addParameters( currentContextPtr_->subcircuitParameters_ );
    N_UTL_Param parameter;
    vector<N_UTL_Param>::iterator paramIter;
    vector<N_UTL_Param>::iterator start =
      asYetUnresolvedSubcircuitParameters.begin();
    vector<N_UTL_Param>::iterator end =
      asYetUnresolvedSubcircuitParameters.end();
    for (paramIter = start; paramIter != end; ++paramIter)
    {
      N_UTL_Param par = *paramIter;
#ifdef Xyce_DEBUG_IO
      cout << " N_IO_CircuitContext::resolve attempting to resolve" << par.uTag()<< endl;
#endif
      if (par.getType() == STR && !par.isNumeric())
      {
        ExtendedString arg ( par.sVal() );
        arg.toUpper();
        if (arg.possibleParam())
          par.setVal(string("{" + arg + "}"));
      }

      if (!resolveParameter(par))
      {
#ifdef Xyce_DEBUG_IO
        cout << "Unable to resolve subcircuit param " << par.uTag() << endl;
#endif
        retryParams.push_back(par);
        somethingLeftToDo=true;
      }
      else
      {
#ifdef Xyce_DEBUG_IO
        cout << "resolveParameter returned true on parameter " << par.uTag() << " after resolution its type is " << par.getType() << "with value " ;
        switch (par.getType()) {
        case STR:
          cout << par.sVal();
          break;
        case DBLE:
          cout << par.dVal();
          break;
        case EXPR:
          cout << par.ePtr()->get_expression();
          break;
        default:
          cout << par.sVal();
        }
        cout << endl;
#endif

        currentContextPtr_->resolvedParams_.addParameter(par);
        resolvedSomethingThisLoop=true;


      }
#ifdef Xyce_DEBUG_IO
      cout << " N_IO_CircuitContext::resolve done attempting to resolve" << par.uTag()<< endl;
#endif
    }
    asYetUnresolvedSubcircuitParameters=retryParams;
    retryParams.clear();

    N_UTL_Param* paramPtr;
    // Reset the subcircuit parameter values with the instance parameter
    // values as needed.
    int i;
    int numInstanceParameters = subcircuitInstanceParams.size();
    for ( i = 0; i < numInstanceParameters; ++i )
    {
#ifdef Xyce_DEBUG_IO
      cout << " N_IO_CircuitContext::resolve resetting subcircuit instance parameter" << subcircuitInstanceParams[i].uTag()<< " with value " << endl;
        switch (subcircuitInstanceParams[i].getType()) {
        case STR:
          cout << subcircuitInstanceParams[i].sVal();
          break;
        case DBLE:
          cout << subcircuitInstanceParams[i].dVal();
          break;
        case EXPR:
          {
            N_UTL_Expression foo(subcircuitInstanceParams[i].eVal());
            cout << "EXPR(" << foo.get_expression() << ")";
            break;
          }
        default:
          cout << subcircuitInstanceParams[i].sVal();
        }
        cout << endl;
#endif

      // Look for the parameter in resolvedParams_, issue
      // a warning and continue if not found. Otherwise, if found
      // set the value.
      parameter.setTag(subcircuitInstanceParams[i].uTag());
      paramPtr = currentContextPtr_->resolvedParams_.findParameter(parameter);

      if ( paramPtr == NULL )
      {
        currentContextPtr_->resolvedParams_.addParameter(subcircuitInstanceParams[i]);
#ifdef Xyce_DEBUG_IO
        cout << " did not find in resolvdParams_, adding parameter " << subcircuitInstanceParams[i].uTag() << endl;
#endif
      }
      else
      {
#ifdef Xyce_DEBUG_IO
        cout << " found in resolvdParams_, setting parameter " << paramPtr->uTag() << endl;
#endif
        paramPtr->setVal(subcircuitInstanceParams[i]);
      }
    }

    // Resolve any .PARAM parameters in the current context and add to
    // the resolvedParams_.
    start = asYetUnresolvedParameters.begin();
    end = asYetUnresolvedParameters.end();
    for (paramIter = start; paramIter != end; ++paramIter)
    {
      parameter = *paramIter;

      if (!resolveParameter(parameter))
      {
        // save it for later, because it might use functions that haven't
        // been resolved yet.
        retryParams.push_back(parameter);
        somethingLeftToDo=true;
      }
      else
      {
        currentContextPtr_->resolvedParams_.addParameter(parameter);
        resolvedSomethingThisLoop=true;
      }
    }
    asYetUnresolvedParameters=retryParams;
    retryParams.clear();

    // Resolve any .GLOBAL_PARAM parameters in the current context and add to
    // the resolvedGlobalParams_.
    start = asYetUnresolvedGlobalParameters.begin();
    end = asYetUnresolvedGlobalParameters.end();
    for (paramIter = start; paramIter != end; ++paramIter)
    {
      parameter = *paramIter;
#ifdef Xyce_DEBUG_IO
      cout << " N_IO_CircuitContext::resolve Attempting to resolve global parameter " << parameter.uTag();
#endif

      if (!resolveParameter(parameter))
      {
        retryParams.push_back(parameter);
        somethingLeftToDo = true;
      }
      else
      {
        if (parameter.getType() == EXPR)
        {
          vector<string> nodes, instances, leads, variables, specials;
          vector<string>::iterator s_i;
          parameter.ePtr()->get_names(XEXP_NODE, nodes);
          parameter.ePtr()->get_names(XEXP_INSTANCE, instances);
          parameter.ePtr()->get_names(XEXP_LEAD, leads);
          parameter.ePtr()->get_names(XEXP_VARIABLE, variables);
          parameter.ePtr()->get_names(XEXP_SPECIAL, specials);

          if (!nodes.empty() || !instances.empty() || !leads.empty())
          {
            // This should be caught earlier, but just in case it is checked here
            string msg("");
            if (!nodes.empty())
            {
              msg += "node(s): ";
              for (s_i=nodes.begin() ; s_i!=nodes.end() ; ++s_i)
              {
                msg += *s_i;
                msg += " ";
              }
            }
            if (!instances.empty())
            {
              msg += "instance(s): ";
              for (s_i=instances.begin() ; s_i!=instances.end() ; ++s_i)
              {
                msg += *s_i;
                msg += " ";
              }
            }
            if (!leads.empty())
            {
              msg += "lead(s): ";
              for (s_i=leads.begin() ; s_i!=leads.end() ; ++s_i)
              {
                msg += *s_i;
                msg += " ";
              }
            }
            msg += "Not allowed in global param expression: ";
            msg += parameter.ePtr()->get_expression();
            N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::DEV_FATAL_0, msg);
          }
          if (!variables.empty())
          {
            // If variables are found, they must be previously defined global params
            for (s_i=variables.begin() ; s_i!=variables.end() ; ++s_i)
            {
              if (currentContextPtr_->resolvedGlobalParams_.findParameter(*s_i) == NULL)
              {
                string msg("Unknown parameter (");
                msg += *s_i;
                msg += ") found in global param expression: ";
                msg += parameter.ePtr()->get_expression();
                N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_FATAL_0, msg);
              }
            }
          }
          if (!specials.empty())
          {
            if (specials.size() > 1 || specials[0] != "TIME")
            {
              string msg("Unknown special var(s):");
              for (s_i=specials.begin() ; s_i!=specials.end() ; ++s_i)
              {
                msg += " ";
                msg += *s_i;
              }
              N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::DEV_FATAL_0, msg);
            }
          }
        }
        currentContextPtr_->resolvedGlobalParams_.addParameter(parameter);
        resolvedSomethingThisLoop=true;
      }
    }
    asYetUnresolvedGlobalParameters=retryParams;
    retryParams.clear();

    // Resolve functions in the current context.
    int numFunctions = asYetUnresolvedFunctions.size();
    for ( i = 0; i < numFunctions; ++i )
    {
      // Define the N_UTL_Param for the function prototype (the function name
      // and arguments) and function body. Pass to the circuit for resolution.
      string functionName(asYetUnresolvedFunctions[i].functionName);
      string functionNameAndArgs(asYetUnresolvedFunctions[i].functionNameAndArgs);
      string functionBody(asYetUnresolvedFunctions[i].functionBody);

      vector<string> functionArgs =
        asYetUnresolvedFunctions[i].functionArgs;
      N_UTL_Param functionParameter(functionNameAndArgs, functionBody);
      // this is only an error if we have no parameters left to resolve
      if (!resolveParameter(functionParameter, functionArgs))
      {
          retryFunctions.push_back(asYetUnresolvedFunctions[i]);
          somethingLeftToDo=true;
      }
      else
      {
        // After resolution, the only strings allowed in the function
        // are the function arguments, check that this holds.
        N_UTL_Expression functionBodyExpression( functionParameter.sVal() );
        vector<string> strings;
        bool canResolveAll=true;

        functionBodyExpression.get_names(XEXP_STRING, strings);

        int numStrings = strings.size();
        for (int j = 0; j < numStrings; ++j)
        {
          // Look for string in functionArgs.
          vector<string>::iterator stringIter =
            find(functionArgs.begin(), functionArgs.end(), strings[j]);
          if (stringIter == functionArgs.end() &&
              currentContextPtr_->resolvedGlobalParams_.findParameter(strings[j]) == NULL)
          {
            retryFunctions.push_back(asYetUnresolvedFunctions[j]);
            somethingLeftToDo=true;
            canResolveAll=false;
          }
        }

        if (canResolveAll)
        {
          // Reset the function body.
          functionBody = functionParameter.sVal();

          // Add the function to resolvedFunctions_.
          currentContextPtr_->resolvedFunctions_[functionName] =
            N_UTL_Param(functionNameAndArgs, functionBody);
          resolvedSomethingThisLoop=true;
        }
      }
    }
    asYetUnresolvedFunctions=retryFunctions;
    retryFunctions.clear();

  }

  if (somethingLeftToDo)
  {
    // we failed to resolve everything, and the last loop did nothing.
    N_UTL_Param parameter;
    vector<N_UTL_Param>::iterator paramIter;
    vector<N_UTL_Param>::iterator start =
      asYetUnresolvedSubcircuitParameters.begin();
    vector<N_UTL_Param>::iterator end =
      asYetUnresolvedSubcircuitParameters.end();
    for (paramIter = start; paramIter != end; ++paramIter)
    {
      string msg("Unable to resolve .subckt parameter ");
      msg += paramIter->uTag();
      msg += " found in .PARAM statement\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
    }

    start=asYetUnresolvedParameters.begin();
    end=asYetUnresolvedParameters.end();
    for (paramIter = start; paramIter != end; ++paramIter)
    {
      string msg("Unable to resolve parameter ");
      msg += paramIter->uTag();
      msg += " found in .PARAM statement\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
    }

    start=asYetUnresolvedGlobalParameters.begin();
    end=asYetUnresolvedGlobalParameters.end();
    for (paramIter = start; paramIter != end; ++paramIter)
    {
      string msg("Unable to resolve global parameter ");
      msg += paramIter->uTag();
      msg += " found in .PARAM statement\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
    }

    vector<N_IO_FunctionBlock>::iterator funcIter;
    vector<N_IO_FunctionBlock>::iterator funcStart =
      asYetUnresolvedFunctions.begin();
    vector<N_IO_FunctionBlock>::iterator funcEnd =
      asYetUnresolvedFunctions.end();
    for (funcIter = funcStart; funcIter != funcEnd; ++funcIter)
    {
      string functionName(funcIter->functionName);
      string functionNameAndArgs(funcIter->functionNameAndArgs);
      string functionBody(funcIter->functionBody);

      vector<string> functionArgs=funcIter->functionArgs;
      string msg("Error resolving function " + functionName);
      msg += ": \n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg );

      N_UTL_Param functionParameter(functionNameAndArgs, functionBody);
      if (!resolveParameter(functionParameter, functionArgs))
      {
        string msg("The function " + functionNameAndArgs);
        msg += " contains an undefined parameter or function.\n";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg );
      }
      else
      {
        N_UTL_Expression functionBodyExpression( functionParameter.sVal() );
        vector<string> strings;
        functionBodyExpression.get_names(XEXP_STRING, strings);

        int numStrings = strings.size();
        for (int i = 0; i < numStrings; ++i)
        {
          // Look for string in functionArgs.
          vector<string>::iterator stringIter =
            find(functionArgs.begin(), functionArgs.end(), strings[i]);
          if (stringIter == functionArgs.end() &&
              currentContextPtr_->resolvedGlobalParams_.findParameter(strings[i]) == NULL)
          {
            string msg("Definition of function ");
            msg += functionNameAndArgs;
            msg += " contains unknown parameter " + strings[i] + "\n";
            N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg );
          }
        }
      }
    }
    string msg("Terminating parsing due to parameter/function resolution issues ");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg );


  }

  return true;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::setContext
// Purpose        : Set the current context to that corresponding to the given
//                  subcircuit name. Save the previous context on the stack for
//                  later retrieval. If the given subcircuit name is not found,
//                  delcare an error and abort.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/07/2003
//----------------------------------------------------------------------------
bool N_IO_CircuitContext::setContext(string const& subcircuitName,
    string const& subcircuitPrefixIn,
    list<string> const& instanceNodes,
    N_IO_CircuitContext * previousContext )
{
  bool success = false;

  map< string, N_IO_CircuitContext* >::iterator ccIter =
   currentContextPtr_->circuitContextTable_.find( subcircuitName );

  if ( ccIter != currentContextPtr_->circuitContextTable_.end() )
  {
      // Save current context and reset the current context pointer.
      if (previousContext == NULL)
        contextList_.push_front(currentContextPtr_);
      else
        contextList_.push_front(previousContext);

      currentContextPtr_ = ccIter->second;

      // If subcircuitPrefix and instanceNodes were given, set the prefix
      // in the new context, and build the node map.
      currentContextPtr_->nodeMap_.clear();
      currentContextPtr_->subcircuitPrefix_ = subcircuitPrefixIn;

      if (!instanceNodes.empty())
      {
        list<string>::const_iterator nodeIter = instanceNodes.begin();
        for (int i = 0; i < currentContextPtr_->nodeList_.size() && nodeIter != instanceNodes.end(); ++i, ++nodeIter)
        {
          if (currentContextPtr_->nodeMap_.find(currentContextPtr_->nodeList_[i]) !=
              currentContextPtr_->nodeMap_.end())
          {
            if (currentContextPtr_->nodeMap_[currentContextPtr_->nodeList_[i]] !=
                *nodeIter)
            {
              string msg = "Duplicate nodes in .subckt point to different nodes in X line invocation";
              N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_INFO_0, msg);
              return false;
            }
          }
          else
          {
            if ((currentContextPtr_->nodeList_[i].size() >= 2 &&
                 currentContextPtr_->nodeList_[i].substr(0,2) == "$G") ||
                 globalNode(currentContextPtr_->nodeList_[i]))
            {
              if (currentContextPtr_->nodeList_[i] != *nodeIter)
              {
                string msg = "Global node in subcircuit invocation must match same name in .subckt";
                N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_INFO_0, msg);
                return false;
              }
            }
            currentContextPtr_->nodeMap_[currentContextPtr_->nodeList_[i]] = *nodeIter;
          }
        }
      }

      return true;
  }

  // The context with the given name was not found in the current context,
  // recursively search parent contexts.
  if (currentContextPtr_->parentContextPtr_ != NULL)
  {
    if (previousContext == NULL)
    {
      previousContext = currentContextPtr_;
    }

    currentContextPtr_ = currentContextPtr_->parentContextPtr_;
    success = setContext(subcircuitName, subcircuitPrefixIn, instanceNodes, previousContext);
  }

  if (!success)
  {
    string msg = "Global node in subcircuit invocation must match same name in .subckt";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_INFO_0, msg);
  }

  return success;
}

void N_IO_CircuitContext::setContext(N_IO_CircuitContext* context)
{
  contextList_.push_front(currentContextPtr_);
  currentContextPtr_ = context;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::restorePreviousContext
// Purpose        : Reset the context the context prior to the last invocation
//                  of setContext.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/07/2003
//----------------------------------------------------------------------------
void N_IO_CircuitContext::restorePreviousContext()
{
  // Restore the previous context from the list of contexts unless
  // the list is empty.
  if (!contextList_.empty())
  {
    currentContextPtr_ = contextList_.front();
    contextList_.pop_front();
  }
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::restorePreviousContext
// Purpose        : Reset the context the context prior to the last invocation
//                  of setContext.
// Special Notes  :
// Scope          :
// Creator        : Dave Shirley, PSSI
// Creation Date  : 10/20/2009
//----------------------------------------------------------------------------
bool N_IO_CircuitContext::globalNode (const string &nodeName)
{
  bool stat;

  if (!(currentContextPtr_->parentContextPtr_))
  {
    if (globalNodes_.find(nodeName) == globalNodes_.end())
      return false;
    else
      return true;
  }

  setContext(currentContextPtr_->parentContextPtr_);
  stat = globalNode ( nodeName );
  restorePreviousContext();

  return stat;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::resolveParameter
// Purpose        : Parameter whose value may be an expression that must be
//                  resolved. If the input parameter has an expression value,
//                  replace replace the parameters and functions in the
//                  expression with their actual values.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/10/2003
//----------------------------------------------------------------------------
bool N_IO_CircuitContext::resolveParameter(N_UTL_Param& parameter,
    vector<string> exceptionStrings)
{
  if (parameter.hasExpressionValue() || parameter.hasExpressionTag() )
  {
#ifdef Xyce_DEBUG_IO
    cout << "N_IO_CircuitContext::resolveParameter parameter " << parameter.uTag()
         << " has expression value ";
#endif
    // Extract the expression from the parameter value by stripping off
    // the enclosing braces.  Only strip if it's there!
    string expressionString;
    if (parameter.sVal()[0] == '{')
      expressionString = parameter.sVal().substr(1, parameter.sVal().size()-2);
    else
      expressionString = parameter.sVal();

#ifdef Xyce_DEBUG_IO
    cout << expressionString << endl;
#endif

    // Parse the expression:
    N_UTL_Expression expression(expressionString);

    // Resolve the strings in the expression. Unresolved strings
    // may be parameters defined in a .PARAM statement or global
    // parameters defined in .GLOBAL_PARAM statement or may
    // be due to function arguments if the expression is the
    // body of a function defined in a .FUNC statement.
    bool stringsResolved = resolveStrings(expression, exceptionStrings);

    // Resolve functions in the expression.
    bool functionsResolved = resolveFunctions(expression);

    // resolve variables in the function body
    if ( expression.get_num(XEXP_STRING) > 0 )
    {
      resolveStrings(expression, exceptionStrings);
    }


    if (expression.get_num( XEXP_LEAD ) > 0)
    {
      parameter.setVal(expression);
      return false;
    }

    if (stringsResolved && functionsResolved)
    {
      // Check the expression for nodes or instance (probably a B-source
      // expression, in which case the parameter value should be
      // expressionString. Also check for "specials", the only special
      // allowed is "time" for time dependent parameters.
      vector<string> nodes, instances, leads, variables, specials;
      expression.get_names(XEXP_NODE, nodes);
      expression.get_names(XEXP_INSTANCE, instances);
      expression.get_names(XEXP_LEAD, leads);
      expression.get_names(XEXP_VARIABLE, variables);
      expression.get_names(XEXP_SPECIAL, specials);

      if (!nodes.empty() || !instances.empty() || !leads.empty() ||
          !variables.empty() || !specials.empty())
      {
#ifdef Xyce_DEBUG_IO
        cout << "N_IO_CircuitContext::resolveParameter:  nodes, instances, leads, variables or specials not empty. " << endl;
        if (!nodes.empty())
        {
          cout << " Nodes: " << endl;
          for (int foo=0; foo<nodes.size(); ++foo)
            cout << foo << " : " << nodes[foo] << endl;
        }
        if (!instances.empty())
        {
          cout << " Instances: " << endl;
          for (int foo=0; foo<instances.size(); ++foo)
            cout << foo << " : " << instances[foo] << endl;
        }
        if (!leads.empty())
        {
          cout << " Leads: " << endl;
          for (int foo=0; foo<leads.size(); ++foo)
            cout << foo << " : " << leads[foo] << endl;
        }
        if (!variables.empty())
        {
          cout << " Variables: " << endl;
          for (int foo=0; foo<variables.size(); ++foo)
            cout << foo << " : " << variables[foo] << endl;
        }
        if (!specials.empty())
        {
          cout << " Specials: " << endl;
          for (int foo=0; foo<specials.size(); ++foo)
            cout << foo << " : " << specials[foo] << endl;
        }
#endif

        parameter.setVal(expression);
#ifdef Xyce_DEBUG_IO
        cout << "N_IO_CircuitContext::resolveParameter: After all expression handling, get_expression returns "
             << expression.get_expression() << endl;
        cout << " after setting the parameter " << parameter.uTag() << ", its type is " << parameter.getType() << endl;
        cout << " and its value is ";
        switch (parameter.getType()) {
        case STR:
          cout << parameter.sVal();
          break;
        case DBLE:
          cout << parameter.dVal();
          break;
        case EXPR:
          cout << parameter.ePtr()->get_expression();
          break;
        default:
          cout << parameter.sVal();
        }
        cout << endl;
#endif
      }
      else
      {
        if (exceptionStrings.empty())
        {
          // Reset the parameter value to the value of the expression.
          double value;
          expression.evaluateFunction ( value );
          parameter.setVal( value );
          // we have resolved the context so set it and the constant value to 
          // make later look ups easier.
          parameter.setSimContextAndData( CONSTANT, value );
#ifdef Xyce_DEBUG_IO
          cout << " N_IO_CircuitContext::resolveParameter --  Resetting parameter value from " << expressionString << " to " << value << " because exceptionStrings empty." << endl;
#endif
        }
        else
        {
          parameter.setVal( expression );
        }
      }
    }
    else
    {
      // Reset the parameter value to the value of the expression with
      // as much resolution as could be achieved.
      parameter.setVal(expression);
    }

#ifdef Xyce_DEBUG_IO
        cout << "N_IO_CircuitContext::resolveParameter: right before returns "
             << endl;
        cout << " after setting the parameter " << parameter.uTag() << ", its type is " << parameter.getType() << endl;
        cout << " and its value is ";
        switch (parameter.getType()) {
        case STR:
          cout << parameter.sVal();
          break;
        case DBLE:
          cout << parameter.dVal();
          break;
        case EXPR:
          cout << parameter.ePtr()->get_expression();
          break;
        default:
          cout << parameter.sVal();
        }
        cout << endl;
#endif

    return stringsResolved && functionsResolved;

  }
  // Handle quoted parameters e.g. "filename" (which get turned into
  // TABLEs)
  resolveQuote(parameter);

  // The parameter does not have an expression value and does not need
  // resolving.
  return true;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::resolveStrings
// Purpose        : Determine if expression has any unresolved strings
//                  and resolve appropriately. Return true if all strings are
//                  resolved otherwise return false.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/11/2003
//----------------------------------------------------------------------------
bool N_IO_CircuitContext::resolveStrings( N_UTL_Expression & expression,
    vector<string> exceptionStrings)
{
  // Strings in the expression must be previously resolved parameters
  // that appear in paramList else there is an error.
  bool unresolvedStrings = false;
  if ( expression.get_num(XEXP_STRING) > 0 )
  {
    // Get the list of strings in the expression.
    vector<string> strings;
    expression.get_names(XEXP_STRING, strings);

    // If the expression is resolvable, each string in the current expression
    // must appear as a resolved parameter in netlistParameters. Get the value
    // if it appears there.
    ExtendedString parameterName("");
    int numStrings = strings.size();
    for (int i = 0; i < numStrings; ++i)
    {
#ifdef Xyce_DEBUG_IO
      cout <<" N_IO_CircuitContext::resolveStrings resolving " << strings[i] << endl;
#endif
      // Skip the current string if it is in exception strings. This prevents
      // a function argument from being improperly resolved when there is
      // a parameter in a .param statement with the same name as the function
      // argument.
      if (!exceptionStrings.empty())
      {
        vector<string>::iterator stringIter = find(exceptionStrings.begin(),
            exceptionStrings.end(),
            strings[i]);
        if (stringIter != exceptionStrings.end())
        {
#ifdef Xyce_DEBUG_IO
          cout <<" N_IO_CircuitContext::resolveStrings skipping exception string " << strings[i] << endl;
#endif
          continue;
        }
      }

      // Look for the string in netlistParameters.
      parameterName = strings[i];
      parameterName.toUpper();

      N_UTL_Param expressionParameter(parameterName, "");
      bool parameterFound = getResolvedParameter(expressionParameter);
      if (parameterFound)
      {
#ifdef Xyce_DEBUG_IO
        cout <<" N_IO_CircuitContext::resolveStrings string " << strings[i] << " is a resolved parameter " << expressionParameter.uTag() << " with type "
             << expressionParameter.getType() << " and value ";
        switch (expressionParameter.getType()) {
        case STR:
          cout << expressionParameter.sVal();
          break;
        case DBLE:
          cout << expressionParameter.dVal();
          break;
        case EXPR:
          cout << "EXPR("<<expressionParameter.ePtr()->get_expression()<< ")";
          break;
        default:
          cout << expressionParameter.sVal();
        }
        cout << endl;

#endif
        if ( expressionParameter.getType() == STR ||
             expressionParameter.getType() == DBLE )
        {
          if (!expression.make_constant(strings[i], expressionParameter.dVal()))
          {
            string msg("Problem converting parameter " + parameterName);
            msg += " to its value.\n";
            N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg );
          }
        }
        else if (expressionParameter.getType() == EXPR)
        {
          string expressionString=expression.get_expression();
          if (expression.replace_var(strings[i], *(expressionParameter.ePtr())) != 0)
          {
            string msg("Problem inserting expression "+
                       expressionParameter.ePtr()->get_expression()+
                       " as substitute for " + parameterName +
                       " in expression " + expressionString);
            N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg );
          }
        }
      }
      else
      {
        parameterFound = getResolvedGlobalParameter(expressionParameter);
#ifdef Xyce_DEBUG_IO
	cout << "N_IO_CircuitContext::resolveStrings attempting to resolve "
              <<  " parameter " << expressionParameter.uTag() << endl;
        if (parameterFound)
        {
            cout << "Found it." << endl;
        }
        else
        {
            cout << " Did not find a resolved global parameter named "
                 << expressionParameter.uTag()	<< endl;
        }
#endif
        if (parameterFound)
        {
          if (!expression.make_var(strings[i]))
          {
            string msg("Problem converting parameter " + parameterName);
            msg += " to its value.\n";
            N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg );
          }
        }
        else
        {
          if (N_UTL::isBool(strings[i]))
          {
            bool stat = false;
            if (N_UTL::Bval(strings[i]))
              stat = expression.make_constant(strings[i], static_cast<double>(1));
            else
              stat = expression.make_constant(strings[i], static_cast<double>(0));
            if (!stat)
            {
              string msg("Problem converting parameter " + parameterName);
              msg += " to its value.\n";
              N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg );
            }
          }
          else
          {
            unresolvedStrings = true;
          }
        }
      }
    }
  }

  return !unresolvedStrings;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::resolveFunctions
// Purpose        : Determine if expression has any unresolved functions
//                  and resolve appropriately. Return true if all functions
//                  are resolved otherwise return false.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/11/2003
//----------------------------------------------------------------------------
bool N_IO_CircuitContext::resolveFunctions(N_UTL_Expression & expression)
{
  // Functions in the expression must be previously defined functions
  // that appear in resolvedFunctions_ else there is an error.
  bool unresolvedFunctions = false;
  if ( expression.get_num(XEXP_FUNCTION) > 0)
  {
    // Get the list of strings in the expression.
    vector<string> functions;
    expression.get_names(XEXP_FUNCTION, functions);

    // If the expression is resolvable, each function in the current expression
    // must appear as a defined function in resolvedFunctions_.
    int numFunctions = functions.size();
    for (int i = 0; i < numFunctions; ++i)
    {
      // Look for the function in resolvedFunctions_.
      N_UTL_Param functionParameter(functions[i], "");
      bool functionfound = getResolvedFunction(functionParameter);
      if (functionfound)
      {
        string functionPrototype(functionParameter.tag());
        string functionBody(functionParameter.sVal());

        // The function prototype is defined here as a string whose
        // value is the  function name together with its parenthese
        // enclosed comma separated argument list. To resolve a
        // function, create an expression from the function prototype
        // and get its ordered list of arguments via get_names, then
        // create an expression from the function definition and
        // order its names from that list. Finally, replace the
        // function in the expression to be resolved.
        N_UTL_Expression prototypeExression(functionPrototype);
        vector<string> arguments;
        prototypeExression.get_names(XEXP_STRING, arguments);
        N_UTL_Expression functionExpression(functionBody);
        functionExpression.order_names(arguments);

        if (expression.replace_func(functions[i], functionExpression,
            static_cast<int>(arguments.size())) < 0)
        {
          string msg("Wrong number of arguments for user defined function: ");
          msg += functionPrototype;
          msg += " in expression: ";
          msg += expression.get_expression();
          N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg);
        }

        // Set the expression value.
#ifdef Xyce_DEBUG_IO
        cout << "N_IO_CircuitContext::resolveFunctions: After all expression handling, get_expression returns "
             << expression.get_expression() << endl;
#endif
      }
      else
      {
        unresolvedFunctions = true;
      }
    }
  }

  return !unresolvedFunctions;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::getResolvedParameter
// Purpose        : Look for a parameter with tag parameterName in the current
//                  context's set of resolved parameters. Check the current
//                  context and recurively check parent contexts.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/11/2003
//----------------------------------------------------------------------------
bool N_IO_CircuitContext::getResolvedParameter(N_UTL_Param & parameter)
{
  bool success = false;

  N_UTL_Param* parameterPtr =
    currentContextPtr_->resolvedParams_.findParameter(parameter);
  if (parameterPtr != NULL)
  {
    // Found a parameter with given name, set the value
    // of parameter and return.
    parameter.setVal(*parameterPtr);
    success = true;
  }
  else if (currentContextPtr_->parentContextPtr_ != NULL)
  {
    setContext(currentContextPtr_->parentContextPtr_);
    success = getResolvedParameter(parameter);
    restorePreviousContext();
  }

  return success;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::getResolvedGlobalParameter
// Purpose        : Look for a parameter with tag parameterName in the current
//                  context's set of resolved parameters. Check the current
//                  context and recurively check parent contexts.
// Special Notes  :
// Scope          : public
// Creator        : Dave Shirley, PSSI
// Creation Date  : 11/17/2005
//----------------------------------------------------------------------------
bool N_IO_CircuitContext::getResolvedGlobalParameter(N_UTL_Param & parameter)
{
  bool success = false;

  N_UTL_Param* parameterPtr =
    currentContextPtr_->resolvedGlobalParams_.findParameter(parameter);
  if (parameterPtr != NULL)
  {
    // Found a parameter with given name, set the value
    // of parameter and return.
    parameter.setVal(*parameterPtr);
    success = true;
  }
  else if (currentContextPtr_->parentContextPtr_ != NULL)
  {
    setContext(currentContextPtr_->parentContextPtr_);
    success = getResolvedGlobalParameter(parameter);
    restorePreviousContext();
  }

  return success;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CircuitContext::getResolvedFunction
// Purpose       : Look for a function with functionName in resolvedFunctions_.
//                 Check current context and recursively check each parent
//                 context.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/27/2001
//-----------------------------------------------------------------------------
bool N_IO_CircuitContext::getResolvedFunction(N_UTL_Param & parameter)
{
  bool success = false;

  string functionToFind(parameter.uTag());

  if (currentContextPtr_->resolvedFunctions_.find(functionToFind) !=
      currentContextPtr_->resolvedFunctions_.end())
  {
    parameter = currentContextPtr_->resolvedFunctions_[functionToFind];
    success = true;
  }
  else if (currentContextPtr_->parentContextPtr_ != NULL)
  {
    setContext(currentContextPtr_->parentContextPtr_);
    success = getResolvedFunction(parameter);
    restorePreviousContext();
  }

  return success;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::findModel
// Purpose        : Search the models in the current context for the model of
//                  the given name. If it is not found, recursively
//                  search each parent context. Return a pointer to
//                  the parameterBlock for the model if it is found,
//                  otherwise return NULL. Also, if the model is found,
//                  construct the appropriate model prefix.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/12/2003
//----------------------------------------------------------------------------
bool N_IO_CircuitContext::findModel(string const& modelName,
    N_IO_ParameterBlock* & modelPtr,
    string& modelPrefix)
{
  bool success = false;

  map<string, N_IO_ParameterBlock*>::iterator modelIter
    = currentContextPtr_->models_.find(modelName);
  if (modelIter != currentContextPtr_->models_.end())
  {
    modelPtr = modelIter->second;
    if (modelPtr->hasExpressionValuedParams())
    {
      modelPrefix = currentContextPtr_->getPrefix();
    }
    else
    {
      string prefix = currentContextPtr_->getCurrentContextName();
      if (prefix == "")
        modelPrefix = "";
      else
        modelPrefix = prefix + ":";
    }

    return true;
  }

  // The model was not found in the current circuit's model list,
  // recursively search the parent circuit's model list.
  if ( currentContextPtr_->parentContextPtr_ != NULL )
  {
    setContext(currentContextPtr_->parentContextPtr_);
    success = findModel( modelName, modelPtr, modelPrefix );
    restorePreviousContext();
  }

  return success;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::findModel
// Purpose        : Overloaded version of findModel for cases when the model
//                  prefix is not needed.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/13/2003
//----------------------------------------------------------------------------
bool N_IO_CircuitContext::findModel(
    string const& modelName,
    N_IO_ParameterBlock* & modelPtr)
{
  bool success;
  string temp;
  success = findModel(modelName, modelPtr, temp);

  return success;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::hasSubcircuitParams
// Purpose        : Check whether a subcircuit context is dependent on
//                  subcircuit parameters. These are parameters on the
//                  .subckt line identified by "params:" keyword. The result
//                  should be true if either the current subcircuit context
//                  or any subcircuit context in the hierarchy containing the
//                  current subcircuit context has subcircuit parameters.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 07/15/2003
//----------------------------------------------------------------------------
bool N_IO_CircuitContext::hasSubcircuitParams()
{
  bool foundSubcircuitParams = false;

  if (!currentContextPtr_->subcircuitParameters_.empty())
  {
    foundSubcircuitParams = true;
  }
  else if (currentContextPtr_->parentContextPtr_ != NULL)
  {
    setContext(currentContextPtr_->parentContextPtr_);
    foundSubcircuitParams = hasSubcircuitParams();
    restorePreviousContext();
  }

  return foundSubcircuitParams;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::getTotalDeviceCount
// Purpose        : Calculate the total number of devices starting at current
//                  context and including all subcircuit instances.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/13/2003
//----------------------------------------------------------------------------
int N_IO_CircuitContext::getTotalDeviceCount()
{
  // Get device count for the current context.
  int count = currentContextPtr_->deviceCount_;

  // Determine the device count associated with each subcircuit instance.
  list<string>::iterator instanceIter;
  list<string>::iterator start = currentContextPtr_->instanceList_.begin();
  list<string>::iterator end = currentContextPtr_->instanceList_.end();

  for (instanceIter = start; instanceIter != end; ++instanceIter)
  {
    bool result = setContext(*instanceIter);
    if (!result)
    {
      int lineNumber = currentContextPtr_->
        instanceErrorInfo_[*instanceIter].first;
      string fileName(currentContextPtr_->instanceErrorInfo_[*instanceIter].second);
      string msg("Problem locating subcircuit: ");
      msg += *instanceIter + "\n";
      msg += "Check for missing or misspelled subcircuit name.\n";
      N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_FATAL_0, msg,
          fileName, lineNumber);
    }
    count += getTotalDeviceCount();
    restorePreviousContext();
  }

  return count;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::totalMutualInductanceCount
// Purpose        : Calculate the total number of MIs starting at current
//                  context and including all subcircuit instances.
// Special Notes  :
// Scope          : public
// Creator        : Rob Hoekstra
// Creation Date  : 08/27/2004
//----------------------------------------------------------------------------
int N_IO_CircuitContext::totalMutualInductanceCount()
{
  // Get device count for the current context.
  int count = currentContextPtr_->mutualInductances_.size();

  // Determine the device count associated with each subcircuit instance.
  list<string>::iterator instanceIter;
  list<string>::iterator start = currentContextPtr_->instanceList_.begin();
  list<string>::iterator end = currentContextPtr_->instanceList_.end();

  for (instanceIter = start; instanceIter != end; ++instanceIter)
  {
    bool result = setContext(*instanceIter);
    if (!result)
    {
      int lineNumber = currentContextPtr_->
        instanceErrorInfo_[*instanceIter].first;
      string fileName(currentContextPtr_->instanceErrorInfo_[*instanceIter].second);
      string msg("Problem locating subcircuit: ");
      msg += *instanceIter + "\n";
      msg += "Check for missing or misspelled subcircuit name.\n";
      N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_FATAL_0, msg,
          fileName, lineNumber);
    }
    count += totalMutualInductanceCount();
    restorePreviousContext();
  }

  return count;
}


//-----------------------------------------------------------------------------
// Function      : N_IO_CircuitContext::instance
// Purpose       : implement packable
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Packable * N_IO_CircuitContext::instance() const
{
  return new N_IO_CircuitContext(metadata_, contextList_ , currentContextPtr_);
}


//-----------------------------------------------------------------------------
// Function      : N_IO_CircuitContext::packedByteCount
// Purpose       : Counts bytes needed to pack block.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int N_IO_CircuitContext::packedByteCount() const
{
  int byteCount = 0;
  int size, i;

  // count name
  byteCount += sizeof( int );
  byteCount += name_.length();

  // count device count
  byteCount += sizeof( int );

  // count models
  map< string, N_IO_ParameterBlock* >::const_iterator it_spbM = models_.begin();
  map< string, N_IO_ParameterBlock* >::const_iterator it_speM = models_.end();
  byteCount += sizeof( int );
  for( ; it_spbM != it_speM; ++it_spbM )
  {
    // ---- count key
    byteCount += sizeof( int );
    byteCount += it_spbM->first.length();

    // ---- count data
    byteCount += it_spbM->second->packedByteCount();
  }

  // count unresolved functions
  size = unresolvedFunctions_.size();
  byteCount += sizeof( int );
  for( i = 0; i < size; ++i )
  {
    byteCount += unresolvedFunctions_[ i ].packedByteCount();
  }

  // count instance list
  list< string >::const_iterator it_stbL = instanceList_.begin();
  list< string >::const_iterator it_steL = instanceList_.end();
  byteCount += sizeof( int );
  for( ; it_stbL != it_steL; ++it_stbL )
  {
    byteCount += sizeof( int );
    byteCount += it_stbL->length();
  }

  // count node list
  size = nodeList_.size();
  byteCount += sizeof( int );
  for( i = 0; i < size; ++i )
  {
    byteCount += sizeof( int );
    byteCount += nodeList_[ i ].length();
  }

  // count subcircuit parameters
  size = subcircuitParameters_.size();
  byteCount += sizeof( int );
  for (i = 0; i < size; ++i)
  {
    byteCount += subcircuitParameters_[ i ].packedByteCount();
  }

  // count unresolved params
  size = unresolvedParams_.size();
  byteCount += sizeof( int );
  for (i = 0; i < size; ++i)
  {
    byteCount += unresolvedParams_[ i ].packedByteCount();
  }

  // pack global node names
  set<string>::const_iterator globalNodes_i, globalNodes_end;
  globalNodes_i = globalNodes_.begin();
  globalNodes_end = globalNodes_.end();
  byteCount += sizeof( int );
  for ( ; globalNodes_i != globalNodes_end ; ++globalNodes_i )
  {
    byteCount += sizeof( int );
    byteCount += globalNodes_i->length();
  }

  // count global params
  size = unresolvedGlobalParams_.size();
  byteCount += sizeof( int );
  for (i = 0; i < size; ++i)
  {
    byteCount += unresolvedGlobalParams_[ i ].packedByteCount();
  }

  //count mutual inductances
  size = mutualInductances_.size();
  byteCount += sizeof( int );
  for( i = 0; i < size; ++i )
  {
    byteCount += mutualInductances_[i].packedByteCount();
  }

  // count subcircuit contexts
  size = circuitContextTable_.size();
  byteCount += sizeof( int );
  map< string, N_IO_CircuitContext * >::const_iterator itsc;
  map< string, N_IO_CircuitContext * >::const_iterator itsc_end = circuitContextTable_.end();
  for ( itsc = circuitContextTable_.begin(); itsc != itsc_end; ++itsc )
  {
    byteCount += sizeof( int );
    byteCount += itsc->first.length();
    byteCount += itsc->second->packedByteCount();
  }

  return byteCount;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CircuitContext::pack
// Purpose       : Packs circuit context into char buffer using MPI_PACK.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_IO_CircuitContext::pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const
{
  int size, length, i;
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+packedByteCount();
#endif

  // pack name
  length = name_.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( name_.c_str(), length, buf, bsize, pos );

  // pack device count
  comm->pack( &deviceCount_, 1, buf, bsize, pos );

  // pack models_
  map< string, N_IO_ParameterBlock* >::const_iterator it_spbM = models_.begin();
  map< string, N_IO_ParameterBlock* >::const_iterator it_speM = models_.end();
  size = models_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( ; it_spbM != it_speM; ++it_spbM )
  {
    // ---- pack key
    length = it_spbM->first.length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( it_spbM->first.c_str(), length, buf, bsize, pos );

    // ---- pack data
    it_spbM->second->pack( buf, bsize, pos, comm );
  }

  // pack unresolved functions
  size = unresolvedFunctions_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( i = 0; i < size; ++i )
  {
    unresolvedFunctions_[ i ].pack( buf, bsize, pos, comm );
  }

  // pack instance list
  list< string >::const_iterator it_stbL = instanceList_.begin();
  list< string >::const_iterator it_steL = instanceList_.end();
  size = instanceList_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( ; it_stbL != it_steL; ++it_stbL )
  {
    length = it_stbL->length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( it_stbL->c_str(), length, buf, bsize, pos );
  }

  // pack node list
  size = nodeList_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( i = 0; i < size; ++i )
  {
    length = nodeList_[ i ].length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( nodeList_[ i ].c_str(), length, buf, bsize, pos );
  }

  // pack subcircuit parameters
  size = subcircuitParameters_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( i = 0; i < size; ++i )
  {
    subcircuitParameters_[ i ].pack( buf, bsize, pos, comm );
  }

  // pack unresolved params
  size = unresolvedParams_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( i = 0; i < size; ++i )
  {
    unresolvedParams_[ i ].pack( buf, bsize, pos, comm );
  }

  // pack global node names
  size = globalNodes_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  set<string>::const_iterator globalNodes_i, globalNodes_end;
  globalNodes_i = globalNodes_.begin();
  globalNodes_end = globalNodes_.end();
  for ( ; globalNodes_i != globalNodes_end ; ++globalNodes_i )
  {
    length = globalNodes_i->length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( globalNodes_i->c_str(), length, buf, bsize, pos );
  }

  // pack global params
  size = unresolvedGlobalParams_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( i = 0; i < size; ++i )
  {
    unresolvedGlobalParams_[ i ].pack( buf, bsize, pos, comm );
  }

  // pack mutual inductances
  size = mutualInductances_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( i = 0; i < size; ++i )
  {
    mutualInductances_[ i ].pack( buf, bsize, pos, comm );
  }

  // pack circuitContextTable_
  size = circuitContextTable_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  map< string, N_IO_CircuitContext * >::const_iterator itsc = circuitContextTable_.begin();
  map< string, N_IO_CircuitContext * >::const_iterator itsc_end = circuitContextTable_.end();
  for ( ; itsc != itsc_end; ++itsc )
  {
    length = itsc->first.length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( itsc->first.c_str(), length, buf, bsize, pos );
    itsc->second->pack( buf, bsize, pos, comm );
  }

#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    string msg("Predicted pos does not match actual pos in N_IO_CircuitContext::pack");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING, msg );
  }
#endif
}


//-----------------------------------------------------------------------------
// Function      : N_IO_CircuitContext::unpack
// Purpose       : Unpacks circuit context from char buffer using MPI_UNPACK.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_IO_CircuitContext::unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm )
{
  int size, length, i;

  // unpack name
  comm->unpack( pB, bsize, pos, &length, 1 );
  name_ = string( ( pB + pos ), length );
  pos += length;

  // unpack device count
  comm->unpack( pB, bsize, pos, &deviceCount_, 1 );

  // unpack models_
  comm->unpack( pB, bsize, pos, &size, 1 );

  for( i = 0; i < size; ++i )
  {

    // ---- unpack key
    comm->unpack( pB, bsize, pos, &length, 1 );
    string aString(string( ( pB + pos ), length ));
    pos += length;

    // ---- unpack data
    models_[ aString ] = new N_IO_ParameterBlock();
    models_[ aString ]->unpack( pB, bsize, pos, comm );
  }

  // unpack unresolved functions (vector of N_IO_FunctionBlock)
  comm->unpack( pB, bsize, pos, &size, 1 );
  for( i = 0; i < size; ++i )
  {
    N_IO_FunctionBlock aFuncBlk;
    aFuncBlk.unpack( pB, bsize, pos, comm );
    unresolvedFunctions_.push_back( aFuncBlk );
  }

  // unpack instance list
  comm->unpack( pB, bsize, pos, &size, 1 );
  for( i = 0;  i < size; ++i )
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    instanceList_.push_back( string( ( pB + pos ), length ) );
    pos += length;
  }

  // unpack node list
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (i = 0; i < size; ++i)
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    nodeList_.push_back( string( ( pB + pos ), length ) );
    pos += length;
  }

  // unpack subcircuit parameters
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (i = 0; i < size; ++i)
  {
    N_UTL_Param aParam;
    aParam.unpack( pB, bsize, pos, comm );
    subcircuitParameters_.push_back( aParam );
  }

  // unpack unresovled params vector
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (i = 0; i < size; ++i)
  {
    N_UTL_Param aParam;
    aParam.unpack( pB, bsize, pos, comm );
    unresolvedParams_.push_back( aParam );
  }

  // unpack global node names
  comm->unpack( pB, bsize, pos, &size, 1 );
  for( i = 0;  i < size; ++i )
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    globalNodes_.insert( string( ( pB + pos ), length ) );
    pos += length;
  }

  // unpack unresovled params vector
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (i = 0; i < size; ++i)
  {
    N_UTL_Param aParam;
    aParam.unpack( pB, bsize, pos, comm );
    unresolvedGlobalParams_.push_back( aParam );
  }

  // unpack mutual inductances vector
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (i = 0; i < size; ++i)
  {
    MutualInductance MI;
    MI.unpack( pB, bsize, pos, comm );
    mutualInductances_.push_back( MI );
  }

  // unpack circuitContextTable_
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (i = 0; i < size; ++i)
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    string tmp( ( pB + pos ), length );
    pos += length;
    pair< map< string, N_IO_CircuitContext *>::iterator, bool > p =
     circuitContextTable_.insert(
      pair< string, N_IO_CircuitContext *>(
       tmp,
       new N_IO_CircuitContext( metadata_, contextList_, currentContextPtr_ ) ) );

    // set the parent context of my children to me
    p.first->second->setParentContextPtr( this );
    p.first->second->unpack( pB, bsize, pos, comm );
  }

  currentContextPtr_ = this;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CircuitContext::MutualInductance::MutualInductance
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra
// Creation Date : 8/27/04
//-----------------------------------------------------------------------------
N_IO_CircuitContext::MutualInductance::MutualInductance( N_IO_DeviceBlock & device )
{
    int numParameters = device.getNumberOfInstanceParameters();
    N_DEV_Param parameter;
    bool first = true;
    for ( int i = 0; i < numParameters; ++i )
    {
      parameter = device.getInstanceParameter(i);

      if ( parameter.tag() != "COUPLING" )
      {
        string inductorName (parameter.uTag());

        if( first )
        {
          firstInductor = inductorName;
          first = false;
        }

        inductors[inductorName] = 0.0;
      }
      else
        coupling = parameter.sVal();
    }

    model = device.getModelName();

    name = device.getName();
    sharedKey = 0;

}

//-----------------------------------------------------------------------------
// Function      : N_IO_CircuitContext::MutualInductance::packedByteCount
// Purpose       : Counts bytes needed to pack block.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra
// Creation Date : 8/27/04
//-----------------------------------------------------------------------------
int N_IO_CircuitContext::MutualInductance::packedByteCount() const
{
  int byteCount = 0;

  // coupling value
  byteCount += sizeof(int);
  byteCount += coupling.length();

  // model name
  byteCount += sizeof(int);
  byteCount += model.length();

  // first inductor name
  byteCount += sizeof(int);
  byteCount += firstInductor.length();

  // inductor info
  byteCount += sizeof(int);
  map<string,double>::const_iterator iterI = inductors.begin();
  map<string,double>::const_iterator  endI = inductors.end();
  for( ; iterI != endI; ++iterI )
  {
    byteCount += sizeof(int);
    byteCount += iterI->first.length();
    byteCount += sizeof(double);
  }

  return byteCount;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CircuitContext::MutualInductance::pack
// Purpose       : Packs MI into char buffer using MPI_PACK.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra
// Creation Date : 08/27/04
//-----------------------------------------------------------------------------
void N_IO_CircuitContext::MutualInductance::pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const
{
  int size, length;
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+packedByteCount();
#endif

  // coupling value
  length = coupling.length();
  comm->pack( &length, 1, buf, bsize, pos );
  if( length ) comm->pack( coupling.c_str(), length, buf, bsize, pos );

  // model name
  length = model.length();
  comm->pack( &length, 1, buf, bsize, pos );
  if( length ) comm->pack( model.c_str(), length, buf, bsize, pos );

  // first inductor name
  length = firstInductor.length();
  comm->pack( &length, 1, buf, bsize, pos );
  if( length ) comm->pack( firstInductor.c_str(), length, buf, bsize, pos );

  // pack inductors
  size = inductors.size();
  comm->pack( &size, 1, buf, bsize, pos );
  map<string,double>::const_iterator iterI = inductors.begin();
  map<string,double>::const_iterator  endI = inductors.end();
  for( ; iterI != endI; ++iterI )
  {
    length = iterI->first.length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( iterI->first.c_str(), length, buf, bsize, pos );
    comm->pack( &(iterI->second), 1, buf, bsize, pos );
  }
#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    string msg("Predicted pos does not match actual pos in MutualInductance::pack");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING, msg );
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CircuitContext::MutualInductance::unpack
// Purpose       : Unpacks MI from char buffer using MPI_UNPACK.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra
// Creation Date : 08/27/04
//-----------------------------------------------------------------------------
void N_IO_CircuitContext::MutualInductance::unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm )
{
  int size, length, i;
  double val;

  // coupling value
  comm->unpack( pB, bsize, pos, &length, 1 );
  if( length )
  {
    coupling = string( ( pB + pos ), length );
    pos += length;
  }

  // model name
  comm->unpack( pB, bsize, pos, &length, 1 );
  if( length )
  {
    model = string( ( pB + pos ), length );
    pos += length;
  }

  // first inductor name
  comm->unpack( pB, bsize, pos, &length, 1 );
  if( length )
  {
    firstInductor = string( ( pB + pos ), length );
    pos += length;
  }

  // unpack inductors
  inductors.clear();
  comm->unpack( pB, bsize, pos, &size, 1 );
  for( i = 0; i < size; ++i )
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    string name(string( ( pB + pos ), length ));
    pos += length;
    comm->unpack( pB, bsize, pos, &val, 1 );
    inductors[name] = val;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CircuitContext::bundleMIs
// Purpose       : Convert all MIs into tokenized device lines
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_IO_CircuitContext::bundleMIs()
{
  string type;
  int i, j, mTableSize, mTableRowSize;
  N_IO_SpiceSeparatedFieldTool::StringToken field;
  vector< string > tmpInductorList, tmpCouplingList;
  vector< N_IO_SpiceSeparatedFieldTool::StringToken > tmpLine;

  // retrieve number of K lines in this (sub)circuit
  mTableSize = currentContextPtr_->allIndexedMIs_.size();

  for( i = 0; i < mTableSize; ++i )
  {
    // reset lists, tables and indices
    tmpInductorList.clear();
    tmpCouplingList.clear();

    // reset tmp vars
    string tmpName("");
    string tmpModel("");
    tmpLine.clear();

    vector< set< string > > & tmpTable =
     currentContextPtr_->getSharedInductorTable();

    mTableRowSize = currentContextPtr_->allIndexedMIs_[i].size();

    for( j = 0; j < mTableRowSize; ++j )
    {
      // select current mutual inductance
      MutualInductance & mutind =
       currentContextPtr_->mutualInductances_[
       currentContextPtr_->allIndexedMIs_[i][j]];

      // collect name segments
      if ( tmpName == "" )
        tmpName = mutind.name;
      else
        tmpName += "_" + mutind.name;

      // flag nonlinear coupling
      if( mutind.model != "" )
      {
        tmpModel = mutind.model;
        type = "N";
      }

      // else linear coupling
      else
      {
        type = "L";
      }

      // assemble inductor list
      map< string, double >::iterator mIter = mutind.inductors.begin();
      for( ; mIter != mutind.inductors.end(); ++mIter )
      {
        // add inductor data if not already present
        if( ( type == "N" ) ||
         ( tmpTable[mutind.sharedKey].find( (*mIter).first ) !=
         tmpTable[mutind.sharedKey].end() ) )
        {
          // retrieve the inductor name
          tmpInductorList.push_back( (*mIter).first );

          //Do some error-checking to make sure that the inductor has been
          //defined elsewhere!

          map<string,vector<string> >::iterator inducIter;
          inducIter = mutind.terminals.find((*mIter).first);
          if (inducIter == mutind.terminals.end())
          {
            string msg("Undefined inductor ");
            msg += (*mIter).first;
            msg += " in mutual inductor ";
            msg += mutind.name;
            msg += " definition.";
            N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0,msg);
          }
          else
          {
            // retrieve the inductor nodes
            tmpInductorList.push_back( mutind.terminals[(*mIter).first][0] );
            tmpInductorList.push_back( mutind.terminals[(*mIter).first][1] );
          }

          // retrieve the inductance value
          stringstream cnvtr;
          cnvtr << (*mIter).second;
          tmpInductorList.push_back( cnvtr.str() );

          // remove inductor from list to prevent duplicates
          if( type == "L" )
          {
            tmpTable[mutind.sharedKey].erase( (*mIter).first );
          }

          /*
          // manually check for nonlin/lin inductor overlap & exit
          else
          {
            // FIXME
          }
          */
        }

        // append to coupling list; names are a bit redundant, but easy to track
        tmpCouplingList.push_back( (*mIter).first );
      }

      // append coupling value to string
      stringstream cnvtr;
      cnvtr << mutind.coupling;
      tmpCouplingList.push_back( cnvtr.str() );
    }

    // fully contruct line from components and type
    if( tmpName != "" )
    {
      // append type
      field.string_ = "YMI" + type;
      tmpLine.push_back( field );

      // append name
      field.string_ = tmpName;
      tmpLine.push_back( field );

      // append list of (inductor, terminals, inductance)+
      for( j = 0; j < tmpInductorList.size(); ++j )
      {
        field.string_ = tmpInductorList[j];
        tmpLine.push_back( field );
      }

      // append coupling list
      for( j = 0; j < tmpCouplingList.size(); ++j )
      {
        field.string_ = tmpCouplingList[j];
        tmpLine.push_back( field );
      }

      // append model name
      if( type == "N" )
      {
        field.string_ = tmpModel;
        tmpLine.push_back( field );
      }

      // store in current context
      currentContextPtr_->kLines_.push_back( tmpLine );
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : N_IO_CircuitContext::getMILine
// Purpose       : Retrieve one tokenized device line
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
vector< N_IO_SpiceSeparatedFieldTool::StringToken > &
 N_IO_CircuitContext::getMILine( int i )
{
  if( ( i < 0 ) || ( i > currentContextPtr_->kLines_.size() ) )
  {
    // bounds checking error exit
    string msg("Request exceeds number of mutual induntances ");
    msg += "in this subcircuit.\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL, msg );
  }

  return currentContextPtr_->kLines_[i];
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::augmentTotalDeviceCount
// Purpose        : Augment total device count after we process k-lines
//
// Special Notes  :
// Scope          : public
// Creator        : Keith Santarelli
// Creation Date  : 09/22/08
//----------------------------------------------------------------------------
void N_IO_CircuitContext::augmentTotalDeviceCount(int kLineCount,
                                                  int coupledICount,
                                                  int YDeviceCount)
{
  // Get device count for the current context.
  int count = currentContextPtr_->deviceCount_;
  count += YDeviceCount - kLineCount - coupledICount;

  if (count < 0)
  {
    string msg("Error in N_IO_CircuitContext::augmentTotalDeviceCount:  ");
    msg += "augmented number of devices is less than 0.";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg );
  }
  else
  {
    currentContextPtr_->deviceCount_ = count;
  }
}

