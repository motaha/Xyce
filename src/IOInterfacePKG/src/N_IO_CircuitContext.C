//-------------------------------------------------------------------------
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
// Revision Date   : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <sstream>

#include <N_UTL_fwd.h>
#include <N_IO_Op.h>

#include <N_DEV_DeviceBlock.h>

#include <N_ERH_ErrorMgr.h>

#include <N_IO_CircuitContext.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_ParameterBlock.h>
#include <N_IO_DeviceBlock.h>
#include <N_IO_Op.h>

#include <N_UTL_Expression.h>
#include <N_UTL_Misc.h>
#include <N_UTL_fwd.h>

#include <N_PDS_Comm.h>
#include <N_IO_CircuitMetadata.h>

namespace Xyce {
namespace IO {

//----------------------------------------------------------------------------
// Function       : CircuitContext::CircuitContext
// Purpose        : Constructor
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 01/21/2003
//----------------------------------------------------------------------------
CircuitContext::CircuitContext( CircuitMetadata & md,
                        std::list<CircuitContext*> & cL,
                        CircuitContext *& ccPtr)
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
// Function       : CircuitContext::~CircuitContext
// Purpose        : Destructor
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 01/21/2003
//----------------------------------------------------------------------------
CircuitContext::~CircuitContext()
{
  // Delete each subcircuit context in this context.
  std::map< std::string, CircuitContext * >::iterator itsc = circuitContextTable_.begin();
  std::map< std::string, CircuitContext * >::const_iterator itsc_end = circuitContextTable_.end();

  for ( ; itsc != itsc_end; ++itsc )
  {
    delete itsc->second;
  }

  circuitContextTable_.clear();


  // We delete all our own stored model pointers.  Because we also
  // delete all our subcircuit contexts, *their* destructors take care of
  // deleting their model pointers.
  ModelMap::iterator modelIter = models_.begin();
  ModelMap::iterator modelIterEnd = models_.end();
  for (; modelIter != modelIterEnd; ++modelIter)
    delete modelIter->second;
  models_.clear();
}

//-----------------------------------------------------------------------------
// Function      : CircuitContext::operator=
// Purpose       : assignment operator
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/18/2006
//-----------------------------------------------------------------------------
CircuitContext & CircuitContext ::operator=(const CircuitContext & right)
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
// Function       : CircuitContext::beginSubcircuitContext
// Purpose        : Add a subcircuit context to the current context.
// Special Notes  : This routine sets the current context to a newly
//                  created CircuitContext object. This context
//                  remains the current context until it is explicitly
//                  terminated with a call to endSubcircuitContext.
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 01/24/2003
//----------------------------------------------------------------------------
bool CircuitContext::beginSubcircuitContext(
    std::string const& netlistFileName,
    std::vector<N_IO_SpiceSeparatedFieldTool::StringToken> & subcircuitLine)

{
#ifdef Xyce_DEBUG_IO
  Xyce::dout() << "CircuitContext::beginSubcircuitContext" << std::endl;
  Xyce::dout() << "*******************************************" << std::endl;
  Xyce::dout() << " subcircuit file name: " << netlistFileName << std::endl;
  Xyce::dout() << "  size of line vector: " << subcircuitLine.size() << std::endl;
#endif

  // remove parens
  std::vector<N_IO_SpiceSeparatedFieldTool::StringToken>::iterator iterLine =
   subcircuitLine.begin();
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
  CircuitContext* subcircuitContextPtr =
    new CircuitContext(metadata_, contextList_, currentContextPtr_);

  // Set the parent context, save the current context and reset it to the
  // newly created context.
  subcircuitContextPtr->parentContextPtr_ = currentContextPtr_;
  contextList_.push_front(currentContextPtr_);
  currentContextPtr_ = subcircuitContextPtr;

  // Extract the subcircuit data from subcircuitLine.
  int numFields = subcircuitLine.size();

  if (numFields < 2)
  {
    Report::UserError0().at(netlistFileName, subcircuitLine[0].lineNumber_)
      << "Subcircuit name required";
    return false;
  }

  if (numFields < 3)
  {
    Report::UserError0().at(netlistFileName, subcircuitLine[0].lineNumber_)
      << "At least one node required for subcircuit " << subcircuitLine[1].string_;
    return false;
  }

#ifdef Xyce_DEBUG_IO
  Xyce::dout() << "  Subcircuit name: " << subcircuitLine[1].string_ << std::endl;
  Xyce::dout() << "  Number of fields on subckt line = " << numFields << std::endl;
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
  Util::Param parameter("","");

  bool result = true;
  ++i; // Advance to start of parameters.
  while (i+2 < numFields)
  {
    ExtendedString fieldES = subcircuitLine[i].string_;
    fieldES.toUpper();
    if (!fieldES.possibleParam())
    {
      Report::UserError0().at(netlistFileName, subcircuitLine[i].lineNumber_)
        << "Parameter name " << subcircuitLine[i].string_ <<  " contains illegal character(s)";
      result = false;
    }
    parameter.setTag( fieldES );

    if ( (parameter.uTag() == "TEMP") || (parameter.uTag() == "VT") ||
         (parameter.uTag() == "GMIN") || (parameter.uTag() == "TIME") )
    {
      Report::UserError0().at(netlistFileName, subcircuitLine[i].lineNumber_)
        << "Parameter name " << parameter.uTag() << " not allowed in subcircuit parameter list for subcircuit "
        << getName();
      result = false;
    }

    if ( subcircuitLine[i+1].string_ != "=" )
    {
      Report::UserError0().at(netlistFileName, subcircuitLine[0].lineNumber_)
        << "Equal sign required between parameter and value in PARAM list for subcircuit " << getName();
      result = false;
    }

    else {
      i+=2; // Advance past "=" sign

      fieldES =  subcircuitLine[i].string_;
      fieldES.toUpper();
      parameter.setVal(fieldES);
      subcircuitContextPtr->subcircuitParameters_.push_back(parameter);
#ifdef Xyce_DEBUG_IO
      Xyce::dout() << " Adding parameter " << parameter.uTag() << " with value " << fieldES << std::endl;
#endif
    }
    ++i;
  }

  // check for truncated params: list
  if ( result && i < numFields )
  {
    Report::UserError0().at(netlistFileName, subcircuitLine[0].lineNumber_)
      << "Parameter list error in subcircuit " << getName();
    result = false;
  }

#ifdef Xyce_DEBUG_IO
  Xyce::dout() << "End of CircuitContext::beginSubcircuitContext" << std::endl;
  Xyce::dout() << "*******************************************" << std::endl;
#endif

  return result;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::endSubcircuitContext
// Purpose        : End the current context, push it onto the previous context's
//                  list of contexts, and reset the current context pointer to
//                  the previous context.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 01/24/2003
//----------------------------------------------------------------------------
void CircuitContext::endSubcircuitContext()
{
#ifdef Xyce_DEBUG_IO
  Xyce::dout() << "CircuitContext::endSubcircuitContext" << std::endl;
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
    Report::UserError0() << "Error encountered while traversing subcircuits.";
  }
#ifdef Xyce_DEBUG_IO
  Xyce::dout() << "End of CircuitContext::endSubcircuitContext" << std::endl;
#endif
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::addModel
// Purpose        : Add a model to the current context.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/03/2003
//----------------------------------------------------------------------------
void CircuitContext::addModel(N_IO_ParameterBlock * modelPtr)
{
  N_IO_ParameterBlock* tmpModel;
  if (findModel(modelPtr->getName(), tmpModel))
  {
    Report::UserWarning0 message;
    message << "Reading model named " << modelPtr->getName() << " in the ";

    if (getCurrentContextName() == "")
      message << "main circuit";
    else
      message << "subcircuit " << getCurrentContextName();

    message << " and found one or more models previously defined in this scope";
  }

  currentContextPtr_->models_[modelPtr->getName()] = modelPtr;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::addParams
// Purpose        : Add a set of .PARAM parameters to the current context.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/03/2003
//----------------------------------------------------------------------------
void CircuitContext::addParams(N_IO_OptionBlock const& param)
{
  int numberOfParams = param.getNumberOfParameters();

  Util::Param parameter;
  for (int i = 0; i < numberOfParams; ++i)
  {
    parameter = param.getParameter(i);
    resolveQuote(parameter);
    currentContextPtr_->unresolvedParams_.push_back(parameter);
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::addGlobalParams
// Purpose        : Add a set of .GLOBAL_PARAM parameters to the current context.
// Special Notes  :
// Scope          :
// Creator        : Dave Shirley, PSSI
// Creation Date  : 11/17/2005
//----------------------------------------------------------------------------
void CircuitContext::addGlobalParams(N_IO_OptionBlock const& param)
{
  int numberOfParams = param.getNumberOfParameters();

  Util::Param parameter;
  for (int i = 0; i < numberOfParams; ++i)
  {
    parameter = param.getParameter(i);
    resolveQuote(parameter);
    currentContextPtr_->unresolvedGlobalParams_.push_back(parameter);
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::addGlobalNode
// Purpose        : Add a global node from a .GLOBAL line
// Special Notes  :
// Scope          :
// Creator        : Dave Shirley, PSSI
// Creation Date  : 10/02/2009
//----------------------------------------------------------------------------
void CircuitContext::addGlobalNode( std::string & gnode)
{
  currentContextPtr_->globalNodes_.insert(gnode);
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::resolveQuote
// Purpose        : Resolve quoted parameters as soon as they are encountered.
// Special Notes  : This avoids every processor in a parallel run trying to
//                  access the same file.  Also exit handling does not work on
//                  error for parallel runs otherwise.
// Scope          :
// Creator        : Dave Shirley, PSSI
// Creation Date  :
//----------------------------------------------------------------------------
void CircuitContext::resolveQuote(Util::Param & parameter)
{
  if (parameter.isQuoted())
  {
    // The parameter is time dependent with its time history defined by the set
    // of time-value pairs in the file given by the value of the parameter.
    // Open and read these values, using them to build a "Table" expression
    // for the value of the parameter.
    std::ifstream paramDataIn;
    std::string parameterFile(parameter.stringValue().substr(1,parameter.stringValue().size()-2));
    paramDataIn.open(parameterFile.c_str(), std::ios::in);
    if ( !paramDataIn.is_open() )
    {
      Report::UserError0() << "Could not find file " << parameterFile;
    }

    std::string table("table(time");
    std::string time;
    std::string value;
    while ( paramDataIn >> time )
    {
      if ( paramDataIn >> value )
      {
        table += "," + time + "," + value;
      }
      else
      {
        Report::UserError0() << "Reached end of file in " << parameterFile << " while expecting another value";
      }
    }

    paramDataIn.close();

    if( table.size() <= 10 ) // the length of "table(time" from above
    {
      Report::UserError0() << "Failed to successfully read contents of " << parameterFile;
    }

    table += ")";

    parameter.setVal( Util::Expression(table) );
    return;
  }
}
//----------------------------------------------------------------------------
// Function       : CircuitContext::addFunction
// Purpose        : Add a .FUNC function to the current context.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/03/2003
//----------------------------------------------------------------------------
void CircuitContext::addFunction(N_IO_FunctionBlock const& function)
{
  currentContextPtr_->unresolvedFunctions_.push_back(function);
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::resolve
// Purpose        : Resolve parameters in the current context.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/04/2003
//----------------------------------------------------------------------------
bool CircuitContext::resolve( std::vector<N_DEV_Param> const& subcircuitInstanceParams )
{
  std::vector<N_UTL_Param> retryParams;
  std::vector<N_IO_FunctionBlock> retryFunctions;

  if (currentContextPtr_->subcircuitParameters_.empty() &&
      subcircuitInstanceParams.empty() &&
      currentContextPtr_->resolved_)
  {
    // If there are no subcircuit parameters and this context has already
    // been resolved, then no work to do.
    return true;
  }

#ifdef Xyce_DEBUG_IO
  Xyce::dout() << " CircuitContext::resolve have something to do..." << std::endl;
#endif

  currentContextPtr_->resolved_ = true;

  // Clear resolved params and funcs.
  currentContextPtr_->resolvedParams_.clearParameters();
  currentContextPtr_->resolvedGlobalParams_.clearParameters();
  currentContextPtr_->resolvedFunctions_.clear();


  std::vector<N_UTL_Param> asYetUnresolvedSubcircuitParameters=currentContextPtr_->subcircuitParameters_;
  std::vector<N_UTL_Param> asYetUnresolvedParameters=currentContextPtr_->unresolvedParams_;
  std::vector<N_UTL_Param> asYetUnresolvedGlobalParameters=currentContextPtr_->unresolvedGlobalParams_;
  std::vector<N_IO_FunctionBlock> asYetUnresolvedFunctions=currentContextPtr_->unresolvedFunctions_;

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
    Util::Param parameter;
    std::vector<N_UTL_Param>::iterator paramIter;
    std::vector<N_UTL_Param>::iterator start =
      asYetUnresolvedSubcircuitParameters.begin();
    std::vector<N_UTL_Param>::iterator end =
      asYetUnresolvedSubcircuitParameters.end();
    for (paramIter = start; paramIter != end; ++paramIter)
    {
      Util::Param par = *paramIter;
#ifdef Xyce_DEBUG_IO
      Xyce::dout() << " CircuitContext::resolve attempting to resolve" << par.uTag()<< std::endl;
#endif
      if (par.getType() == Xyce::Util::STR && !par.isNumeric())
      {
        ExtendedString arg ( par.stringValue() );
        arg.toUpper();
        if (arg.possibleParam())
          par.setVal(std::string("{" + arg + "}"));
      }

      if (!resolveParameter(par))
      {
#ifdef Xyce_DEBUG_IO
        Xyce::dout() << "Unable to resolve subcircuit param " << par.uTag() << std::endl;
#endif
        retryParams.push_back(par);
        somethingLeftToDo=true;
      }
      else
      {
#ifdef Xyce_DEBUG_IO
        Xyce::dout() << "resolveParameter returned true on parameter " << par.uTag() << " after resolution its type is " << par.getType() << "with value " ;
        switch (par.getType()) {
        case Xyce::Util::STR:
          Xyce::dout() << par.stringValue();
          break;
        case Xyce::Util::DBLE:
          Xyce::dout() << par.getImmutableValue<double>();
          break;
        case Xyce::Util::EXPR:
          Xyce::dout() << par.getValue<Util::Expression>().get_expression();
          break;
        default:
          Xyce::dout() << par.stringValue();
        }
        Xyce::dout() << std::endl;
#endif

        currentContextPtr_->resolvedParams_.addParameter(par);
        resolvedSomethingThisLoop=true;


      }
#ifdef Xyce_DEBUG_IO
      Xyce::dout() << " CircuitContext::resolve done attempting to resolve" << par.uTag()<< std::endl;
#endif
    }
    asYetUnresolvedSubcircuitParameters=retryParams;
    retryParams.clear();

    Util::Param* paramPtr;
    // Reset the subcircuit parameter values with the instance parameter
    // values as needed.
    int i;
    int numInstanceParameters = subcircuitInstanceParams.size();
    for ( i = 0; i < numInstanceParameters; ++i )
    {
#ifdef Xyce_DEBUG_IO
      Xyce::dout() << " CircuitContext::resolve resetting subcircuit instance parameter" << subcircuitInstanceParams[i].uTag()<< " with value " << std::endl;
        switch (subcircuitInstanceParams[i].getType()) {
        case Xyce::Util::STR:
          Xyce::dout() << subcircuitInstanceParams[i].stringValue();
          break;
        case Xyce::Util::DBLE:
          Xyce::dout() << subcircuitInstanceParams[i].getImmutableValue<double>();
          break;
        case Xyce::Util::EXPR:
          {
            Util::Expression foo(subcircuitInstanceParams[i].getValue<Util::Expression>());
            Xyce::dout() << "EXPR(" << foo.get_expression() << ")";
            break;
          }
        default:
          Xyce::dout() << subcircuitInstanceParams[i].stringValue();
        }
        Xyce::dout() << std::endl;
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
        Xyce::dout() << " did not find in resolvdParams_, adding parameter " << subcircuitInstanceParams[i].uTag() << std::endl;
#endif
      }
      else
      {
#ifdef Xyce_DEBUG_IO
        Xyce::dout() << " found in resolvdParams_, setting parameter " << paramPtr->uTag() << std::endl;
#endif
        paramPtr->setVal(static_cast<const Util::Param &>(subcircuitInstanceParams[i]));
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
      Xyce::dout() << " CircuitContext::resolve Attempting to resolve global parameter " << parameter.uTag();
#endif

      if (!resolveParameter(parameter))
      {
        retryParams.push_back(parameter);
        somethingLeftToDo = true;
      }
      else
      {
        if (parameter.getType() == Xyce::Util::EXPR)
        {
          std::vector<std::string> nodes, instances, leads, variables, specials;

          parameter.getValue<Util::Expression>().get_names(XEXP_NODE, nodes);
          parameter.getValue<Util::Expression>().get_names(XEXP_INSTANCE, instances);
          parameter.getValue<Util::Expression>().get_names(XEXP_LEAD, leads);
          parameter.getValue<Util::Expression>().get_names(XEXP_VARIABLE, variables);
          parameter.getValue<Util::Expression>().get_names(XEXP_SPECIAL, specials);

          if (!nodes.empty() || !instances.empty() || !leads.empty())
          {
            Report::UserError0 message;
            message << "The following are not allowed in global param expression: " << parameter.getValue<Util::Expression>().get_expression();
            
            // This should be caught earlier, but just in case it is checked here
            if (!nodes.empty())
            {
              message << std::endl << "node(s):";
              for (std::vector<std::string>::iterator s_i=nodes.begin() ; s_i!=nodes.end() ; ++s_i)
              {
                message << " " << *s_i;
              }
            }
            if (!instances.empty())
            {
              message << std::endl << "instance(s): ";
              for (std::vector<std::string>::iterator s_i=instances.begin() ; s_i!=instances.end() ; ++s_i)
              {
                message << " " << *s_i;
              }
            }
            if (!leads.empty())
            {
              message << std::endl << "lead(s): ";
              for (std::vector<std::string>::iterator s_i=leads.begin() ; s_i!=leads.end() ; ++s_i)
              {
                message << " " << *s_i;
              }
            }
          }

          if (!variables.empty())
          {
            // If variables are found, they must be previously defined global params
            for (std::vector<std::string>::iterator s_i=variables.begin() ; s_i!=variables.end() ; ++s_i)
            {
              if (currentContextPtr_->resolvedGlobalParams_.findParameter(*s_i) == NULL)
              {
                Report::UserError0() << "Unknown parameter (" << *s_i << ") found in global param expression: "
                                     << parameter.getValue<Util::Expression>().get_expression();
              }
            }
          }
          if (!specials.empty())
          {
            if (specials.size() > 1 || specials[0] != "TIME")
            {
              Report::UserError0 message;
              message << "Unknown special var(s):";
              for (std::vector<std::string>::iterator s_i=specials.begin() ; s_i!=specials.end() ; ++s_i)
              {
                message << " " << *s_i;
              }
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
      // Define the Util::Param for the function prototype (the function name
      // and arguments) and function body. Pass to the circuit for resolution.
      std::string functionName(asYetUnresolvedFunctions[i].functionName);
      std::string functionNameAndArgs(asYetUnresolvedFunctions[i].functionNameAndArgs);
      std::string functionBody(asYetUnresolvedFunctions[i].functionBody);

      std::vector<std::string> functionArgs =
        asYetUnresolvedFunctions[i].functionArgs;
      Util::Param functionParameter(functionNameAndArgs, functionBody);
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
        Util::Expression functionBodyExpression( functionParameter.stringValue() );
        std::vector<std::string> strings;
        bool canResolveAll=true;

        functionBodyExpression.get_names(XEXP_STRING, strings);

        int numStrings = strings.size();
        for (int j = 0; j < numStrings; ++j)
        {
          // Look for string in functionArgs.
          std::vector<std::string>::iterator stringIter =
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
          functionBody = functionParameter.stringValue();

          // Add the function to resolvedFunctions_.
          currentContextPtr_->resolvedFunctions_[functionName] =
            Util::Param(functionNameAndArgs, functionBody);
          resolvedSomethingThisLoop=true;
        }
      }
    }
    asYetUnresolvedFunctions=retryFunctions;
    retryFunctions.clear();
  }

  if (somethingLeftToDo)    // we failed to resolve everything, and the last loop did nothing.
  {
    for (std::vector<N_UTL_Param>::iterator it = asYetUnresolvedSubcircuitParameters.begin(); it != asYetUnresolvedSubcircuitParameters.end(); ++it)
    {
      Report::UserError0() << "Unable to resolve .subckt parameter " << (*it).uTag() << " found in .PARAM statement";
    }

    for (std::vector<N_UTL_Param>::iterator it = asYetUnresolvedParameters.begin(); it != asYetUnresolvedParameters.end(); ++it)
    {
      Report::UserError0() << "Unable to resolve parameter " << (*it).uTag() << " found in .PARAM statement";
    }

    for (std::vector<N_UTL_Param>::iterator it = asYetUnresolvedGlobalParameters.begin(); it != asYetUnresolvedGlobalParameters.end(); ++it)
    {
      Report::UserError0() << "Unable to resolve global parameter " << (*it).uTag() << " found in .PARAM statement";
    }

    for (std::vector<N_IO_FunctionBlock>::iterator func_it = asYetUnresolvedFunctions.begin(); func_it != asYetUnresolvedFunctions.end(); ++func_it)
    {
      const std::string &functionName = (*func_it).functionName;
      const std::string &functionNameAndArgs = (*func_it).functionNameAndArgs;
      const std::string &functionBody = (*func_it).functionBody;
      const std::vector<std::string> &functionArgs = (*func_it).functionArgs;

      Util::Param functionParameter(functionNameAndArgs, functionBody);
      if (!resolveParameter(functionParameter, (*func_it).functionArgs))
      {
        Report::UserError0() << functionNameAndArgs << " contains an undefined parameter or function.";
      }
      else
      {
        std::vector<std::string> strings;

        Util::Expression functionBodyExpression(functionParameter.stringValue());
        functionBodyExpression.get_names(XEXP_STRING, strings);
        for (std::vector<std::string>::const_iterator it = strings.begin(); it != strings.end(); ++it)
        {
          if (find(functionArgs.begin(), functionArgs.end(), *it) == functionArgs.end() && currentContextPtr_->resolvedGlobalParams_.findParameter(*it) == NULL)
          {
            Report::UserError0() << functionNameAndArgs << " contains unknown parameter " << (*it);
          }
        }
      }
    }
  }

  return true;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::setContext
// Purpose        : Set the current context to that corresponding to the given
//                  subcircuit name. Save the previous context on the stack for
//                  later retrieval. If the given subcircuit name is not found,
//                  declare an error and abort.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/07/2003
//----------------------------------------------------------------------------
bool CircuitContext::setContext(
  const std::string &                   subcircuitName,
  const std::string &                   subcircuitPrefixIn,
  const std::list<std::string> &        instanceNodes,
  CircuitContext *                      previousContext ) const
{
  bool success = false;

  std::map< std::string, CircuitContext* >::iterator ccIter =
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
      std::list<std::string>::const_iterator nodeIter = instanceNodes.begin();
      for (int i = 0; i < currentContextPtr_->nodeList_.size() && nodeIter != instanceNodes.end(); ++i, ++nodeIter)
      {
        if (currentContextPtr_->nodeMap_.find(currentContextPtr_->nodeList_[i]) !=
            currentContextPtr_->nodeMap_.end())
        {
          if (currentContextPtr_->nodeMap_[currentContextPtr_->nodeList_[i]] !=
              *nodeIter)
          {
            Report::UserError0()
              <<  "Duplicate nodes in .subckt " << subcircuitName << " point to different nodes in X line invocation";
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
              Report::UserError0()
                << "Global node in subcircuit invocation must match same name in .subckt";
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

  // if (!success)
  // {
  //   Report::UserError0() << "Global node in subcircuit invocation must match same name in .subckt";
  // }

  return success;
}

void CircuitContext::setContext(CircuitContext* context) const
{
  contextList_.push_front(currentContextPtr_);
  currentContextPtr_ = context;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::restorePreviousContext
// Purpose        : Reset the context the context prior to the last invocation
//                  of setContext.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/07/2003
//----------------------------------------------------------------------------
void CircuitContext::restorePreviousContext() const
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
// Function       : CircuitContext::restorePreviousContext
// Purpose        : Reset the context the context prior to the last invocation
//                  of setContext.
// Special Notes  :
// Scope          :
// Creator        : Dave Shirley, PSSI
// Creation Date  : 10/20/2009
//----------------------------------------------------------------------------
bool CircuitContext::globalNode (const std::string &nodeName) const
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
// Function       : CircuitContext::resolveParameter
// Purpose        : Parameter whose value may be an expression that must be
//                  resolved. If the input parameter has an expression value,
//                  replace replace the parameters and functions in the
//                  expression with their actual values.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/10/2003
//----------------------------------------------------------------------------
bool CircuitContext::resolveParameter(Util::Param& parameter,
    std::vector<std::string> exceptionStrings)
{
  if (hasExpressionTag(parameter) || parameter.hasExpressionValue() )
  {
#ifdef Xyce_DEBUG_IO
    Xyce::dout() << "CircuitContext::resolveParameter parameter " << parameter.uTag()
         << " has expression value ";
#endif
    // Extract the expression from the parameter value by stripping off
    // the enclosing braces.  Only strip if it's there!
    std::string expressionString;
    if (parameter.stringValue()[0] == '{')
      expressionString = parameter.stringValue().substr(1, parameter.stringValue().size()-2);
    else
      expressionString = parameter.stringValue();

#ifdef Xyce_DEBUG_IO
    Xyce::dout() << expressionString << std::endl;
#endif

    // Parse the expression:
    Util::Expression expression(expressionString);

    if (!expression.parsed())
      return false;

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
      std::vector<std::string> nodes, instances, leads, variables, specials, nodecomps;
      expression.get_names(XEXP_NODE, nodes);
      expression.get_names(XEXP_INSTANCE, instances);
      expression.get_names(XEXP_LEAD, leads);
      expression.get_names(XEXP_VARIABLE, variables);
      expression.get_names(XEXP_SPECIAL, specials);
      expression.get_names(XEXP_NODAL_COMPUTATION, nodecomps);

      if (!nodes.empty() || !instances.empty() || !leads.empty() ||
          !variables.empty() || !specials.empty() || !nodecomps.empty())
      {
#ifdef Xyce_DEBUG_IO
        Xyce::dout() << "CircuitContext::resolveParameter:  nodes, instances, leads, variables or specials not empty. " << std::endl;
        if (!nodes.empty())
        {
          Xyce::dout() << " Nodes: " << std::endl;
          for (int foo=0; foo<nodes.size(); ++foo)
            Xyce::dout() << foo << " : " << nodes[foo] << std::endl;
        }
        if (!instances.empty())
        {
          Xyce::dout() << " Instances: " << std::endl;
          for (int foo=0; foo<instances.size(); ++foo)
            Xyce::dout() << foo << " : " << instances[foo] << std::endl;
        }
        if (!leads.empty())
        {
          Xyce::dout() << " Leads: " << std::endl;
          for (int foo=0; foo<leads.size(); ++foo)
            Xyce::dout() << foo << " : " << leads[foo] << std::endl;
        }
        if (!variables.empty())
        {
          Xyce::dout() << " Variables: " << std::endl;
          for (int foo=0; foo<variables.size(); ++foo)
            Xyce::dout() << foo << " : " << variables[foo] << std::endl;
        }
        if (!specials.empty())
        {
          Xyce::dout() << " Specials: " << std::endl;
          for (int foo=0; foo<specials.size(); ++foo)
            Xyce::dout() << foo << " : " << specials[foo] << std::endl;
        }
#endif

        parameter.setVal(expression);
#ifdef Xyce_DEBUG_IO
        Xyce::dout() << "CircuitContext::resolveParameter: After all expression handling, get_expression returns "
             << expression.get_expression() << std::endl;
        Xyce::dout() << " after setting the parameter " << parameter.uTag() << ", its type is " << parameter.getType() << std::endl;
        Xyce::dout() << " and its value is ";
        switch (parameter.getType()) {
        case Xyce::Util::STR:
          Xyce::dout() << parameter.stringValue();
          break;
        case Xyce::Util::DBLE:
          Xyce::dout() << parameter.getImmutableValue<double>();
          break;
        case Xyce::Util::EXPR:
          Xyce::dout() << parameter.getValue<Util::Expression>().get_expression();
          break;
        default:
          Xyce::dout() << parameter.stringValue();
        }
        Xyce::dout() << std::endl;
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
          // parameter.addOp(Util::CONSTANT, new IO::ConstantOp(parameter.tag(), value));
#ifdef Xyce_DEBUG_IO
          Xyce::dout() << " CircuitContext::resolveParameter --  Resetting parameter value from " << expressionString << " to " << value << " because exceptionStrings empty." << std::endl;
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
        Xyce::dout() << "CircuitContext::resolveParameter: right before returns "
             << std::endl;
        Xyce::dout() << " after setting the parameter " << parameter.uTag() << ", its type is " << parameter.getType() << std::endl;
        Xyce::dout() << " and its value is ";
        switch (parameter.getType()) {
        case Xyce::Util::STR:
          Xyce::dout() << parameter.stringValue();
          break;
        case Xyce::Util::DBLE:
          Xyce::dout() << parameter.getImmutableValue<double>();
          break;
        case Xyce::Util::EXPR:
          Xyce::dout() << parameter.getValue<Util::Expression>().get_expression();
          break;
        default:
          Xyce::dout() << parameter.stringValue();
        }
        Xyce::dout() << std::endl;
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
// Function       : CircuitContext::resolveStrings
// Purpose        : Determine if expression has any unresolved strings
//                  and resolve appropriately. Return true if all strings are
//                  resolved otherwise return false.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/11/2003
//----------------------------------------------------------------------------
bool CircuitContext::resolveStrings( Util::Expression & expression,
    std::vector<std::string> exceptionStrings)
{
  // Strings in the expression must be previously resolved parameters
  // that appear in paramList else there is an error.
  bool unresolvedStrings = false;
  if ( expression.get_num(XEXP_STRING) > 0 )
  {
    // Get the list of strings in the expression.
    std::vector<std::string> strings;
    expression.get_names(XEXP_STRING, strings);

    // If the expression is resolvable, each string in the current expression
    // must appear as a resolved parameter in netlistParameters. Get the value
    // if it appears there.
    ExtendedString parameterName("");
    int numStrings = strings.size();
    for (int i = 0; i < numStrings; ++i)
    {
#ifdef Xyce_DEBUG_IO
      Xyce::dout() <<" CircuitContext::resolveStrings resolving " << strings[i] << std::endl;
#endif
      // Skip the current string if it is in exception strings. This prevents
      // a function argument from being improperly resolved when there is
      // a parameter in a .param statement with the same name as the function
      // argument.
      if (!exceptionStrings.empty())
      {
        std::vector<std::string>::iterator stringIter = find(exceptionStrings.begin(),
            exceptionStrings.end(),
            strings[i]);
        if (stringIter != exceptionStrings.end())
        {
#ifdef Xyce_DEBUG_IO
          Xyce::dout() <<" CircuitContext::resolveStrings skipping exception string " << strings[i] << std::endl;
#endif
          continue;
        }
      }

      // Look for the string in netlistParameters.
      parameterName = strings[i];
      parameterName.toUpper();

      Util::Param expressionParameter(parameterName, "");
      bool parameterFound = getResolvedParameter(expressionParameter);
      if (parameterFound)
      {
#ifdef Xyce_DEBUG_IO
        Xyce::dout() <<" CircuitContext::resolveStrings string " << strings[i] << " is a resolved parameter " << expressionParameter.uTag() << " with type "
             << expressionParameter.getType() << " and value ";
        switch (expressionParameter.getType()) {
        case Xyce::Util::STR:
          Xyce::dout() << expressionParameter.stringValue();
          break;
        case Xyce::Util::DBLE:
          Xyce::dout() << expressionParameter.getImmutableValue<double>();
          break;
        case Xyce::Util::EXPR:
          Xyce::dout() << "EXPR("<<expressionParameter.getValue<Util::Expression>().get_expression()<< ")";
          break;
        default:
          Xyce::dout() << expressionParameter.stringValue();
        }
        Xyce::dout() << std::endl;

#endif
        if ( expressionParameter.getType() == Xyce::Util::STR ||
             expressionParameter.getType() == Xyce::Util::DBLE )
        {
          if (!expression.make_constant(strings[i], expressionParameter.getImmutableValue<double>()))
          {
            Report::UserWarning0() << "Problem converting parameter " << parameterName << " to its value.";
          }
        }
        else if (expressionParameter.getType() == Xyce::Util::EXPR)
        {
          std::string expressionString=expression.get_expression();
          if (expression.replace_var(strings[i], expressionParameter.getValue<Util::Expression>()) != 0)
          {
            Report::UserWarning0() << "Problem inserting expression " << expressionParameter.getValue<Util::Expression>().get_expression()
                                   << " as substitute for " << parameterName << " in expression " << expressionString;
          }
        }
      }
      else
      {
        parameterFound = getResolvedGlobalParameter(expressionParameter);
#ifdef Xyce_DEBUG_IO
        Xyce::dout() << "CircuitContext::resolveStrings attempting to resolve "
              <<  " parameter " << expressionParameter.uTag() << std::endl;
        if (parameterFound)
        {
            Xyce::dout() << "Found it." << std::endl;
        }
        else
        {
            Xyce::dout() << " Did not find a resolved global parameter named "
                 << expressionParameter.uTag()	<< std::endl;
        }
#endif
        if (parameterFound)
        {
          if (!expression.make_var(strings[i]))
          {
            Report::UserWarning0() << "Problem converting parameter " << parameterName <<" to its value";
          }
        }
        else
        {
          if (Util::isBool(strings[i]))
          {
            bool stat = false;
            if (Util::Bval(strings[i]))
              stat = expression.make_constant(strings[i], static_cast<double>(1));
            else
              stat = expression.make_constant(strings[i], static_cast<double>(0));
            if (!stat)
            {
              Report::UserWarning0() << "Problem converting parameter " << parameterName << " to its value";
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
// Function       : CircuitContext::resolveFunctions
// Purpose        : Determine if expression has any unresolved functions
//                  and resolve appropriately. Return true if all functions
//                  are resolved otherwise return false.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/11/2003
//----------------------------------------------------------------------------
bool CircuitContext::resolveFunctions(Util::Expression & expression)
{
  // Functions in the expression must be previously defined functions
  // that appear in resolvedFunctions_ else there is an error.
  bool unresolvedFunctions = false;
  if ( expression.get_num(XEXP_FUNCTION) > 0)
  {
    // Get the list of strings in the expression.
    std::vector<std::string> functions;
    expression.get_names(XEXP_FUNCTION, functions);

    // If the expression is resolvable, each function in the current expression
    // must appear as a defined function in resolvedFunctions_.
    int numFunctions = functions.size();
    for (int i = 0; i < numFunctions; ++i)
    {
      // Look for the function in resolvedFunctions_.
      Util::Param functionParameter(functions[i], "");
      bool functionfound = getResolvedFunction(functionParameter);
      if (functionfound)
      {
        std::string functionPrototype(functionParameter.tag());
        std::string functionBody(functionParameter.stringValue());

        // The function prototype is defined here as a string whose
        // value is the  function name together with its parenthese
        // enclosed comma separated argument list. To resolve a
        // function, create an expression from the function prototype
        // and get its ordered list of arguments via get_names, then
        // create an expression from the function definition and
        // order its names from that list. Finally, replace the
        // function in the expression to be resolved.
        Util::Expression prototypeExression(functionPrototype);
        std::vector<std::string> arguments;
        prototypeExression.get_names(XEXP_STRING, arguments);
        Util::Expression functionExpression(functionBody);
        functionExpression.order_names(arguments);

        if (expression.replace_func(functions[i], functionExpression,
            static_cast<int>(arguments.size())) < 0)
        {
          Report::UserError0() << "Wrong number of arguments for user defined function " << functionPrototype << " in expression " << expression.get_expression();
        }

        // Set the expression value.
#ifdef Xyce_DEBUG_IO
        Xyce::dout() << "CircuitContext::resolveFunctions: After all expression handling, get_expression returns "
                     << expression.get_expression() << std::endl;
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
// Function       : CircuitContext::getResolvedParameter
// Purpose        : Look for a parameter with tag parameterName in the current
//                  context's set of resolved parameters. Check the current
//                  context and recurively check parent contexts.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/11/2003
//----------------------------------------------------------------------------
bool CircuitContext::getResolvedParameter(Util::Param & parameter)
{
  bool success = false;

  Util::Param* parameterPtr =
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
// Function       : CircuitContext::getResolvedGlobalParameter
// Purpose        : Look for a parameter with tag parameterName in the current
//                  context's set of resolved parameters. Check the current
//                  context and recurively check parent contexts.
// Special Notes  :
// Scope          : public
// Creator        : Dave Shirley, PSSI
// Creation Date  : 11/17/2005
//----------------------------------------------------------------------------
bool CircuitContext::getResolvedGlobalParameter(Util::Param & parameter)
{
  bool success = false;

  Util::Param* parameterPtr =
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
// Function      : CircuitContext::getResolvedFunction
// Purpose       : Look for a function with functionName in resolvedFunctions_.
//                 Check current context and recursively check each parent
//                 context.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/27/2001
//-----------------------------------------------------------------------------
bool CircuitContext::getResolvedFunction(Util::Param & parameter)
{
  bool success = false;

  std::string functionToFind(parameter.uTag());

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
// Function       : CircuitContext::findModel
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
bool CircuitContext::findModel(
  const std::string &           modelName,
  N_IO_ParameterBlock* &        modelPtr,
  std::string &                 modelPrefix) const
{
  bool success = false;

  ModelMap::const_iterator modelIter = currentContextPtr_->models_.find(modelName);
  if (modelIter != currentContextPtr_->models_.end())
  {
    modelPtr = modelIter->second;
    if (modelPtr->hasExpressionValuedParams())
    {
      modelPrefix = currentContextPtr_->getPrefix();
    }
    else
    {
      std::string prefix = currentContextPtr_->getCurrentContextName();
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
// Function       : CircuitContext::findModel
// Purpose        : Overloaded version of findModel for cases when the model
//                  prefix is not needed.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/13/2003
//----------------------------------------------------------------------------
bool CircuitContext::findModel(
  const std::string &           modelName,
  N_IO_ParameterBlock* &        modelPtr) const
{
  bool success;
  std::string temp;
  success = findModel(modelName, modelPtr, temp);

  return success;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::hasSubcircuitParams
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
bool CircuitContext::hasSubcircuitParams()
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
// Function       : CircuitContext::getTotalDeviceCount
// Purpose        : Calculate the total number of devices starting at current
//                  context and including all subcircuit instances.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/13/2003
//----------------------------------------------------------------------------
int CircuitContext::getTotalDeviceCount()
{
  // Get device count for the current context.
  int count = currentContextPtr_->deviceCount_;

  // Determine the device count associated with each subcircuit instance.
  std::list<std::string>::iterator instanceIter;
  std::list<std::string>::iterator start = currentContextPtr_->instanceList_.begin();
  std::list<std::string>::iterator end = currentContextPtr_->instanceList_.end();

  for (instanceIter = start; instanceIter != end; ++instanceIter)
  {
    bool result = setContext(*instanceIter);
    if (result)
    {
      count += getTotalDeviceCount();
    }
    // else
    // {
    //   Report::UserError0().at(currentContextPtr_->instanceErrorInfo_[*instanceIter])
    //     << "Subcircuit " << *instanceIter << " has not been defined";
    // }

    restorePreviousContext();
  }

  return count;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::totalMutualInductanceCount
// Purpose        : Calculate the total number of MIs starting at current
//                  context and including all subcircuit instances.
// Special Notes  :
// Scope          : public
// Creator        : Rob Hoekstra
// Creation Date  : 08/27/2004
//----------------------------------------------------------------------------
int CircuitContext::totalMutualInductanceCount()
{
  // Get device count for the current context.
  int count = currentContextPtr_->mutualInductances_.size();

  // Determine the device count associated with each subcircuit instance.
  std::list<std::string>::iterator instanceIter;
  std::list<std::string>::iterator start = currentContextPtr_->instanceList_.begin();
  std::list<std::string>::iterator end = currentContextPtr_->instanceList_.end();

  for (instanceIter = start; instanceIter != end; ++instanceIter)
  {
    bool result = setContext(*instanceIter);
    if (result)
    {
      count += totalMutualInductanceCount();
    }
    // else
    // {
    //   Report::UserError0().at(currentContextPtr_->instanceErrorInfo_[*instanceIter])
    //     << "Subcircuit " <<  *instanceIter << " has not been defined";
    // }

    restorePreviousContext();
  }

  return count;
}


//-----------------------------------------------------------------------------
// Function      : CircuitContext::instance
// Purpose       : implement packable
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Packable * CircuitContext::instance() const
{
  return new CircuitContext(metadata_, contextList_ , currentContextPtr_);
}


//-----------------------------------------------------------------------------
// Function      : CircuitContext::packedByteCount
// Purpose       : Counts bytes needed to pack block.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int CircuitContext::packedByteCount() const
{
  int byteCount = 0;
  int size, i;

  // count name
  byteCount += sizeof( int );
  byteCount += name_.length();

  // count device count
  byteCount += sizeof( int );

  // count models
  ModelMap::const_iterator it_spbM = models_.begin();
  ModelMap::const_iterator it_speM = models_.end();
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
  std::list< std::string >::const_iterator it_stbL = instanceList_.begin();
  std::list< std::string >::const_iterator it_steL = instanceList_.end();
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
  std::set<std::string>::const_iterator globalNodes_i, globalNodes_end;
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
  std::map< std::string, CircuitContext * >::const_iterator itsc;
  std::map< std::string, CircuitContext * >::const_iterator itsc_end = circuitContextTable_.end();
  for ( itsc = circuitContextTable_.begin(); itsc != itsc_end; ++itsc )
  {
    byteCount += sizeof( int );
    byteCount += itsc->first.length();
    byteCount += itsc->second->packedByteCount();
  }

  return byteCount;
}

//-----------------------------------------------------------------------------
// Function      : CircuitContext::pack
// Purpose       : Packs circuit context into char buffer using MPI_PACK.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void CircuitContext::pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const
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
  ModelMap::const_iterator it_spbM = models_.begin();
  ModelMap::const_iterator it_speM = models_.end();
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
  std::list< std::string >::const_iterator it_stbL = instanceList_.begin();
  std::list< std::string >::const_iterator it_steL = instanceList_.end();
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
  std::set<std::string>::const_iterator globalNodes_i, globalNodes_end;
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
  std::map< std::string, CircuitContext * >::const_iterator itsc = circuitContextTable_.begin();
  std::map< std::string, CircuitContext * >::const_iterator itsc_end = circuitContextTable_.end();
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
    DevelFatal(*this, "CircuitContext::pack") << "Predicted pos does not match actual pos";
  }
#endif
}


//-----------------------------------------------------------------------------
// Function      : CircuitContext::unpack
// Purpose       : Unpacks circuit context from char buffer using MPI_UNPACK.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void CircuitContext::unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm )
{
  int size, length, i;

  // unpack name
  comm->unpack( pB, bsize, pos, &length, 1 );
  name_ = std::string( ( pB + pos ), length );
  pos += length;

  // unpack device count
  comm->unpack( pB, bsize, pos, &deviceCount_, 1 );

  // unpack models_
  comm->unpack( pB, bsize, pos, &size, 1 );

  for( i = 0; i < size; ++i )
  {

    // ---- unpack key
    comm->unpack( pB, bsize, pos, &length, 1 );
    std::string aString(std::string( ( pB + pos ), length ));
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
    instanceList_.push_back( std::string( ( pB + pos ), length ) );
    pos += length;
  }

  // unpack node list
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (i = 0; i < size; ++i)
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    nodeList_.push_back( std::string( ( pB + pos ), length ) );
    pos += length;
  }

  // unpack subcircuit parameters
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (i = 0; i < size; ++i)
  {
    Util::Param aParam;
    aParam.unpack( pB, bsize, pos, comm );
    subcircuitParameters_.push_back( aParam );
  }

  // unpack unresovled params vector
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (i = 0; i < size; ++i)
  {
    Util::Param aParam;
    aParam.unpack( pB, bsize, pos, comm );
    unresolvedParams_.push_back( aParam );
  }

  // unpack global node names
  comm->unpack( pB, bsize, pos, &size, 1 );
  for( i = 0;  i < size; ++i )
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    globalNodes_.insert( std::string( ( pB + pos ), length ) );
    pos += length;
  }

  // unpack unresovled params vector
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (i = 0; i < size; ++i)
  {
    Util::Param aParam;
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
    std::string tmp( ( pB + pos ), length );
    pos += length;
    std::pair< std::map< std::string, CircuitContext *>::iterator, bool > p =
     circuitContextTable_.insert(
      std::pair< std::string, CircuitContext *>(
       tmp,
       new CircuitContext( metadata_, contextList_, currentContextPtr_ ) ) );

    // set the parent context of my children to me
    p.first->second->setParentContextPtr( this );
    p.first->second->unpack( pB, bsize, pos, comm );
  }

  currentContextPtr_ = this;
}

//-----------------------------------------------------------------------------
// Function      : CircuitContext::MutualInductance::MutualInductance
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra
// Creation Date : 8/27/04
//-----------------------------------------------------------------------------
CircuitContext::MutualInductance::MutualInductance( N_IO_DeviceBlock & device )
{
    int numParameters = device.getNumberOfInstanceParameters();
    Device::Param parameter;
    bool first = true;
    for ( int i = 0; i < numParameters; ++i )
    {
      parameter = device.getInstanceParameter(i);

      if ( parameter.tag() != "COUPLING" )
      {
        std::string inductorName (parameter.uTag());

        if( first )
        {
          firstInductor = inductorName;
          first = false;
        }

        inductors[inductorName] = 0.0;
      }
      else
        coupling = parameter.stringValue();
    }

    model = device.getModelName();

    name = device.getName();
    sharedKey = 0;

}

//-----------------------------------------------------------------------------
// Function      : CircuitContext::MutualInductance::packedByteCount
// Purpose       : Counts bytes needed to pack block.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra
// Creation Date : 8/27/04
//-----------------------------------------------------------------------------
int CircuitContext::MutualInductance::packedByteCount() const
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
  std::map<std::string,double>::const_iterator iterI = inductors.begin();
  std::map<std::string,double>::const_iterator  endI = inductors.end();
  for( ; iterI != endI; ++iterI )
  {
    byteCount += sizeof(int);
    byteCount += iterI->first.length();
    byteCount += sizeof(double);
  }

  return byteCount;
}

//-----------------------------------------------------------------------------
// Function      : CircuitContext::MutualInductance::pack
// Purpose       : Packs MI into char buffer using MPI_PACK.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra
// Creation Date : 08/27/04
//-----------------------------------------------------------------------------
void CircuitContext::MutualInductance::pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const
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
  std::map<std::string,double>::const_iterator iterI = inductors.begin();
  std::map<std::string,double>::const_iterator  endI = inductors.end();
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
    DevelFatal(*this, "CircuitContext::MutualInductance::pack") << "Predicted pos does not match actual pos";
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : CircuitContext::MutualInductance::unpack
// Purpose       : Unpacks MI from char buffer using MPI_UNPACK.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra
// Creation Date : 08/27/04
//-----------------------------------------------------------------------------
void CircuitContext::MutualInductance::unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm )
{
  int size, length, i;
  double val;

  // coupling value
  comm->unpack( pB, bsize, pos, &length, 1 );
  if( length )
  {
    coupling = std::string( ( pB + pos ), length );
    pos += length;
  }

  // model name
  comm->unpack( pB, bsize, pos, &length, 1 );
  if( length )
  {
    model = std::string( ( pB + pos ), length );
    pos += length;
  }

  // first inductor name
  comm->unpack( pB, bsize, pos, &length, 1 );
  if( length )
  {
    firstInductor = std::string( ( pB + pos ), length );
    pos += length;
  }

  // unpack inductors
  inductors.clear();
  comm->unpack( pB, bsize, pos, &size, 1 );
  for( i = 0; i < size; ++i )
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    std::string name(std::string( ( pB + pos ), length ));
    pos += length;
    comm->unpack( pB, bsize, pos, &val, 1 );
    inductors[name] = val;
  }
}

//-----------------------------------------------------------------------------
// Function      : CircuitContext::bundleMIs
// Purpose       : Convert all MIs into tokenized device lines
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void CircuitContext::bundleMIs()
{
  std::string type;
  int i, j, mTableSize, mTableRowSize;
  N_IO_SpiceSeparatedFieldTool::StringToken field;
  std::vector< std::string > tmpInductorList, tmpCouplingList;
  std::vector< N_IO_SpiceSeparatedFieldTool::StringToken > tmpLine;

  // retrieve number of K lines in this (sub)circuit
  mTableSize = currentContextPtr_->allIndexedMIs_.size();

  for( i = 0; i < mTableSize; ++i )
  {
    // reset lists, tables and indices
    tmpInductorList.clear();
    tmpCouplingList.clear();

    // reset tmp vars
    std::string tmpName("");
    std::string tmpModel("");
    tmpLine.clear();

    std::vector< std::set< std::string > > & tmpTable =
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
      std::map< std::string, double >::iterator mIter = mutind.inductors.begin();
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

          std::map<std::string,std::vector<std::string> >::iterator inducIter;
          inducIter = mutind.terminals.find((*mIter).first);
          if (inducIter == mutind.terminals.end())
          {
            Report::UserError0() << "Undefined inductor " << (*mIter).first << " in mutual inductor " << mutind.name << " definition.";
          }
          else
          {
            // retrieve the inductor nodes
            tmpInductorList.push_back( mutind.terminals[(*mIter).first][0] );
            tmpInductorList.push_back( mutind.terminals[(*mIter).first][1] );
          }

          // retrieve the inductance value
          std::stringstream cnvtr;
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
      std::stringstream cnvtr;
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
// Function      : CircuitContext::getMILine
// Purpose       : Retrieve one tokenized device line
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
std::vector< N_IO_SpiceSeparatedFieldTool::StringToken > &
 CircuitContext::getMILine( int i )
{
  if( ( i < 0 ) || ( i > currentContextPtr_->kLines_.size() ) )
  {
    // bounds checking error exit
    Report::UserError() << "Request exceeds number of mutual inductances in this subcircuit";
  }

  return currentContextPtr_->kLines_[i];
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::augmentTotalDeviceCount
// Purpose        : Augment total device count after we process k-lines
//
// Special Notes  :
// Scope          : public
// Creator        : Keith Santarelli
// Creation Date  : 09/22/08
//----------------------------------------------------------------------------
void CircuitContext::augmentTotalDeviceCount(int kLineCount,
                                                  int coupledICount,
                                                  int YDeviceCount)
{
  // Get device count for the current context.
  int count = currentContextPtr_->deviceCount_;
  count += YDeviceCount - kLineCount - coupledICount;

  if (count < 0)
  {
    Report::DevelFatal() << "Augmented number of devices is less than 0.";
  }
  else
  {
    currentContextPtr_->deviceCount_ = count;
  }
}

} // namespace IO
} // namespace Xyce
