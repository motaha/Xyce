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
// Filename       : $RCSfile: N_IO_ParameterBlock.C,v $
//
// Purpose        : Define the N_IO_ParameterBlock class instantiations of
//                  which are associated with netlist .MODEL lines.
//
// Special Notes  : ERK.  It seems that the name "N_IO_ModelBlock" would have been
//                  more appropriate and less confusing, or possibly
//                  N_IO_ModelParametersBlock.  Calling it "parameters block"
//                  makes it sound much more general than it really is.
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 09/10/2001
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.72.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:43 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>
#include <iostream>
#include <algorithm>

#include <sstream>

// ----------   Xyce Includes   ----------
#include <N_IO_CircuitMetadata.h>
#include <N_IO_ParameterBlock.h>
#include <N_DEV_Param.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Misc.h>
#include <N_UTL_Expression.h>
#include <N_PDS_Comm.h>

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::N_IO_ParameterBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
N_IO_ParameterBlock::N_IO_ParameterBlock(
    string const& fileName,
    vector<N_IO_SpiceSeparatedFieldTool::StringToken> const& parsedInputLine)
: netlistFileName_(fileName),
  defaultApplied_(false),
  parsedLine(parsedInputLine)
{
  setLevel("1");
  lineNumber_ = parsedInputLine[0].lineNumber_;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::print
// Purpose       : Output the details of a device block to standard out.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
void N_IO_ParameterBlock::print()
{
  cout << endl;
  cout << "Parameter Block Information" << endl;
  cout << "---------------------------" << endl;

  cout << "netlist line:" << endl;
  int lineSize = parsedLine.size();
  for ( int i = 0; i < lineSize; ++i )
  {
    cout << "  " << parsedLine[i].string_;
  }
  cout << endl;
  cout << "  name : " << getName() << endl;
  cout << "  type : " << getType() << endl;
  cout << "  level: " << getLevel() << endl;

  cout << "  parameters: " << endl;
  int numParameters = getNumberOfParameters();
  for (int i = 0; i < numParameters; ++i)
  {
    cout << "  " << getParameter(i).uTag() << " : ";
    cout << getParameter(i).sVal();
    if ( getParameter(i).isTimeDependent() )
    {
      cout << "  time dependent";
    }
    cout << endl;
  }
  cout << endl;

  cout << "------------------" << endl;
  cout << endl;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::extractModelData
// Purpose       : Extract model data from parsedLine using model metadata.
//
// Special Notes : ERK - this function doesn't appear to actually use the
//                 metadata object.
//
//                 Unfortunately, this model data gets extracted before
//                 the models defaults get determined.  That makes it much
//                 harder to deal with unusual cases like VECTOR and
//                 VECTOR-COMPOSITE.  Defaults are set up later.
//
//                 Instances are handled in the opposite order.
//
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
bool N_IO_ParameterBlock::extractModelData( N_IO_CircuitMetadata & metadata )
{
  int numFields = parsedLine.size();

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 3 )
  {
    string msg(".model line has too few fields \n");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
        netlistFileName_, parsedLine[0].lineNumber_);
  }

  // Extract the model name and type from parsedLine.
  ExtendedString field ( parsedLine[1].string_ );
  field.toUpper();
  setName( field );
  field = parsedLine[2].string_;
  field.toUpper();
  setType( field );
  setLineNumber(netlistFileName_, parsedLine[0].lineNumber_);

  // Set the start and end position of the model parameters in parsedLine
  // accounting for optional enclosing parentheses.
  int paramStartPos;
  int offsetEq;
  bool levelSet = false;
  vector<N_DEV_Param> inputParameters;
  if ( numFields > 3 )
  {
    int paramEndPos;
    bool beginParenFound = false;

    if ( parsedLine[3].string_ == "(" )
    {
      paramStartPos = 4;
      beginParenFound = true;
    }
    else
    {
      paramStartPos = 3;
    }

    if ( (beginParenFound && parsedLine[numFields - 1].string_ != ")") ||
         (!beginParenFound && parsedLine[numFields - 1].string_ == ")") )
    {
      string msg("Unmatched parentheses for .model " + getName() + "\n");
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
          netlistFileName_, parsedLine[0].lineNumber_);
    }
    else if ( parsedLine[numFields - 1].string_ == ")" )
    {
      paramEndPos = numFields - 2;
    }
    else
    {
      paramEndPos = numFields - 1;
    }

    // flag equal sign usage in parameter list
    if ( ( paramStartPos + 1 >= numFields ) ||
         ( parsedLine[paramStartPos + 1].string_ == "=" ) )
    {
      offsetEq = 0;
    }
    else
    {
      offsetEq = -1;
    }

    // Extract the parameters from parsedLine.
    N_DEV_Param parameter("", "");

    // new way (ERK), which handles traditional name=value format params,
    // VECTOR params, and VECTOR-COMPOSITE params.
    for (int i = paramStartPos; i <= paramEndPos-2-offsetEq; i+=3+offsetEq )
    {
      bool paramDone(false);
      if ( parsedLine[i+1].string_ == "=" && parsedLine[i+2].string_ == "{" )
      {
        // This might be a VECTOR or a VECTOR-COMPOSITE
        //   1. Find the matching final brace.
        //   2. check for "=" signs.  If present between curly braces,
        //   then this is a vector-composite.
        int unmatchedLeftBraces = 1;
        int tmpIndex=i+3;
        int numEqualSigns=0;
        int numCommas=0;
        int finalLeftBraceIndex=i+3;

        for (; tmpIndex <= paramEndPos;++tmpIndex)
        {
          if ( parsedLine[tmpIndex].string_  == "}")
          {
            unmatchedLeftBraces--;
            if (unmatchedLeftBraces==0)
            {
              finalLeftBraceIndex=tmpIndex;
              break;
            }
          }

          if ( parsedLine[tmpIndex].string_  == "{" )
          {
            ++unmatchedLeftBraces;
          }

          if ( parsedLine[tmpIndex].string_ == "=" )
          {
            ++numEqualSigns;
          }

          if ( parsedLine[tmpIndex].string_ == "," )
          {
            ++numCommas;
          }
        }

        if (unmatchedLeftBraces > 0)
        {
          string msg("Unmatched curly braces \n");
          msg += "in .model " + getName() + "\n";
          N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
              netlistFileName_, parsedLine[i+1].lineNumber_);
        }

        bool composite(numEqualSigns > 0);

        if (composite)
        {
#ifdef Xyce_DEBUG_IO
          cout << " N_IO_ParameterBlock::extractModelData: YES! I found a VECTOR-COMPOSITE!!!!" << endl;
#endif
          int linePosition=i+3;
          int savedLinePosition=linePosition;
          int numBlocks = 0;
          int numComponents = 0;
          ExtendedString paramBase ( parsedLine[i].string_ );
          paramBase.toUpper();
/////

// count the number of blocks, get stawman structure
          vector<N_DEV_Param> tmpParamVec;
          while (parsedLine[linePosition].string_ != "}")
          {
            ExtendedString component ( parsedLine[linePosition].string_ );
            component.toUpper();
            ++numComponents;

            parameter.set(component, 0.0);
            tmpParamVec.push_back(parameter);
#ifdef Xyce_DEBUG_IO
            cout << "parameter = " << parameter << endl;
#endif

            linePosition += 2;

            int blockIndex=0;
            bool startVC = true;
            while (startVC || parsedLine[linePosition].string_ == ",")
            {
              if (!startVC) ++linePosition;
              startVC = false;

              ++linePosition;
              ++blockIndex;
            } // end of while loop.

            numBlocks = blockIndex;
          } // end of while loop.
#ifdef Xyce_DEBUG_IO
          cout << "numBlocks = " << numBlocks << endl;
#endif

          inputCompositeParamVecMap[paramBase].resize(numBlocks);
          int iblock=0;
          for (iblock=0;iblock<numBlocks;++iblock)
          {
            inputCompositeParamVecMap[paramBase][iblock] = tmpParamVec;
#ifdef Xyce_DEBUG_IO
            cout << "paramBase = " << paramBase;
            cout << "  iblock = " << iblock << endl;
            for (int ieric=0;ieric<tmpParamVec.size();++ieric)
            {
              cout << "par["<<ieric<<"]="<<tmpParamVec[ieric]<<endl;
            }
#endif
          }
/////
          linePosition = savedLinePosition;

// fill in the allocated data structure
          int componentIndex=-1;
          while (parsedLine[linePosition].string_ != "}")
          {
            ExtendedString component ( parsedLine[linePosition].string_ );
            component.toUpper();
            ++componentIndex;

            // In the instance version of vector-composite processing, there is
            // an error check at this point to see if the specified component
            // exists or not.  We can't do that here, because this function
            // is called prior to metadata being set up.

            linePosition += 2;

            // Mark the component as given in components. Later we will
            // add all of the components and their defaults that were
            // not given.
            //paramIter->setGiven(true);

            int blockIndex=0;
            bool startVC = true;
            while (startVC || parsedLine[linePosition].string_ == ",")
            {
              if (!startVC) ++linePosition;
              startVC = false;

              ExtendedString value ( parsedLine[linePosition].string_ );
#ifdef Xyce_DEBUG_IO
              cout << " N_IO_ParameterBlock::extractModelData: value = " << value << endl;
#endif
              //
              // This line commented  out by TVR on 4 Feb 2009.  It is not
              // clear that it is necessary in any case to upper-case the
              // values of parameters as a matter of course, and it definitely
              // breaks things when the parameter is a file name and the
              // platform has a file system that honors case (e.g. UNIX as
              // opposed to Losedows or Mac).
              // It is probably better to let the code that *uses* a parameter
              // decide whether to force it to upper.
              //
              // value.toUpper();

              inputCompositeParamVecMap[paramBase][blockIndex][componentIndex].setVal(value);
              inputCompositeParamVecMap[paramBase][blockIndex][componentIndex].setGiven( true );

              ++linePosition;
              ++blockIndex;
            } // end of while loop.

            //numComponents = blockIndex;
          } // end of while loop.
/////
          paramDone=true;
          i=finalLeftBraceIndex-2;

#ifdef Xyce_DEBUG_IO

          cout << "User-input .MODEL statement vector-composite for param = "<< paramBase << endl;
          vector <vector<N_DEV_Param> > & tmpParamVecVec = inputCompositeParamVecMap[paramBase];
          for (int ivec=0;ivec<tmpParamVecVec.size();++ivec)
          {

            cout << "Composite for column (block) "<< ivec<< ":" <<endl;
            vector<N_DEV_Param> & tmpVec = tmpParamVecVec[ivec];
            int vecSize = tmpVec.size();

            for (int ip=0; ip<vecSize; ++ip)
            {
              cout << "param["<<ip<<"] = " << tmpVec[ip];
            }
            cout << "-----------" << endl;
          }

#endif

          continue;
        }
        else
        {
#ifdef Xyce_DEBUG_IO
          cout << " N_IO_ParameterBlock::extractModelData: YES! I found a VECTOR!!!!" << endl;
          for (int k = i; k<=finalLeftBraceIndex; k++)
            cout << "\"" << parsedLine[k].string_ << "\" ";
          cout << std::endl;
#endif
          // Extract parameter name and value from parsedLine and add to
          // parameter list, check for LEVEL parameter treat as special case.
          field = parsedLine[i].string_;
          field.toUpper();
          vector<string> values;
          for( int k=i+3; k<finalLeftBraceIndex; k+=2 )
          {
            //std::cout << "adding \"" << parsedLine[k].string_ << "\"" << std::endl;
            values.push_back( parsedLine[k].string_ );
          }
          parameter.set( field, values );
          parameter.setGiven( true );
          inputParameters.push_back( parameter );
          paramDone=true;
          i=finalLeftBraceIndex-2;
          continue;
        }
      }

      // If paramDone=false, then it was not a VECTOR or a VECTOR-COMPOSITE.
      if ( ( parsedLine[i+1].string_ == "=" && !paramDone ) ||
           ( offsetEq == -1 && !paramDone ) )
      {
        // Extract parameter name and value from parsedLine and add to
        // parameter list, check for LEVEL parameter treat as special case.
        field = parsedLine[i].string_;
        field.toUpper();
        parameter.set( field, parsedLine[i+2+offsetEq].string_ );
        if ( parameter.uTag() != "LEVEL" )
        {
          parameter.setGiven( true );
          inputParameters.push_back( parameter );
        }
        else
        {
          setLevel( parameter.sVal() );
          levelSet = true;
        }
        paramDone=true;
      }

      if (!paramDone)
      {
        string msg("Parameter not formatted correctly \n");
        msg += "in .model " + getName() + "\n";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
            netlistFileName_, parsedLine[i+1].lineNumber_);
      }
    }
  }

  if (!levelSet)
  {
    setLevel("1");
  }

#ifdef Xyce_DEBUG_IO
  cout << " N_IO_ParameterBlock::extractModelData.  inputParameters: " << endl;

  //vector<N_DEV_Param> inputParameters;
  int iIPSize = inputParameters.size();
  for (int ieric=0;ieric<iIPSize;++ieric)
  {
    cout << inputParameters[ieric];
  }
#endif

  // For now, accept any parameter.  We will have to wait until the model metadata is
  // generated before checking that the parameters are valid.  This will not happen
  // until pass two.

  // Insert the vector of input parameters into the official parameter vector
  // for this model.
  addParameters(inputParameters);

  // ERK: 12/14/09 note.  at this stage the input parameters were conditionally added
  // to the expressionValuedParms object.  However, in order to include parameters
  // from vector-composites, this process is now done at the end of the
  // addDefaultModelParameters function (at which point vector-composites have been
  // added in to modelData.params object).

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::addDefaultModelParameters
// Purpose       : Add the default parameters for a model from metadata.
// Special Notes : While not required, it is generally expected that
//                 the parameters given in the input netlist will have
//                 already been extracted via extractModelData().
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/17/2001
// Last Modified : 05/28/2002
//-----------------------------------------------------------------------------
void N_IO_ParameterBlock::addDefaultModelParameters(N_IO_CircuitMetadata & metadata )
{
  vector<N_DEV_Param>::iterator paramIter, pparamIter;
  int i(0);
  map<string,bool> paramMetadataExistMap;
  map<string,bool>::iterator pmap_i;

  if (defaultApplied_)
    return;

  // get model metadata
  vector<N_DEV_Param> * modelParameterPtr = metadata.getPtrToModelParameters(getType(), getLevel());

  // check for null ptr (bad model name)
  if( modelParameterPtr == NULL )
  {
    string msg("The model " /* + getName() + " of */"type " + getType() +
     " is not supported in Xyce.\n");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg);
  }

  // Process the metadata.
#ifdef Xyce_DEBUG_IO
  cout << "BEFORE adding final params:"<<endl;
  int itmp3=0;
  for (itmp3=0;itmp3<modelData.params.size();++itmp3)
  {
    cout << modelData.params[itmp3];
  }
  cout << "End of BEFORE parameter List"<<endl;
#endif

  for (i=0; i< modelParameterPtr->size(); ++i)
  {
    N_DEV_Param & param = (*modelParameterPtr)[i];
    string paramSVal(param.sVal());

    paramMetadataExistMap[param.tag()] = true;

    // ERK:  If this is a vector-composite param, we have to do some special stuff.
    // The metadata is stored in a different structure, and the user-specified
    // stuff is in a different structure as well.  Ultimately, however, the
    // vector-composite components has to be converted to the param structure.
    if (paramSVal == "VECTOR-COMPOSITE")
    {
      modelData.params.push_back(param);
      addDefaultCompositeModelParameters_ ( metadata, param, paramMetadataExistMap );
    }

    if (paramSVal == "VECTOR")
    {
      std::cout << "N_IO_ParameterBlock::addDefaultModelParameters() Parameter type is VECTOR" << std::endl;
    }
  }

#ifdef Xyce_DEBUG_IO
  cout << " N_IO_ParameterBlock::addDefaultModelParameters: metadata modelParameterPtr[] (list of defaults):"<<endl;
  for (i=0; i< modelParameterPtr->size(); ++i)
  {
      cout << (*modelParameterPtr)[i];
  }
#endif

  // loop over user-specified params.  At this point the vector composites (if they
  // exist), including defaults, have already been added the params vector.
  int pos = 0;
  map<string,vector<N_DEV_Param>::iterator> paramDefaultProcessed;
  vector<N_DEV_Param>::iterator paramEnd = modelData.params.end();
  for (paramIter=modelData.params.begin() ; paramIter!=paramEnd ; )
  {
    N_DEV_Param & param = (*paramIter);

    if (param.tag() == "INDEPENDENT;PARAM")
    {
      vector<N_DEV_Param> defaults;
      for (pparamIter=modelParameterPtr->begin() ; pparamIter != modelParameterPtr->end() ; ++pparamIter)
      {
        if (paramDefaultProcessed.find((*pparamIter).tag()) == paramDefaultProcessed.end())
        {
          defaults.push_back(*pparamIter);
        }
      }
      paramDefaultProcessed.clear();

      modelData.params.insert(paramIter, defaults.begin(), defaults.end());
      pos += defaults.size()+1;
      paramIter = modelData.params.begin() + pos;
      paramEnd = modelData.params.end();
    }
    else if (paramMetadataExistMap.find(param.tag()) == paramMetadataExistMap.end())
    {
      string msg("No model parameter " + param.tag());
      msg += " found for model " + getName();
      msg += " of type " + getType() + ".\n";
      msg += " This parameter ignored.\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg,
          netlistFileName_, parsedLine[0].lineNumber_);
      paramIter=modelData.params.erase(paramIter);
      paramEnd = modelData.params.end();
    }
    else
    {
      paramDefaultProcessed[(*paramIter).tag()] = paramIter;
      ++paramIter;
      ++pos;
    }
  }

  // Any parameters, that have not been pushed back on the params vector, but
  // which do exist in metadata, should be pushed back at this point.  Vector-composite
  // defaults are already here, but other defaults may not be.
  for (paramIter=modelParameterPtr->begin() ; paramIter != modelParameterPtr->end() ; ++paramIter)
  {
    if (paramDefaultProcessed.find((*paramIter).tag()) == paramDefaultProcessed.end())
    {
      modelData.params.push_back(*paramIter);
    }
  }

  // Add to the vector of expression parameter pointers.
  // ERK. 12/14/09 note:  this procedure used to be done
  // at the end of the extractModelData function.  I moved
  // it here to include vector-composite params.
  int numParameters = modelData.params.size();
  for ( int i = 0; i < numParameters; ++i )
  {
    N_DEV_Param & parameter = modelData.params[i];
    if (parameter.hasExpressionValue() || parameter.isQuoted())
    {
      expressionValuedParams_.push_back(parameter);
    }
  }

#ifdef Xyce_DEBUG_IO
  cout << "AFTER adding final default params:"<<endl;
  int itmp2=0;
  for (itmp2=0;itmp2<modelData.params.size();++itmp2)
  {
    cout << modelData.params[itmp2];
  }
  cout << "End of AFTER parameter List"<<endl;
#endif


  defaultApplied_ = true;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::addDefaultCompositeModelParameters_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 05/15/2008
//-----------------------------------------------------------------------------
void N_IO_ParameterBlock::addDefaultCompositeModelParameters_
  ( N_IO_CircuitMetadata & metadata,
    N_DEV_Param & baseParam ,
    map<string,bool> & paramMetadataExistMap)
{
  int icomp(0), ic(0), iNumCols(0);
  string typ(getType());
  int lev(getLevel());

  // Get the defaultComponents from metadata.
  vector<N_DEV_Param> defaultComponents;
  string baseTag(baseParam.uTag());
  //metadata.getModelCompositeComponents(string(""),
  metadata.getModelCompositeComponents(getType(), baseTag, lev, defaultComponents);

  int iDefaultCompSize=defaultComponents.size();

#ifdef Xyce_DEBUG_IO
  cout << "FOUND A COMPOSITE! ";
  cout << "N_IO_ParameterBlock::addDefaultCompositeModelParameters_: " << baseParam;

  for (icomp=0;icomp<iDefaultCompSize;++icomp)
  {
    cout << "component["<<icomp<<"] : ";
    cout << defaultComponents[icomp];
  }
#endif

  if (inputCompositeParamVecMap.find(baseTag) == inputCompositeParamVecMap.end())
  {
#if 0
    // did the user specify this composite?  If not, then set up a single
    // "column" of the composite.
    inputCompositeParamVecMap[baseTag].resize(1);
    for (icomp=0;icomp<iDefaultCompSize;++icomp)
    {
      N_DEV_Param & paramDefault = defaultComponents[icomp];
      inputCompositeParamVecMap[baseTag][0].push_back(paramDefault);
    }
#else

    // do nothing.  If user did not specify, then don't include it in params.

#endif
  }
  else
  {

#if 0
    // Error checks:  (flesh these out later)


    // double check that the input composite has correct defaultComponents in it.
    iNumCols = inputCompositeParamVecMap[baseTag].size();
    for (ic=0;ic<iNumCols;++ic)
    {

    }

    // double check that all the input composite parameters exist in metadata

    // double check that the input composite is a matrix (each column same size) .

#endif

    // push defaults into the inputCompositeParamVecMap, as needed.
    for (icomp=0;icomp<iDefaultCompSize;++icomp)
    {
      N_DEV_Param & paramDefault = defaultComponents[icomp];
      string utag1(defaultComponents[icomp].uTag());

      bool found=false;
      int itmp=0;
      int iInputCompSize =  inputCompositeParamVecMap[baseTag][0].size();
      for (itmp=0;itmp<iInputCompSize;++itmp)
      {
        string utag2(inputCompositeParamVecMap[baseTag][0][itmp].uTag());
        if (utag1 == utag2)
        {
          found=true;
          break;
        }
      }

      // if not found, then user did not specify.
      // Vector-composite is always specified as a matrix.  So, if the parameter
      // wasn't specified for one column, it wasn't specified for any of them.
      if (!found)
      {
        iNumCols = inputCompositeParamVecMap[baseTag].size();
        for (ic=0;ic<iNumCols;++ic)
        {
          inputCompositeParamVecMap[baseTag][ic].push_back(paramDefault);
        }
      }
    }

    int iInputCompSize =  inputCompositeParamVecMap[baseTag][0].size();
#ifdef Xyce_DEBUG_IO
    cout << "N_IO_ParameterBlock::addDefaultCompositeModelParameters_: input composite for " << baseTag << ": " <<endl;
    iNumCols = inputCompositeParamVecMap[baseTag].size();
    cout << "Column: " ;
    for (ic=0;ic<iNumCols;++ic)
    {
      cout << "\t\t" << ic;
    }
    cout << endl;
    for (icomp=0;icomp<iInputCompSize;++icomp)
    {
      cout.width(14);
      cout << inputCompositeParamVecMap[baseTag][0][icomp].tag();
      for (ic=0;ic<iNumCols;++ic)
      {
        if (icomp >= inputCompositeParamVecMap[baseTag][ic].size())
        {
          string msg("N_IO_ParameterBlock::addDefaultCompositeModelParameters_: inputCompositeParamVecMap overrun!");
          N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg );
        }

        N_DEV_Param & tmpPar = inputCompositeParamVecMap[baseTag][ic][icomp];
        cout << "\t";
        cout.width(20);
        cout << tmpPar.sVal();
        cout << "(given=";
        if (tmpPar.given()) cout << "TRUE) ";
        else                cout << "FALSE)";
      }
      cout << endl;
    }
    cout << endl;
#endif

    // convert from composite components to conventional parameters.
    ostringstream paramName;
    ExtendedString paramBase ( baseTag );
    paramBase.toUpper();
    vector<N_DEV_Param> convertedCompositeParams;

    iInputCompSize =  inputCompositeParamVecMap[baseTag][0].size();
    iNumCols = inputCompositeParamVecMap[baseTag].size();
    for (icomp=0;icomp<iInputCompSize;++icomp)
    {
      for (ic=0;ic<iNumCols;++ic)
      {
        N_DEV_Param newParam = inputCompositeParamVecMap[baseTag][ic][icomp];
        paramName.str("");
        paramName << paramBase << ic << "." << newParam.uTag();
        newParam.setTag(paramName.str());

        convertedCompositeParams.push_back(newParam);
        paramMetadataExistMap[paramName.str()] = true;
      }
    }

    modelData.params.insert
    (modelData.params.end(),
    convertedCompositeParams.begin(),
    convertedCompositeParams.end());

#ifdef Xyce_DEBUG_IO
    cout << "After adding composite params:"<<endl;
    int itmp2=0;
    for (itmp2=0;itmp2<modelData.params.size();++itmp2)
    {
      cout << modelData.params[itmp2];
    }
#endif

  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::setParameterValues
// Purpose       : Look for expression valued parameters in the parameter
//                 list, evaluate expressions found and reset the parameter
//                 value accordingly.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/19/2001
//-----------------------------------------------------------------------------
void N_IO_ParameterBlock::setParameterValues(N_IO_CircuitContext* contextPtr)
{
  N_DEV_Param* parameterPtr;
  int numParameters = expressionValuedParams_.size();
  int i;

  for ( i = 0; i < numParameters; ++i )
  {
    parameterPtr = findParameter(expressionValuedParams_[i]);
#ifdef Xyce_DEBUG_IO
    cout << " setParameterValues, processing expressionValuedParams[" << i << "]" << endl;
    if (parameterPtr == NULL)
    {
      cout << "N_IO_ParameterBlock::setParmaeterValues.  expressionValuedParams_["<<i<<"].uTag = " <<
      expressionValuedParams_[i].uTag() << " is NOT found in model data!" << endl;
      exit(0);
    }

    cout << "   Tag is " << parameterPtr->uTag() << endl;
    cout << "   value is " << parameterPtr->sVal() << endl;
#endif
    if (parameterPtr != NULL)
    {
      if (!contextPtr->resolveParameter(*parameterPtr))
      {
        N_UTL_Expression & expr = *(parameterPtr->ePtr());

        string msg("Parameter " + parameterPtr->uTag() + " for model ");
        msg += getName();
        if (expr.get_num(XEXP_LEAD) > 0)
        {
          msg += " contains illegal use of lead current: ";
        }
        else
        {
          msg += " contains unrecognized symbols: \n";
        }
        msg += expr.get_expression();
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg,
            netlistFileName_, parsedLine[0].lineNumber_);
      }
    }
#ifdef Xyce_DEBUG_IO
    cout << " after resolution, "  << endl;
    cout << "   Tag is " << parameterPtr->uTag() << endl;
    if (parameterPtr->hasExpressionValue())
    {
      cout << " Expression did not fully resolve, still has value "
           << parameterPtr->ePtr()->get_expression() << endl;
    }
    else
    {
      cout << "   value is " << parameterPtr->dVal() << endl;
    }
#endif
  }

  vector<N_DEV_Param> & modelDataParams = modelData.params;
  int modelDataSize = modelDataParams.size();
  for (int iparam=0;iparam<modelDataSize;++iparam)
  {
    N_DEV_Param & param = modelDataParams[iparam];

    if (param.getType() == STR && !param.isNumeric())
    {
      ExtendedString paramNameOrig(param.sVal());

      // since we're going to upcase paramNameOrig to try resolving it as an
      // expression, save the original just in case we *can't* resolve it
      // and need to restore it.  This is necessary because the string value
      //  might later be used as a file name, for example, and upcasing the
      // value is Just Wrong.
      ExtendedString paramNameSave(paramNameOrig);

      paramNameOrig.toUpper();
      if (paramNameOrig.possibleParam())
      {
        param.setVal(string("{" + paramNameOrig + "}"));

        // try to resolve this string as a simple parameter.  If it can't
        // resolve, restore it to its original value (it's real original
        // value, not its upcased value)
        if (!contextPtr->resolveParameter(param))
          param.setVal(paramNameSave);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::findParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
N_DEV_Param* N_IO_ParameterBlock::findParameter( N_DEV_Param const& parameter )
{
  vector<N_DEV_Param>::iterator paramIter;
  paramIter = find( modelData.params.begin(),
                    modelData.params.end(),
                    parameter );
  if ( paramIter != modelData.params.end() )
  {
    return &(*paramIter);
  }
  else
  {
    return NULL;
  }
}



//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::instance
// Purpose       : implement packable
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Packable * N_IO_ParameterBlock::instance() const
{
  return new N_IO_ParameterBlock();
}


//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::packedByteCount
// Purpose       : Counts bytes needed to pack block.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int N_IO_ParameterBlock::packedByteCount() const
{
  int byteCount = 0;
  int size, length, j;

  // count pardedLine
  size = parsedLine.size();
  byteCount += sizeof( int );
  for( j = 0; j < size; ++j)
  {
    byteCount += parsedLine[ j ].packedByteCount();
  }

  // count modelData
  byteCount += modelData.packedByteCount();

  // count netlistfileName_
  byteCount += sizeof( int );
  byteCount += netlistFileName_.length();

  // count lineNumber_
  byteCount += sizeof( int );

  // count defaultApplied_
  byteCount += sizeof( int );

  // count expressionValuedParams_
  size = expressionValuedParams_.size();
  byteCount += sizeof( int );
  for( j = 0; j < size; ++j )
  {
    byteCount += expressionValuedParams_[ j ].packedByteCount();
  }

  // count the input composite size
  byteCount += sizeof( int ); // size

  // count the internals of the composite, if not empty.
  if ( !(inputCompositeParamVecMap.empty()) )
  {
    map < string, vector <vector<N_DEV_Param> > >::const_iterator  iter;
    map < string, vector <vector<N_DEV_Param> > >::const_iterator  begin = inputCompositeParamVecMap.begin();
    map < string, vector <vector<N_DEV_Param> > >::const_iterator  end = inputCompositeParamVecMap.end();
    for (iter=begin;iter!=end;++iter)
    {
      // count paramBase
      string paramBase = iter->first;
      byteCount += sizeof( int ); // size
      byteCount += paramBase.length();

      // count size of paramVecVec.
      const vector <vector<N_DEV_Param> > & tmpParamVecVec = iter->second;
      int sizeTmpParamVec = tmpParamVecVec.size();
      byteCount += sizeof( int ); // size

      for (int ivec=0;ivec< sizeTmpParamVec; ++ivec)
      {
        const vector<N_DEV_Param> & tmpVec = tmpParamVecVec[ivec];
        int vecSize = tmpVec.size();
        byteCount += sizeof( int ); // size

        for (int ip=0; ip<vecSize; ++ip)
        {
          byteCount += tmpVec[ip].packedByteCount();
        }
      }
    }
  }

  return byteCount;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::pack
// Purpose       : Packs parameter block into char buffer using MPI_PACK.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_IO_ParameterBlock::pack( char * buf, int bsize, int & pos,
 N_PDS_Comm * comm ) const
{
  int size, length, j;
  int def;
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+packedByteCount();
#endif

  // pack parsedLine
  size = parsedLine.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( j = 0; j < size; ++j)
  {
    parsedLine[ j ].pack( buf, bsize, pos, comm );
  }

  // pack modelData;
  modelData.pack( buf, bsize, pos, comm );

  // pack netlistFileName_
  length = netlistFileName_.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( netlistFileName_.c_str(), length, buf, bsize, pos );

  // pack line number
  comm->pack( &lineNumber_, 1, buf,  bsize, pos );

  if (defaultApplied_)
    def = 1;
  else
    def = 0;
  comm->pack( &def, 1, buf,  bsize, pos );

  // pack expressionValuedParams_
  size = expressionValuedParams_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( j = 0; j < size; ++j )
  {
    expressionValuedParams_[ j ].pack( buf, bsize, pos, comm );
  }

  // pack inputCompositeParamVecMap
  int sizeComposite=inputCompositeParamVecMap.size();
  comm->pack( &sizeComposite, 1, buf, bsize, pos );

  if ( !(inputCompositeParamVecMap.empty()) )
  {
    map < string, vector <vector<N_DEV_Param> > >::const_iterator  iter;
    map < string, vector <vector<N_DEV_Param> > >::const_iterator  begin = inputCompositeParamVecMap.begin();
    map < string, vector <vector<N_DEV_Param> > >::const_iterator  end = inputCompositeParamVecMap.end();
    for (iter=begin;iter!=end;++iter)
    {
      // pack paramBase name.
      string paramBase = iter->first;
      int baseLength = paramBase.length();
      comm->pack( &baseLength, 1, buf, bsize, pos );
      comm->pack( paramBase.c_str(), baseLength, buf, bsize, pos );

      // pack size of tmpParamVecVec
      const vector <vector<N_DEV_Param> > & tmpParamVecVec = iter->second;
      int sizeTmpParamVec = tmpParamVecVec.size();
      comm->pack( &sizeTmpParamVec, 1, buf, bsize, pos );

      // pack contents of tmpParamVecVec
      for (int ivec=0;ivec< sizeTmpParamVec; ++ivec)
      {
        const vector<N_DEV_Param> & tmpVec = tmpParamVecVec[ivec];
        int vecSize = tmpVec.size();
        comm->pack( &vecSize, 1, buf, bsize, pos );

        for (int ip=0; ip<vecSize; ++ip)
        {
          tmpVec[ip].pack(buf, bsize, pos, comm);
        }
      }
    }
  }

#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    string msg("Predicted pos does not match actual pos in N_IO_ParameterBlock::pack");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING, msg );
  }
#endif
}


//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::unpack
// Purpose       : Unpacks parameter block from char buffer using MPI_UNPACK.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_IO_ParameterBlock::unpack( char * pB, int bsize, int & pos,
 N_PDS_Comm * comm)
{
  int size, length, j;
  int def;

  // unpack parsedLine
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (j = 0; j < size; ++j)
  {
    N_IO_SpiceSeparatedFieldTool::StringToken aStringToken;
    aStringToken.unpack( pB, bsize, pos, comm );
    parsedLine.push_back( aStringToken );
  }

  // unpack optionData
  modelData.unpack( pB, bsize, pos, comm );

  // unpack netlistFileName_
  comm->unpack( pB, bsize, pos, &length, 1 );
  netlistFileName_ = string( ( pB + pos ), length );
  pos += length;

  comm->unpack( pB, bsize, pos, &lineNumber_, 1 );

  comm->unpack( pB, bsize, pos, &def, 1 );
  if (def != 0)
    defaultApplied_=true;
  else
    defaultApplied_=false;

  // unpack expressionValuedParams_
  comm->unpack( pB, bsize, pos, &size, 1 );
  for( j = 0; j < size; ++j )
  {
    expressionValuedParams_.push_back( N_DEV_Param() );
    expressionValuedParams_[ j ].unpack( pB, bsize, pos, comm );
  }


//#ifdef Xyce_FIX_COMPOSITE_PACK
  // unpack inputCompositeParamVecMap
  int sizeComposite=0;
  comm->unpack( pB, bsize, pos, &size, 1 );
  sizeComposite=size;

  if ( sizeComposite>0 )
  {
    for (int ic=0;ic<sizeComposite;++ic)
    {
      // unpack paramBase name.
      comm->unpack( pB, bsize, pos, &length, 1 );
      string paramBase = string( ( pB + pos ), length );
      pos += length;

      // unpack size of paramVecVec
      int sizeTmpParamVec(0);
      comm->unpack( pB, bsize, pos, &sizeTmpParamVec, 1 );
      vector <vector<N_DEV_Param> > tmpParamVecVec;
      tmpParamVecVec.resize(sizeTmpParamVec);

      for (int ivec=0;ivec< sizeTmpParamVec; ++ivec)
      {
        // unpack size of paramVec
        comm->unpack( pB, bsize, pos, &length, 1 );
        int vecSize = length;
        tmpParamVecVec[ivec].resize(vecSize);

        for (int ip=0; ip<vecSize; ++ip)
        {
          tmpParamVecVec[ivec][ ip ].unpack( pB, bsize, pos, comm );
        }
      }
      inputCompositeParamVecMap[paramBase] = tmpParamVecVec;
    }
  }
//#endif

}

