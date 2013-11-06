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
// Filename      : $RCSfile: N_IO_OutputPrn.C,v $
// Purpose       : Base class measure functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 12/06/2012
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.2.2.2 $
// Revision Date  : $Date: 2013/10/03 17:23:42 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include<N_IO_OutputPrn.h>
#include <N_ERH_ErrorMgr.h>

//-----------------------------------------------------------------------------
// Function      : N_IO_OutputPrn::N_IO_OutputPrn() 
// Purpose       : Constructor for PRN file IO
// Special Notes :
// Creator       : Rich Schiek
// Creation Date : 4/16/13
//-----------------------------------------------------------------------------
N_IO_OutputPrn::N_IO_OutputPrn() 
{
}


//-----------------------------------------------------------------------------
// Function      : N_IO_OutputPrn::~N_IO_OutputPrn() 
// Purpose       : Destructor for PRN file IO
// Special Notes :
// Creator       : Rich Schiek
// Creation Date : 4/16/13
//-----------------------------------------------------------------------------
N_IO_OutputPrn::~N_IO_OutputPrn() 
{
}


// Note The following two functions should be made more general 
// to extract strings given some set of deliminators chars 
// and the same for doubles.  Then the could be moved to the base class
// and re-used for other text based formats.

//-----------------------------------------------------------------------------
// Function      : N_IO_OutputPrn::openFileForRead() 
// Purpose       : Get header line of PRN file for var names.
// Special Notes : return false if we cannot read from file.
// Creator       : Rich Schiek
// Creation Date : 4/16/13
//-----------------------------------------------------------------------------
bool N_IO_OutputPrn::getOutputVarNames( vector< string > & varNames )
{
  bool retVal = true;
  stringstream extractedVarName;
  bool doneWithReadLine = false;
  bool withinWord=false;
  while( !doneWithReadLine )
  {
    char characterRead=0;
    istreamPtr_->get( characterRead );
    if( characterRead == '\n' || characterRead == '\r' )
    {
      doneWithReadLine=true;
    }
    if( characterRead != ' ' && characterRead != '\t' && 
        characterRead != '\r' && characterRead != '\n' )
    {
      withinWord = true;
      extractedVarName.put(characterRead);
    }
    else
    {
      if( withinWord )
      {
        // just found first white space after word so 
        // transfer extractedVarName to varNames array 
        string name;
        extractedVarName >> name;
        varNames.push_back( name );
        extractedVarName.clear();
        withinWord=false;
      }
    }
  }
  
  if( varNames.size() == 0 )
  {
    // nothing read so return false 
    retVal = false;
  }  
  return retVal;
}


//-----------------------------------------------------------------------------
// Function      : N_IO_OutputPrn::openFileForRead() 
// Purpose       : Get a line of data from PRN file.  
// Special Notes : Return false if we hit the end of the file.
// Creator       : Rich Schiek
// Creation Date : 4/16/13
//-----------------------------------------------------------------------------
bool N_IO_OutputPrn::getOutputNextVarValues( N_LAS_Vector * varValues )
{
  bool retVal = true;
  stringstream extractedVarValue;
  bool doneWithReadLine = false;
  bool withinWord=false;
  const string validNumberChars="0123456789";
  // assume we start filling the array at index 0.
  int varIndex = 0;
  while( !doneWithReadLine )
  {
    char characterRead=0;
    istreamPtr_->get( characterRead );
    if( istreamPtr_->eof() )
    {
      // hit end of file. return false 
      doneWithReadLine = true;
      retVal = false;
    }
    if( characterRead == '\n' || characterRead == '\r' )
    {
      doneWithReadLine=true;
    }
    if( characterRead != ' ' && characterRead != '\t' && 
        characterRead != '\r' && characterRead != '\n' )
    {
      // need to ensure that characters read are valid numbers 
      // 0123456789
      if( withinWord || (validNumberChars.find( characterRead ) != string::npos ) )
      {
        withinWord = true;
        extractedVarValue.put(characterRead);
      }
    }
    else
    {
      if( withinWord )
      {
        // just found first white space after word so 
        // transfer extractedVarName to varValues array 
        double value;
        extractedVarValue >> value;
        (*varValues)[varIndex] = value ;
        varIndex++;
        extractedVarValue.clear();
        withinWord=false;
      }
    }
  }
  
  if( varIndex == 0 )
  {
    // nothing read so return false 
    retVal = false;
  }  
  return retVal;
}


