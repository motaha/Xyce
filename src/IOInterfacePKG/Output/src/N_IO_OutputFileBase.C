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
// Filename      : $RCSfile: N_IO_OutputFileBase.C,v $
// Purpose       : Base class measure functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 12/5/2012
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.10 $
// Revision Date  : $Date: 2014/02/24 23:49:20 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_IO_OutputFileBase.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Vector.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : OutputFileBase::OutputFileBase
// Purpose       : Constructor
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------
OutputFileBase::OutputFileBase() :
  ostreamPtr_(NULL),
  outputFileBaseName_(""),
  fileSuffix_(""),
  simulationSuffix_(""),
  fileFormatName_("FormatUndefined"),
  appendOutputFlag_(false)    
{};


//-----------------------------------------------------------------------------
// Function      : OutputFileBase::~OutputFileBase
// Purpose       : Destructor
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------
OutputFileBase::~OutputFileBase()
{
  if ( (ostreamPtr_ != &std::cout) && (ostreamPtr_ != NULL) ) 
  {
    Report::DevelFatal().in("OutputFileBase::~OutputFileBase()") << "Non-null ostreamPtr_ from " << fileFormatName_ << " derived class.";
  }
  
  if( !(istreamPtr_ != & std::cin) && (istreamPtr_ != 0)  )
  { 
    delete istreamPtr_;
    istreamPtr_ = 0;
  }
}



//-----------------------------------------------------------------------------
// Function      : OutputFileBase::openFile
// Purpose       : open output file for writing.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------
void OutputFileBase::openFile( std::string basename, std::string simulationSuffix )
{
  simulationSuffix_ = simulationSuffix;
  outputFileBaseName_ = basename;
  fullFileName_ = outputFileBaseName_ + simulationSuffix_ + fileSuffix_;
  
  if( ostreamPtr_ == 0 )
  {
    // only try to open a new ostream if the pointer is NULL
    if (fullFileName_ == "CONSOLE")
    {
      ostreamPtr_ = &std::cout;
    }
    else
    {
      if( appendOutputFlag_ == true ) 
      {
        ostreamPtr_ = new std::ofstream(fullFileName_.c_str(), std::ios_base::app );
      }
      else
      {
        ostreamPtr_ = new std::ofstream(fullFileName_.c_str());
      }
    }
    if( !ostreamPtr_ )
    {
      std::string msg = "Could not open file, \""  + fullFileName_ + "\" for output.";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, msg);
    }
  }
  else
  {
    Report::DevelFatal().in("void OutputFileBase::openFile( std::string basename, std::string simulationSuffix )")
      << "ostreamPtr_ was not NULL.";
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputFileBase::closeOutput
// Purpose       : close the output file.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------   
void OutputFileBase::closeOutput()
{
  // Try just closing the output file.  Report any errors.
  if ( ostreamPtr_ != &std::cout && ostreamPtr_ ) 
  {
    delete ostreamPtr_;
    ostreamPtr_ = 0;
  }
}


//-----------------------------------------------------------------------------
// Function      : OutputFileBase::openFileForRead
// Purpose       : opens an existing output file for reading.
// Special Notes : returns false if open fails
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/15/13
//----------------------------------------------------------------------------- 
bool OutputFileBase::openFileForRead( std::string filename )
{
  bool retVal=true;
  istreamPtr_ = new std::ifstream(filename.c_str());
  retVal = istreamPtr_->is_open();
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : OutputFileBase::getOutputVarNames
// Purpose       : Read the header for far names
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/15/13
//----------------------------------------------------------------------------- 
bool OutputFileBase::getOutputVarNames(std::vector<std::string> & varNames)
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : OutputFileBase::getOutputNextVarValues
// Purpose       : Get one row of simulation data
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/15/13
//----------------------------------------------------------------------------- 
bool OutputFileBase::getOutputNextVarValues(N_LAS_Vector * varValues)
{
  return false;
}


//-----------------------------------------------------------------------------
// Function      : OutputFileBase::closeFileForRead
// Purpose       : close the old output file.
// Special Notes : returns false if close fails
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/15/13
//----------------------------------------------------------------------------- 
bool OutputFileBase::closeFileForRead()
{
  bool retVal=true;
  if( istreamPtr_ != 0 )
  {
    istreamPtr_->close();
    retVal = !(istreamPtr_->is_open());
  }
  if( !(istreamPtr_ != & std::cin) && (istreamPtr_ != 0)  )
  { 
    delete istreamPtr_;
    istreamPtr_ = 0;
  }
  return retVal;
}


//-----------------------------------------------------------------------------
// Function      : OutputFileBase::convertOutputNamesToSolVarNames
// Purpose       : Converts names as they appear in an output header to
//                 variable names as they would appear in Xyce's solution 
//                 vector. E.g.: v(a) to "a", I(a) to "a_dev" and Ix(a) to "a_devx"
// Special Notes : Voltage difference v(a,b) and expressions cannot be reduced 
//                 to solution var names so they are left as is.  Also, lead 
//                 currents through devices are of the format deviceName_dev but
//                 currents through voltage sources are deviceName_branch 
//                 This routine ignores that and just tacks on the _dev ending.  
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/17/13
//----------------------------------------------------------------------------- 
void OutputFileBase::convertOutputNamesToSolVarNames( std::vector< std::string > & varNames )
{
  int numNames = varNames.size();
  for( int i=0; i<numNames; i++ )
  {
    // now I need to make the var names from the header match what would normal 
    // be found in the solution/state/store variable map.  
    // Essentially: V(name) becomes name,  I(name) becomes name_dev or name_branch
    // and Ix(name) becomes name_devx.  Expression on the print line 
    // such ax {V(x)-V(y)} can't be reduced further than the expression they are; 
    // a complication we will have to deal with.
    
    if( varNames[i].length() > 1 )
    {
      std::string::size_type begloc = varNames[i].find_first_of('(');  
      std::string::size_type endloc = varNames[i].find_first_of(')');
      if( begloc != std::string::npos &&  endloc != std::string::npos )
      {
        // this takes care of v(a) and i(a)
        // but not v(a,b) as this just becomes a,b -- can't do more with that.
        // stil need to handle i(a) --> a_dev or a_branch
        // and ix(a) becomes a_devx
        std::string nodeName;
        nodeName.assign(varNames[i], begloc+1, endloc-begloc-1);
        if (varNames[i][1] == 'I' || varNames[i][1] == 'i' )
        {
          nodeName.append("_dev");
          // if begloc != 1 then there is an additional lead suffix 
          // on the I(a) as in Ib(a) 
          if( begloc != 1 )
          {
            nodeName.append(varNames[i], 1, 1);
          }
        }
        varNames[i].assign( nodeName );
      }
      else
      {
        // couldn't find '(' and ')' characters as part of V( x ) or I( y )
        // so just leave the name as it is.
      }
    }
    else
    {
      // unexpectedly short name, just leave it as is.
    }
  }

}

} // namespace IO
} // namespace Xyce
