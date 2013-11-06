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
// Filename      : $RCSfile: N_IO_OutputFileBase.C,v $
// Purpose       : Base class measure functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 12/5/2012
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.3.2.2 $
// Revision Date  : $Date: 2013/10/03 17:23:42 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_IO_OutputFileBase.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Vector.h>

//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::N_IO_OutputFileBase
// Purpose       : Constructor
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------
N_IO_OutputFileBase::N_IO_OutputFileBase() :
  ostreamPtr_(NULL),
  outputFileBaseName_(""),
  fileSuffix_(""),
  simulationSuffix_(""),
  fileFormatName_("FormatUndefined"),
  appendOutputFlag_(false)    
{};


//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::~N_IO_OutputFileBase
// Purpose       : Destructor
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------
N_IO_OutputFileBase::~N_IO_OutputFileBase()
{
  if ( (ostreamPtr_ != &cout) && (ostreamPtr_ != NULL) ) 
  {
    string msg = "N_IO_OutputFileBase destructor called with non-null ostreamPtr_ from " 
      + fileFormatName_ + " derived class.";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg);
  }
  
  if( !(istreamPtr_ != & std::cin) && (istreamPtr_ != 0)  )
  { 
    delete istreamPtr_;
    istreamPtr_ = 0;
  }
}



//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::openFile
// Purpose       : open output file for writing.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------
void N_IO_OutputFileBase::openFile( string basename, string simulationSuffix )
{
  simulationSuffix_ = simulationSuffix;
  outputFileBaseName_ = basename;
  fullFileName_ = outputFileBaseName_ + simulationSuffix_ + fileSuffix_;
  
  if( ostreamPtr_ == 0 )
  {
    // only try to open a new ostream if the pointer is NULL
    if (fullFileName_ == "CONSOLE")
    {
      ostreamPtr_ = &(cout);
    }
    else
    {
      if( appendOutputFlag_ == true ) 
      {
        ostreamPtr_ = new ofstream(fullFileName_.c_str(), ios_base::app );
      }
      else
      {
        ostreamPtr_ = new ofstream(fullFileName_.c_str());
      }
    }
    if( !ostreamPtr_ )
    {
      string msg = "Could not open file, \""  + fullFileName_ + "\" for output.";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, msg);
    }
  }
  else
  {
    string msg = "N_IO_OutputFileBase::openFile(filenmae, suffix) called when ostreamPtr_ was not NULL.";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg);
  }
}


//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::outputHeader
// Purpose       : output the appropriate header info for this file type.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//----------------------------------------------------------------------------- 
void N_IO_OutputFileBase::outputHeader()
{
  string msg = "N_IO_OutputFileBase::outputHeader called in base class from "
    + fileFormatName_ + " derived class.";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg);
}


//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::outputHOMOTOPYHeader
// Purpose       : output the homotopy header info for this file type.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//----------------------------------------------------------------------------- 
void N_IO_OutputFileBase::outputHOMOTOPYHeader( const vector<string> & paramNames)
{
  string msg = "N_IO_OutputFileBase::outputHOMOTOPYHeader called in base class from "
    + fileFormatName_ + " derived class.";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg);
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::outputDC
// Purpose       : output DC results for this file type.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------   
void N_IO_OutputFileBase::outputDC( 
    const int dcNumber, 
    const int maxDC, 
    const vector<N_ANP_SweepParam> & dcParamVec1,
    N_LAS_Vector * solnVecPtr,
    N_LAS_Vector * stateVecPtr,
    N_LAS_Vector * storeVecPtr )
{
  string msg = "N_IO_OutputFileBase::outputDC called when "
    + fileFormatName_ + " format does not support this analysis type.";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg);
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::outputTran
// Purpose       : output transient results for this file type.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------       
void N_IO_OutputFileBase::outputTran(
    const double & time, 
    N_LAS_Vector * solnVecPtr,
    N_LAS_Vector * stateVecPtr,
    N_LAS_Vector * storeVecPtr )
{
  string msg = "N_IO_OutputFileBase::outputTran called when "
    + fileFormatName_ + " format does not support this analysis type.";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg);
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::outputStep
// Purpose       : output Step results for this file type.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------     
void N_IO_OutputFileBase::outputStep(
    const int stepNumber, 
    const int maxStep, 
    const vector<N_ANP_SweepParam> & stepParamVec1,
    N_LAS_Vector * solnVecPtr,
    N_LAS_Vector * stateVecPtr,
    N_LAS_Vector * storeVecPtr )
{
  string msg = "N_IO_OutputFileBase::outputStep called when "
    + fileFormatName_ + " format does not support this analysis type.";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg);
}      


//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::outputAC
// Purpose       : output AC results for this file type.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------       
void N_IO_OutputFileBase::outputAC( 
    const double & freq,
    N_LAS_Vector * freqDomainSolnVecReal,
    N_LAS_Vector * freqDomainSolnVecImaginary)
{
  string msg = "N_IO_OutputFileBase::outputAC called when "
    + fileFormatName_ + " format does not support this analysis type.";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg);
}      


//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::outputHOMOTOPY
// Purpose       : output Homotopy results for this file type.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------       
void N_IO_OutputFileBase::outputHOMOTOPY( 
    const vector<string> & paramNames,
    const vector<double> & paramVals,
    N_LAS_Vector * solnVecPtr )
{
  string msg = "N_IO_OutputFileBase::outputHOMOTOPY called when "
    + fileFormatName_ + " format does not support this analysis type.";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg);
}      


//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::outputMPDE
// Purpose       : output MPDE results for this file type.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------   
void N_IO_OutputFileBase::outputMPDE(const double & time, N_LAS_Vector * solnVecPtr )
{
  string msg = "N_IO_OutputFileBase::outputMPDE called when "
    + fileFormatName_ + " format does not support this analysis type.";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg);
}      

 
//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::outputHB
// Purpose       : output MPDE results for this file type.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------  
void N_IO_OutputFileBase::outputHB(
      const N_LAS_BlockVector & timeDomainSolnVec, 
      const N_LAS_BlockVector & freqDomainSolnVecReal,
      const N_LAS_BlockVector & freqDomainSolnVecImaginary)
{
  string msg = "N_IO_OutputFileBase::outputHB called when "
    + fileFormatName_ + " format does not support this analysis type.";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg);
}  
  
//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::outputMOR()
// Purpose       : output MOR results for this file type.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//----------------------------------------------------------------------------- 
void N_IO_OutputFileBase::outputMOR()
{
  string msg = "N_IO_OutputFileBase::outputMOR called when "
    + fileFormatName_ + " format does not support this analysis type.";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg);
}  

//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::finishOutput
// Purpose       : Output any closing lines.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------   
void N_IO_OutputFileBase::finishOutput()
{
  string msg = "N_IO_OutputFileBase::finishOutput called in base class from "
    + fileFormatName_ + " derived class.";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg);
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::closeOutput
// Purpose       : close the output file.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------   
void N_IO_OutputFileBase::closeOutput()
{
  // Try just closing the output file.  Report any errors.
  if ( ostreamPtr_ != &cout && ostreamPtr_ ) 
  {
    delete ostreamPtr_;
    ostreamPtr_ = 0;
  }
}


//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::openFileForRead
// Purpose       : opens an existing output file for reading.
// Special Notes : returns false if open fails
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/15/13
//----------------------------------------------------------------------------- 
bool N_IO_OutputFileBase::openFileForRead( string filename )
{
  bool retVal=true;
  istreamPtr_ = new ifstream(filename.c_str());
  retVal = istreamPtr_->is_open();
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::getOutputVarNames
// Purpose       : Read the header for far names
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/15/13
//----------------------------------------------------------------------------- 
bool N_IO_OutputFileBase::getOutputVarNames(vector<string> & varNames)
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::getOutputNextVarValues
// Purpose       : Get one row of simulation data
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/15/13
//----------------------------------------------------------------------------- 
bool N_IO_OutputFileBase::getOutputNextVarValues(N_LAS_Vector * varValues)
{
  return false;
}


//-----------------------------------------------------------------------------
// Function      : N_IO_OutputFileBase::closeFileForRead
// Purpose       : close the old output file.
// Special Notes : returns false if close fails
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/15/13
//----------------------------------------------------------------------------- 
bool N_IO_OutputFileBase::closeFileForRead()
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
// Function      : N_IO_OutputFileBase::convertOutputNamesToSolVarNames
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
void N_IO_OutputFileBase::convertOutputNamesToSolVarNames( vector< string > & varNames )
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
        string nodeName;
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

