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
// Filename      : $RCSfile: N_IO_MeasureRelativeError.C,v $
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.13.2.1 $
// Revision Date  : $Date: 2014/03/06 19:22:47 $
// Current Owner  : $Author: rlschie $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_MeasureRelativeError.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Misc.h>
#include <Epetra_SerialDenseVector.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : RelativeError::RelativeError()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
RelativeError::RelativeError( const Util::OptionBlock & measureBlock, IO::OutputMgr &outputMgr ):
  Base(measureBlock, outputMgr)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;
  if( comparisonFunctionName_.empty() )
    comparisonFunctionName_ = "L2NORM";

  if( !dataFileName_.empty() ) 
  {
    // load the data file
    std::ifstream * dataFileStream( openStreamFromFileName( dataFileName_ ) );
    bool result = readData( *dataFileStream, varNames, indepVarValues, indep2VarValues, dataValues );
    if( !result ) 
    {
      Report::UserError0() << "Could not open file in the measure function \"" << name_ << "\". Failed filename was = \"" << dataFileName_ << "\"";
    }
    closeStream( dataFileStream  );
    indepVarValues.resize( dataValues.size() );
    // copy over data into column vector for independent var
   
    // Xyce::dout() << " Extracting independent and dependent vars  from columns : " << independentVarColumn_ << ", " << dependentVarColumn_ << std::endl;
    if( independentVarColumn_ > dataValues[0].size() )
    {
      Report::UserError0() << "In measure function \"" << name_ 
        << "\". Requested column for independent variable, " 
	<< independentVarColumn_ << ", did not exist in the data file (last column number was)" << dataValues[0].size();
    }

    for(int i=0; i< dataValues.size(); i++ )
    {
      indepVarValues[i] = dataValues[i][independentVarColumn_];
      //Xyce::dout() << indepVarValues[i] << std::endl;
    }
 
    // results from the simulation to compare to a column in dataValues; 
    simulationDataVals_.resize( dataValues.size() );
    simulationDataValsFound_.resize( dataValues.size() );
    for( int i=0; i< dataValues.size(); i++ )
    {
      simulationDataVals_[i] = 0.0;
      simulationDataValsFound_[i] = 0;
    }
    lastIndepIndex_=0;
 
  }
}

void RelativeError::prepareOutputVariables()
{
  // this measurement should have only one dependent variable.
  // Error for now if it doesn't
  numOutVars_ = outputVars_.size();
  
  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for relative error measure, \"" + name_ + "\" Exiting.";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL, msg);
  }

  outVarValues_.resize( numOutVars_, 0.0 );
  lastDepVar_ = 0;
  lastOutVarValues_.resize( numOutVars_, 0.0);
 
}


//-----------------------------------------------------------------------------
// Function      : RelativeError::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void RelativeError::updateTran( const double circuitTime, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{  

  //Xyce::dout() << "in RelativeError updateTran " << std::endl;
  if( !calculationDone_ && withinFromToWindow( circuitTime ) )
  {
    // we're in the transient window, now we need to calculate the value of this
    // measure and see if it triggers any specified rise, fall, cross windowing.

    // update our outVarValues_ vector
    for( int i=0; i< numOutVars_; i++ )
    {
      outVarValues_[i] = getOutputValue(outputVars_[i], solnVec, stateVec, storeVec, 0);
      //Xyce::dout() << " outVarValues_[" << i << "]= " << outVarValues_[i] << std::endl;
    }

    // scan the external data supplied to see if we are either right on a data point 
    // or need to interpolate around it. 
     
    for( int i=lastIndepIndex_; i< indepVarValues.size(); i++ )
    {
      double currentDiff = indepVarValues[i] - circuitTime;
      double lastDiff = indepVarValues[i] - lastDepVar_;

      if( std::fabs( currentDiff ) < minval_ )
      {
        // close enough that we can call this as being on the last data point
        simulationDataVals_[i] = outVarValues_[0];
        simulationDataValsFound_[i] = 1;
        lastIndepIndex_=i;
        //Xyce::dout() << " found matching value: " << circuitTime << ", " << currentDiff << ", " << outVarValues_[0] 
        //  << ", : " << dataValues[i][ dependentVarColumn_ ] << std::endl;
        break; 
      }
      else if( ((currentDiff > 0.0) && (lastDiff < 0.0)) || ((currentDiff < 0.0) && (lastDiff > 0.0)) ) 
      {
        // interpolate between these two points 
        double slope = (outVarValues_[0] - lastOutVarValues_[0]) / (circuitTime-lastDepVar_);
        simulationDataVals_[i] = lastOutVarValues_[0] + slope*(indepVarValues[i] - lastDepVar_);
        simulationDataValsFound_[i] = 1;
        lastIndepIndex_=i;
        //Xyce::dout() << " found boundings values: " << circuitTime << ", " << currentDiff << ", " << lastDiff << ", " << outVarValues_[0]
        //  << "intep= " << simulationDataVals_[i] << ", : " << dataValues[i][ dependentVarColumn_ ] << std::endl;
        break;
      }
    }
     
    // store any values needed for interpolation to given datapoints 
    for( int i=0; i<outVarValues_.size(); i++)
      lastOutVarValues_[i] = outVarValues_[i];

    lastDepVar_=circuitTime; 

  }

}


//-----------------------------------------------------------------------------
// Function      : RelativeError::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void RelativeError::updateDC( const std::vector<N_ANP_SweepParam> & dcParamsVec, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{
  //Xyce::dout() << "in RelativeError updateDC " << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : RelativeError::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2014
//-----------------------------------------------------------------------------
void RelativeError::updateAC( const double frequency, const N_LAS_Vector * solnVec, const N_LAS_Vector *imaginaryVec)
{
  //Xyce::dout() << "in RelativeError updateAC " << std::endl;
  //if( !calculationDone_ && withinFromToWindow( circuitTime ) )
  {
    // we're in the transient window, now we need to calculate the value of this
    // measure and see if it triggers any specified rise, fall, cross windowing.

    // update our outVarValues_ vector
    for( int i=0; i< numOutVars_; i++ )
    {
      outVarValues_[i] = getOutputValue(outputVars_[i], solnVec, 0, 0, imaginaryVec);
      //Xyce::dout() << " outVarValues_[" << i << "]= " << outVarValues_[i] << std::endl;
    }

    // scan the external data supplied to see if we are either right on a data point 
    // or need to interpolate around it. 
     
    for( int i=lastIndepIndex_; i< indepVarValues.size(); i++ )
    {
      double currentDiff = indepVarValues[i] - frequency;
      double lastDiff = indepVarValues[i] - lastDepVar_;

      if( std::fabs( currentDiff ) < minval_ )
      {
        // close enough that we can call this as being on the last data point
        simulationDataVals_[i] = outVarValues_[0];
        simulationDataValsFound_[i] = 1;
        lastIndepIndex_=i;
        //Xyce::dout() << " found matching value: " << frequency << ", " << currentDiff << ", " << outVarValues_[0] 
        //  << ", : " << dataValues[i][ dependentVarColumn_ ] << std::endl;
        break; 
      }
      else if( ((currentDiff > 0.0) && (lastDiff < 0.0)) || ((currentDiff < 0.0) && (lastDiff > 0.0)) ) 
      {
        // interpolate between these two points 
        double slope = (outVarValues_[0] - lastOutVarValues_[0]) / (frequency-lastDepVar_);
        simulationDataVals_[i] = lastOutVarValues_[0] + slope*(indepVarValues[i] - lastDepVar_);
        simulationDataValsFound_[i] = 1;
        lastIndepIndex_=i;
        //Xyce::dout() << " found boundings values: " << frequency << ", " << currentDiff << ", " << lastDiff << ", " << outVarValues_[0]
        //  << "intep= " << simulationDataVals_[i] << ", : " << dataValues[i][ dependentVarColumn_ ] << std::endl;
        break;
      }
    }
     
    // store any values needed for interpolation to given datapoints 
    for( int i=0; i<outVarValues_.size(); i++)
      lastOutVarValues_[i] = outVarValues_[i];

    lastDepVar_=frequency; 

  }


}

//-----------------------------------------------------------------------------
// Function      : Average::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double RelativeError::getMeasureResult()
{

  // make an epetra vector to hold difference values
  Epetra_SerialDenseVector differenceVector(simulationDataVals_.size());
  int numberPointsFound=0;
  for( int i=0; i<simulationDataValsFound_.size(); i++ )
  {
    if( simulationDataValsFound_[i] == 1 )
      numberPointsFound++;
  }
  //Xyce::dout() << "Found " << numberPointsFound << " out of a possible " << dataValues.size() << std::endl;
  
  // load the difference vector
  for (int i=0 ; i<simulationDataVals_.size(); i++)
  {
    //differenceVector[i] = simulationDataValsFound_[i]*(simulationDataVals_[i]-dataValues[i][dependentVarColumn_]);
    differenceVector[i] = (simulationDataVals_[i]-dataValues[i][dependentVarColumn_]);
  }

  if( comparisonFunctionName_ == "L1NORM" )
  {
    calculationResult_ = differenceVector.Norm1();
  }
  else if( comparisonFunctionName_ == "L2NORM" )
  {
    calculationResult_ = differenceVector.Norm2();
  }
  else
  {
    calculationResult_ = differenceVector.NormInf();
  }

  return calculationResult_;
}




//
// utility functions for opening and reading in a file 
//

std::ifstream * openStreamFromFileName( const std::string fileName )
{
  return new std::ifstream( fileName.c_str() );
}

bool readData( std::istream & inputStream, 
               std::vector<std::string> & varNames, 
               std::vector<double> & indepVarValues,
               std::vector<double> & indep2VarValues,
               std::vector< std::vector< double > > & dataValues )
{
  // read in the first line and try to devine the format of the file.
  // some possibliities are:
  //
  // (1) white space deliminated table with the first line optionally being 
  //     the name of each column and potentially a line of text 
  //     at the end.  This is the same as Xyce's prn format.
  //     Example: 
  //
  //     [ text  text  text ]
  //     number number number 
  //     ... 
  //     number number number 
  //     [ optional end of simulation text ] 
  //
  // (2) a compressed DC sweep format where the first line is the 
  //     name of the first independent variable and the second line
  //     gives the name of the second independent variable and 
  //     the values of the seconde indepentent variable for each column 
  //
  //     Example:
  //     text
  //     text value1 ... valueN
  //     value0 ... valueN
  //     value0 ... valueN
  //     value0 ... valueN
  //
  //     real example
  //     VD
  //     VG   -1.5 -0.5 0.0 0.5 1.0
  //     -0.9 0.00  0.1 0.2 0.1 0.1
  //     -0.8 0.00  0.1 0.2 0.1 0.1
  //
  // (3) probe or "CSD" formatted.  Starts with "#H" for header and "#N" for data


  // current line being parsed
  std::string workingLine;

  // try to read the first line
  std::getline( inputStream, workingLine );
  if( !(inputStream.good()) )
  {
    // throw error that read failed.
  }

  // a stringstream used in parsing up the line that was just read
  std::stringstream streamForWorkingLine;
  int dataLinesRead=0;
  int headerLinesRead=0;
  int linecount=0;
  // store up data read for each line so it can 
  // later be saved in the vector< vector< double > >
  std::vector<double> aDataLine;

  
  while( inputStream.good() )
  { 
    // if fist line is "#H" then we know we have a probe formatted file
    if( workingLine.find( "#H" ) == 0 )
    {
      // have a probe formatted file so read that 
    }
  
    // look for one or more text headers (there may be none)
    // trying to decern "text" from "number" at this stage 
    //streamForWorkingLine.flush();
    streamForWorkingLine.clear();
    streamForWorkingLine << workingLine;
    std::string  subString;
    // break up line by whitespace 
    while( streamForWorkingLine >> std::ws >> subString )
    {
      if( Util::isValue( subString ) )
      {
        // found a numeric value 
        aDataLine.push_back( Util::Value(subString) );
        //dataValues[linecount-1].push_back( Util::Value(subString) );
      }
      else 
      {
        varNames.push_back(subString);
        headerLinesRead++;
      }
    }
    
    if( !aDataLine.empty() )
    {
      dataValues.push_back( aDataLine );
      aDataLine.clear();
    }
    // 
    std::getline( inputStream, workingLine );
    linecount++;
  } 

/*
  // for debugging purposes output the data read 
  Xyce::dout() << "Data read from stream " << std::endl;
  for( int i=0; i<varNames.size(); i++ )
    Xyce::dout() << " varNames[" << i << "]=" << varNames[i] << std::endl;

  Xyce::dout() << "Data values: " << std::endl;
  for( int i=0; i<dataValues.size(); i++ )
  {
    for( int j=0; j<dataValues[i].size(); j++ )
      Xyce::dout() << dataValues[i][j] << "\t";
    Xyce::dout() << std::endl;
  }
*/


  return true;
}

void closeStream( std::ifstream * streamToClose )
{
  streamToClose->close();
  delete streamToClose;
}



} // namespace Measure
} // namespace IO
} // namespace Xyce
