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
// Filename       : $RCSfile: N_IO_MeasureRelativeError.h,v $
//
// Purpose        : Measure relative error of an output variable
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 03/10/2009
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12 $
// Revision Date  : $Date: 2014/02/24 23:49:20 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureRelativeError_h
#define Xyce_N_IO_MeasureRelativeError_h

#include <N_IO_MeasureBase.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : RelativeError
// Purpose       : Measure relative error of an output variable
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-------------------------------------------------------------------------
class RelativeError : public Base
{
  public:
    RelativeError( const Util::OptionBlock & measureBlock, IO::OutputMgr &outputMgr );
    ~RelativeError() {};
    
    void prepareOutputVariables();
    void updateTran( const double circuitTime, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec);
    void updateDC( const std::vector<N_ANP_SweepParam> & dcParamsVec, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec);
    void updateAC( const double frequency, const N_LAS_Vector *solnVec, const N_LAS_Vector *imaginaryVec);
    double getMeasureResult();

  private:
    // these are used to hold data from an external file
    std::vector< std::string > varNames;
    std::vector< double > indepVarValues;
    std::vector< double > indep2VarValues;
    std::vector< std::vector< double > > dataValues;
    
    int numOutVars_;
    std::vector<double> outVarValues_;

    // useful for interpolation
    double lastDepVar_; 
    std::vector<double> lastOutVarValues_;
    int lastIndepIndex_;

    // results from the simulation to compare to a column in dataValues; 
    std::vector<double> simulationDataVals_;
    std::vector<int> simulationDataValsFound_;
     
};

// non-class utility functions

std::ifstream * openStreamFromFileName( const std::string fileName );
bool readData( std::istream & inputStream, 
               std::vector<std::string> & varNames, 
               std::vector<double> & indepVarValues,
               std::vector<double> & indep2VarValues,
               std::vector< std::vector< double > > & dataValues );
void closeStream( std::ifstream * streamToClose );

} // namespace Measure
} // namespace IO
} // namespace Xyce

typedef Xyce::IO::Measure::RelativeError N_IO_MeasureRelativeError;

#endif
