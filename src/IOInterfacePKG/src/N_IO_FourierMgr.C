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
// Filename      : $RCSfile: N_IO_FourierMgr.C,v $
// Purpose       : This file contains the functions to manage fourier objects
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.1.2.10 $
// Revision Date  : $Date: 2013/10/03 17:23:43 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_IO_FourierMgr.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputMgr.h>

//-----------------------------------------------------------------------------
// Function      : N_IO_FourierMgr::N_IO_FourierMgr
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
N_IO_FourierMgr::N_IO_FourierMgr( N_IO_OutputMgr &outputManager )
  : outputManager_(outputManager),
    numFreq_(10),
    gridSize_(200),
    calculated_(false)
{
  // Initialize first output variable pointer.
  outputVarsPtr_.push_back( 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_IO_FourierMgr::N_IO_FourierMgr
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
N_IO_FourierMgr::~N_IO_FourierMgr()
{

}

//-----------------------------------------------------------------------------
// Function      : N_IO_FourierMgr::addFourierAnalysis
// Purpose       : Entry point when .four lines are pased in the netlist
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
bool N_IO_FourierMgr::addFourierAnalysis(const N_UTL_OptionBlock & fourierBlock)
{
  // based on what's in the option block passed in, we
  // create the needed fourier instance
#ifdef Xyce_DEBUG_IO
  std::cout << "In N_IO_FourierMgr::addFourier" << std::endl;
  std::cout << ".FOUR line passed was: " << std::endl << fourierBlock << std::endl;
#endif

  if( !fourierBlock.tagExists( "FREQ" ) )
  {
    // this shouldn't happen, but catch it if does
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, "Missing FREQ in .FOUR line!");
  }

  int numSolVars = 0;

  list<N_UTL_Param>::iterator currentParamIt = const_cast<N_UTL_OptionBlock &>(fourierBlock).begin();
  list<N_UTL_Param>::iterator endParamIt = const_cast<N_UTL_OptionBlock &>(fourierBlock).end();
  while( currentParamIt != endParamIt )
  {
    std::string name = "";
    string tag = currentParamIt->tag();

    if( tag == "FREQ" )
    {
      freqVector_.push_back( currentParamIt->dVal() );
    }
    else if( (tag == "V") || (tag == "I") )
    {
      int nodes = currentParamIt->iVal();
      N_UTL_Param aParam;
      aParam.set( tag, nodes );
      name += tag;
      name += "(";

      // here we just store the needed parts of V(a) or v(a,b) or I(device).
      // only the v(a,b) case will need an extra node in the outputVars_ array.
      numSolVars++;
      // use list::insert() to keep an iterator pointing to this spot
      list<N_UTL_Param>::iterator newDepSolVarItr = outputVars_.insert( outputVars_.end(), aParam);
      depSolVarIterVector_.push_back( newDepSolVarItr );
      for( int i=0; i<nodes; i++ )
      {
        currentParamIt++;
        aParam.set( currentParamIt->tag(), currentParamIt->dVal() );
        name += currentParamIt->tag(); 
        if (i != nodes-1 && nodes > 1) name += ",";
        outputVars_.push_back( aParam );
      }
      name += ")";
    }
    currentParamIt++;

    // Save voltage or current variable name.
    if (name != "")
    {
      names_.push_back( name );
    }
  }

  // Store number of solution variables for this .four line.
  int oVPsize = outputVarsPtr_.size();
  outputVarsPtr_.push_back( outputVarsPtr_[oVPsize-1] + numSolVars );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_FourierMgr::updateFourierData
// Purpose       : Called during the simulation to update the fourier objects
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
void N_IO_FourierMgr::updateFourierData( const double circuitTime, RCP< N_LAS_Vector > solnVecRCP)
{
  // Save the time.
  time_.push_back(circuitTime);

  for (unsigned int i=0; i<depSolVarIterVector_.size(); ++i)
  {
     if( depSolVarIterVector_[i]->getSimContext() == UNDEFINED )
     {
       // call set param context
       outputManager_.setParamContextType_( depSolVarIterVector_[i] );
     }
     double retVal = outputManager_.getPrintValue( depSolVarIterVector_[i], solnVecRCP.getRawPtr() );

     outputVarsValues_.push_back(retVal);
  }
}


//-----------------------------------------------------------------------------
// Function      : N_IO_FourierMgr::outputResults
// Purpose       : Output fourier results at end of simulation
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
void N_IO_FourierMgr::outputResults(std::ostream& outputStream)
{
  // Only calculate something if a .four line was encountered and transient data was collected.
  int numOutVars = depSolVarIterVector_.size();

  if ( numOutVars && !time_.empty() && !calculated_ )
  {
    // Calculated the fourier coefficients for the given nodes.
    getLastPeriod_();
    interpolateData_();
    calculateFT_();
    calculated_ = true;
  }

  // Output the information to the outputStream 
  printResult_( outputStream );

}


//-----------------------------------------------------------------------------
// Function      : N_IO_FourierMgr::getLastPeriod_()
// Purpose       : finds the indices to access the last period of simulation
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
void N_IO_FourierMgr::getLastPeriod_()
{
  // We want to do the analysis on only the last period of the transient window. So here we find the indices
  // to access the endpoints of that interval.
  int numPoints = time_.size();
  int prdEnd = numPoints - 1;
  double endTime = time_[prdEnd];

  int nFreq = freqVector_.size();
  lastPrdStart_.resize( nFreq );
  prdStart_.resize( nFreq );

  for (int i=0; i<nFreq; ++i)
  {
    double period = 1.0/freqVector_[i];
    lastPrdStart_[i] = endTime - period;
  
    if (lastPrdStart_[i] > 0)
    {
      // Initialize prdStart_ to be the index of the last element in time_.
      // Then scan until time_[i] <= endTime_ - period_.
      prdStart_[i] = prdEnd;
      while (time_[prdStart_[i]] > lastPrdStart_[i])
      {
        prdStart_[i] -= 1;
      }
    }
    else if (lastPrdStart_[i] == 0)
    {
      prdStart_[i] = 0;
    } 
    else
    {
      string msg = "Error: The period is greater than the length of the time simulation. Exiting.";
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL, msg);
    }
  }
} 

//-----------------------------------------------------------------------------
// Function      : N_IO_FourierMgr::interpolateData_()
// Purpose       : evaluates interpolating polynomial at equidistant time pts
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
bool N_IO_FourierMgr::interpolateData_()
{
  double A, B, C;

  int numOutVars = depSolVarIterVector_.size();
  int numPoints = time_.size();
  int nFreq = freqVector_.size();

  newTime_.resize(nFreq * gridSize_,0.0);
  newValues_.resize(numOutVars * gridSize_,0.0);

  for (int j=0; j<nFreq; ++j)
  {
    // Get the number of solution variables associated with this frequency
    int numCurrVars = outputVarsPtr_[j+1] - outputVarsPtr_[j];

    // Get the number of time points in the final period
    int nData = numPoints-prdStart_[j];
    double period = 1.0/freqVector_[j];
    vector<double> h(nData-1, 0.0);
    vector<double> b(nData-1, 0.0);
    vector<double> u(nData-1, 0.0);
    vector<double> v(nData-1, 0.0);
    vector<double> z(nData, 0.0);

    // Compute new, equally spaced time points.
    double step = period/gridSize_;

    newTime_[gridSize_*j] = lastPrdStart_[j];
    for (int i = 1; i < gridSize_; i++) 
    {
      newTime_[gridSize_*j + i] = newTime_[gridSize_*j + i-1] + step;
    }
  
    // Loop over all the data from the output variables associated with this frequency 
    for (int k=0; k<numCurrVars; k++)
    { 
      int offset = outputVarsPtr_[j] + k;

      // Cubic spline interpolation. We first need to find the z's.
      for (int i = 0; i < nData-1; i++)
      {
        h[i] = time_[i+1+prdStart_[j]]-time_[i+prdStart_[j]];
        b[i] = (6/h[i])*(outputVarsValues_[(i+1+prdStart_[j])*numOutVars + offset]
                         - outputVarsValues_[(i+prdStart_[j])*numOutVars + offset]);
      }

      u[1] = 2*(h[0]+h[1]);
      v[1] = b[1]-b[0];

      for (int i=2; i < nData-1; i++)
      {
        u[i] = 2*(h[i]+h[i-1])-((h[i-1])*(h[i-1]))/u[i-1];
        v[i] = b[i]-b[i-1]-(h[i-1]*v[i-1])/u[i-1];
      }

      z[nData-1] = 0;
      for (int i=nData-2; i > 0; i--)
      {
        z[i] = (v[i]-h[i]*z[i+1])/u[i];
      }
      z[0] = 0;

      // Calculate the new values at the new time points.
      for (int i = 0; i < gridSize_; i++)
      { 
        int idx = nData-1;
        while ((newTime_[gridSize_*j + i]-time_[idx+prdStart_[j]]) < 0)
        {
          idx--;
        }
        A = (z[idx+1]-z[idx])/(6*h[idx]);
        B = z[idx]/2;
        C = -(h[idx]/6)*z[idx+1]-(h[idx]/3)*z[idx]+(outputVarsValues_[(idx+1+prdStart_[j])*numOutVars + offset]
                                                    - outputVarsValues_[(idx+prdStart_[j])*numOutVars + offset])/h[idx];
        newValues_[gridSize_*offset + i] = outputVarsValues_[(idx+prdStart_[j])*numOutVars + offset] + (newTime_[gridSize_*j + i]-time_[idx+prdStart_[j]])*(C+(newTime_[gridSize_*j + i]-time_[idx+prdStart_[j]])*(B+(newTime_[gridSize_*j + i]-time_[idx+prdStart_[j]])*A));
      }
    } 
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_FourierMgr::calculateFT_()
// Purpose       : performs fourier analysis
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
void N_IO_FourierMgr::calculateFT_()
{ 
  int numOutVars = depSolVarIterVector_.size();
  int nFreq = freqVector_.size();

  mag_.resize(numFreq_*numOutVars, 0.0);
  phase_.resize(numFreq_*numOutVars, 0.0);
  nmag_.resize(numFreq_*numOutVars, 0.0);
  nphase_.resize(numFreq_*numOutVars, 0.0);
  freq_.resize(numFreq_*nFreq, 0.0);
  thd_.resize(numOutVars, 0.0);

  // Compute frequencies 
  for (int j=0; j<nFreq; ++j)
  {
    for (int i=0; i < numFreq_; i++)
    {
      freq_[j*numFreq_ + i] = i * freqVector_[j];
    }
  }
   
  // Perform Fourier analysis for all the output variables. 
  for (int k=0; k < numOutVars; k++)
  { 
    for (int i=0; i < gridSize_; i++)
    { 
      for (int j=0; j < numFreq_; j++)
      {
        mag_[numFreq_*k + j] += (newValues_[gridSize_*k+i]*sin(j*2.0*M_PI*i/((double) gridSize_)));
        phase_[numFreq_*k + j] += (newValues_[gridSize_*k+i]*cos(j*2.0*M_PI*i/((double) gridSize_)));
      }
    }
  
    mag_[numFreq_*k] = phase_[numFreq_*k]/gridSize_;
    phase_[numFreq_*k] = 0;
    thd_[k] = 0; 

    double convRadDeg = 180.0/M_PI;

    for(int i = 1; i < numFreq_ ; i++)
    { 
      double tmp = mag_[numFreq_*k+i]*2.0 /gridSize_;
      phase_[numFreq_*k+i] *= 2.0/gridSize_;
      mag_[numFreq_*k+i] = sqrt(tmp*tmp+phase_[numFreq_*k+i]*phase_[numFreq_*k+i]);
      phase_[numFreq_*k+i] = atan2(phase_[numFreq_*k+i],tmp)*convRadDeg;
      nmag_[numFreq_*k+i] = mag_[numFreq_*k+i]/mag_[numFreq_*k+1];
      nphase_[numFreq_*k+i] = phase_[numFreq_*k+i]-phase_[numFreq_*k+1];
      if(i>1) thd_[k] += nmag_[numFreq_*k+i]*nmag_[numFreq_*k+i];
    }
    thd_[k] = 100*sqrt(thd_[k]);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_IO_FourierMgr::printResult_( std::ostream& os )
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
std::ostream& N_IO_FourierMgr::printResult_( std::ostream& os )
{ 
  // Compute measure.
  if (calculated_)
  { 
  
    int nFreq = freqVector_.size();
    for (unsigned int i=0; i<nFreq; ++i)
    {
      for (int j=outputVarsPtr_[i]; j<outputVarsPtr_[i+1]; j++)
      {
        int colWidth = 12, colWidth2 = 5;
        os << "Fourier analysis for " << names_[j] << ":" << std::endl;
        os << "  No. Harmonics: " << numFreq_ << ", THD: " << thd_[j] << ", Gridsize: " << gridSize_
           << ", Interpolation Type: Cubic Spline" << std::endl;
  
        os << std::setw(colWidth) << "Harmonic" << std::setw(colWidth) << "Frequency"
           << std::setw(colWidth) << "Magnitude" << std::setw(colWidth) << "Phase"
           << std::setw(colWidth) << "Norm. Mag" << std::setw(colWidth) << "Norm. Phase" << std::endl;
        for (int k = 0; k < numFreq_; k++)
        {
           os << std::setw(colWidth) << k << std::setw(colWidth) << freq_[numFreq_*i+k]
              << std::setw(colWidth) << std::setprecision(colWidth2) << mag_[numFreq_*j+k]
              << std::setw(colWidth) << std::setprecision(colWidth2) << phase_[numFreq_*j+k]
              << std::setw(colWidth) << std::setprecision(colWidth2) << nmag_[numFreq_*j+k]
              << std::setw(colWidth) << std::setprecision(colWidth2) << nphase_[numFreq_*j+k] << std::endl;
        }
        os << std::endl;
      }  
    }
  }
  return os;
}

