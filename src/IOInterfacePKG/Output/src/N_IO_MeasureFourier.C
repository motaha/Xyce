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
// Filename      : $RCSfile: N_IO_MeasureFourier.C,v $
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/05/2013
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.12 $
// Revision Date  : $Date: 2014/02/24 23:49:20 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <iomanip>

#include <N_IO_MeasureFourier.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : Fourier::Fourier()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/05/2013
//-----------------------------------------------------------------------------
Fourier::Fourier( const Util::OptionBlock & measureBlock, N_IO_OutputMgr &outputMgr )
  : Base(measureBlock, outputMgr),
  initialized_(false),
  calculated_(false)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;
}

void Fourier::prepareOutputVariables() 
{
  // this measurement should have only one dependent variable.
  // Error for now if it doesn't
  numOutVars_ = outputVars_.size();

  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for statistical measure, \"" + name_ + "\" Exiting.";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL, msg);
  }
}


//-----------------------------------------------------------------------------
// Function      : Fourier::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/05/2013
//-----------------------------------------------------------------------------
void Fourier::updateTran( const double circuitTime, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{
  if( !calculationDone_ && withinTransientWindow( circuitTime ) )
  {
     // we're in the transient window, now we need to store the time and output values

    if( !initialized_  )
    {
      initialized_ = true;
    }

    time_.push_back(circuitTime);
    outVarValues_.push_back(getOutputValue(outputVars_[0], solnVec, stateVec, storeVec, 0));

  }
}

//-----------------------------------------------------------------------------
// Function      : Fourier::getLastPeriod_()
// Purpose       : finds the indices to access the last period of simulation
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
void Fourier::getLastPeriod_()
{
  // We want to do the analysis on only the last period of the transient window. So here we find the indices
  // to access the endpoints of that interval.
  period_ = 1/at_;
  int numPoints = time_.size();
  int prdEnd = numPoints-1;
  double endTime = time_[prdEnd];
  lastPrdStart_ = endTime-period_;

  if (lastPrdStart_ > 0)
  {
    // Initialize prdStart_ to be the index of the last element in time_.
    // Then scan until time_[i] <= endTime - period_.
    prdStart_ = prdEnd;
    while (time_[prdStart_] > lastPrdStart_)
    {
      prdStart_--;
    }
  }
  else if (lastPrdStart_ == 0)
  {
    prdStart_ = 0;
  }
  else
  {
    std::string msg = "Error: The period is greater than the length of the time simulation. Exiting.";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL, msg);
  }

}

//-----------------------------------------------------------------------------
// Function      : Fourier::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
void Fourier::updateDC( const std::vector<N_ANP_SweepParam> & dcParamsVec, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{
}

//-----------------------------------------------------------------------------
// Function      : Fourier::interpolateData_()
// Purpose       : evaluates interpolating polynomial at equidistant time pts
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
bool Fourier::interpolateData_()
{
  double A, B, C;
  int nData, j;

  // nData is the number of data points in the time_ and outVarValues_ vectors over the last period.
  int numPoints = time_.size();
  nData = numPoints-prdStart_;
  std::vector<double> h(nData-1, 0.0);
  std::vector<double> b(nData-1, 0.0);
  std::vector<double> u(nData-1, 0.0);
  std::vector<double> v(nData-1, 0.0);
  std::vector<double> z(nData, 0.0);

  // Cubic spline interpolation. We first need to find the z's.
  for (int i = 0; i < nData-1; i++)
  {
    h[i] = time_[i+1+prdStart_]-time_[i+prdStart_];
    b[i] = (6/h[i])*(outVarValues_[i+1+prdStart_]-outVarValues_[i+prdStart_]);
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

  // Compute new, equally spaced time points.
  newTime_.resize(gridSize_,0.0);
  newValues_.resize(gridSize_,0.0);
  double step = period_/gridSize_;

  newTime_[0] = lastPrdStart_;
  for (int i = 1; i < gridSize_; i++)
  {
    newTime_[i] = newTime_[i-1] + step;
  }

  // Calculate the new values at the new time points.
  for (int i = 0; i < gridSize_; i++)
  {
    j = nData-1;
    while ((newTime_[i]-time_[j+prdStart_]) < 0)
    {
      j--;
    }
    A = (z[j+1]-z[j])/(6*h[j]);
    B = z[j]/2;
    C = -(h[j]/6)*z[j+1]-(h[j]/3)*z[j]+(outVarValues_[j+1+prdStart_]-outVarValues_[j+prdStart_])/h[j];
    newValues_[i] = outVarValues_[j+prdStart_] + (newTime_[i]-time_[j+prdStart_])*(C+(newTime_[i]-time_[j+prdStart_])*(B+(newTime_[i]-time_[j+prdStart_])*A));
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Fourier::calculateFT_()
// Purpose       : performs fourier analysis
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
void Fourier::calculateFT_()
{
  double tmp;
  mag_.resize(numFreq_,0.0);
  phase_.resize(numFreq_,0.0);
  nmag_.resize(numFreq_,0.0);
  nphase_.resize(numFreq_,0.0);
  freq_.resize(numFreq_,0.0);

  for (int i=0; i < gridSize_; i++)
  {
    for (int j=0; j < numFreq_; j++)
    {
      mag_[j] += (newValues_[i]*sin(j*2.0*M_PI*i/((double) gridSize_)));
      phase_[j] += (newValues_[i]*cos(j*2.0*M_PI*i/((double) gridSize_)));
    }
  }

  double convRadDeg = 180.0/M_PI;

  mag_[0] = phase_[0]/gridSize_;
  phase_[0] = 0.0;
  thd_ = 0.0;
  for(int i = 1; i < numFreq_ ; i++)
  {
    tmp = mag_[i]*2.0 /gridSize_;
    phase_[i] *= 2.0/gridSize_;
    freq_[i] = i * at_;
    mag_[i] = sqrt(tmp*tmp+phase_[i]*phase_[i]);
    phase_[i] = atan2(phase_[i],tmp)*convRadDeg;
    nmag_[i] = mag_[i]/mag_[1];
    nphase_[i] = phase_[i]-phase_[1];
    if(i>1) thd_ += nmag_[i]*nmag_[i];
  }
  thd_ = 100*sqrt(thd_);

}

//-----------------------------------------------------------------------------
// Function      : Fourier::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
double Fourier::getMeasureResult()
{
  // Only output results if transient data is available to analyze.
  if( initialized_ && !time_.empty() )
  {
    getLastPeriod_();

    interpolateData_();

    calculateFT_();

    calculated_ = true;
  }
  // Total harmonic distortion will be the calculation result.
  // printMeasureResult will be overloaded to print out all the harmonic information
  return calculationResult_ = thd_;
}

//-----------------------------------------------------------------------------
// Function      : Fourier::printMeasureResult( std::ostream& os )
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
std::ostream& Fourier::printMeasureResult( std::ostream& os )
{

  // Only output results if transient data is available to analyze.
  if ( !time_.empty() )
  {
    // Compute measure.
    if (!calculated_)
    {
      this->getMeasureResult();
    }

    int colWidth = 12, colWidth2 = 5;
    os << name_ << ":  No. Harmonics: " << numFreq_ << ", THD: " << thd_ << ", Gridsize: " << gridSize_
       << ", Interpolation Type: Cubic Spline" << std::endl;

    os << std::setw(colWidth) << "Harmonic" << std::setw(colWidth) << "Frequency"
       << std::setw(colWidth) << "Magnitude" << std::setw(colWidth) << "Phase"
       << std::setw(colWidth) << "Norm. Mag" << std::setw(colWidth) << "Norm. Phase" << std::endl;
    for (int i = 0; i < numFreq_; i++)
    {
       os << std::setw(colWidth) << i << std::setw(colWidth) << freq_[i]
          << std::setw(colWidth) << std::setprecision(colWidth2) << mag_[i]
          << std::setw(colWidth) << std::setprecision(colWidth2) << phase_[i]
          << std::setw(colWidth) << std::setprecision(colWidth2) << nmag_[i]
          << std::setw(colWidth) << std::setprecision(colWidth2) << nphase_[i] << std::endl;
    }
  }
  
  return os;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce
