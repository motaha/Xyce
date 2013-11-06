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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_LOCA_StepSizeControl.h,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : 
//
// Creation Date  : 
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:47 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------





#ifndef N_NLS_LOCA_STEPSIZECONTROL_H
#define N_NLS_LOCA_STEPSIZECONTROL_H

#include "LOCA_StepSize_Generic.H"  // base class

namespace N_NLS_LOCA {

  //! %Adaptive step size control strategy
  /*!
    This class implements an adaptive step size control strategy derived
    from the strategy implemented in the LOCA::StepSize::Constant class.
    If the previous step was unsucessful, the step size is cut in half as
    in the constant strategy, but if the step was sucessful this strategy
    increases the step size based on the number of nonlinear solver 
    iterations required in the previous step.  In particular, the new
    step size \f$\Delta s_{new}\f$ is given by
    \f[
    \Delta s_{new} = \Delta s_{old}\left(1 + a\left(\frac{N_{max} - N}{N_{max}}\right)^2\right)
    \f]
    where \f$a\in[0,1]\f$ is an aggressiveness factor, \f$N\f$ is the 
    number of nonlinear solver iterations in the previous step, and
    \f$N_{max}\f$ is the maximum number of nonlinear solver iterations.
    
    The parameters used by this class supplied in the constructor or reset
    method are the same as used by the Constant class in addition to:
    <ul>
    <li> "Aggressiveness" - Aggressiveness factor \f$a\f$ (Default 0.0)
    </ul>
  */
  class StepSizeControl : public LOCA::StepSize::Generic {
    
  public:
    
    //! Constructor. 
    StepSizeControl();
    
    //! Destructor
    virtual ~StepSizeControl();
    
    virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& params);
    
    virtual NOX::Abstract::Group::ReturnType 
      compute(LOCA::Continuation::ExtendedGroup& curGroup,
	      const LOCA::Continuation::ExtendedVector& predictor,
	      const NOX::Solver::Generic& solver,
	      const LOCA::Abstract::Iterator::StepStatus& stepStatus,
	      const LOCA::Stepper& stepper,
	      double& stepSize);
    
    virtual NOX::Abstract::Group::ReturnType 
      compute(LOCA::MultiContinuation::AbstractStrategy& curGroup,
	      const LOCA::MultiContinuation::ExtendedVector& predictor,
	      const NOX::Solver::Generic& solver,
	      const LOCA::Abstract::Iterator::StepStatus& stepStatus,
	      const LOCA::NewStepper& stepper,
	      double& stepSize);

    virtual double getPrevStepSize() const;

    virtual double getStartStepSize() const;
    
  protected:
    
    virtual NOX::Abstract::Group::ReturnType
      clipStepSize(double& stepSize);
    
  protected:
    //! Maximum step size
    double maxStepSize;
    
    //! Minimum step size
    double minStepSize;
    
    //! Initial step size
    double startStepSize;
    
    //! Factor by which step size is reduced after a failed step
    double failedFactor;
    
    //! Factor by which step size is increased after a successful step
    double successFactor;
    
    //! Previous step size
    double prevStepSize;

    //! Flag indicating if this is the first step
    bool isFirstStep;
    
    //! Stores the aggressiveness factor \f$a\f$
    double agrValue;
  }; 
} 

#endif
