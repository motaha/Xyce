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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: Charon_CircuitInterface.h,v $
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
// Revision Number: $Revision: 1.9 $
//
// Revision Date  : $Date: 2014/02/24 23:49:15 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------


#ifndef _CHARON_CIRCUIT_INTERFACE_H_
#define _CHARON_CIRCUIT_INTERFACE_H_

#include <map>
#include <vector>
#include "Teuchos_RefCountPtr.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace charon {

  namespace sc {

    /**
     * @brief - Provides the interface for external circuit codes to call charon to take a transient step.
     *
     * This is a singleton so we can avoid punching layers into xyce
     * to pass this in.
     */
    class CircuitInterface {
      
    public:
      
      //! Destructor.
      ~CircuitInterface();

      //! Returns an instance of this object. 
      static CircuitInterface& getInstance();

      //! Calls the CharonClient::Run_Step_CircuitSimulation routine to start a charon solve.
      bool takeStep(
            const Teuchos::RefCountPtr<Teuchos::ParameterList>& inputList,
	    const std::map<std::string, double>& inputMap,
	    const Teuchos::RefCountPtr<Teuchos::ParameterList>& outputList,
	    std::vector<double>& outputVector,
	    std::vector< std::vector<double> >& outputJacobian);

      //! Tells Charon to accept the time step and write output.
      void acceptTimeStep(bool& is_active);

    private:
      
      /*! Tells Charon to reset the time integrator and that we will be making another solve at the current time step.  This should be called after every solve of charon unless we accept the step, otherwise, the transient history will be corrupted.  When Charon attempts a step, it pushes the new solution and time step size onto the stored history stack.  If the step fails we need to remove it from the stack.
       */
      void declineTimeStep();

      //! Singleton - disallow constructor.
      CircuitInterface();

      //! Singleton - disallow copy constructor.
      CircuitInterface(const CircuitInterface& source) {};

    protected:

      //! Used to determine if the next step should be reset. 
      bool accepted_step_;

    };
    
  } // END "namespace sc"
  
} // END "namespace charon"

#endif 
