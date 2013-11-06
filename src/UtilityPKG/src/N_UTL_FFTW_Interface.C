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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_FFTW_Interface.C,v $
//
// Purpose        : This file contains specializations for the FFTW interface
//                  for various vector types.
//
// Special Notes  : 
//
// Creator        : Heidi Thornquist
//
// Creation Date  : 5/27/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.4.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:52 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
// ---------- Standard Includes ----------

#include <Xyce_config.h>

#include <N_LAS_BlockVector.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>

#include <N_UTL_FFTW_Interface.hpp>

// ----------   Other Includes   ----------

#include <iostream>
#include <vector>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// FFTW Interface specialization for std::vector
//-----------------------------------------------------------------------------

  template<>
  void N_UTL_FFTW_Interface<std::vector<double> >::calculateFFT()
  {
    // If the plan needs to be constructed do that first, then execute the plan
    if (firstForwardFFT_)
    {
      // Although we used a const on input to show that we aren't changing fftInData_
      // we need to cast that away as the FFT library takes non-const pointers.
      std::vector<double>::const_iterator inDataItr = (this->fftInData_)->begin();
      double * inDataPtr = const_cast< double * >( &(*inDataItr) );

      // We need to create a temp vector for converting the storage format and
      // set it as extra data on the fftOutData_ RCP.
      Teuchos::RCP<std::vector<double> > outResultTmp = 
        Teuchos::rcp( new std::vector<double>((this->fftOutData_)->size(),0.0) );
      Teuchos::set_extra_data( outResultTmp, "outResultTmp", inOutArg(this->fftOutData_) );  

      forwardPlan_ = fftw_plan_r2r_1d(signalLength_, inDataPtr, &(*outResultTmp)[0], 
                                     FFTW_R2HC, FFTW_ESTIMATE );
      firstForwardFFT_ = false;
    }

    // Execute the FFT.
    fftw_execute(forwardPlan_);
 
    // Now get the vector back that was used to store the results from FFTW and copy
    // the data into the appropriate place in fftOutData_. 
    Teuchos::RCP<std::vector<double> > outResultTmp = 
      Teuchos::get_extra_data<Teuchos::RCP<std::vector<double> > >(this->fftOutData_, "outResultTmp"); 
    int n2 = (int)(signalLength_/2);
    (*(this->fftOutData_))[0] = (*outResultTmp)[0];
    (*(this->fftOutData_))[1] = 0.0;
    for(int i=1; i<=n2; ++i)
    { 
      (*(this->fftOutData_))[2*i] = (*outResultTmp)[i];
      (*(this->fftOutData_))[2*i+1] = (*outResultTmp)[signalLength_-i];
    }
  }

  // Calculate IFT with the vectors that have been registered.
  template<>
  void N_UTL_FFTW_Interface<std::vector<double> >::calculateIFT()
  {
    // If the plan needs to be constructed do that first, then execute the plan
    if (firstInverseFFT_)
    {
      // We need to create a temp vector for converting the storage format and
      // set it as extra data on the iftOutData_ RCP.
      Teuchos::RCP<std::vector<double> > inDataTmp =
        Teuchos::rcp( new std::vector<double>((this->iftInData_)->size(),0.0) );
      Teuchos::set_extra_data( inDataTmp, "inDataTmp", inOutArg(this->iftInData_) );

      inversePlan_ = fftw_plan_r2r_1d(signalLength_, &(*inDataTmp)[0], &(*this->iftOutData_)[0], 
                                     FFTW_HC2R, FFTW_ESTIMATE );
      firstInverseFFT_ = false;
    }

    // Now get the vector back that was used to store the results from FFTW and copy
    // the data into the appropriate place in iftInData_. 
    Teuchos::RCP<std::vector<double> > inDataTmp = 
      Teuchos::get_extra_data<Teuchos::RCP<std::vector<double> > >(this->iftInData_, "inDataTmp"); 

    int n2 = (int)(signalLength_/2);
    (*inDataTmp)[0] = (*(this->iftInData_))[0];
    for(int i=1; i<=n2; ++i)
    {
      (*inDataTmp)[i] = (*(this->iftInData_))[2*i];
      (*inDataTmp)[signalLength_-i] = (*(this->iftInData_))[2*i+1];
    }

    // Execute the IFT.
    fftw_execute(inversePlan_);

    // Scale the output by "n"
    for (int i=0; i<this->iftOutData_->size(); ++i)
    {
      (*(this->iftOutData_))[i] /= signalLength_;
    }
  }
