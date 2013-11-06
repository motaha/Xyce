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
// Filename       : $RCSfile: N_UTL_IntelFFT_Interface.C,v $
//
// Purpose        : This file contains specializations for the Intel FFT interface
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
// Revision Number: $Revision: 1.2.6.2 $
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

#include <N_UTL_IntelFFT_Interface.hpp>

// ----------   Other Includes   ----------

#include <iostream>
#include <vector>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Intel FFT Interface specialization for std::vector
//-----------------------------------------------------------------------------

  template<>
  void N_UTL_IntelFFT_Interface<std::vector<double> >::calculateFFT()
  {
    // Although we used a const on input to show that we aren't changing fftInData_
    // we need to cast that away as the FFT library takes non-const pointers.
    std::vector<double>::const_iterator inDataItr = (this->fftInData_)->begin();
    double * inDataPtr = const_cast< double * >( &(*inDataItr) );
    std::vector<double>::iterator outResultItr = (this->fftOutData_)->begin();
    double * outResultPtr = &(*outResultItr);

    long status = DftiComputeForward( fftDescriptor, inDataPtr, outResultPtr);
    checkAndTrapErrors( status );
  }

  // Calculate IFT with the vectors that have been registered.
  template<>
  void N_UTL_IntelFFT_Interface<std::vector<double> >::calculateIFT()
  {
    // Although we used a const on input to show that we aren't changing iftInData_
    // we need to cast that away as the FFT library takes non-const pointers.
    std::vector<double>::const_iterator inDataItr = (this->iftInData_)->begin();
    double * inDataPtr = const_cast< double * >( &(*inDataItr) );
    std::vector<double>::iterator outResultItr = (this->iftOutData_)->begin();
    double * outResultPtr = &(*outResultItr);

    long status = DftiComputeBackward( fftDescriptor, inDataPtr, outResultPtr);
    checkAndTrapErrors( status );
  }
