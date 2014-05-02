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
// Filename       : $RCSfile: N_MPDE_WarpedPhaseCondition.C,v $
//
// Purpose        : Specification of Warped MPDE phase condition (graph & load)
//
// Special Notes  :
//
// Creator        : Todd Coffey, SNL, 1414
//
// Creation Date  : 12/7/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.17 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


#include <N_MPDE_WarpedPhaseCondition.h>
#include <Epetra_MultiVector.h>
#include <N_PDS_Comm.h>

Teuchos::RefCountPtr<std::vector<int> > N_MPDE_WarpedPhaseCondition::getPhaseGraph()
{
  if ((warpMPDEOSCOUT_ == -1) && (warpPhase_ != 0))
  {
    std::string msg;
    msg = "N_MPDE_WarpedPhaseCondition::getPhaseGraph";
    msg += " No value for oscout which is required by specified phase equation\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  Teuchos::RefCountPtr<std::vector<int> > phaseGraph = Teuchos::rcp(new std::vector<int>);
  if ( warpPhase_ == 0 )
  {
    // Add graph for phase condition equation for MPDE:  omega-1=0
    phaseGraph->push_back(omegaGID_);
  }
  else if ( warpPhase_ == 1 )
  {
    // Add graph for phase condition 1 for WaMPDE:\hat{x}_w(t_1,0) - alpha = 0
    phaseGraph->push_back(warpMPDEOSCOUT_);
  }
  else if ( warpPhase_ == 2 )
  {
    // Add graph for phase condition 2 for WaMPDE:
    // \omega*((\hat{x}_w(t_1,0)-hat{x}_w(t_1,-h2))/h2 - alpha) = 0
    phaseGraph->push_back( warpMPDEOSCOUT_ );
    phaseGraph->push_back( warpMPDEOSCOUT_ + (size_-1)*offset_ );
    phaseGraph->push_back( omegaGID_ );
  }
  else if ( warpPhase_ == 3 )
  {
    // Add graph for phase condition 2 for WaMPDE:
    // \omega*((\hat{x}_w(t_1,h2)-hat{x}_w(t_1,-h2))/(2*h2) - alpha) = 0
    phaseGraph->push_back( warpMPDEOSCOUT_ + (   1   )*offset_ );
    phaseGraph->push_back( warpMPDEOSCOUT_ + (size_-1)*offset_ );
    phaseGraph->push_back( omegaGID_ );
  }
  else
  {
    std::string msg;
    msg = "N_MPDE_WarpedPhaseCondition::getPhaseGraph";
    msg += " Unrecognized value for WaMPDE Phase option\n";
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_WARNING, msg);
  }
  return(phaseGraph);
}


double N_MPDE_WarpedPhaseCondition::getPhaseCondition(
    N_LAS_BlockVector* bXptr,
    std::vector<double>& fastTimes
    )
{
  N_LAS_BlockVector& bX = *bXptr;
  int BlockCount = bX.blockCount();
  double phaseValue = 0.0, tmpPhaseValue = 0.0;
  int omegaLID = bX.pmap()->globalToLocalIndex(omegaGID_);
 
  // Get omega to all processors since it is not certain that the omega node is
  // on the same processor as the warpMPDEOSCOUT node.
  double omega = 0.0, tmpOmega = 0.0;
  if (omegaLID >= 0)
  {
    tmpOmega = bX[omegaLID];
  }
  bX.pmap()->pdsComm()->sumAll( &tmpOmega, &omega, 1 );

  // Get the local ID for the warpMPDEOSCOUT node, assume all block shifts are also
  // on the same processor.
  int warpMPDEOSCOUTLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_);
  if (warpMPDEOSCOUTLID >= 0)
  {
    if ( warpPhase_ == 0 )
    {
      tmpPhaseValue = omega - 1.0;
    }
    else if ( warpPhase_ == 1 )
    {
      // Add graph for phase condition 1 for WaMPDE:  \hat{x}_w(t_1,0) - alpha = 0
      tmpPhaseValue = bX[warpMPDEOSCOUTLID] - warpPhaseCoeff_;
    }
    else if ( warpPhase_ == 2 )
    {
      // Add graph for phase condition 2 for WaMPDE:
      // \omega*((\hat{x}_w(t_1,0)-hat{x}_w(t_1,-h2))/h2 - alpha) = 0
      int shift = (BlockCount-1)*offset_;
      int xwmLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shift);
      double xwp = bX[warpMPDEOSCOUTLID];
      double xwm = bX[xwmLID];
      double invh2 = fastTimes[BlockCount] - fastTimes[BlockCount-1];
      tmpPhaseValue = omega*(xwp - xwm)*invh2 - warpPhaseCoeff_ ;
    }
    else if ( warpPhase_ == 3 )
    {
      // Add graph for phase condition 2 for WaMPDE:
      // \omega*((\hat{x}_w(t_1,h2)-hat{x}_w(t_1,-h2))/(2*h2) - alpha) = 0
      int shiftp = offset_;
      int shiftm = (BlockCount-1)*offset_;
      int xwpLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shiftp);
      int xwmLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shiftm);
      double xwp = bX[xwpLID];
      double xwm = bX[xwmLID];
      double invh2 = (fastTimes[BlockCount] - fastTimes[BlockCount-1]) + (fastTimes[1] - fastTimes[0]);
      tmpPhaseValue = omega*(xwp - xwm)*(0.5*invh2) - warpPhaseCoeff_ ;
    }
    else
    {
      std::string msg;
      msg = "N_MPDE_WarpedPhaseCondition::getPhaseCondition";
      msg += " Unrecognized value for WaMPDE Phase option\n";
      N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_WARNING, msg);
    }
  }
 
  // Communicate to all processors.
  bX.pmap()->pdsComm()->sumAll( &tmpPhaseValue, &phaseValue, 1 );
   
  return(phaseValue);
}

void N_MPDE_WarpedPhaseCondition::getPhaseConditionDerivative(
    N_LAS_BlockVector* bXptr,
    std::vector<double>& fastTimes,
    std::vector<int>* colIndicesPtr,
    std::vector<double>* coeffsPtr
    )
{
  N_LAS_BlockVector & bX = *bXptr;
  std::vector<int> & colIndices = *colIndicesPtr;
  std::vector<double> & coeffs = *coeffsPtr;
  int BlockCount = bX.blockCount();
  
  // Get omega to all processors since it is not certain that the omega node is
  // on the same processor as the warpMPDEOSCOUT node.
  double omega = 0.0, tmpOmega = 0.0;
  int omegaLID = bX.pmap()->globalToLocalIndex(omegaGID_);
  if (omegaLID >= 0)
  {
    tmpOmega = bX[omegaLID];
  }
  bX.pmap()->pdsComm()->sumAll( &tmpOmega, &omega, 1 );

  int warpMPDEOSCOUTLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_);
  if ( warpMPDEOSCOUTLID >= 0 )
  {
    if ( warpPhase_ == 0 )
    {
      // put derivative wrt x_v of (omega-1) = 0 into omegaGID_,x_v locations
      coeffs.resize(1);
      coeffs[0] = 1.0;
      colIndices.resize(1);
      colIndices[0] = omegaGID_;
    }
    else if ( warpPhase_ == 1 )
    {
      // put derivative wrt x_v of (x_v-warpMPDECoeff) = 0 into omegaGID_,x_v location
      coeffs.resize(1);
      coeffs[0] = 1.0;
      colIndices.resize(1);
      colIndices[0] = warpMPDEOSCOUT_;
    }
    else if ( warpPhase_ == 2 )
    {
      // put derivative wrt x_v, xv+shift, and omega of
      // \omega*((\hat{x}_w(t_1,0)-\hat{x}_w(t_1,-h2))/h2 - alpha)
      //   = omega/h2 into  omegaGID_,x_v location
      //   = -omega/h2 into omegaGID_,x_v+shift location
      //   = (\hat{x}_w(t_1,0)-\hat{x}_w(t_1,-h2))/h2 into omegaGID_,omegaGID_ location
      int shift = (BlockCount-1)*offset_;
      int xwmLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shift);
      double xwp = bX[warpMPDEOSCOUTLID];
      double xwm = bX[xwmLID];
      double invh2 = fastTimes[BlockCount] - fastTimes[BlockCount-1];
      coeffs.resize(3);
      coeffs[0] =  omega*invh2;
      coeffs[1] = -omega*invh2;
      coeffs[2] = (xwp - xwm)*invh2;
      colIndices.resize(3);
      colIndices[0] = warpMPDEOSCOUT_;
      colIndices[1] = warpMPDEOSCOUT_+shift;
      colIndices[2] = omegaGID_;
    }
    else if ( warpPhase_ == 3 )
    {
      // put derivative wrt x_v, xv+shift, and omega of
      // \omega*((\hat{x}_w(t_1,h2)-hat{x}_w(t_1,-h2))/(2*h2) - alpha)
      //   = omega/(2*h2) into  omegaGID_,x_v+shiftp location
      //   = -omega/(2*h2) into omegaGID_,x_v+shiftm location
      //   = (\hat{x}_w(t_1,h2)-\hat{x}_w(t_1,-h2))/(2*h2) into omegaGID_,omegaGID_ location
      int shiftp = offset_;
      int shiftm = (BlockCount-1)*offset_;
      int xwpLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shiftp);
      int xwmLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shiftm);
      double xwp = bX[xwpLID];
      double xwm = bX[xwmLID];
      double invh2 = (fastTimes[BlockCount] - fastTimes[BlockCount-1]) + (fastTimes[1] - fastTimes[0]);
      coeffs.resize(3);
      coeffs[0] = omega*0.5*invh2;
      coeffs[1] = -omega*0.5*invh2;
      coeffs[2] = (xwp-xwm)*0.5*invh2;
      colIndices.resize(3);
      colIndices[0] = warpMPDEOSCOUT_+shiftp;
      colIndices[1] = warpMPDEOSCOUT_+shiftm;
      colIndices[2] = omegaGID_;
    }
    else
    {
      std::string msg;
      msg = "N_MPDE_WarpedPhaseCondition::getPhaseConditionDerivative";
      msg += " Unrecognized value for WaMPDE Phase option\n";
      N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_WARNING, msg);
    }
  }
  else 
  {
    coeffs.resize(0);
  }

#ifdef Xyce_PARALLEL_MPI
  // Now get this information to the rest of the processors.
  int numCoeffs = 0, tmpNumCoeffs = coeffs.size();
  bX.pmap()->pdsComm()->maxAll( &tmpNumCoeffs, &numCoeffs, 1 );
  std::vector<int> tmpColIndices(numCoeffs, -1);
  std::vector<double> tmpCoeffs(numCoeffs, 0.0);
  if ( warpMPDEOSCOUTLID < 0 )
  {
    coeffs.resize(numCoeffs);
    colIndices.resize(numCoeffs);
  }
  else
  {
    tmpColIndices = colIndices;
    tmpCoeffs = coeffs;
  }
  bX.pmap()->pdsComm()->sumAll( &tmpCoeffs[0], &coeffs[0], numCoeffs );
  bX.pmap()->pdsComm()->maxAll( &tmpColIndices[0], &colIndices[0], numCoeffs );
#endif 
}
