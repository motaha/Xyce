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
// Revision Number: $Revision: 1.11.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:47 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


#include <N_MPDE_WarpedPhaseCondition.h>
#include <Epetra_MultiVector.h>

Teuchos::RefCountPtr<vector<int> > N_MPDE_WarpedPhaseCondition::getPhaseGraph()
{
  if ((warpMPDEOSCOUT_ == -1) && (warpPhase_ != 0))
  {
    string msg;
    msg = "N_MPDE_WarpedPhaseCondition::getPhaseGraph";
    msg += " No value for oscout which is required by specified phase equation\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  Teuchos::RefCountPtr<vector<int> > phaseGraph = Teuchos::rcp(new vector<int>);
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
    string msg;
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
  double phaseValue = 0.0;
  int xwpLID;
  int xwmLID;
  int omegaLID = bX.pmap()->globalToLocalIndex(omegaGID_);
  double omega = bX[omegaLID]; 
  if ( warpPhase_ == 0 )
  {
    phaseValue = bX[omegaLID] - 1.0;
  }
  else if ( warpPhase_ == 1 )
  {
    // Add graph for phase condition 1 for WaMPDE:  \hat{x}_w(t_1,0) - alpha = 0
    int warpMPDEOSCOUTLID_ = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_);
    phaseValue = bX[warpMPDEOSCOUTLID_] - warpPhaseCoeff_;
  }
  else if ( warpPhase_ == 2 )
  {
    // Add graph for phase condition 2 for WaMPDE:
    // \omega*((\hat{x}_w(t_1,0)-hat{x}_w(t_1,-h2))/h2 - alpha) = 0
    int shift = (BlockCount-1)*offset_;
    xwpLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_);
    xwmLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shift);
    double xwp = bX[xwpLID];
    double xwm = bX[xwmLID];
    double invh2 = fastTimes[BlockCount] - fastTimes[BlockCount-1];
    phaseValue = omega*(xwp - xwm)*invh2 - warpPhaseCoeff_ ;
  }
  else if ( warpPhase_ == 3 )
  {
    // Add graph for phase condition 2 for WaMPDE:
    // \omega*((\hat{x}_w(t_1,h2)-hat{x}_w(t_1,-h2))/(2*h2) - alpha) = 0
    int shiftp = offset_;
    int shiftm = (BlockCount-1)*offset_;
    xwpLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shiftp);
    xwmLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shiftm);
    double xwp = bX[xwpLID];
    double xwm = bX[xwmLID];
    double invh2 = (fastTimes[BlockCount] - fastTimes[BlockCount-1]) + (fastTimes[1] - fastTimes[0]);
    phaseValue = omega*(xwp - xwm)*(0.5*invh2) - warpPhaseCoeff_ ;
  }
  else
  {
    string msg;
    msg = "N_MPDE_WarpedPhaseCondition::getPhaseCondition";
    msg += " Unrecognized value for WaMPDE Phase option\n";
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_WARNING, msg);
  }
  return(phaseValue);
}

void N_MPDE_WarpedPhaseCondition::getPhaseConditionDerivative(
    N_LAS_BlockVector* bXptr, 
    std::vector<double>& fastTimes,
    vector<int>* colIndicesPtr,
    vector<double>* coeffsPtr
    )
{
  N_LAS_BlockVector & bX = *bXptr;
  vector<int> & colIndices = *colIndicesPtr;
  vector<double> & coeffs = *coeffsPtr;
  int BlockCount = bX.blockCount();
  int xwpLID;
  int xwmLID;
  int omegaLID = bX.pmap()->globalToLocalIndex(omegaGID_);
  double omega = bX[omegaLID]; 
  int warpMPDEOSCOUTLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_);
  if ( warpPhase_ == 0 )
  {
    // put derivative wrt x_v of (omega-1) = 0 into omegaGID_,x_v locations
    coeffs.resize(1);
    coeffs[0] = 1.0;
    colIndices.resize(1);
    colIndices[0] = omegaLID;
  }
  else if ( warpPhase_ == 1 )
  {
    // put derivative wrt x_v of (x_v-warpMPDECoeff) = 0 into omegaGID_,x_v location
    coeffs.resize(1);
    coeffs[0] = 1.0;
    colIndices.resize(1);
    colIndices[0] = warpMPDEOSCOUTLID;
  }
  else if ( warpPhase_ == 2 )
  {
    // put derivative wrt x_v, xv+shift, and omega of 
    // \omega*((\hat{x}_w(t_1,0)-\hat{x}_w(t_1,-h2))/h2 - alpha) 
    //   = omega/h2 into  omegaGID_,x_v location
    //   = -omega/h2 into omegaGID_,x_v+shift location
    //   = (\hat{x}_w(t_1,0)-\hat{x}_w(t_1,-h2))/h2 into omegaGID_,omegaGID_ location
    int shift = (BlockCount-1)*offset_;
    xwpLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_);
    xwmLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shift);
    double xwp = bX[xwpLID];
    double xwm = bX[xwmLID];
    double invh2 = fastTimes[BlockCount] - fastTimes[BlockCount-1];
    coeffs.resize(3);
    coeffs[0] =  omega*invh2;
    coeffs[1] = -omega*invh2;
    coeffs[2] = (xwp - xwm)*invh2;
    colIndices.resize(3);
    colIndices[0] = xwpLID;
    colIndices[1] = xwmLID;
    colIndices[2] = omegaLID;
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
    xwpLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shiftp);
    xwmLID = bX.pmap()->globalToLocalIndex(warpMPDEOSCOUT_+shiftm);
    double xwp = bX[xwpLID];
    double xwm = bX[xwmLID];
    double invh2 = (fastTimes[BlockCount] - fastTimes[BlockCount-1]) + (fastTimes[1] - fastTimes[0]);
    coeffs.resize(3);
    coeffs[0] = omega*0.5*invh2;
    coeffs[1] = -omega*0.5*invh2;
    coeffs[2] = (xwp-xwm)*0.5*invh2;
    colIndices.resize(3);
    colIndices[0] = xwpLID;
    colIndices[1] = xwmLID;
    colIndices[2] = omegaLID;
  }
  else
  {
    string msg;
    msg = "N_MPDE_WarpedPhaseCondition::getPhaseConditionDerivative";
    msg += " Unrecognized value for WaMPDE Phase option\n";
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_WARNING, msg);
  }
}
