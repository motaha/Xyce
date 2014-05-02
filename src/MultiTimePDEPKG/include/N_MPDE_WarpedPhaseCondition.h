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
// Filename       : $RCSfile: N_MPDE_WarpedPhaseCondition.h,v $
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
// Revision Number: $Revision: 1.10 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef  Xyce_MPDE_WARPED_PHASE_CONDITION_H
#define  Xyce_MPDE_WARPED_PHASE_CONDITION_H

// ---------- Standard Includes ----------

#include <string>
#include <vector>

// ----------   Trilinos Includes   ----------

#include <Teuchos_RefCountPtr.hpp>

// ----------  Xyce Includes   ----------

#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockMatrix.h>
#include <N_ERH_ErrorMgr.h>
#include <N_PDS_ParMap.h>

//-----------------------------------------------------------------------------
// Class         : N_MPDE_WarpedPhaseCondition
// Purpose       :
// Special Notes :
// Creator       : Todd Coffey, SNL, 1414
// Watcher of Creator:  Heidi Thornquist, SNL, 1437
// Creation Date : 12/7/06
//-----------------------------------------------------------------------------
class N_MPDE_WarpedPhaseCondition 
{

 public:

  // Default Constructor
  N_MPDE_WarpedPhaseCondition( const int warpPhase,
                               const double warpPhaseCoeff,
                               const int warpMPDEOSCOUT,
                               const int omegaGID,
                               const int offset,
                               const int size,
                               const int phiGID,
                               const int augProcID)
  : warpPhase_(warpPhase),
    warpPhaseCoeff_(warpPhaseCoeff),
    warpMPDEOSCOUT_(warpMPDEOSCOUT),
    omegaGID_(omegaGID),
    offset_(offset),
    size_(size),
    phiGID_(phiGID),
    augProcID_(augProcID)
  {}

  // Destructor
  ~N_MPDE_WarpedPhaseCondition() {};

  void setWarpPhase( const int warpPhase )
  { 
    warpPhase_=warpPhase;
  }

  void setWarpPhaseCoeff( const int warpPhaseCoeff )
  { 
    warpPhaseCoeff_=warpPhaseCoeff;
  }

  void setWarpMPDEOSCOUT( const int warpMPDEOSCOUT )
  { 
    warpMPDEOSCOUT_=warpMPDEOSCOUT;
  }

  void setOmegaGID( const int omegaGID )
  { 
    omegaGID_=omegaGID;
  }

  void setOffset( const int offset )
  { 
    offset_=offset;
  }

  void setSize( const int size )
  { 
    size_=size;
  }

  void setPhiGID( const int phiGID )
  { 
    phiGID_=phiGID;
  }

  int getWarpPhase( )
  { return warpPhase_; }
  
  double getWarpPhaseCoeff( )
  { return warpPhaseCoeff_; }

  int getWarpMPDEOSCOUT( )
  { return warpMPDEOSCOUT_; }

  int getOmegaGID( )
  { return omegaGID_; }

  int getOffset( )
  { return offset_; }

  int getSize( )
  { return size_; }

  int getPhiGID( )
  { return phiGID_; }

  // Graph and Load functions
    Teuchos::RefCountPtr<std::vector<int> > getPhaseGraph();
  
  double getPhaseCondition(N_LAS_BlockVector* bX, std::vector<double>& fastTimes);

  // The column indices returned by this method are global ids.
  void getPhaseConditionDerivative(
      N_LAS_BlockVector* bX, 
      std::vector<double>& fastTimes,
      std::vector<int>* colIndicesPtr,
      std::vector<double>* coeffsPtr
      );

private:

  int warpPhase_;
  double warpPhaseCoeff_;
  int warpMPDEOSCOUT_;
  int omegaGID_;
  int offset_;
  int size_;
  int phiGID_;
  int augProcID_;
};

#endif // Xyce_MPDE_WARPED_PHASE_CONDITION_H


