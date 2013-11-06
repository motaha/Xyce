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
// Filename      : $RCSfile: N_TIA_Assembler.h,v $
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/07
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.9.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:49 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ASSEMBER_H
#define Xyce_N_ASSEMBER_H

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>
#include <N_UTL_Xyce.h>

// ---------- Forward Declarations ----------
class N_LOA_Loader;
class N_TIA_TimeIntegrationMethod;
class N_TIA_DataStore;
class N_TIA_WorkingIntegrationMethod;
class N_LAS_Vector;
                       
class N_PDS_Manager;

class N_UTL_Timer;

//-----------------------------------------------------------------------------
// Class         : N_TIA_Assembler
// Purpose       : Base class
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 1/30/07
//-----------------------------------------------------------------------------
class N_TIA_Assembler
{
public:

  // Constructor.
  N_TIA_Assembler( N_TIA_DataStore & ds,
                   //N_ANP_AnalysisManager & ca,
                   N_LOA_Loader & loader,
                   N_TIA_WorkingIntegrationMethod & wim,
                   N_PDS_Manager & pds,
                   bool daeStateDerivFlag
                   ) ;

  // Destructor
  virtual ~N_TIA_Assembler();

  virtual bool loadRHS () = 0;
  virtual bool loadJacobian () = 0;
  virtual bool applyJacobian (const N_LAS_Vector& input, N_LAS_Vector& result) 
  { return false; };

protected:
  N_UTL_Timer * residualTimerPtr_;
  N_UTL_Timer * jacobianTimerPtr_;
  double        residualTime_;
  double        jacobianTime_;

  bool daeStateDerivFlag_;

  N_TIA_DataStore & ds_;
  N_LOA_Loader & loader_;
  N_TIA_WorkingIntegrationMethod & wim_;

  N_PDS_Manager & pdsMgr;
  
  friend class N_ANP_AnalysisInterface;
};

//-----------------------------------------------------------------------------
// Class         : N_TIA_DAE_Assembler
// Purpose       : Derrived class
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 1/30/07
//-----------------------------------------------------------------------------
class N_TIA_DAE_Assembler : public N_TIA_Assembler
{
public:

  // Default constructor.
  N_TIA_DAE_Assembler( N_TIA_DataStore & ds,
                       //N_ANP_AnalysisManager & ca,
                       N_LOA_Loader & loader,
                       N_TIA_WorkingIntegrationMethod & wim,
                       N_PDS_Manager & pds,
                       bool daeStateDerivFlag
                       );

  // Destructor
  virtual ~N_TIA_DAE_Assembler();

  // load functions, new-DAE form:
  virtual bool loadRHS ();
  virtual bool loadJacobian ();
  virtual bool applyJacobian (const N_LAS_Vector& input, N_LAS_Vector& result);

  friend class N_ANP_AnalysisInterface;

};


#endif
