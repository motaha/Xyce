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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_ParamMgr.h,v $
//
// Purpose        : 
//
// Special Notes  :
//
// Creator        : Eric R. Keiter
//
// Creation Date  : 10/25/02
//
// Revision Information:
// ----------------------
//
// Revision Number: $Revision: 1.14 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_ParamMgr_h
#define Xyce_N_NLS_ParamMgr_h

#include <vector>

#include <N_IO_fwd.h>

#include <N_NLS_NLParams.h>


//-----------------------------------------------------------------------------
// Class         : N_NLS_ParamMgr
//
// Purpose       : It is not unusual for a Xyce simulation to have several
//                 different parameter sets.  Each parameter set (for dcop,
//                 transient, etc.) is stored in a N_NLS_NLParams data
//                 structure.
//
//                 The management of which parameter set is currently being
//                 used is handled by this class.
//
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/26/02
//-----------------------------------------------------------------------------
class N_NLS_ParamMgr
{
public:
  N_NLS_ParamMgr(N_IO_CmdParse & cp);
  ~N_NLS_ParamMgr();

  bool addParameterSet (AnalysisMode mode, N_NLS_NLParams & nlp);

  bool getParams (AnalysisMode mode, N_NLS_NLParams & nlp); 

  bool getCurrentParams (N_NLS_NLParams & nlp); 
  inline bool setAnalysisMode (AnalysisMode mode);
  inline void resetAnalysisMode();
  inline AnalysisMode getAnalysisMode() const;

protected:
private:

public:

protected:

private:
  std::vector<N_NLS_NLParams> paramVector_;
  AnalysisMode   currentMode_;

  bool modeToggled_;
  bool gcp_calledBefore_;
  bool paramsChanged_;
};

//-----------------------------------------------------------------------------
// Function      : N_NLS_ParamMgr::setAnalysisMode
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/26/02
//-----------------------------------------------------------------------------
inline bool N_NLS_ParamMgr::setAnalysisMode (AnalysisMode mode)
{
  if (currentMode_ != mode) 
  {
    paramsChanged_ = true;
    currentMode_   = mode;
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ParamMgr::resetAnalysisMode
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/26/02
//-----------------------------------------------------------------------------
inline void N_NLS_ParamMgr::resetAnalysisMode()
{
  if (currentMode_ != DC_OP) 
  {
    paramsChanged_ = true;
    currentMode_   = DC_OP;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ParamMgr::getAnalysisMode
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/26/02
//-----------------------------------------------------------------------------
inline AnalysisMode N_NLS_ParamMgr::getAnalysisMode() const
{
  return currentMode_;
}

#endif

