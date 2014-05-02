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
// Filename       : $RCSfile: N_TIA_Dummy.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 5/26/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7 $
//
// Revision Date  : $Date: 2014/02/24 23:49:26 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef  _TIA_DUMMY_H
#define  _TIA_DUMMY_H

enum Vector_Tag
{
  _SOLUTION_VECTOR, // 0
  _STATE_VECTOR // 1
};

//-----------------------------------------------------------------------------
// Class         : N_TIA_TimeIntegrationAlgorithm
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------

class N_TIA_TimeIntegrationAlgorithm
{
public:

  // Constructor
  N_TIA_TimeIntegrationAlgorithm() : dt(1.0e-6), time(0.0) { };

  // Destructor
  ~N_TIA_TimeIntegrationAlgorithm() { };

  double partialTimeDerivative(const int itag)
  {
    return (1.0 / dt);
  }

  double timeDerivative(const int index, const int itag)
  {
    return (1.0 / dt);
  }

  double getTime()
  {
    return time;
  }

protected:
private :

public :
  int iIntegrationScheme;
  double dt;
  double time;

};

#endif
