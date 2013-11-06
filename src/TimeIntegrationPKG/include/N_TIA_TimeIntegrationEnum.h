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
// Filename      : $RCSfile: N_TIA_TimeIntegrationEnum.h,v $
//
// Purpose       : This file defines the classes for the time integration
//                 methods --- the "interface base class" along with the
//                 accompanying integration methods classes which can be
//                 used in the time integration algorithm.
//
// Special Notes :
//
// Creator       : Eric Keiter
//
// Creation Date : 2/11/07
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.4.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:49 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_TIME_INTEG_enum_H
#define Xyce_N_TIA_TIME_INTEG_enum_H

enum Integration_Method_Mode
{ // Time Integrator Method Mode
  TIAMethod_NONE, // 0
  TIAMethod_BACKWARD_EULER, // 1
  TIAMethod_BACKWARD_DIFFERENTIATION_2, // 2
  TIAMethod_TRAPEZOIDAL, // 3
  TIAMethod_VARIABLE_THETA, // 4
  TIAMethod_A_CONTRACTIVE_2, // 5
  TIAMethod_BACKWARD_DIFFERENTIATION_15,// 6
  TIAMethod_ONESTEP, // 7
  TIAMethod_GEAR_12 // 8
};

#endif // Xyce_N_TIA_TIME_INTEG_enum_H

