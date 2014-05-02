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
// Filename       : $RCSfile: N_UTL_Interface_Enum_Types.h,v $
//
// Purpose        : This file contains some enumerated types which are needed 
//                  both inside and outside the expression interface.
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
// Revision Number: $$
//
// Revision Date  : $$
//
// Current Owner  : $$
//-----------------------------------------------------------------------------

#ifndef N_UTL_Interface_Enum_Types_H
#define N_UTL_Interface_Enum_Types_H

// ---------- Enumerations ----------

enum XEXP_TYPES
{
  XEXP_ALL,            //  0
  XEXP_NODE,           //  1
  XEXP_INSTANCE,       //  2
  XEXP_LEAD,           //  3
  XEXP_STRING,         //  4
  XEXP_SPECIAL,        //  5
  XEXP_VARIABLE,       //  6
  XEXP_FUNCTION,       //  7
  XEXP_NODAL_COMPUTATION, //8
  XEXP_COUNT           //  9
};

#endif // N_UTL_EXPRESSION_H
