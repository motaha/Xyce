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
// Filename       : $RCSfile: N_IO_OutputRaw.h,v $
//
// Purpose        : Base class for handling file output of simulation results.
//
// Special Notes  :
//
// Creator        : Richard Schiek, Electrical Systems Modeling, Sandia National Laboratories
//
// Creation Date  : 12/06/12
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputRaw_h
#define Xyce_N_IO_OutputRaw_h

// ----------   Standard Includes   ----------

// ----------   Xyce Includes   ----------
#include<N_IO_OutputFileBase.h>

// ---------- Forward Declarations ----------

class N_IO_OutputRaw : public N_IO_OutputFileBase
{
  public:
  N_IO_OutputRaw();
  ~N_IO_OutputRaw();

};

#endif // Xyce_N_IO_OutputRaw_h
