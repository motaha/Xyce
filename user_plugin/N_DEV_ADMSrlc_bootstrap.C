//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_ADMSrlc_bootstrap.C,v $
//
// Purpose        : 
// A device must be registered with the Xyce device subsystem in order
// for models and classes of the devie to be created.  This class
// implements a call to the device registration function on
// construction and the static object is created upon shareable object
// load.
//
// Special Notes  : 
//
// Creator        : Tom Russo
//
// Creation Date  : 1/11/2012
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision $
//
// Revision Date  : $Date: 2014/03/11 18:50:03 $
//
// Current Owner  : $Author $
//-------------------------------------------------------------------------
#include <N_DEV_ADMSrlc.h>

struct Bootstrap 
{
  Bootstrap() 
  {
    Xyce::Device::ADMSrlc::registerDevice();
  }
};

Bootstrap s_bootstrap;



