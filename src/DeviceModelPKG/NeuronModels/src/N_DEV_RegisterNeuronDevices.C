//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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
// Filename       : $RCSfile: N_DEV_RegisterNeuronDevices.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : David Baur
//
// Creation Date  : 3/15/2013
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1 $
//
// Revision Date  : $Date: 2014/01/29 16:49:36 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_RegisterDevices.h>

#include <N_DEV_Neuron.h>
#include <N_DEV_Neuron2.h>
#include <N_DEV_Neuron3.h>
#include <N_DEV_Neuron4.h>
#include <N_DEV_Neuron5.h>
#include <N_DEV_Neuron6.h>
#include <N_DEV_Neuron7.h>
#include <N_DEV_Neuron8.h>
#include <N_DEV_Neuron9.h>
#include <N_DEV_NeuronPop1.h>
#include <N_DEV_Synapse.h>
#include <N_DEV_Synapse2.h>
#include <N_DEV_Synapse3.h>
#include <N_DEV_Synapse4.h>

namespace Xyce {
namespace Device {

void
registerNeuronDevices()
{
  Neuron::registerDevice();
  Neuron2::registerDevice();
  Neuron3::registerDevice();
  Neuron4::registerDevice();
  Neuron5::registerDevice();
  Neuron6::registerDevice();
  Neuron7::registerDevice();
  Neuron8::registerDevice();
  Neuron9::registerDevice();
  NeuronPop1::registerDevice();
  Synapse::registerDevice();
  Synapse2::registerDevice();
  Synapse3::registerDevice();
  Synapse4::registerDevice();
}

} // namespace Device
} // namespace Xyce
