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
// Filename       : $RCSfile: N_DEV_RegisterOpenDevices.C,v $
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
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2014/02/11 23:01:33 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#ifdef HAVE_DLFCN_H
#include <dlfcn.h>
#endif

#include <N_DEV_RegisterDevices.h>

#include <N_DEV_ACC.h>
#include <N_DEV_ADC.h>
#include <N_DEV_BJT.h>
#include <N_DEV_Bsrc.h>
#include <N_DEV_Capacitor.h>
#include <N_DEV_DAC.h>
#include <N_DEV_Digital.h>
#include <N_DEV_Diode.h>
#include <N_DEV_Inductor.h>
#include <N_DEV_ISRC.h>
#include <N_DEV_JFET.h>
#include <N_DEV_LTRA.h>
#include <N_DEV_MESFET.h>
#include <N_DEV_MOSFET1.h>
#include <N_DEV_MOSFET2.h>
#include <N_DEV_MOSFET3.h>
#include <N_DEV_MOSFET6.h>
#include <N_DEV_MOSFET_B3.h>
#include <N_DEV_MOSFET_B3SOI.h>
#include <N_DEV_MOSFET_B4.h>
#include <N_DEV_MutIndLin.h>
#include <N_DEV_MutIndNonLin2.h>
#include <N_DEV_MutIndNonLin.h>
#include <N_DEV_OpAmp.h>
#include <N_DEV_Resistor3.h>
#include <N_DEV_Resistor.h>
#include <N_DEV_ThermalResistor.h>
#include <N_DEV_ROM.h>
#include <N_DEV_RxnSet.h>
#include <N_DEV_SW.h>
#include <N_DEV_TRA.h>
#include <N_DEV_VCCS.h>
#include <N_DEV_Vcvs.h>
#include <N_DEV_VDMOS.h>
#include <N_DEV_Vsrc.h>
#include <N_DEV_Xygra.h>
#include <N_DEV_TransLine.h>

namespace Xyce {
namespace Device {

void
registerOpenDevices()
{
  Resistor::registerDevice();
  ThermalResistor::registerDevice();
  Resistor3::registerDevice();
  Capacitor::registerDevice();
  Inductor::registerDevice();
  Diode::registerDevice();
  BJT::registerDevice();
  JFET::registerDevice();
  MESFET::registerDevice();
  MOSFET1::registerDevice();
  MOSFET2::registerDevice();
  MOSFET3::registerDevice();
  MOSFET6::registerDevice();
  MOSFET_B3::registerDevice();
  MOSFET_B3SOI::registerDevice();
  MOSFET_B4::registerDevice();
  VDMOS::registerDevice();
  ISRC::registerDevice();
  Vcvs::registerDevice();
  Bsrc::registerDevice();
  VCCS::registerDevice();
  Vsrc::registerDevice();
  LTRA::registerDevice();
  TRA::registerDevice();
  SW::registerDevice();
  ADC::registerDevice();
  DAC::registerDevice();
  MutIndLin::registerDevice();
  MutIndNonLin::registerDevice();
  MutIndNonLin2::registerDevice();
  OpAmp::registerDevice();
  Digital::registerDevice();
  ACC::registerDevice();
  Xygra::registerDevice();
  ROM::registerDevice();
  RxnSet::registerDevice();
  TransLine::registerDevice();
}

} // namespace Device
} // namespace Xyce
