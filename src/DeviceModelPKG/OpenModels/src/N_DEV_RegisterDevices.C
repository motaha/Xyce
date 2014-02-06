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
// Filename       : $RCSfile: N_DEV_RegisterDevices.C,v $
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
// Revision Number: $Revision: 1.6.2.4 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>
#include <N_DEV_Factory.h>

#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_Param.h> 
#include <N_DEV_Const.h>

#ifdef Xyce_RAD_MODELS
#include <N_DEV_ExtendedModelTypes.h>
#else
#include <N_DEV_ModelTypes.h>
#endif

#include <N_DEV_2DPDE.h>
#include <N_DEV_ACC.h>
#include <N_DEV_ADC.h>
#include <N_DEV_BJT.h>
#include <N_DEV_Bsrc.h>
#include <N_DEV_Capacitor.h>
#include <N_DEV_DAC.h>
#include <N_DEV_Digital.h>
#include <N_DEV_Diode.h>
#include <N_DEV_DiodePDE.h>
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
#include <N_DEV_NewDevice.h>
#include <N_DEV_OpAmp.h>
#include <N_DEV_Resistor3.h>
#include <N_DEV_Resistor.h>
#include <N_DEV_ROM.h>
#include <N_DEV_RxnSet.h>
#include <N_DEV_SW.h>
#include <N_DEV_ThermalResistor.h>
#include <N_DEV_TRA.h>
#include <N_DEV_VCCS.h>
#include <N_DEV_Vcvs.h>
#include <N_DEV_VDMOS.h>
#include <N_DEV_Vsrc.h>
#include <N_DEV_Xygra.h>

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

#include <N_DEV_ExternDevice.h>

#include <N_DEV_ADMSHBT_X.h>
#include <N_DEV_ADMSPSP103VA.h>
#include <N_DEV_ADMSvbic.h>

namespace Xyce {
namespace Device {

/**
 * These functions define the factories used to create devices.  The
 * DeviceBuilder::createDeviceByNetlistDeviceType() and DeviceBuilder::createDeviceByModelType()
 * functions create device models, device instance and device masters
 * 
 */

Device *ACCfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Accelerated Object Device");
  const string className("ACC");
  const string defaultModelName("ACC level 1");

  Device *devptr =
    new DeviceTemplate<ACC::Model, ACC::Instance>
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *ADCfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("ADC");
  const string className("ADC");
  const string defaultModelName("YADC level 1 (Analog to Digital Interface)");

  Device *devptr =
    new ADC::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *BJTfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Bipolar Junction Transistor");
  const string className("BJT");
  string defaultModelName ("Q level 1");

  Device *devptr =
    new BJT::Master
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Bsrcfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name ("Expression Based Voltage or Current Source");
  const string className("Bsrc");
  const string defaultModelName ("B level 1");

  pair <string,double> p(string("I"), 0.0);

  Device *devptr =
    new Bsrc::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Capacitorfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Capacitor");
  const string className("Capacitor");
  const string defaultModelName("C level 1");

  Device *devptr =
    new Capacitor::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *DACfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("DAC");
  const string className("DAC");
  const string defaultModelName("YDAC level 1 (Digital to Analog Interface)");

  Device *devptr =
    new DAC::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Digitalfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name ("Behavioral Digital");
  const string className ("Digital");
  const string defaultModelName ("Digital level 1");

  Device *devptr =
    new DeviceTemplate<Digital::Model, Digital::Instance>
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Diodefactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Diode");
  const string className("Diode");
  const string defaultModelName("D level 1,2");

  Device *devptr =
    new Diode::Master
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Inductorfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Inductor");
  const string className("Inductor");
  const string defaultModelName("L level 1");

  Device *devptr =
    new Inductor::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *ISRCfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name ("Independent Current Source");
  const string className("ISRC");
  const string defaultModelName ("I level 1");

  Device *devptr =
    new ISRC::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *JFETfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("JFET");
  const string className("JFET");
  const string defaultModelName("J leve1 1,2");

  Device *devptr =
    new JFET::Master
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *LTRAfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Lossy Transmission Line");
  const string className("LTRA");
  const string defaultModelName ("O level 1");

  Device *devptr =
    new LTRA::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *MESFETfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("MESFET");
  const string className("MESFET");
  const string defaultModelName("Z level 1");

  Device *devptr =
    new MESFET::Master
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *MOSFET1factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("MOSFET level 1");
  const string className("MOSFET1");
  const string defaultModelName("M level 1");

  Device *devptr =
    new MOSFET1::Master
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *MOSFET2factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("MOSFET level 2");
  const string className("MOSFET2");
  const string defaultModelName("M level 2");

  Device *devptr =
    new MOSFET2::Master
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *MOSFET3factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("MOSFET level 3");
  const string className("MOSFET3");
  const string defaultModelName("M level 3");

  Device *devptr =
    new MOSFET3::Master
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *MOSFET6factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("MOSFET level 6");
  const string className("MOSFET6");
  const string defaultModelName("M level 6");

  Device *devptr =
    new MOSFET6::Master
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *MOSFET_B3factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name( "BSIM3" );
  const string className("MOSFET_B3");
  const string defaultModelName("M level 9");

  Device *devptr =
    new MOSFET_B3::Master
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *MOSFET_B3SOIfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name ("BSIM3 SOI");
  const string className("MOSFET_B3SOI");
  const string defaultModelName ("M level 10");

  Device *devptr =
    new MOSFET_B3SOI::Master
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options,
      4);       // number of external variables - matters for soi.      );

  return devptr;
}

Device *MOSFET_B4factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name( "BSIM4" );
  const string className("MOSFET_B4");
  const string defaultModelName("M level 14");

  Device *devptr =
    new MOSFET_B4::Master
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *MutIndLinfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Linear Mutual Inductor");
  const string className("MutIndLin");
  const string defaultModelName("K level 1");

  Device *devptr =
    new MutIndLin::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *MutIndNonLin2factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Nonlinear Mutual Inductor");
  const string className("MutIndNonLin2");
  const string defaultModelName("K level 2");

  Device *devptr =
    new DeviceTemplate<MutIndNonLin2::Model, MutIndNonLin2::Instance>
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *MutIndNonLinfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Nonlinear Mutual Inductor");
  const string className("MutIndNonLin");
  const string defaultModelName("K level 1");

  Device *devptr =
    new MutIndNonLin::Master
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *NewDevicefactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("NewDevice");
  const string className("NewDevice");
  const string defaultModelName("NewDevice level 1");

  Device *devptr =
    new DeviceTemplate<NewDevice::Model, NewDevice::Instance>
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *OpAmpfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Operational Amplifier");
  const string className("OpAmp");
  const string defaultModelName("OpAmp level 1");

  Device *devptr =
    new DeviceTemplate<OpAmp::Model, OpAmp::Instance>
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Resistor3factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name ("Resistor");
  const string className("Resistor3");
  const string defaultModelName("R level 3");

  Device *devptr =
    new Resistor3::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Resistorfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Resistor");
  const string className("N_DEV_Resistor");
  const string defaultModelName("R level 1");

  Device *devptr =
    new Resistor::Master(name, className, defaultModelName, Device::LINEAR_DEVICE, solver_state, device_options);

  return devptr;
}

Device *ROMfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("ROM");
  const string className("ROM");
  const string defaultModelName("ROM level 1");

  Device *devptr =
    new ROM::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *RxnSetfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Rxn Effects Device");
  const string className("RxnSet");
  const string defaultModelName ("YRXN level 1 (Rxn Device)");

  Device *devptr =
    new DeviceTemplate<RxnSet::Model, RxnSet::Instance>
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *SWfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Controlled Switch");
  const string className("SW");
  const string defaultModelName("S level 1");

  Device *devptr =
    new SW::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *ThermalResistorfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Resistor");
  const string className("ThermalResistor");
  const string defaultModelName("R level 2");

  Device *devptr =
    new ThermalResistor::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *TRAfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Ideal Transmission Line");
  const string className("TRA");
  const string defaultModelName ("T level 1");

  Device *devptr =
    new DeviceTemplate<TRA::Model, TRA::Instance>
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *VCCSfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name ("Linear Voltage Controlled Current Source");
  const string className("VCCS");
  const string defaultModelName ("G level 1");

  Device *devptr =
    new VCCS::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Vcvsfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name ("Linear Voltage Controlled Voltage Source");
  const string className("Vcvs");
  const string defaultModelName ("E level 1");

  Device *devptr =
    new Vcvs::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *VDMOSfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name( "Power MOSFET");
  const string className("VDMOS");
  const string defaultModelName("M level 18");

  Device *devptr =
    new VDMOS::Master
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Vsrcfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name ("Independent Voltage Source");
  const string className("Vsrc");
  const string defaultModelName("V level 1");

  Device *devptr =
    new Vsrc::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Xygrafactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Xygra");
  const string className("Xygra");
  const string defaultModelName("Xygra level 1");

  Device *devptr =
    new DeviceTemplate<Xygra::Model, Xygra::Instance>
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *DiodePDEfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("1D PDE Device");
  const string className("N_DEV_DiodePDE");
  const string defaultModelName ("PDE level 3");

  N_DEV_Device *devptr =
    new Xyce::Device::DeviceTemplate<DiodePDE::Model, DiodePDE::Instance>
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *TwoDPDEfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("2D PDE Device");
  const string className("N_DEV_2DPDE");
  const string defaultModelName ("PDE level 2");

  N_DEV_Device *devptr =
    new DeviceTemplate<TwoDPDE::Model, TwoDPDE::Instance>
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Neuron2factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Neuron");
  const string className("Neuron2");
  const string defaultModelName("YNEURON level 2");

  Device *devptr =
    new Xyce::Device::DeviceTemplate<Neuron2::Model, Neuron2::Instance>
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Neuron3factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Neuron");
  const string className("Neuron3");
  const string defaultModelName("YNEURON level 3");

  Device *devptr =
    new Xyce::Device::DeviceTemplate<Neuron3::Model, Neuron3::Instance>
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Neuron4factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Neuron");
  const string className("Neuron4");
  const string defaultModelName("YNEURON level 4");

  Device *devptr =
    new Xyce::Device::DeviceTemplate<Neuron4::Model, Neuron4::Instance>
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Neuron5factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Neuron");
  const string className("Neuron5");
  const string defaultModelName("YNEURON level 5");

  Device *devptr =
    new Xyce::Device::DeviceTemplate<Neuron5::Model, Neuron5::Instance>
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Neuron6factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Neuron");
  const string className("Neuron6");
  const string defaultModelName("YNEURON level 6");

  Device *devptr =
    new Xyce::Device::DeviceTemplate<Neuron6::Model, Neuron6::Instance>
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Neuron7factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Neuron");
  const string className("Neuron7");
  const string defaultModelName("YNEURON level 7");

  Device *devptr =
    new Xyce::Device::DeviceTemplate<Neuron7::Model, Neuron7::Instance>
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Neuron8factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Neuron");
  const string className("Neuron8");
  const string defaultModelName("YNEURON level 8");

  Device *devptr =
    new Xyce::Device::DeviceTemplate<Neuron8::Model, Neuron8::Instance>
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Neuron9factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Neuron");
  const string className("Neuron9");
  const string defaultModelName("YNEURON level 9");

  Device *devptr =
    new Neuron9::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Neuronfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Neuron");
  const string className("Neuron");
  const string defaultModelName("YNEURON level 1");

  Device *devptr =
    new Neuron::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *NeuronPop1factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("NeuronPopulation");
  const string className("NeuronPop1");
  const string defaultModelName("YNEURONPOP level 1");

  Device *devptr =
    new Xyce::Device::DeviceTemplate<NeuronPop1::Model, NeuronPop1::Instance>
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Synapse2factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Synapse");
  const string className("Synapse2");
  const string defaultModelName("YSYNAPSE level 2");

  Device *devptr =
    new Synapse2::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Synapse3factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Synapse, Clopath");
  const string className("Synapse3");
  const string defaultModelName ("YSYNAPSE level 3");

  Device *devptr =
    new Xyce::Device::DeviceTemplate<Synapse3::Model, Synapse3::Instance>
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Synapse4factory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Synapse");
  const string className("Synapse4");
  const string defaultModelName("YSYNAPSE level 4");

  Device *devptr =
    new Synapse4::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *Synapsefactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("Synapse");
  const string className("Synapse");
  const string defaultModelName("YSYNAPSE level 1");

  Device *devptr =
    new Synapse::Master
    ( name,
      className,
      defaultModelName,
      Device::LINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *ExternDevicefactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("External Device");
  const string className("N_DEV_ExternDevice");
  const string defaultModelName ("YEXT level 1 (External Device)");

  N_DEV_Device *devptr =
    new Xyce::Device::DeviceTemplate<ExternDevice::Model, ExternDevice::Instance>
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *ADMSHBT_Xfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("ADMS HBT_X");
  const string className("ADMSHBT_X");
  const string defaultModelName("Q level 23");

  Device *devptr =
    new DeviceTemplate<ADMSHBT_X::Model, ADMSHBT_X::Instance>
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *ADMSPSP103VAfactory(SolverState & solver_state, DeviceOptions & device_options)
{
  const string name("PSP103VA MOSFET");
  const string className("ADMSPSP103VA");
  const string defaultModelName("M level 103");

  Device *devptr =
    new DeviceTemplate<ADMSPSP103VA::Model, ADMSPSP103VA::Instance>
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

Device *ADMSvbicfactory(SolverState & solver_state,DeviceOptions & device_options)
{
  const string name("ADMS vbic");
  const string className("ADMSvbic");
  const string defaultModelName("Q level 10");

  Device *devptr =
    new DeviceTemplate<ADMSvbic::Model, ADMSvbic::Instance>
    ( name,
      className,
      defaultModelName,
      Device::NONLINEAR_DEVICE,
      solver_state,device_options);

  return devptr;
}

void
registerDevices()
{
// Register Model parameters
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("r", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Resistor::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("r", 2), reinterpret_cast<ParametricData<void> &(*)()>(&ThermalResistor::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("r", 3), reinterpret_cast<ParametricData<void> &(*)()>(&Resistor3::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("c", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Capacitor::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("l", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Inductor::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("d", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Diode::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("d", 2), reinterpret_cast<ParametricData<void> &(*)()>(&Diode::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("q", 1), reinterpret_cast<ParametricData<void> &(*)()>(&BJT::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("q", 10), reinterpret_cast<ParametricData<void> &(*)()>(&ADMSvbic::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("q", 23), reinterpret_cast<ParametricData<void> &(*)()>(&ADMSHBT_X::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("j", 1), reinterpret_cast<ParametricData<void> &(*)()>(&JFET::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("j", 2), reinterpret_cast<ParametricData<void> &(*)()>(&JFET::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("z", 1), reinterpret_cast<ParametricData<void> &(*)()>(&MESFET::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("m", 1), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET1::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("m", 2), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET2::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("m", 3), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET3::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("m", 6), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET6::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("m", 9), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET_B3::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("m", 49), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET_B3::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("m", 10), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET_B3SOI::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("m", 57), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET_B3SOI::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("m", 14), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET_B4::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("m", 54), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET_B4::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("m", 18), reinterpret_cast<ParametricData<void> &(*)()>(&VDMOS::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("m", 103), reinterpret_cast<ParametricData<void> &(*)()>(&ADMSPSP103VA::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("i", 1), reinterpret_cast<ParametricData<void> &(*)()>(&ISRC::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("e", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Vcvs::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("f", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Bsrc::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("g", 1), reinterpret_cast<ParametricData<void> &(*)()>(&VCCS::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("h", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Bsrc::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("v", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Vsrc::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("b", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Bsrc::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("o", 1), reinterpret_cast<ParametricData<void> &(*)()>(&LTRA::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("t", 1), reinterpret_cast<ParametricData<void> &(*)()>(&TRA::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("s", 1), reinterpret_cast<ParametricData<void> &(*)()>(&SW::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("pde", 1), reinterpret_cast<ParametricData<void> &(*)()>(&DiodePDE::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("pde", 2), reinterpret_cast<ParametricData<void> &(*)()>(&TwoDPDE::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("ext", 1), reinterpret_cast<ParametricData<void> &(*)()>(&ExternDevice::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("adc", 1), reinterpret_cast<ParametricData<void> &(*)()>(&ADC::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("dac", 1), reinterpret_cast<ParametricData<void> &(*)()>(&DAC::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("mil", 1), reinterpret_cast<ParametricData<void> &(*)()>(&MutIndLin::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("min", 1), reinterpret_cast<ParametricData<void> &(*)()>(&MutIndNonLin::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("min", 2), reinterpret_cast<ParametricData<void> &(*)()>(&MutIndNonLin2::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("opamp", 1), reinterpret_cast<ParametricData<void> &(*)()>(&OpAmp::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("not", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("and", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("nand", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("or", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("nor", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("add", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("xor", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("nxor", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("dff", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("acc", 1), reinterpret_cast<ParametricData<void> &(*)()>(&ACC::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("neuron", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("neuron", 2), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron2::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("neuron", 3), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron3::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("neuron", 4), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron4::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("neuron", 5), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron5::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("neuron", 6), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron6::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("neuron", 7), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron7::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("neuron", 8), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron8::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("neuron", 9), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron9::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("synapse", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Synapse::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("synapse", 2), reinterpret_cast<ParametricData<void> &(*)()>(&Synapse2::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("synapse", 3), reinterpret_cast<ParametricData<void> &(*)()>(&Synapse3::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("synapse", 4), reinterpret_cast<ParametricData<void> &(*)()>(&Synapse4::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("neuronpop", 1), reinterpret_cast<ParametricData<void> &(*)()>(&NeuronPop1::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("newd", 1), reinterpret_cast<ParametricData<void> &(*)()>(&NewDevice::Model::getParametricData));

  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("xygra", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Xygra::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("rom", 1), reinterpret_cast<ParametricData<void> &(*)()>(&ROM::Model::getParametricData));
  Tom::Register(getXyceModelRegistry(), DeviceLevelKey("rxn", 1), reinterpret_cast<ParametricData<void> &(*)()>(&RxnSet::Model::getParametricData));

// Register Instance parameters
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("r", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Resistor::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("r", 2), reinterpret_cast<ParametricData<void> &(*)()>(&ThermalResistor::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("r", 3), reinterpret_cast<ParametricData<void> &(*)()>(&Resistor3::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("c", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Capacitor::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("l", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Inductor::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("d", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Diode::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("d", 2), reinterpret_cast<ParametricData<void> &(*)()>(&Diode::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("q", 1), reinterpret_cast<ParametricData<void> &(*)()>(&BJT::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("q", 10), reinterpret_cast<ParametricData<void> &(*)()>(&ADMSvbic::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("q", 23), reinterpret_cast<ParametricData<void> &(*)()>(&ADMSHBT_X::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("j", 1), reinterpret_cast<ParametricData<void> &(*)()>(&JFET::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("j", 2), reinterpret_cast<ParametricData<void> &(*)()>(&JFET::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("z", 1), reinterpret_cast<ParametricData<void> &(*)()>(&MESFET::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("m", 1), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET1::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("m", 2), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET2::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("m", 3), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET3::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("m", 6), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET6::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("m", 9), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET_B3::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("m", 49), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET_B3::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("m", 10), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET_B3SOI::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("m", 57), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET_B3SOI::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("m", 14), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET_B4::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("m", 54), reinterpret_cast<ParametricData<void> &(*)()>(&MOSFET_B4::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("m", 18), reinterpret_cast<ParametricData<void> &(*)()>(&VDMOS::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("m", 103), reinterpret_cast<ParametricData<void> &(*)()>(&ADMSPSP103VA::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("i", 1), reinterpret_cast<ParametricData<void> &(*)()>(&ISRC::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("e", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Vcvs::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("f", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Bsrc::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("g", 1), reinterpret_cast<ParametricData<void> &(*)()>(&VCCS::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("h", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Bsrc::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("v", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Vsrc::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("b", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Bsrc::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("o", 1), reinterpret_cast<ParametricData<void> &(*)()>(&LTRA::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("t", 1), reinterpret_cast<ParametricData<void> &(*)()>(&TRA::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("s", 1), reinterpret_cast<ParametricData<void> &(*)()>(&SW::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("pde", 1), reinterpret_cast<ParametricData<void> &(*)()>(&DiodePDE::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("pde", 2), reinterpret_cast<ParametricData<void> &(*)()>(&TwoDPDE::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("ext", 1), reinterpret_cast<ParametricData<void> &(*)()>(&ExternDevice::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("adc", 1), reinterpret_cast<ParametricData<void> &(*)()>(&ADC::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("dac", 1), reinterpret_cast<ParametricData<void> &(*)()>(&DAC::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("mil", 1), reinterpret_cast<ParametricData<void> &(*)()>(&MutIndLin::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("min", 1), reinterpret_cast<ParametricData<void> &(*)()>(&MutIndNonLin::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("min", 2), reinterpret_cast<ParametricData<void> &(*)()>(&MutIndNonLin2::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("opamp", 1), reinterpret_cast<ParametricData<void> &(*)()>(&OpAmp::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("not", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("and", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("nand", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("or", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("nor", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("add", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("xor", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("nxor", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("dff", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Digital::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("acc", 1), reinterpret_cast<ParametricData<void> &(*)()>(&ACC::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("neuron", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("neuron", 2), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron2::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("neuron", 3), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron3::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("neuron", 4), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron4::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("neuron", 5), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron5::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("neuron", 6), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron6::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("neuron", 7), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron7::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("neuron", 8), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron8::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("neuron", 9), reinterpret_cast<ParametricData<void> &(*)()>(&Neuron9::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("synapse", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Synapse::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("synapse", 2), reinterpret_cast<ParametricData<void> &(*)()>(&Synapse2::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("synapse", 3), reinterpret_cast<ParametricData<void> &(*)()>(&Synapse3::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("synapse", 4), reinterpret_cast<ParametricData<void> &(*)()>(&Synapse4::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("neuronpop", 1), reinterpret_cast<ParametricData<void> &(*)()>(&NeuronPop1::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("newd", 1), reinterpret_cast<ParametricData<void> &(*)()>(&NewDevice::Instance::getParametricData));

  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("xygra", 1), reinterpret_cast<ParametricData<void> &(*)()>(&Xygra::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("rom", 1), reinterpret_cast<ParametricData<void> &(*)()>(&ROM::Instance::getParametricData));
  Tom::Register(getXyceInstanceRegistry(), DeviceLevelKey("rxn", 1), reinterpret_cast<ParametricData<void> &(*)()>(&RxnSet::Instance::getParametricData));


  Bob::Register<Resistor::Master>(getXyceRegistry2(), DeviceLevelKey("r", 1), Resistorfactory);
  Bob::Register<ThermalResistor::Master>(getXyceRegistry2(), DeviceLevelKey("r", 2), ThermalResistorfactory);
  Bob::Register<Resistor3::Master>(getXyceRegistry2(), DeviceLevelKey("r", 3), Resistor3factory);

  Bob::Register<Capacitor::Master>(getXyceRegistry2(), DeviceLevelKey("c", 1), Capacitorfactory);

  Bob::Register<Inductor::Master>(getXyceRegistry2(), DeviceLevelKey("l", 1), Inductorfactory);

  Bob::Register<Diode::Master>(getXyceRegistry2(), DeviceLevelKey("d", 1), Diodefactory);
  Bob::Register<Diode::Master>(getXyceRegistry2(), DeviceLevelKey("d", 2), Diodefactory);

  Bob::Register<BJT::Master>(getXyceRegistry2(), DeviceLevelKey("q", 1), BJTfactory);
  Bob::Register<DeviceTemplate<ADMSvbic::Model, ADMSvbic::Instance> >(getXyceRegistry2(), DeviceLevelKey("q", 10), ADMSvbicfactory);
  Bob::Register<DeviceTemplate<ADMSHBT_X::Model, ADMSHBT_X::Instance> >(getXyceRegistry2(), DeviceLevelKey("q", 23), ADMSHBT_Xfactory);

  Bob::Register<JFET::Master>(getXyceRegistry2(), DeviceLevelKey("j", 1), JFETfactory);
  Bob::Register<JFET::Master>(getXyceRegistry2(), DeviceLevelKey("j", 2), JFETfactory);

  Bob::Register<MESFET::Master>(getXyceRegistry2(), DeviceLevelKey("z", 1), MESFETfactory);

  Bob::Register<MOSFET1::Master>(getXyceRegistry2(), DeviceLevelKey("m", 1), MOSFET1factory);
  Bob::Register<MOSFET2::Master>(getXyceRegistry2(), DeviceLevelKey("m", 2), MOSFET2factory);
  Bob::Register<MOSFET3::Master>(getXyceRegistry2(), DeviceLevelKey("m", 3), MOSFET3factory);
  Bob::Register<MOSFET6::Master>(getXyceRegistry2(), DeviceLevelKey("m", 6), MOSFET6factory);
  Bob::Register<MOSFET_B3::Master>(getXyceRegistry2(), DeviceLevelKey("m", 9), MOSFET_B3factory);
  Bob::Register<MOSFET_B3::Master>(getXyceRegistry2(), DeviceLevelKey("m", 49), MOSFET_B3factory);
  Bob::Register<MOSFET_B3SOI::Master>(getXyceRegistry2(), DeviceLevelKey("m", 10), MOSFET_B3SOIfactory);
  Bob::Register<MOSFET_B3SOI::Master>(getXyceRegistry2(), DeviceLevelKey("m", 57), MOSFET_B3SOIfactory);
  Bob::Register<MOSFET_B4::Master>(getXyceRegistry2(), DeviceLevelKey("m", 14), MOSFET_B4factory);
  Bob::Register<MOSFET_B4::Master>(getXyceRegistry2(), DeviceLevelKey("m", 54), MOSFET_B4factory);
  Bob::Register<VDMOS::Master>(getXyceRegistry2(), DeviceLevelKey("m", 18), VDMOSfactory);
  Bob::Register<DeviceTemplate<ADMSPSP103VA::Model, ADMSPSP103VA::Instance> >(getXyceRegistry2(), DeviceLevelKey("m", 103), ADMSPSP103VAfactory);

  Bob::Register<ISRC::Master>(getXyceRegistry2(), DeviceLevelKey("i", 1), ISRCfactory);

  Bob::Register<Vcvs::Master>(getXyceRegistry2(), DeviceLevelKey("e", 1), Vcvsfactory);

  Bob::Register<Bsrc::Master>(getXyceRegistry2(), DeviceLevelKey("f", 1), Bsrcfactory);

  Bob::Register<VCCS::Master>(getXyceRegistry2(), DeviceLevelKey("g", 1), VCCSfactory);

  Bob::Register<Bsrc::Master>(getXyceRegistry2(), DeviceLevelKey("h", 1), Bsrcfactory);

  Bob::Register<Vsrc::Master>(getXyceRegistry2(), DeviceLevelKey("v", 1), Vsrcfactory);

  Bob::Register<Bsrc::Master>(getXyceRegistry2(), DeviceLevelKey("b", 1), Bsrcfactory);

  Bob::Register<LTRA::Master>(getXyceRegistry2(), DeviceLevelKey("o", 1), LTRAfactory);

  Bob::Register<DeviceTemplate<TRA::Model, TRA::Instance> >(getXyceRegistry2(), DeviceLevelKey("t", 1), TRAfactory);

  Bob::Register<SW::Master>(getXyceRegistry2(), DeviceLevelKey("s", 1), SWfactory);

  Bob::Register<DeviceTemplate<DiodePDE::Model, DiodePDE::Instance> >(getXyceRegistry2(), DeviceLevelKey("pde", 1), DiodePDEfactory);
  Bob::Register<DeviceTemplate<TwoDPDE::Model, TwoDPDE::Instance> >(getXyceRegistry2(), DeviceLevelKey("pde", 2), TwoDPDEfactory);

  Bob::Register<DeviceTemplate<ExternDevice::Model, ExternDevice::Instance> >(getXyceRegistry2(), DeviceLevelKey("ext", 1), ExternDevicefactory);

  Bob::Register<ADC::Master>(getXyceRegistry2(), DeviceLevelKey("adc", 1), ADCfactory);

  Bob::Register<DAC::Master>(getXyceRegistry2(), DeviceLevelKey("dac", 1), DACfactory);

  Bob::Register<MutIndLin::Master>(getXyceRegistry2(), DeviceLevelKey("mil", 1), MutIndLinfactory);

  Bob::Register<MutIndNonLin::Master>(getXyceRegistry2(), DeviceLevelKey("min", 1), MutIndNonLinfactory);
  Bob::Register<DeviceTemplate<MutIndNonLin2::Model, MutIndNonLin2::Instance> >(getXyceRegistry2(), DeviceLevelKey("min", 2), MutIndNonLin2factory);

  Bob::Register<DeviceTemplate<OpAmp::Model, OpAmp::Instance> >(getXyceRegistry2(), DeviceLevelKey("opamp", 1), OpAmpfactory);

  Bob::Register<DeviceTemplate<Digital::Model, Digital::Instance> >(getXyceRegistry2(), DeviceLevelKey("not", 1), Digitalfactory);
  Bob::Register<DeviceTemplate<Digital::Model, Digital::Instance> >(getXyceRegistry2(), DeviceLevelKey("and", 1), Digitalfactory);
  Bob::Register<DeviceTemplate<Digital::Model, Digital::Instance> >(getXyceRegistry2(), DeviceLevelKey("nand", 1), Digitalfactory);
  Bob::Register<DeviceTemplate<Digital::Model, Digital::Instance> >(getXyceRegistry2(), DeviceLevelKey("or", 1), Digitalfactory);
  Bob::Register<DeviceTemplate<Digital::Model, Digital::Instance> >(getXyceRegistry2(), DeviceLevelKey("nor", 1), Digitalfactory);
  Bob::Register<DeviceTemplate<Digital::Model, Digital::Instance> >(getXyceRegistry2(), DeviceLevelKey("add", 1), Digitalfactory);
  Bob::Register<DeviceTemplate<Digital::Model, Digital::Instance> >(getXyceRegistry2(), DeviceLevelKey("xor", 1), Digitalfactory);
  Bob::Register<DeviceTemplate<Digital::Model, Digital::Instance> >(getXyceRegistry2(), DeviceLevelKey("nxor", 1), Digitalfactory);
  Bob::Register<DeviceTemplate<Digital::Model, Digital::Instance> >(getXyceRegistry2(), DeviceLevelKey("dff", 1), Digitalfactory);

  Bob::Register<DeviceTemplate<ACC::Model, ACC::Instance> >(getXyceRegistry2(), DeviceLevelKey("acc", 1), ACCfactory);

  Bob::Register<Neuron::Master>(getXyceRegistry2(), DeviceLevelKey("neuron", 1), Neuronfactory);
  Bob::Register<DeviceTemplate<Neuron2::Model, Neuron2::Instance> >(getXyceRegistry2(), DeviceLevelKey("neuron", 2), Neuron2factory);
  Bob::Register<DeviceTemplate<Neuron3::Model, Neuron3::Instance> >(getXyceRegistry2(), DeviceLevelKey("neuron", 3), Neuron3factory);
  Bob::Register<DeviceTemplate<Neuron4::Model, Neuron4::Instance> >(getXyceRegistry2(), DeviceLevelKey("neuron", 4), Neuron4factory);
  Bob::Register<DeviceTemplate<Neuron5::Model, Neuron5::Instance> >(getXyceRegistry2(), DeviceLevelKey("neuron", 5), Neuron5factory);
  Bob::Register<DeviceTemplate<Neuron6::Model, Neuron6::Instance> >(getXyceRegistry2(), DeviceLevelKey("neuron", 6), Neuron6factory);
  Bob::Register<DeviceTemplate<Neuron7::Model, Neuron7::Instance> >(getXyceRegistry2(), DeviceLevelKey("neuron", 7), Neuron7factory);
  Bob::Register<DeviceTemplate<Neuron8::Model, Neuron8::Instance> >(getXyceRegistry2(), DeviceLevelKey("neuron", 8), Neuron8factory);
  Bob::Register<Neuron9::Master>(getXyceRegistry2(), DeviceLevelKey("neuron", 9), Neuron9factory);

  Bob::Register<Synapse::Master>(getXyceRegistry2(), DeviceLevelKey("synapse", 1), Synapsefactory);
  Bob::Register<Synapse2::Master>(getXyceRegistry2(), DeviceLevelKey("synapse", 2), Synapse2factory);
  Bob::Register<DeviceTemplate<Synapse3::Model, Synapse3::Instance> >(getXyceRegistry2(), DeviceLevelKey("synapse", 3), Synapse3factory);
  Bob::Register<Synapse4::Master>(getXyceRegistry2(), DeviceLevelKey("synapse", 4), Synapse4factory);

  Bob::Register<DeviceTemplate<NeuronPop1::Model, NeuronPop1::Instance> >(getXyceRegistry2(), DeviceLevelKey("neuronpop", 1), NeuronPop1factory);

  Bob::Register<DeviceTemplate<NewDevice::Model, NewDevice::Instance> >(getXyceRegistry2(), DeviceLevelKey("newd", 1), NewDevicefactory);

  Bob::Register<DeviceTemplate<Xygra::Model, Xygra::Instance> >(getXyceRegistry2(), DeviceLevelKey("xygra", 1), Xygrafactory);
  Bob::Register<ROM::Master>(getXyceRegistry2(), DeviceLevelKey("rom", 1), ROMfactory);
  Bob::Register<DeviceTemplate<RxnSet::Model, RxnSet::Instance> >(getXyceRegistry2(), DeviceLevelKey("rxn", 1), RxnSetfactory);

  Fred::Register<Resistor::Master>(getXyceRegistry(), (int) ModelType::RESISTOR, Resistorfactory);
  Fred::Register<ThermalResistor::Master>(getXyceRegistry(), (int) ModelType::THERMAL_RESISTOR, ThermalResistorfactory);
  Fred::Register<Resistor3::Master>(getXyceRegistry(), (int) ModelType::RESISTOR3, Resistor3factory);
  Fred::Register<JFET::Master>(getXyceRegistry(), (int) ModelType::JFET, JFETfactory);
  Fred::Register<Diode::Master>(getXyceRegistry(), (int) ModelType::DIODE, Diodefactory);
  Fred::Register<BJT::Master>(getXyceRegistry(), (int) ModelType::BJT, BJTfactory);
  Fred::Register<DeviceTemplate<ADMSvbic::Model, ADMSvbic::Instance> >(getXyceRegistry(), (int) ModelType::ADMS_VBIC, ADMSvbicfactory);
  Fred::Register<DeviceTemplate<ADMSHBT_X::Model, ADMSHBT_X::Instance> >(getXyceRegistry(), (int) ModelType::ADMS_HBT_X, ADMSHBT_Xfactory);
  Fred::Register<DeviceTemplate<DiodePDE::Model, DiodePDE::Instance> >(getXyceRegistry(), (int) ModelType::DIODE_PDE, DiodePDEfactory);
  Fred::Register<DeviceTemplate<TwoDPDE::Model, TwoDPDE::Instance> >(getXyceRegistry(), (int) ModelType::TWO_D_PDE, TwoDPDEfactory);
  Fred::Register<MutIndNonLin::Master>(getXyceRegistry(), (int) ModelType::MUTUAL_INDUCTOR_NONLINEAR, MutIndNonLinfactory);
  Fred::Register<DeviceTemplate<MutIndNonLin2::Model, MutIndNonLin2::Instance> >(getXyceRegistry(), (int) ModelType::MUTUAL_INDUCTOR_NONLINEAR_2, MutIndNonLin2factory);
  Fred::Register<Neuron::Master>(getXyceRegistry(), (int) ModelType::NEURON, Neuronfactory);
  Fred::Register<DeviceTemplate<Neuron2::Model, Neuron2::Instance> >(getXyceRegistry(), (int) ModelType::NEURON_2, Neuron2factory);
  Fred::Register<DeviceTemplate<Neuron3::Model, Neuron3::Instance> >(getXyceRegistry(), (int) ModelType::NEURON_3, Neuron3factory);
  Fred::Register<DeviceTemplate<Neuron4::Model, Neuron4::Instance> >(getXyceRegistry(), (int) ModelType::NEURON_4, Neuron4factory);
  Fred::Register<DeviceTemplate<Neuron5::Model, Neuron5::Instance> >(getXyceRegistry(), (int) ModelType::NEURON_5, Neuron5factory);
  Fred::Register<DeviceTemplate<Neuron6::Model, Neuron6::Instance> >(getXyceRegistry(), (int) ModelType::NEURON_6, Neuron6factory);
  Fred::Register<DeviceTemplate<Neuron7::Model, Neuron7::Instance> >(getXyceRegistry(), (int) ModelType::NEURON_7, Neuron7factory);
  Fred::Register<DeviceTemplate<Neuron8::Model, Neuron8::Instance> >(getXyceRegistry(), (int) ModelType::NEURON_8, Neuron8factory);
  Fred::Register<Neuron9::Master>(getXyceRegistry(), (int) ModelType::NEURON_9, Neuron9factory);
  Fred::Register<DeviceTemplate<NeuronPop1::Model, NeuronPop1::Instance> >(getXyceRegistry(), (int) ModelType::NEURONPOP, NeuronPop1factory);
  Fred::Register<Synapse::Master>(getXyceRegistry(), (int) ModelType::SYNAPSE, Synapsefactory);
  Fred::Register<Synapse2::Master>(getXyceRegistry(), (int) ModelType::SYNAPSE_2, Synapse2factory);
  Fred::Register<DeviceTemplate<Synapse3::Model, Synapse3::Instance> >(getXyceRegistry(), (int) ModelType::SYNAPSE_3, Synapse3factory);
  Fred::Register<Synapse4::Master>(getXyceRegistry(), (int) ModelType::SYNAPSE_4, Synapse4factory);
  Fred::Register<Capacitor::Master>(getXyceRegistry(), (int) ModelType::CAPACITOR, Capacitorfactory);
  Fred::Register<Inductor::Master>(getXyceRegistry(), (int) ModelType::INDUCTOR, Inductorfactory);
  Fred::Register<MESFET::Master>(getXyceRegistry(), (int) ModelType::MESFET, MESFETfactory);
  Fred::Register<MOSFET1::Master>(getXyceRegistry(), (int) ModelType::MOSFET1, MOSFET1factory);
  Fred::Register<MOSFET2::Master>(getXyceRegistry(), (int) ModelType::MOSFET2, MOSFET2factory);
  Fred::Register<MOSFET3::Master>(getXyceRegistry(), (int) ModelType::MOSFET3, MOSFET3factory);
  Fred::Register<MOSFET6::Master>(getXyceRegistry(), (int) ModelType::MOSFET6, MOSFET6factory);
  Fred::Register<MOSFET_B3::Master>(getXyceRegistry(), (int) ModelType::MOSFET_B3, MOSFET_B3factory);
  Fred::Register<MOSFET_B4::Master>(getXyceRegistry(), (int) ModelType::MOSFET_B4, MOSFET_B4factory);
  Fred::Register<MOSFET_B3SOI::Master>(getXyceRegistry(), (int) ModelType::MOSFET_B3SOI, MOSFET_B3SOIfactory);
  Fred::Register<VDMOS::Master>(getXyceRegistry(), (int) ModelType::VDMOS, VDMOSfactory);
  Fred::Register<DeviceTemplate<ADMSPSP103VA::Model, ADMSPSP103VA::Instance> >(getXyceRegistry(), (int) ModelType::ADMS_PSP103, ADMSPSP103VAfactory);
  Fred::Register<ISRC::Master>(getXyceRegistry(), (int) ModelType::ISRC, ISRCfactory);
  Fred::Register<Vcvs::Master>(getXyceRegistry(), (int) ModelType::VCVS, Vcvsfactory);
  Fred::Register<VCCS::Master>(getXyceRegistry(), (int) ModelType::VCCS, VCCSfactory);
  Fred::Register<Vsrc::Master>(getXyceRegistry(), (int) ModelType::VSRC, Vsrcfactory);
  Fred::Register<Bsrc::Master>(getXyceRegistry(), (int) ModelType::BSRC, Bsrcfactory);
  Fred::Register<LTRA::Master>(getXyceRegistry(), (int) ModelType::LTRA, LTRAfactory);
  Fred::Register<DeviceTemplate<TRA::Model, TRA::Instance> >(getXyceRegistry(), (int) ModelType::TRA, TRAfactory);
  Fred::Register<SW::Master>(getXyceRegistry(), (int) ModelType::SW, SWfactory);
#ifdef Xyce_EXTDEV
  Fred::Register<DeviceTemplate<ExternDevice::Model, ExternDevice::Instance> >(getXyceRegistry(), (int) ModelType::EXTERN_DEVICE, ExternDevicefactory);
#endif
  Fred::Register<ADC::Master>(getXyceRegistry(), (int) ModelType::ADC, ADCfactory);
  Fred::Register<DAC::Master>(getXyceRegistry(), (int) ModelType::DAC, DACfactory);
  Fred::Register<MutIndLin::Master>(getXyceRegistry(), (int) ModelType::MUTUAL_INDUCTOR_LINEAR, MutIndLinfactory);
  Fred::Register<DeviceTemplate<OpAmp::Model, OpAmp::Instance> >(getXyceRegistry(), (int) ModelType::OPAMP, OpAmpfactory);
  Fred::Register<DeviceTemplate<Digital::Model, Digital::Instance> >(getXyceRegistry(), (int) ModelType::DIGITAL, Digitalfactory);
  Fred::Register<DeviceTemplate<ACC::Model, ACC::Instance> >(getXyceRegistry(), (int) ModelType::ACC, ACCfactory);
  Fred::Register<DeviceTemplate<NewDevice::Model, NewDevice::Instance> >(getXyceRegistry(), (int) ModelType::NEW_DEVICE, NewDevicefactory);
  Fred::Register<DeviceTemplate<Xygra::Model, Xygra::Instance> >(getXyceRegistry(), (int) ModelType::XYGRA, Xygrafactory);
  Fred::Register<ROM::Master>(getXyceRegistry(), (int) ModelType::ROM, ROMfactory);
  Fred::Register<DeviceTemplate<RxnSet::Model, RxnSet::Instance> >(getXyceRegistry(), (int) ModelType::RXNSET, RxnSetfactory);

#ifdef Xyce_RAD_MODELS
  registerSandiaDevices();
#endif

#ifdef Xyce_NONFREE_MODELS
  registerNonFreeDevices();
#endif
}

} // namespace Device
} // namespace Xyce
