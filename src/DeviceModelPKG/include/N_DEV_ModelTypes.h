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
// Filename       : $RCSfile: N_DEV_ModelTypes.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/26/07
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.36.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_ModelTypes_h
#define Xyce_N_DEV_ModelTypes_h

namespace Xyce {
namespace Device {

struct ModelType {

// ---------- Enum Definitions ----------
// device indices:
enum DevIndex {
  DUMMY,                        ///< This is a placeholder.
  RESISTOR,                     ///<
  THERMAL_RESISTOR,             ///< thermal resistor model
  RESISTOR3,                    ///< special case of resistor with zero resistance
  CAPACITOR,                    ///<
  INDUCTOR,                     ///<
  MUTUAL_INDUCTOR_LINEAR,       ///< a linear mutual inductor
  MUTUAL_INDUCTOR_NONLINEAR,    ///< a nonliner mutal inductor
  MUTUAL_INDUCTOR_NONLINEAR_2,  ///< a nonliner mutal inductor lv 2
  DIODE,                        ///<
  BJT,                          ///<
  JFET,                         ///<
  MESFET,                       ///<
  MOSFET1,                      ///< level = 1 mosfet.  (implemented)
  MOSFET2,                      ///< level = 2 mosfet.  (implemented)
  MOSFET3,                      ///< level = 3 mosfet.  (implemented)
  MOSFET6,                      ///< level = 6 mosfet.  (implemented)
  MOSFET_B3,                    ///< BSIM3   mosfet.  (implemented)
  MOSFET_B4,                    ///< BSIM4   mosfet.
  MOSFET_B3SOI,                 ///< BSIM3 silicon on insulator mosfet.  (implemented)
  VDMOS,                        ///< UCCM    mosfet.  (implemented)
  ISRC,                         ///< independent current source
  VCCS,                         ///< voltage controlled current source
  VCVS,                         ///< voltage controlled voltage source
  VSRC,                         ///< independent voltage source
  BSRC,                         ///< expression dependent source
  SW,                           ///< controlled switch
  TRA,                          ///< transmission line
  LTRA,                         ///< lossy transmission line
  ADC,                          ///< Analog to digital converter
  DAC,                          ///< Digital to analog converter
  DIODE_PDE,                    ///< PDE-based diode  (1D PDE problem)
  TWO_D_PDE,                    ///< 2d PDE based device
#ifdef Xyce_EXTDEV
  EXTERN_DEVICE,                ///< External device, could be PDE, or a subcircuit
#endif
  OPAMP,                        ///< Ideal Operational Amplifier
  DIGITAL,                      ///< General purpose Digital Device with Analog I/O
  ACC,                          ///< Accelerated mass device
  NEURON,                       ///< Neuron model device (Hodgkin-Huxley model)
  NEURON_2,                     ///< Neuron model device (Connor Stevens model)
  NEURON_3,                     ///< Cable equation neuron model device (Hodgkin-Huxley model)
  NEURON_4,                     ///< Cable equation neuron model device (Connor Stevens model)
  NEURON_5,                     ///< Generalized Membrane patch (dynamic ion-channel equations)
  NEURON_6,                     ///< Generalized Cable equation (dynamic ion-channel equations)
  NEURON_7,                     ///< non-linear point neuron model Isikavish
  NEURON_8,                     ///< non-linear point neuron model Mihalas & Neibur
  NEURON_9,                     ///< point neuron model described in Brette et al 07 benchmark
  SYNAPSE,                      ///< a synapse joining neurons, calculates neurotransmitter concentration
  SYNAPSE_2,                    ///< a synapse joining neurons, gated by presynaptic voltage
  SYNAPSE_3,                    ///< a synapse joining neurons, simple delay synapse
  SYNAPSE_4,                    ///< a synapse joining neurons, based on NEURON's Exp2Syn
  NEURONPOP,                    ///< a population of neurons
  NEW_DEVICE,                   ///< template for a new device
  XYGRA,                        ///< Xyce/Alegra coupling device
  ADMS_VBIC,                    ///< VBIC model converted from Verilog by ADMS
  ADMS_PSP103,                  ///< PSP 103.1 model converted from Verilog by ADMS
  ADMS_HBT_X,                   ///< FBH HBT v2.3 model converted from Verilog by ADMS
#ifdef Xyce_NONFREE_MODELS
  ADMS_EKV,                     ///< EKV  301.02 model converted from Verilog by ADMS
#endif
  ROM,                          ///< reduced order model
  RXNSET,                       ///< reaction-parser example model
  NUMDEV                        ///< total number of device types
};

};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::ModelType N_DEV_ModelType;

#endif

