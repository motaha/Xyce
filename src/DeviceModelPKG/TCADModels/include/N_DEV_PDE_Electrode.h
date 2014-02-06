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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_PDE_Electrode.h,v $
//
// Purpose        : This is the class for mesh processing/ownership.
//                  of two dimensional meshes.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/21/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_PDE_Electrode__h
#define Xyce_N_DEV_PDE_Electrode__h

#include <N_DEV_CompositeParam.h>


namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : PDE_Electrode
//
// Purpose       : This class contains the base class for user
//                 specifications of electrodes.
//
// Special Notes : In general, this class will ONLY contain info that came
//                 from netlist user-specified vector-composites.  It
//                 does not contain much else.
//
//                 For both the 1D and 2D device, there are other classes
//                 (such as bcData for 1D and deviceInterfaceNode for 2D),
//                 which contain a lot of other information, such as the
//                 names of the circuit nodes, and all the indexing
//                 information.
//
// Creator       : Eric Keiter
// Creation Date : 04/17/03
//-----------------------------------------------------------------------------
class PDE_Electrode : public CompositeParam
{
public:
  PDE_Electrode ()
    : CompositeParam (),
      name   ("ANODE"),// got to call it something...
      nodeName("node1"),
      bcName("bc1"),
      material ("neutral"),
      materialGiven(false),
      oxideBndryFlag(false),
      oxthick(0.0),
      oxcharge(0.0)
  {};


  // PDE_Electrode (const PDE_Electrode & right)
  //   : CompositeParam (right),

  //     name(right.name),
  //     nodeName(right.nodeName),
  //     bcName(right.bcName),
  //     material(right.material),
  //     materialGiven(right.materialGiven),
  //     oxideBndryFlag(right.oxideBndryFlag),
  //     oxthick(right.oxthick),
  //     oxcharge(right.oxcharge)

  // {};
  virtual ~PDE_Electrode () {}

private:
  PDE_Electrode(const PDE_Electrode &);

public:
  virtual void processParams () {}

public:
  string name;     // name of the electrode.
  string nodeName; // name of the ckt node.
  string bcName;   // name of the bc.
  string material;
  bool   materialGiven;
  bool   oxideBndryFlag;
  double oxthick;
  double oxcharge;

};

//-----------------------------------------------------------------------------
// Class         : PDE_1DElectrode
// Purpose       : This class contains user specification of a 1D electrode.
//
// Special Notes :
//
// Creator       : Eric Keiter
// Creation Date : 04/17/03
//-----------------------------------------------------------------------------
class PDE_1DElectrode : public PDE_Electrode
{
  friend class ParametricData<PDE_1DElectrode>;

public:
  static ParametricData<PDE_1DElectrode> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

  PDE_1DElectrode ();

  // PDE_1DElectrode (const PDE_1DElectrode & right)
  //   : PDE_Electrode (right),
  //     area(right.area),
  //     location(right.location),
  //     sideGiven(right.sideGiven),
  //     side(right.side)
  // {}

  virtual ~PDE_1DElectrode () {}
  virtual void processParams ();

private:
  PDE_1DElectrode(const PDE_1DElectrode &);

public:
  double area;
  bool areaGiven;
  double location;
  bool sideGiven;
  string side;    // this class implicitly assumes that the device is
  // one dimensional.  The options for side in this
  // class are: left (x=0), middle (0<x<xmax), right (x=xmax)

};

//-----------------------------------------------------------------------------
// Class         : PDE_2DElectrode
// Purpose       : This class contains user specification of a 2D electrode.
//
// Special Notes :
//
// Creator       : Eric Keiter
// Creation Date : 04/17/03
//-----------------------------------------------------------------------------
class PDE_2DElectrode : public PDE_Electrode
{
  friend class ParametricData<PDE_2DElectrode>;

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

public:
  static ParametricData<PDE_2DElectrode> &getParametricData();

  PDE_2DElectrode ();
  // PDE_2DElectrode (const PDE_2DElectrode & right)
  //   : PDE_Electrode (right),

  //     start(right.start),
  //     end(right.end),
  //     startGiven(right.startGiven),
  //     endGiven(right.endGiven),
  //     sideGiven(right.sideGiven),
  //     side(right.side),
  //     iA (right.iA),
  //     iB(right.iB),
  //     uLabel(right.uLabel)

  // {}

  virtual ~PDE_2DElectrode () {}
  virtual void processParams ();

private:
  PDE_2DElectrode(const PDE_2DElectrode &);

public:
  double start;     // beginning location.
  double end;       // ending location.

  bool startGiven;  // beginning location.
  bool endGiven;    // ending location.

  bool sideGiven;
  string side;    // this class implicitly assumes that
  // the device is a 4-sided parallelogram.  Any
  // electrode, therefore, is on the top, bottom, left,
  // or right side.

  int iA, iB;
  int uLabel;   // label index
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::PDE_Electrode N_DEV_PDE_Electrode;
typedef Xyce::Device::PDE_1DElectrode N_DEV_PDE_1DElectrode;
typedef Xyce::Device::PDE_2DElectrode N_DEV_PDE_2DElectrode;

#endif
