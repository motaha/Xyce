//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright(C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
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
// Filename       : $RCSfile: N_DEV_TransportHelper.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Thomas V. Russo, SNL, Component Information and Models
//
// Creation Date  : 08/19/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_TransportHelper_h
#define Xyce_N_DEV_TransportHelper_h

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : TransportHelper
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/06
//-----------------------------------------------------------------------------
class TransportHelper
{
public:
  TransportHelper () :
    D_specie  (3.6e-11),
    transportFlag(false),
    flux_bc1(0.0),
    flux_bc2(0.0),
    bcScale1(1.0),
    bcScale2(1.0)
  {};

  TransportHelper (double D, std::string & n) :
    D_specie  (D),
    name(n),
    transportFlag(false),
    flux_bc1(0.0),
    flux_bc2(0.0),
    bcScale1(1.0),
    bcScale2(1.0)
  {
    if (D!= 0.0) transportFlag = true;
  };

  std::string name;
  std::vector<int> regSubIndexVec;
  std::vector<double> fluxVec;

  double flux_bc1;
  double flux_bc2;
  double bcScale1;
  double bcScale2;

  std::vector<int> specie_id;
  double D_specie; // diffusion constant.

  bool transportFlag;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::TransportHelper N_DEV_TransportHelper;

#endif

