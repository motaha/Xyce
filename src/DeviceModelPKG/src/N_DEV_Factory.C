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
// Filename       : $RCSfile: N_DEV_Factory.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
//
// Creation Date  : 2013/04/18 18:01:27
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.3.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>
#include <N_DEV_Factory.h>

namespace Xyce {
namespace Device {


/** 
 * getXyceRegistry 
 *
 *
 * @date   Mon Aug  5 10:11:59 2013
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * 
 *
 * @return 
 */
Registry &getXyceRegistry() {
  static Registry theRegistry;

  return theRegistry;
}

/** 
 * getXyceRegistry2 
 *
 *
 * @date   Mon Aug  5 10:12:06 2013
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * 
 *
 * @return 
 */
Registry2 &getXyceRegistry2() {
  static Registry2 theRegistry;

  return theRegistry;
}

/** 
 * getXyceInstanceRegistry 
 *
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Mon Aug  5 10:12:10 2013
 * 
 *
 * @return 
 */
Registry2 &getXyceInstanceRegistry() {
  static Registry2 theRegistry;

  return theRegistry;
}

/** 
 * getXyceModelRegistry 
 *
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Mon Aug  5 10:12:14 2013
 * 
 *
 * @return 
 */
Registry2 &getXyceModelRegistry() {
  static Registry2 theRegistry;

  return theRegistry;
}

} // namespace Device
} // namespace Xyce
