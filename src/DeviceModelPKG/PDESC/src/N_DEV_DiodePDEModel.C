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
//
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_DiodePDEModel.C,v $
//
// Purpose        : One dimensional PDE device, model class implementation.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/06/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.13.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------  Standard Includes ----------
#include <N_UTL_Misc.h>
#ifdef Xyce_DEBUG_DEVICE
#include <iostream>
#endif

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_DiodePDE.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<DiodePDE::Model>::ParametricData()
{}


namespace DiodePDE {


ParametricData<Model> &Model::getParametricData() {
  static ParametricData<Model> parMap;

  return parMap;
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Model::processParams (string param)
{
  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirely, PSSI
// Creation Date : 03/23/06
//----------------------------------------------------------------------------
bool Model::processInstanceParams(string param)
{

  vector<Instance*>::iterator iter;
  vector<Instance*>::iterator first = instanceContainer.begin();
  vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    (*iter)->processParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/00
//-----------------------------------------------------------------------------

Model::Model(
  const ModelBlock & MB,
  SolverState &      ss1,
  DeviceOptions &    do1)
  : DevicePDEModel(MB,ss1,do1)
{
  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/00
//-----------------------------------------------------------------------------
Model::~Model()
{
  vector<Instance*>::iterator iter;
  vector<Instance*>::iterator first = instanceContainer.begin();
  vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    delete (*iter);
  }

}

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << endl;
  os << "    name     getModelName()  Parameters" << endl;
  for (i = 0, iter = first; iter != last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
    os << (*iter)->getModelName();
    //os << "   C = " << (*iter)-> C;
    //os << "  IC = " << (*iter)->IC;
    os << endl;
  }
  os << endl;
}

} // namespace DiodePDE
} // namespace Device
} // namespace Xyce
