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
// Filename       : $RCSfile: N_DEV_MOSFET_B2.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MOSFET_B2_h
#define Xyce_N_DEV_MOSFET_B2_h


// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_Device.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : MOSFET_B2Instance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

class MOSFET_B2Instance : public DeviceInstance
{
  public:
    MOSFET_B2Instance(const MOSFET_B2Instance &right);
    ~MOSFET_B2Instance();

  protected:
  private:

  private:
    friend class MOSFET_B2;
    friend class MOSFET_B2Model;
};

//-----------------------------------------------------------------------------
// Class         : MOSFET_B2Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class MOSFET_B2Model : public DeviceModel
{
  public:
    MOSFET_B2Model (const MOSFET_B2Model &right);
    ~MOSFET_B2Model ();

    const MOSFET_B2Instance * get_the_MOSFET_B2Instance () const;
    void set_the_MOSFET_B2Instance (MOSFET_B2Instance * value);

  protected:

  private:

  private:
    MOSFET_B2Instance *the_MOSFET_B2Instance;

    friend class MOSFET_B2;
    friend class MOSFET_B2Instance;
};


//-----------------------------------------------------------------------------
// Class         : MOSFET_B2
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class MOSFET_B2 : public Device
{
  public:
//    MOSFET_B2(const MOSFET_B2 &right);
    ~MOSFET_B2();

    const MOSFET_B2Model * get_the_MOSFET_B2Model () const;
    void set_the_MOSFET_B2Model (MOSFET_B2Model * value);

//    static Device * factory(SolverState & ss1, DeviceOptions & do1);

  protected:

  private:
    MOSFET_B2(SolverState & ss1,
                          DeviceOptions & do1);

  private:
    MOSFET_B2Model *the_MOSFET_B2Model;

    friend class MOSFET_B2Model;
    friend class MOSFET_B2Instance;
};

// Class MOSFET_B2Instance

// Class MOSFET_B2Model

//-----------------------------------------------------------------------------
// Function      : MOSFET2Model::get_the_MOSFET_B2Instance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline const MOSFET_B2Instance * MOSFET_B2Model::get_the_MOSFET_B2Instance () const
{
  return the_MOSFET_B2Instance;
}

//-----------------------------------------------------------------------------
// Function      : MOSFET2Model::set_the_MOSFET_B2Instance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline void MOSFET_B2Model::set_the_MOSFET_B2Instance (MOSFET_B2Instance * value)
{
  the_MOSFET_B2Instance = value;
}

// Class MOSFET_B2

//-----------------------------------------------------------------------------
// Function      : MOSFET2::get_the_MOSFET_B2Model
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline const MOSFET_B2Model * MOSFET_B2::get_the_MOSFET_B2Model () const
{
  return the_MOSFET_B2Model;
}

//-----------------------------------------------------------------------------
// Function      : MOSFET2::set_the_MOSFET_B2Model
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline void MOSFET_B2::set_the_MOSFET_B2Model (MOSFET_B2Model * value)
{
  the_MOSFET_B2Model = value;
}

//-----------------------------------------------------------------------------
// Function      : MOSFET2::factory
// Purpose       : Factory function for the MOSFET_B2 class.
//
// Special Notes : ERK.  10/16/2005.  This used to be a singleton (ie a
//                 static pointer was returned) but had to be changed
//                 so that the library version of Xyce would work
//                 correctly.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
// inline Device * MOSFET_B2::factory(SolverState & ss1,
//                           DeviceOptions & do1)
// {
//    Device * devptr = new MOSFET_B2(ss1,do1);
//    return devptr;
// }

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::MOSFET_B2Instance N_DEV_MOSFET_B2Instance;
typedef Xyce::Device::MOSFET_B2Model N_DEV_MOSFET_B2Model;

#endif
