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

#include <N_DEV_Message.h>

#include <N_DEV_RegisterDevices.h>
#include <N_DEV_RegisterOpenDevices.h>
#include <N_DEV_RegisterTCADDevices.h>
#include <N_DEV_RegisterNeuronDevices.h>
#include <N_DEV_RegisterADMSDevices.h>

#ifdef Xyce_EXTDEV
#include <N_DEV_ExternDevice.h>
#endif

#ifdef Xyce_NONFREE_MODELS
#include <N_DEV_RegisterNonFreeDevices.h>
#endif

#ifdef Xyce_RAD_MODELS
#include <N_DEV_RegisterSandiaDevices.h>
#endif

namespace Xyce {
namespace Device {

void
registerDevices()
{
  static bool initialized = false;

  if (!initialized) {
    initialized = true;
    
    registerOpenDevices();
    registerNeuronDevices();
    registerADMSDevices();
    registerTCADDevices();

#ifdef Xyce_EXTDEV
    ExternDevice::registerDevice();
#endif

#ifdef Xyce_RAD_MODELS
    registerSandiaDevices();
#endif

#ifdef Xyce_NONFREE_MODELS
    registerNonFreeDevices();
#endif
  }
}

void registerDL(const char *so_path, const char *function_key = 0);

void
registerPlugin(const char *name) 
{
  registerDL(name);
}

typedef void (*dl_register_t)();

void registerDL(const char *so_path, const char *function_key) {
#ifdef HAVE_DLFCN_H
  void *dl = dlopen(so_path, RTLD_NOW);
  if (!dl) {
    const char *error = dlerror();
    Report::UserError0() << "Failed to load plugin " << error;
  }
  else {
    if (function_key) {
      std::string s = std::strlen(function_key) ? function_key : "dl_register";

      dl_register_t f = (dl_register_t) dlsym(dl, s.c_str());
      if (!f) {
        f = (dl_register_t) dlsym(dl, s.c_str());
      }

      if (f) {
        (*f)();
      }
      else {
        if (std::strlen(function_key)) {
          Report::UserError0() << "Executing dynamic library " << so_path << " function " << s << "()";
        }
      }
    }
  }
#endif // HAVE_DLFCN_H
}

} // namespace Device
} // namespace Xyce
