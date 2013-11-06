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
// Filename       : $RCSfile: N_DEV_fwd.h,v $
//
// Purpose        : Forward declarations
//
// Special Notes  : Forward declaring everything as a class breaks if the implementation of the type changes (like during
//                  templatization)
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
//
// Creation Date  : 2013/04/18 18:01:27
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_fwd_h
#define Xyce_N_DEV_fwd_h

#include <Xyce_config.h>
#include <map>
#include <string>
#include <utility>

namespace Xyce {
namespace Device {

class CompositeParam;
class Depend;
class Device;
class DeviceBuilder;
class DeviceEntity;
class DeviceInstance;
class DeviceInterface;
class DeviceMgr;
class DeviceModel;
class DeviceOptions;
class DeviceState;
class DeviceSupport;
class ExternalSimulationData;
class ExtendedDeviceBuilder;
class ExternCodeInterface;
class ExternData;
class InstanceBlock;
class MatrixLoadData;
class ModelBlock;
class NumericalJacobian;
class Param;
class Region;
class RegionData;
class RxnRegion;
class RxnRegion2;
class RxnRegionData;

class SolverState;
class SourceInstance;
class XyceInterface;

class ACData;
class ConstData;
class ExpData;
class PWLinData;
class PulseData;
class SFFMData;
class SinData;
class SmoothData;
class SmoothPulseData;
class SourceData;

class DevicePDEInstance;
class DevicePDEModel;

class PDE_Electrode;
class PDE_1DElectrode;
class PDE_2DElectrode;

class XygraCoilData;

class SpecieSource;

class ScalingVars;

namespace Xygra {
class Instance;
class Model;
}

namespace ExternDevice {
class Instance;
class Model;
}

namespace Vsrc {
class Instance;
class Model;
}

typedef std::map<std::string, CompositeParam *> CompositeMap;

struct DeviceLevelKey;

void registerDevices();

#ifdef Xyce_RAD_MODELS
void registerSandiaDevices();
#endif

#ifdef Xyce_NONFREE_MODELS
void registerNonFreeDevices();
#endif

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::CompositeParam N_DEV_CompositeParam;
typedef Xyce::Device::Depend sDepend;
typedef Xyce::Device::Device N_DEV_Device;
typedef Xyce::Device::DeviceBuilder N_DEV_DeviceBuilder;
typedef Xyce::Device::DeviceEntity N_DEV_DeviceEntity;
typedef Xyce::Device::DeviceInstance N_DEV_DeviceInstance;
typedef Xyce::Device::DeviceInterface N_DEV_DeviceInterface;
typedef Xyce::Device::DeviceMgr N_DEV_DeviceMgr;
typedef Xyce::Device::DeviceModel N_DEV_DeviceModel;
typedef Xyce::Device::DeviceOptions N_DEV_DeviceOptions;
typedef Xyce::Device::DeviceState N_DEV_DeviceState;
typedef Xyce::Device::DeviceSupport N_DEV_DeviceSupport;
typedef Xyce::Device::ExternalSimulationData N_DEV_ExternalSimulationData;
typedef Xyce::Device::ExtendedDeviceBuilder N_DEV_ExtendedDeviceBuilder;
typedef Xyce::Device::ExternCodeInterface N_DEV_ExternCodeInterface;
typedef Xyce::Device::ExternData N_DEV_ExternData;
typedef Xyce::Device::InstanceBlock N_DEV_InstanceBlock;
typedef Xyce::Device::MatrixLoadData N_DEV_MatrixLoadData;
typedef Xyce::Device::ModelBlock N_DEV_ModelBlock;
typedef Xyce::Device::NumericalJacobian N_DEV_NumericalJacobian;
typedef Xyce::Device::Param N_DEV_Param;
typedef Xyce::Device::Region N_DEV_Region;
typedef Xyce::Device::RegionData N_DEV_RegionData;
typedef Xyce::Device::RxnRegion N_DEV_RxnRegion;
typedef Xyce::Device::RxnRegion2 N_DEV_RxnRegion2;
typedef Xyce::Device::RxnRegionData N_DEV_RxnRegionData;
typedef Xyce::Device::SolverState N_DEV_SolverState;
typedef Xyce::Device::SourceInstance N_DEV_SourceInstance;
typedef Xyce::Device::XyceInterface N_DEV_XyceInterface;

typedef Xyce::Device::SourceData N_DEV_SourceData;
typedef Xyce::Device::SmoothData N_DEV_SmoothData;
typedef Xyce::Device::SinData N_DEV_SinData;
typedef Xyce::Device::ExpData N_DEV_ExpData;
typedef Xyce::Device::ACData N_DEV_ACData;
typedef Xyce::Device::PulseData N_DEV_PulseData;
typedef Xyce::Device::PWLinData N_DEV_PWLinData;
typedef Xyce::Device::SFFMData N_DEV_SFFMData;
typedef Xyce::Device::ConstData N_DEV_ConstData;
typedef Xyce::Device::SmoothPulseData N_DEV_SmoothPulseData;

typedef Xyce::Device::DevicePDEInstance N_DEV_DevicePDEInstance;
typedef Xyce::Device::DevicePDEModel N_DEV_DevicePDEModel;

typedef Xyce::Device::PDE_Electrode N_DEV_PDE_Electrode;
typedef Xyce::Device::PDE_1DElectrode N_DEV_PDE_1DElectrode;
typedef Xyce::Device::PDE_2DElectrode N_DEV_PDE_2DElectrode;

typedef Xyce::Device::Xygra::Instance N_DEV_XygraInstance;
typedef Xyce::Device::Xygra::Model N_DEV_XygraModel;

typedef Xyce::Device::XygraCoilData N_DEV_XygraCoilData;

typedef Xyce::Device::ScalingVars N_DEV_ScalingVars;

typedef Xyce::Device::SpecieSource N_DEV_SpecieSource;

typedef Xyce::Device::ExternDevice::Instance N_DEV_ExternDeviceInstance;
typedef Xyce::Device::ExternDevice::Model N_DEV_ExternDeviceModel;


#endif // Xyce_N_DEV_fwd_h

