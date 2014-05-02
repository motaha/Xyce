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
// Filename       : $RCSfile: N_IO_Measure_fwd.h,v $
//
// Purpose        : Forward declarations
//
// Special Notes  : Forward declaring everything as a class breaks if the implementation of the type changes (like during
//                  templatization)
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1.2.1 $
//
// Revision Date  : $Date: 2014/03/03 18:29:28 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_IO_Measure_fwd_h
#define Xyce_N_IO_Measure_fwd_h

namespace Xyce {
namespace IO {
namespace Measure {

class Base;
class Average;
class DerivativeEvaluation;
class Duty;
class EquationEvaluation;
class FindWhen;
class Fourier;
class Frequency;
class IntegralEvaluation;
class Max;
class Min;
class OffTime;
class OnTime;
class PeakToPeak;
class RMS;
class RelativeError;
class RiseFallDelay;

} // namespace Measure
} // namespace IO
} // namespace Xyce

typedef Xyce::IO::Measure::Base N_IO_MeasureBase;
typedef Xyce::IO::Measure::Average N_IO_MeasureAverage;
typedef Xyce::IO::Measure::DerivativeEvaluation N_IO_MeasureDerivativeEvaluation;
typedef Xyce::IO::Measure::Duty N_IO_MeasureDuty;
typedef Xyce::IO::Measure::EquationEvaluation N_IO_MeasureEquationEvaluation;
typedef Xyce::IO::Measure::FindWhen N_IO_MeasureFindWhen;
typedef Xyce::IO::Measure::Fourier N_IO_MeasureFourier;
typedef Xyce::IO::Measure::Frequency N_IO_MeasureFrequency;
typedef Xyce::IO::Measure::IntegralEvaluation N_IO_MeasureIntegralEvaluation;
typedef Xyce::IO::Measure::Max N_IO_MeasureMax;
typedef Xyce::IO::Measure::Min N_IO_MeasureMin;
typedef Xyce::IO::Measure::OffTime N_IO_MeasureOffTime;
typedef Xyce::IO::Measure::OnTime N_IO_MeasureOnTime;
typedef Xyce::IO::Measure::PeakToPeak N_IO_MeasurePeakToPeak;
typedef Xyce::IO::Measure::RMS N_IO_MeasureRMS;
typedef Xyce::IO::Measure::RelativeError N_IO_MeasureRelativeError;
typedef Xyce::IO::Measure::RiseFallDelay N_IO_MeasureRiseFallDelay;

#endif // Xyce_N_IO_Measure_fwd_h
