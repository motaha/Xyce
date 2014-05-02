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
// Filename       : $RCSfile: N_UTL_fwd.h,v $
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
// Revision Number: $Revision: 1.12.2.1 $
//
// Revision Date  : $Date: 2014/03/03 18:29:29 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_fwd_h
#define Xyce_N_UTL_fwd_h

#include <complex>
#include <iosfwd>
#include <map>
#include <list>
#include <string>
#include <utility>

#include <N_UTL_NoCase.h>

namespace Xyce {

std::ostream &lout();
std::ostream &dout();
std::ostream &pout();

extern const char *section_divider;             // Defined in LogStream
extern const char *subsection_divider;          // Defined in LogStream

typedef std::map<std::string, std::pair<int, double>, LessNoCase > NodeNamePairMap;

template <class T>

struct DataTypeTrait;

typedef std::complex<double> complex;

namespace Util {

class BreakPoint;
class Expression;
class ExpressionData;
class ExpressionInternals;
class MachineDependentParams;
class OptionBlock;
class Param;
class Timer;

template<class T, class R, class E = T>
struct Op;

template<class Ch, class Tr>
std::basic_ostream<Ch, Tr> &push(std::basic_ostream<Ch, Tr> &os);

template<class Ch, class Tr>
std::basic_ostream<Ch, Tr> &pop(std::basic_ostream<Ch, Tr> &os);

typedef std::list<Util::Param> ParameterList;

} // namespace Util
} // namespace Xyce

typedef Xyce::NodeNamePairMap N_UTL_NodePairMap;

typedef Xyce::Util::MachineDependentParams N_UTL_MachineDependentParams;
typedef Xyce::Util::OptionBlock N_UTL_OptionBlock;
typedef Xyce::Util::Param N_UTL_Param;
typedef Xyce::Util::Timer N_UTL_Timer;
typedef Xyce::Util::Expression N_UTL_Expression;
typedef Xyce::Util::ExpressionInternals N_UTL_ExpressionInternals;
typedef Xyce::Util::ExpressionData N_UTL_ExpressionData;
typedef Xyce::Util::BreakPoint N_UTL_BreakPoint;

#endif // Xyce_N_UTL_fwd_h
