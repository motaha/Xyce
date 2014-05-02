//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
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
// Filename      : $RCSfile: N_ANP_Report.C,v $
//
// Purpose       : Contains the class definition for the N_ERH_ErrorMgr
//                 class.  This class is the error handler for Xyce.
//
// Special Notes : 
//
// Creator       : Eric Keiter, SNL,  Parallel Computational Sciences
//
// Creation Date : 3/15/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.3.2.1 $
//
// Revision Date  : $Date: 2014/02/27 00:52:16 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <N_ANP_Report.h>
#include <N_ANP_AnalysisBase.h>

namespace Xyce {
namespace Analysis {

class messageHeader {
public:
  messageHeader(const AnalysisBase &analysis_base, const std::string &function_name = "")
    : analysisBase_(analysis_base),
      functionName_(function_name)
  {}

  const AnalysisBase &        analysisBase_;
  const std::string &           functionName_;
};


std::ostream &operator<<(std::ostream &os, const messageHeader &x)
  {
    if (!x.functionName_.empty())
      os << " in function " << x.functionName_;

    return os;
  }


UserWarning::UserWarning(const AnalysisBase &analysis_base)
{
  os() << messageHeader(analysis_base) << ": ";
}

UserWarning0::UserWarning0(const AnalysisBase &analysis_base)
{
  os() << messageHeader(analysis_base) << ": ";
}

UserFatal::UserFatal(const AnalysisBase &analysis_base)
{
  os() << messageHeader(analysis_base) << ": ";
}

UserFatal0::UserFatal0(const AnalysisBase &analysis_base)
{
  os() << messageHeader(analysis_base) << ": ";
}

DevelFatal::DevelFatal(const AnalysisBase &analysis_base, const std::string &function_name)
{
  os() << messageHeader(analysis_base, function_name) << ": ";
}

DevelFatal0::DevelFatal0(const AnalysisBase &analysis_base, const std::string &function_name)
{
  os() << messageHeader(analysis_base, function_name) << ": ";
}

DevelFatal0::DevelFatal0(const std::string &function_name)
{
  if (!function_name.empty())
    os() << " in function " << function_name;
  os() << ": ";
}

} // namespace Analysis
} // namespace Xyce
