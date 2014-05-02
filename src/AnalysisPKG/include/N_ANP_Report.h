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
// Filename      : $RCSfile: N_ANP_Report.h,v $
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

#ifndef Xyce_N_ANP_Report_H
#define Xyce_N_ANP_Report_H

#include <N_ANP_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_UTL_LogStream.h>

namespace Xyce {
namespace Analysis {

struct UserWarning : public Report::UserWarning
{
    UserWarning(const AnalysisBase &analysis_base);
};

struct UserWarning0 : public Report::UserWarning0
{
    UserWarning0(const AnalysisBase &analysis_base);
};

struct UserFatal : public Report::UserFatal
{
    UserFatal(const AnalysisBase &analysis_base);
};

struct UserFatal0 : public Report::UserFatal0
{
    UserFatal0(const AnalysisBase &analysis_base);
};

struct DevelFatal : public Report::DevelFatal
{
    DevelFatal(const AnalysisBase &analysis_base, const std::string &function_name = "");
};

struct DevelFatal0 : public Report::DevelFatal0
{
    DevelFatal0(const AnalysisBase &analysis_base, const std::string &function_name = "");
    DevelFatal0(const std::string &function_name = "");
};

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_Report_H
