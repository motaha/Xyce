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

/**
 * @file   N_UTL_NoCase.h
 * @author David G. Baur  KTech Corp.  Sandia National Laboratories 9143 
 * @date   Mon Apr 22 08:48:30 2013
 *
 * @brief Case insensitive functions and functors.
 *
 * Revision Number: $Revision: 1.4.2.1 $
 *
 * Revision Date  : $Date: 2014/03/03 18:29:29 $
 */

#include <N_UTL_NoCase.h>

namespace Xyce {

int compare_nocase(const char *s0, const char *s1) {
  while (tolower(*s0) == tolower(*s1)) {
    if (*s0++ == '\0')
      return 0;
    s1++;
  }
  return tolower(*s0) - tolower(*s1);
}

} // namepace Xyce
