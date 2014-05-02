//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_LTRA_Faddeeva.h,v $
//
// Purpose        : 
//
// Special Notes : This file is, at present, only required by Windows
//                 builds because of the Intel compilers on that
//                 platform not putting erfc() in the standard math
//                 library.
//
// Creator        : Gary Hennigan
//
// Creation Date  : 12/7/2012
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.9.2.1 $
//
// Revision Date  : $Date: 2014/03/04 23:50:54 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

// Copyright (c) 2012, 2013 Massachusetts Institute of Technology
// 
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
// 
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// 

//  Available at: http://ab-initio.mit.edu/Faddeeva

//  Header file for Faddeeva.cc; see that file for more information.

#ifndef N_DEV_LTRA_Faddeeva_h
#define N_DEV_LTRA_Faddeeva_h

#include <complex>

namespace Xyce {
namespace Device {
namespace Faddeeva {

// compute erfcx(z) = exp(z^2) erfc(z)
extern double erfcx(double x); // special case for real x

// compute erf(z), the error function of complex arguments
extern double erf(double x); // special case for real x

// compute erfc(z) = 1 - erf(z), the complementary error function
extern double erfc(double x); // special case for real x

} // namespace Faddeeva
} // namespace Device
} // namespace Xyce

#endif // N_DEV_LTRA_Faddeeva_h
