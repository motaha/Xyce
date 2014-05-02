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
// Filename      : $RCSfile: N_MPDE_Discretization.h,v $
//
// Purpose       : Tool for Discretization Properties
//
// Special Notes : 
//
// Creator       : Robert Hoekstra, 9233, Computational Sciences
//
// Creation Date : 5/6/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_MPDE_DISCRETIZATION_H
#define Xyce_MPDE_DISCRETIZATION_H

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>

#include <vector>

// ---------- Forward Declarations ----------

// ---------- Enum Definitions ----------

//-----------------------------------------------------------------------------
// Class         : N_MPDE_Discretization
// Purpose       : 
// Special Notes :
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 5/6/04
//-----------------------------------------------------------------------------
class N_MPDE_Discretization
{
 public:

  enum Type{ Backward, Centered, Forward };

  //Constructor
  N_MPDE_Discretization( Type type, int Order );

  int Start() const { return start_; }
  int Width() const { return width_; }
  int Order() const { return order_; }

  const std::vector<double> & Coeffs() const { return coeffs_; }

 private:

  Type type_;
  int order_;

  int start_;
  int width_;

  std::vector<double> coeffs_;

  void GenerateCoeffs_( Type type, int Order, std::vector<double> & coeffs );
};

#endif //Xyce_MPDE_DISCRETIZATION_H
