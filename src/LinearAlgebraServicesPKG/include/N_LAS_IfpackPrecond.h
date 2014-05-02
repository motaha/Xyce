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
// Filename       : $RCSfile: N_LAS_IfpackPrecond.h,v $
//
// Purpose        : interface to preconditioner
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 09/27/07
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.16 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_IfpackPrecond_h
#define Xyce_N_LAS_IfpackPrecond_h

#include <N_UTL_fwd.h>

#include <N_LAS_Preconditioner.h>
#include <N_LAS_Problem.h>

// ----------  Fwd Declares     ----------

class N_LAS_MultiVector;

class Epetra_Operator;

class Ifpack_Preconditioner;
class Ifpack_CrsRiluk;
class Ifpack_IlukGraph;

//-----------------------------------------------------------------------------
// Class         : N_LAS_IfpackPrecond
// Purpose       : interface to ifpack preconditioner
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
class N_LAS_IfpackPrecond : public N_LAS_Preconditioner
{

public:
  // Constructors
  N_LAS_IfpackPrecond();

  // Destructor
  virtual ~N_LAS_IfpackPrecond() {}

  // Set the preconditioner options
  bool setOptions(const N_UTL_OptionBlock & OB);
  bool setDefaultOptions();
  bool setDefaultOption( const std::string & option );

  // Set individual preconditioner options
  bool setParam( const N_UTL_Param & param );

  // Set or reset the matrix pattern for the preconditioner using problem.
  // \note The preconditioner will be recreated each time this is called.
  bool initGraph( const Teuchos::RCP<N_LAS_Problem> & problem );

  // Set the matrix values for the preconditioner
  bool initValues( const Teuchos::RCP<N_LAS_Problem> & problem );

  // Compute the preconditioner using the current matrix values.
  bool compute();

  // Apply the preconditioner; y = M*x.
  int apply( N_LAS_MultiVector & x, N_LAS_MultiVector & y );

  // Return the preconditioner as an Epetra_Operator object.
  Teuchos::RCP<Epetra_Operator> epetraObj() { return epetraPrec_; }

private:

  // Indicate whether the old or new Ifpack factory interface is being used.
  bool useFactory_;

  // If the factory is used, which preconditioner should be used.
  std::string ifpackType_;

  // Diagonal perturbation, if necessary, to compute better conditioned preconditioner.
  double diagPerturb_;

  // Overlap for Additive Schwarz preconditioner
  int overlap_;
  // Drop tolerance for LU or ILUT preconditioners
  double dropTol_;
  // Fill-factor for ILUT preconditioner
  double ilutFill_;
  // Diagonal shifting relative threshold
  double rThresh_;
  // Diagonal shifting absolute threshold
  double aThresh_;

  // Default preconditioner values
  static const bool useFactory_default_;
  static const std::string ifpackType_default_;
  static const double diagPerturb_default_;
  static const int overlap_default_;
  static const double dropTol_default_;
  static const double ilutFill_default_;
  static const double rThresh_default_;
  static const double aThresh_default_;

  // Old Ifpack preconditioning classes.
  Teuchos::RCP<Ifpack_IlukGraph> ilukGraph_;
  Teuchos::RCP<Ifpack_CrsRiluk>  rILUK_;

  // New Ifpack abstract interface to preconditioners.
  Teuchos::RCP<Ifpack_Preconditioner> ifpackPrecond_;

  // Current matrix being preconditioned.
  Teuchos::RCP<N_LAS_Problem> problem_;

  // Preconditioner as an Epetra_Operator object.
  Teuchos::RCP<Epetra_Operator> epetraPrec_;

  // Options block.
  Teuchos::RCP<N_UTL_OptionBlock> options_;

  // No copying
  N_LAS_IfpackPrecond(const N_LAS_IfpackPrecond & right);
  N_LAS_IfpackPrecond & operator=(const N_LAS_IfpackPrecond & right);

  // No comparison
  bool operator==(const N_LAS_IfpackPrecond & right) const;
  bool operator!=(const N_LAS_IfpackPrecond & right) const;

};

#endif
