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
// Filename       : $RCSfile: N_LAS_HBFDJacobianPrecond.h,v $
//
// Purpose        : block jacobi preconditioner designed for harmonic balance 
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 11/11/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_HBFDJacobianPrecond_h
#define Xyce_N_LAS_HBFDJacobianPrecond_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_LAS_Preconditioner.h>
#include <N_LAS_Problem.h>
#include <N_DEV_DeviceInterface.h>
#include <N_UTL_OptionBlock.h>

// ----------  Fwd Declares     ----------

class N_LAS_Problem;
class N_LAS_System;
class N_LAS_MultiVector;
class Epetra_Operator;
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_LinearProblem;
class Amesos_BaseSolver;
class N_LAS_HBBuilder;
class N_LOA_HBLoader;
class N_LAS_Builder;

//-----------------------------------------------------------------------------
// Class         : N_LAS_HBFDJacobianPrecond
// Purpose       : interface to block jacobi preconditioner
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
class N_LAS_HBFDJacobianPrecond : public N_LAS_Preconditioner
{

public:
  // Constructors
  N_LAS_HBFDJacobianPrecond();

  // Destructor
  virtual ~N_LAS_HBFDJacobianPrecond() {}                                                                                          
	
  // Set the preconditioner options
  bool setOptions(const N_UTL_OptionBlock & OB);
  bool setDefaultOptions();
  bool setDefaultOption( const std::string & option );

  // Set individual preconditioner options
  bool setParam( const N_UTL_Param & param );

  // Set the time steps being used in the HB analysis.
  void setTimeSteps( const std::vector<double> & timeSteps )
    { timeSteps_ = timeSteps; }

  // Register the application system builder
  void registerAppBuilder( const Teuchos::RCP<N_LAS_Builder>& appBuilderPtr )
    { appBuilderPtr_ = appBuilderPtr; }

  // Register the HB loader
  void registerHBLoader( const Teuchos::RCP<N_LOA_HBLoader> & hbLoaderPtr )
    { hbLoaderPtr_ = hbLoaderPtr; }

  // Register the HB builder
  void registerHBBuilder( const Teuchos::RCP<N_LAS_HBBuilder> & hbBuilderPtr )
    { hbBuilderPtr_ = hbBuilderPtr; }

  // Register the linear system pointer
  void registerLinearSystem( const Teuchos::RCP< N_LAS_System >& lasSysPtr )
    { lasSysPtr_ = lasSysPtr; }
                                    
  // Set or reset the matrix pattern for the preconditioner using problem.
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

  // Fourier information.
  // N_ is the number of time points.
  int N_;   

  // Time steps.
  std::vector<double> timeSteps_;

  // Harmonic Balance loader.
  Teuchos::RCP<N_LOA_HBLoader> hbLoaderPtr_;

  // Harmonic Balance builder.
  Teuchos::RCP<N_LAS_HBBuilder> hbBuilderPtr_;

  // Application builder.
  Teuchos::RCP<N_LAS_Builder> appBuilderPtr_; 

  // Linear system.
  Teuchos::RCP<N_LAS_System> lasSysPtr_;

  // Amesos interface to Klu solver.
  std::vector<Teuchos::RCP<Amesos_BaseSolver> > amesosPtr_;                                         

  // N_LAS_Matrix for each matrix, the underlying Epetra_CrsMatrix is passed to the direct solver.
  std::vector<Teuchos::RCP<N_LAS_Matrix> > diagMatrix_;

  // N_LAS_Vector for each vector, the underlying Epetra_MultiVector is passed to the direct solver.
  Teuchos::RCP<N_LAS_Vector> diagRHS_, diagSoln_;
                                               
  // Current problems being preconditioned.
  std::vector<Teuchos::RCP<Epetra_LinearProblem> > epetraProblem_;
                                         
  // Preconditioner as an Epetra_Operator object.                                                 
  Teuchos::RCP<Epetra_Operator> epetraPrec_;
					
  // Preconditioning matrix.  
  Teuchos::RCP<N_LAS_Matrix> matrix_;

  // Options block.
  Teuchos::RCP<const N_UTL_OptionBlock> options_;

  // No copying
  N_LAS_HBFDJacobianPrecond(const N_LAS_HBFDJacobianPrecond & right);
  N_LAS_HBFDJacobianPrecond & operator=(const N_LAS_HBFDJacobianPrecond & right);

  // No comparison
  bool operator==(const N_LAS_HBFDJacobianPrecond & right) const;
  bool operator!=(const N_LAS_HBFDJacobianPrecond & right) const;
                                                      
};
                                                                                          
#endif
