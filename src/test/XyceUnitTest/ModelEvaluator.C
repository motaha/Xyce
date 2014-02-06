//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: ModelEvaluator.C,v $
// Purpose       : This file contains unit tests for the Model Evaluator interface
// Special Notes :
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 05/26/09
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.3 $
// Revision Date  : $Date: 2011/11/29 22:46:38 $
// Current Owner  : $Author: hkthorn $
//-----------------------------------------------------------------------------


#ifdef Xyce_TRILINOS_DEV

// ---------- Trilinos Includes ----------
#include <Teuchos_UnitTestHarness.hpp>
#include <EpetraExt_ModelEvaluator.h>
#include <Epetra_LinearProblem.h>
#include <Amesos_Klu.h>
// ---------- Trilinos::Rythmos Includes ---------
#ifdef Xyce_TRILINOS_DEV_RYTHMOS
#include <Rythmos_ConfigDefs.h>
#include <Rythmos_TimeStepNonlinearSolver.hpp>
#include <Thyra_NonlinearSolverBase.hpp>
#include <Rythmos_BackwardEulerStepper.hpp>
#include <Rythmos_ImplicitRKStepper.hpp>
#include <Rythmos_RKButcherTableauBuilder.hpp>
#include <Rythmos_RKButcherTableau.hpp>
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_EpetraThyraWrappers.hpp>
#include <Thyra_ModelEvaluator.hpp>
#include <Thyra_EpetraModelEvaluator.hpp>
#endif // Xyce_TRILINOS_DEV_RYTHMOS
// ---------- Standard Includes ----------
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <EpetraExt_MatrixMatrix.h>
#include <map>
// ---------- Special Includes ----------
#include <N_ANP_ModelEvaluator.h>
#include <N_ANP_ModelEvaluator_Stateless.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockVector.h>

using Teuchos::RCP;
using Teuchos::rcp;
 
std::vector<int> getVariableIndices(
    const std::vector<std::string> & inNames,
    const std::vector<std::string> & outNames
    )
{
  TEUCHOS_ASSERT( inNames.size() == outNames.size() );
  std::vector<int> indices;
  for (int i = 0 ; i < outNames.size() ; ++i) {
    std::vector<std::string>::const_iterator result = 
      std::find(inNames.begin(),inNames.end(),outNames[i]);
    int ind = result-inNames.begin();
    indices.push_back(ind);
  }
  return indices;
}

TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, create ) {
  RCP<N_ANP_ModelEvaluator> model = rcp(new N_ANP_ModelEvaluator());
  TEST_EQUALITY_CONST( Teuchos::is_null(model), false );
}
TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, getVariableNames ) {
  int iargs = 2;
  char *cargs[iargs];
  cargs[0] = "Xyce";
  cargs[1] = "modelEvaluatorTest.cir";
  // Use simple RC circuit netlist as input
  // Set up t,x,xdot,s,sdot,f,W
  RCP<N_ANP_ModelEvaluator> model = rcp(new N_ANP_ModelEvaluator());
  model->initialize(iargs,cargs);
  std::vector<std::string> names = model->getVariableNames();
  TEST_ASSERT( names.size() == 3);
  std::vector<std::string> myNames;
  {
    myNames.push_back("V(2)");
    myNames.push_back("I(VSRC)");
    myNames.push_back("V(1)");
  }
  std::vector<int> ind = getVariableIndices(names,myNames);
  TEST_EQUALITY_CONST( names[ind[0]], "V(2)" );
  TEST_EQUALITY_CONST( names[ind[1]], "I(VSRC)" );
  TEST_EQUALITY_CONST( names[ind[2]], "V(1)" );
}

// In all the tests below, we've assumed the indices are as follows:
// x[0] = V(2)
// x[1] = I(VSRC)
// x[2] = V(1)
// f[0] = (v2-v1)/R;
// f[1] = v1-vsrc;
// f[2] = (v1-v2)/R+i_vsrc;
// f[3] = C*v2;
// f[4] = 0.0;
// f[5] = 0.0;
//void determineIndices(
//    Teuchos::Array<int>& xi,
//    Teuchos::Array<int>& fi
//    )
//{
//  int iargs = 2;
//  char *cargs[iargs];
//  cargs[0] = "Xyce";
//  cargs[1] = "modelEvaluatorTest.cir";
//  // Use simple RC circuit netlist as input
//  // Set up t,x,xdot,s,sdot,f,W
//  RCP<N_ANP_ModelEvaluator> model = rcp(new N_ANP_ModelEvaluator());
//  model->initialize(iargs,cargs);
//  // The solution variables can be figured out through Xyce:
//  //std::vector<std::string> names = model->getVariableNames();
//  //int L = Teuchos::as<int>(names.size());
//  //std::vector<std::string> myNames;
//  //{
//  //  myNames.push_back("V(2)");
//  //  myNames.push_back("I(VSRC)");
//  //  myNames.push_back("V(1)");
//  //}
//  //xi.clear();
//  //xi = getVariableIndices(names,myNames);
//  //for (int i=0 ; i<L ; ++i) {
//  //  xi.push_back(xi[i]+L);
//  //}
//  
//  // To figure out the F variables, we need to work a little harder:
//  // Set V(2) = 1, V(1) = 0, I(VSRC) = 0, then 
//  // F(0) = 1/R, F(1) = 0, F(2) = -1/R
//  // F(3) = C, F(4) = 0, F(5) = 0
//  // This should tell us what order everything is inside F
//  RCP<const Epetra_Map> x_map = model->get_x_map();
//  RCP<const Epetra_Map> s_map = model->get_p_map(0);
//  RCP<const Epetra_Map> sdot_map = model->get_p_map(1);
//  RCP<const Epetra_Map> f_map = model->get_f_map();
//  RCP<Epetra_Vector> xbar = rcp(new Epetra_Vector(*x_map)); 
//  RCP<Epetra_Vector> xbardot = rcp(new Epetra_Vector(*x_map)); 
//  RCP<Epetra_Vector> s = rcp(new Epetra_Vector(*s_map)); 
//  RCP<Epetra_Vector> sdot = rcp(new Epetra_Vector(*sdot_map)); 
//  RCP<Epetra_Vector> fbar = rcp(new Epetra_Vector(*f_map)); 
//  RCP<Epetra_Vector> voltLimF = rcp(new Epetra_Vector(*f_map)); 
//  RCP<Epetra_Vector> voltLimQ = rcp(new Epetra_Vector(*f_map)); 
//  RCP<Epetra_Operator> W = model->create_W();
//  int L = 3;
//  double C = 1.0e-7;
//  double R = 1.0e6;
//  double t = 0.0;
//  double alpha = 0.0;
//  double beta = 0.0;
//  Epetra_Vector& xb = *xbar;
//  xb[0] = 1.0;
//  xb[1] = 2.0;
//  xb[2] = 3.0;
//  xb[3] = 0.0;
//  xb[4] = 0.0;
//  xb[5] = 0.0;
//  for (int i=0 ; i<2*L ; ++i) {
//    std::cout << "xb["<<i<<"] = " << xb[i] << std::endl;
//  }
//  std::cout << std::endl;
//  {
//    EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
//    inArgs.set_t(t);
//    inArgs.set_x(xbar);
//    inArgs.set_x_dot(xbardot);
//    EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
//		outArgs.set_g(0,s);
//    model->evalModel(inArgs,outArgs);
//  }
//  {
//    EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
//    inArgs.set_t(t);
//    inArgs.set_x(xbar);
//    inArgs.set_x_dot(xbardot);
//    inArgs.set_p(0,s);
//    inArgs.set_p(1,sdot);
//    inArgs.set_alpha(alpha);
//    inArgs.set_beta(beta);
//    EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
//    fbar->PutScalar(0.0);
//    outArgs.set_f(fbar);
//    outArgs.set_W(W);
//    outArgs.set_g(1,voltLimQ);
//    outArgs.set_g(2,voltLimF);
//    model->evalModel(inArgs,outArgs);
//  }
//  Epetra_Vector& fb = *fbar;
//  for (int i=0 ; i<2*L ; ++i) {
//    std::cout << "fb[" << i << "] = " << fb[i] << std::endl;
//  }
//  std::cout << std::endl;
//  // Assign fi its values:
//  fi.clear();
//  fi.resize(L);
//
//  // Enumerate the possible combinations for x:
//  Array<Array<int> > xpos;
//  xpos.push_back(tuple<int>( 0, 1, 2 ));
//  xpos.push_back(tuple<int>( 0, 2, 1 ));
//  xpos.push_back(tuple<int>( 1, 0, 2 ));
//  xpos.push_back(tuple<int>( 1, 2, 0 ));
//  xpos.push_back(tuple<int>( 2, 0, 1 ));
//  xpos.push_back(tuple<int>( 2, 1, 0 ));
//  // Enumerate the possible combinations for f:
//  Array<Array<double> > fpos;
//  fpos.push_back(tuple<double>( -2.0/R, 3.0, +2.0/R+2.0 ));
//  fpos.push_back(tuple<double>( -1.0/R, 2.0, +1.0/R+3.0 ));
//  fpos.push_back(tuple<double>( -1.0/R, 3.0, +1.0/R+1.0 ));
//  fpos.push_back(tuple<double>( +1.0/R, 1.0, -1.0/R+3.0 ));
//  fpos.push_back(tuple<double>( +1.0/R, 2.0, -1.0/R+1.0 ));
//  fpos.push_back(tuple<double>( +2.0/R, 1.0, -2.0/R+2.0 ));
//
//  for (int x_ind=0 ; x_ind<3 ; ++x_ind) {
//    for (int f_ind=0 ; f_ind<6 ; ++f_ind) {
//      if (fb[x_ind] == fpos[f_ind][2]) {
//        xi = xpos[f_ind];
//        fi[x_ind] = 2;
//        // Now we need to determine the other two fi values
//        int x_indp  = (x_ind+1) % 3;
//        int x_indpp = (x_ind+2) % 3;
//        if (fb[x_indp] == fpos[f_ind][1]) {
//          fi[x_indp] = 1;
//          fi[x_indpp] = 0;
//        }
//        else {
//          fi[x_indp] = 0;
//          fi[x_indpp] = 1;
//        }
//      }
//      break;
//    }
//  }
//  // Test that we got it right.
//  //TEUCHOS_ASSERT( fb[fi[0]] == C );
//  for (int i=0 ; i<L ; ++i) {
//    xi.push_back(xi[i]+L);
//    fi.push_back(fi[i]+L);
//  }
//  for (int i = 0 ; i<2*L ; ++i) {
//    std::cout << "xi[" << i << "] = " << xi[i] << std::endl;
//  }
//  std::cout << std::endl;
//  for (int i = 0 ; i<2*L ; ++i) {
//    std::cout << "fi[" << i << "] = " << fi[i] << std::endl;
//  }
//}

TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, local_getVariableIndices ) {
  std::vector<std::string> inNames;
  std::vector<std::string> outNames;
  inNames.push_back("Hello");
  inNames.push_back("World");
  inNames.push_back("foo");

  outNames.push_back("World");
  outNames.push_back("foo");
  outNames.push_back("Hello");

  std::vector<int> indices = getVariableIndices(inNames,outNames);
  TEST_ASSERT( indices.size() == 3 );
  TEST_EQUALITY_CONST( indices[0], 1 );
  TEST_EQUALITY_CONST( indices[1], 2 );
  TEST_EQUALITY_CONST( indices[2], 0 );
}

#if 1 
TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, RCLoad ) {
  int iargs = 2;
  char *cargs[iargs];
  cargs[0] = "Xyce";
  cargs[1] = "modelEvaluatorTest.cir";
  // Use simple RC circuit netlist as input
  // Set up t,x,xdot,s,sdot,f,W
  RCP<N_ANP_ModelEvaluator> model = rcp(new N_ANP_ModelEvaluator());
  model->initialize(iargs,cargs);
  RCP<const Epetra_Map> x_map = model->get_x_map();
  RCP<const Epetra_Map> f_map = model->get_f_map();
  RCP<const Epetra_Map> s_map = model->get_p_map(0);
  RCP<const Epetra_Map> sdot_map = model->get_p_map(1);
  RCP<Epetra_Operator> W = model->create_W();
  RCP<Epetra_Operator> W_exact = model->create_W();

  Array<int> xi,fi;
  //determineIndices(xi,fi);
  // Get indexing correct for the solution vector and the F vector:
  //std::vector<std::string> names = model->getVariableNames();
  //int L = Teuchos::as<int>(names.size());
  //std::vector<std::string> myNames;
  //{
  //  myNames.push_back("V(2)");
  //  myNames.push_back("I(VSRC)");
  //  myNames.push_back("V(1)");
  //}
  //std::vector<int> xi = getVariableIndices(names,myNames);
  //for (int i=0 ; i<L ; ++i) {
  //  xi.push_back(xi[i]+L);
  //}
  //std::vector<int> fi;
  //fi = xi;
  //{
  //  for (int i=0; i<2*L ; ++i) {
  //    fi.push_back(i);
  //  }
  //}

  /*
	std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "W_exact->epetraObj().Graph() = " << std::endl;
  W_exact->epetraObj().Graph().Print(std::cout);
  std::cout << "-----------------------------------------------" << std::endl;
  */
  double alpha,beta;
  double t;
  RCP<Epetra_Vector> xbar = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> xbardot = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> fbar = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> s = rcp(new Epetra_Vector(*s_map)); 
  RCP<Epetra_Vector> sdot = rcp(new Epetra_Vector(*sdot_map)); 
  RCP<Epetra_Vector> voltLimF = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> voltLimQ = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> fbar_exact = rcp(new Epetra_Vector(*f_map));
  RCP<Epetra_Vector> fq_exact = rcp(new Epetra_Vector(*f_map));
  {
    Epetra_Vector& xb = *xbar;
    Epetra_Vector& xbdot = *xbardot;
    Epetra_Vector& sb = *s;
    Epetra_Vector& sbdot = *sdot;
    Epetra_Vector& fb_exact = *fbar_exact;
    Epetra_Vector& fq_e = *fq_exact;
    t = 0.0001; // 0.00025;
    //two blocks of three elements each.
    xb[0] = 2.0; // v_2
    xb[1] = 3.0; // i_vsrc
    xb[2] = 4.0; // v_1
    double v2 = xb[0]; // 2.0
    double i_vsrc = xb[1]; // 3.0
    double v1 = xb[2]; // 4.0

    xb[3] = 5.0; // z
    xb[4] = 6.0; 
    xb[5] = 7.0; 

    xbdot[0] = 8.0; // v_2 dot
    xbdot[1] = 9.0; // i_vsrc dot
    xbdot[2] = 10.0; // v_1 dot

    xbdot[3] = 11.0; // z dot
    xbdot[4] = 12.0; 
    xbdot[5] = 13.0; 

    sb[0] = 14.0; // s

    sbdot[0] = 17.0; // sdot

    alpha = 20;
    beta = 30;

    double C = 1.0e-7;
    double R = 1.0e6;
    double sinfreq = 1000;
    double pi = 4.0*atan(1.0);
    double vsrc = 12.0*sin(2*pi*sinfreq*t);
    // f:
    fq_e[0] = (v2-v1)/R;
    fq_e[1] = v1-vsrc;
    fq_e[2] = (v1-v2)/R+i_vsrc;
    // q:
    fq_e[3] = C*v2;
    fq_e[4] = 0.0;
    fq_e[5] = 0.0;

    fb_exact[0] = xbdot[3] + fq_e[0];
    fb_exact[1] = xbdot[4] + fq_e[1];
    fb_exact[2] = xbdot[5] + fq_e[2];
    fb_exact[3] = fq_e[3] - xb[3];
    fb_exact[4] = fq_e[4] - xb[4];
    fb_exact[5] = fq_e[5] - xb[5];

    // df/dx:
    // df0/dv2, df0/di, df0/dv1 = [ 1/R, 0, -1/R ]
    // df1/dv2, df1/di, df1/dv1 = [ 0,   0,  1   ]
    // df2/dv2, df2/di, df2/dv1 = [-1/R, 1,  1/R ]

    // dq/dx:
    // dq0/dv2, dq0/di, dq0/dv1 = [ C, 0, 0 ]
    // dq1/dv2, dq1/di, dq1/dv1 = [ 0, 0, 0 ]
    // dq2/dv2, dq2/di, dq2/dv1 = [ 0, 0, 0 ]

    // W = [ beta*df/dx , alpha*I ]
    //     [ beta*dq/dx , -beta*I  ]
    std::string label = "N_LAS_BlockMatrix";
    RCP<N_LAS_BlockMatrix> bW_exact = Teuchos::get_extra_data<RCP<N_LAS_BlockMatrix> >(W_exact,label);

    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "bW_exact->epetraObj().Graph() = " << std::endl;
    //bW_exact->epetraObj().Graph().Print(std::cout);
    //std::cout << "-----------------------------------------------" << std::endl;
    bW_exact->put(0.0);
    {
      N_LAS_Matrix & ul = bW_exact->block(0,0);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ul.epetraObj().Graph() = " << std::endl;
      //ul.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
     
      /* Added by Prateek while playing with the code	*/
      //N_LAS_Matrix & u2 = bW_exact->block(1,0);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "u2.epetraObj().Graph() = " << std::endl;
      //u2.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;

      int row;
      int length;
      double coeffs[3];
      int colIndices[3];
      row = 0;
      length=2;
      coeffs[0] = beta/R;
      coeffs[1] = -beta/R;
      colIndices[0] = 0;
      colIndices[1] = 2;
      ul.putRow(row,length,coeffs,colIndices);
      row = 1;
      length = 1;
      coeffs[0] = beta;
      colIndices[0] = 2;
      ul.putRow(row,length,coeffs,colIndices);
      row = 2;
      length=3;
      coeffs[0] = -beta/R;
      coeffs[1] = beta;
      coeffs[2] = beta/R;
      colIndices[0] = 0;
      colIndices[1] = 1;
      colIndices[2] = 2;
      ul.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ul = " << std::endl;
      //ul.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    {
      N_LAS_Matrix & ur = bW_exact->block(0,1);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ur.epetraObj().Graph() = " << std::endl;
      //ur.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
      int row;
      int length = 1;
      double coeffs[1];
      int colIndices[1];
      row = 0;
      coeffs[0] = alpha;
      colIndices[0] = 0;
      ur.putRow(row,length,coeffs,colIndices);
      row = 1;
      colIndices[0] = 1;
      ur.putRow(row,length,coeffs,colIndices);
      row = 2;
      colIndices[0] = 2;
      ur.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ur = " << std::endl;
      //ur.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    {
      N_LAS_Matrix & ll = bW_exact->block(1,0);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ll.epetraObj().Graph() = " << std::endl;
      //ll.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
      int row;
      int length;
      double coeffs[3];
      int colIndices[3];
      row = 0;
      length = 1;
      coeffs[0] = beta*C;
      colIndices[0] = 0;
      ll.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ll = " << std::endl;
      //ll.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    {
      N_LAS_Matrix & lr = bW_exact->block(1,1);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "lr.epetraObj().Graph() = " << std::endl;
      //lr.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
      int row;
      int length = 1;
      double coeffs[1];
      int colIndices[1];
      row = 0;
      coeffs[0] = -beta;
      colIndices[0] = 0;
      lr.putRow(row,length,coeffs,colIndices);
      row = 1;
      colIndices[0] = 1;
      lr.putRow(row,length,coeffs,colIndices);
      row = 2;
      colIndices[0] = 2;
      lr.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "lr = " << std::endl;
      //lr.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    bW_exact->fillComplete();
    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "bW_exact->epetraObj().Graph() = " << std::endl;
    //bW_exact->epetraObj().Graph().Print(std::cout);
    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "bW_exact = " << std::endl;
    //bW_exact->printPetraObject();
    //std::cout << "-----------------------------------------------" << std::endl;

  }

  /* Why does Todd do two evamodels???? */
  // First we call the model to load the state vector s.
  {
    EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
    inArgs.set_t(t);
    inArgs.set_x(xbar);
    inArgs.set_x_dot(xbardot);
    //inArgs.set_p(0,s);
    //inArgs.set_p(1,sdot);
    EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
    
		/* Question: what is g?? */
		outArgs.set_g(0,s);
    //fbar->PutScalar(0.0);
    //outArgs.set_f(fbar);

    // evalModel
    model->evalModel(inArgs,outArgs);
  }
  // Verify model loaded s correctly:
  // TODO
  // Second we "calculate" sdot
  // Third we call the model to load f
  {
    EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
    inArgs.set_t(t);
    inArgs.set_x(xbar);
    inArgs.set_x_dot(xbardot);
    inArgs.set_p(0,s);
    inArgs.set_p(1,sdot);
    inArgs.set_alpha(alpha);
    inArgs.set_beta(beta);
    EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
    fbar->PutScalar(0.0);
    outArgs.set_f(fbar);
    outArgs.set_W(W);
    outArgs.set_g(1,voltLimQ);
    outArgs.set_g(2,voltLimF);

    // evalModel
    model->evalModel(inArgs,outArgs);
  }

  //std::cout << "xbar = " << *xbar << std::endl;
  //std::cout << "xbardot = " << *xbardot << std::endl;
  //std::cout << "s = " << *s << std::endl;
  //std::cout << "sdot = " << *sdot << std::endl;
  //std::cout << "fbar = " << *fbar << std::endl;
  //std::cout << "fq_exact = " << *fq_exact << std::endl;
  //std::cout << "fbar_exact = " << *fbar_exact << std::endl;
  
  // Verify model loaded F correctly
  double tol = 1.0e-10;
  TEST_FLOATING_EQUALITY( (*fbar)[0], (*fbar_exact)[0], tol );
  TEST_FLOATING_EQUALITY( (*fbar)[1], (*fbar_exact)[1], tol );
  TEST_FLOATING_EQUALITY( (*fbar)[2], (*fbar_exact)[2], tol );

  TEST_FLOATING_EQUALITY( (*fbar)[3], (*fbar_exact)[3], tol );
  TEST_FLOATING_EQUALITY( (*fbar)[4], (*fbar_exact)[4], tol );
  TEST_FLOATING_EQUALITY( (*fbar)[5], (*fbar_exact)[5], tol );

  // Verify model loaded W correctly
  RCP<Epetra_CrsMatrix> W_matrix_exact = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_exact,true);
  //std::cout << "W_matrix_exact = " << std::endl;
  //W_matrix_exact->Print(std::cout);
  //std::cout << "-----------------------------------------------" << std::endl;
  RCP<Epetra_CrsMatrix> W_matrix = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W,true);
  //std::cout << "W_matrix = " << std::endl;
  //W_matrix->Print(std::cout);
  //std::cout << "-----------------------------------------------" << std::endl;

  for (int i=0 ; i<6 ; ++i) {
    int GlobalRow = i;
    if (i>2) {
      GlobalRow = (i-3)+10;
    }
    int Length = 6;
    int NumEntries;
    std::vector<double> Values; Values.resize(Length);
    std::vector<int> Indices ; Indices.resize(Length);
    W_matrix->ExtractGlobalRowCopy(GlobalRow, Length, NumEntries, &Values[0], &Indices[0]);

    int NumEntries_exact;
    std::vector<double> Values_exact; Values_exact.resize(Length);
    std::vector<int> Indices_exact ; Indices_exact.resize(Length);
    W_matrix_exact->ExtractGlobalRowCopy(GlobalRow, Length, NumEntries_exact, &Values_exact[0], &Indices_exact[0]);

    TEST_EQUALITY( NumEntries, NumEntries_exact );
    for (int j=0 ; j<NumEntries ; ++j) {
      TEST_FLOATING_EQUALITY( Values[j], Values_exact[j], tol );
    }
  }
}
#endif

// Verify model loaded W correctly at differrent time instants t_1, t_2, t_3
// This copy paste of above test we change the time below to get W at different
// timee
#if 1 
TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, RCLoad_diff_times_same_x0 ) {
  int iargs = 2;
  char *cargs[iargs];
  cargs[0] = "Xyce";
  cargs[1] = "modelEvaluatorTest.cir";
  // Use simple RC circuit netlist as input
  // Set up t,x,xdot,s,sdot,f,W
  RCP<N_ANP_ModelEvaluator> model = rcp(new N_ANP_ModelEvaluator());
  model->initialize(iargs,cargs);
  RCP<const Epetra_Map> x_map = model->get_x_map();
  RCP<const Epetra_Map> f_map = model->get_f_map();
  RCP<const Epetra_Map> s_map = model->get_p_map(0);
  RCP<const Epetra_Map> sdot_map = model->get_p_map(1);
  RCP<Epetra_Operator> W = model->create_W();
  RCP<Epetra_Operator> W_exact = model->create_W();
  
	double alpha,beta;
  double t;
  RCP<Epetra_Vector> xbar = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> xbardot = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> fbar = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> s = rcp(new Epetra_Vector(*s_map)); 
  RCP<Epetra_Vector> sdot = rcp(new Epetra_Vector(*sdot_map)); 
  RCP<Epetra_Vector> voltLimF = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> voltLimQ = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> fbar_exact = rcp(new Epetra_Vector(*f_map));
  RCP<Epetra_Vector> fq_exact = rcp(new Epetra_Vector(*f_map));
  {
    Epetra_Vector& xb = *xbar;
    Epetra_Vector& xbdot = *xbardot;
    Epetra_Vector& sb = *s;
    Epetra_Vector& sbdot = *sdot;
    Epetra_Vector& fb_exact = *fbar_exact;
    Epetra_Vector& fq_e = *fq_exact;
    t = 0.0001; // 0.00025;
    //two blocks of three elements each.
    xb[0] = 2.0; // v_2
    xb[1] = 3.0; // i_vsrc
    xb[2] = 4.0; // v_1
    double v2 = xb[0]; // 2.0
    double i_vsrc = xb[1]; // 3.0
    double v1 = xb[2]; // 4.0

    xb[3] = 5.0; // z
    xb[4] = 6.0; 
    xb[5] = 7.0; 

    xbdot[0] = 8.0; // v_2 dot
    xbdot[1] = 9.0; // i_vsrc dot
    xbdot[2] = 10.0; // v_1 dot

    xbdot[3] = 11.0; // z dot
    xbdot[4] = 12.0; 
    xbdot[5] = 13.0; 

    sb[0] = 14.0; // s

    sbdot[0] = 17.0; // sdot

    alpha = 20;
    beta = 30;

    double C = 1.0e-7;
    double R = 1.0e6;
    double sinfreq = 1000;
    double pi = 4.0*atan(1.0);
    double vsrc = 12.0*sin(2*pi*sinfreq*t);
    // f:
    fq_e[0] = (v2-v1)/R;
    fq_e[1] = v1-vsrc;
    fq_e[2] = (v1-v2)/R+i_vsrc;
    // q:
    fq_e[3] = C*v2;
    fq_e[4] = 0.0;
    fq_e[5] = 0.0;

    fb_exact[0] = xbdot[3] + fq_e[0];
    fb_exact[1] = xbdot[4] + fq_e[1];
    fb_exact[2] = xbdot[5] + fq_e[2];
    fb_exact[3] = fq_e[3] - xb[3];
    fb_exact[4] = fq_e[4] - xb[4];
    fb_exact[5] = fq_e[5] - xb[5];

    // df/dx:
    // df0/dv2, df0/di, df0/dv1 = [ 1/R, 0, -1/R ]
    // df1/dv2, df1/di, df1/dv1 = [ 0,   0,  1   ]
    // df2/dv2, df2/di, df2/dv1 = [-1/R, 1,  1/R ]

    // dq/dx:
    // dq0/dv2, dq0/di, dq0/dv1 = [ C, 0, 0 ]
    // dq1/dv2, dq1/di, dq1/dv1 = [ 0, 0, 0 ]
    // dq2/dv2, dq2/di, dq2/dv1 = [ 0, 0, 0 ]

    // W = [ beta*df/dx , alpha*I ]
    //     [ beta*dq/dx , -beta*I  ]
    std::string label = "N_LAS_BlockMatrix";
    RCP<N_LAS_BlockMatrix> bW_exact = Teuchos::get_extra_data<RCP<N_LAS_BlockMatrix> >(W_exact,label);

    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "bW_exact->epetraObj().Graph() = " << std::endl;
    //bW_exact->epetraObj().Graph().Print(std::cout);
    //std::cout << "-----------------------------------------------" << std::endl;
    bW_exact->put(0.0);
    {
      N_LAS_Matrix & ul = bW_exact->block(0,0);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ul.epetraObj().Graph() = " << std::endl;
      //ul.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
     
      /* Added by Prateek while playing with the code	*/
      N_LAS_Matrix & u2 = bW_exact->block(1,0);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "u2.epetraObj().Graph() = " << std::endl;
      //u2.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;

      int row;
      int length;
      double coeffs[3];
      int colIndices[3];
      row = 0;
      length=2;
      coeffs[0] = beta/R;
      coeffs[1] = -beta/R;
      colIndices[0] = 0;
      colIndices[1] = 2;
      ul.putRow(row,length,coeffs,colIndices);
      row = 1;
      length = 1;
      coeffs[0] = beta;
      colIndices[0] = 2;
      ul.putRow(row,length,coeffs,colIndices);
      row = 2;
      length=3;
      coeffs[0] = -beta/R;
      coeffs[1] = beta;
      coeffs[2] = beta/R;
      colIndices[0] = 0;
      colIndices[1] = 1;
      colIndices[2] = 2;
      ul.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ul = " << std::endl;
      //ul.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    {
      N_LAS_Matrix & ur = bW_exact->block(0,1);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ur.epetraObj().Graph() = " << std::endl;
      //ur.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
      int row;
      int length = 1;
      double coeffs[1];
      int colIndices[1];
      row = 0;
      coeffs[0] = alpha;
      colIndices[0] = 0;
      ur.putRow(row,length,coeffs,colIndices);
      row = 1;
      colIndices[0] = 1;
      ur.putRow(row,length,coeffs,colIndices);
      row = 2;
      colIndices[0] = 2;
      ur.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ur = " << std::endl;
      //ur.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    {
      N_LAS_Matrix & ll = bW_exact->block(1,0);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ll.epetraObj().Graph() = " << std::endl;
      //ll.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
      int row;
      int length;
      double coeffs[3];
      int colIndices[3];
      row = 0;
      length = 1;
      coeffs[0] = beta*C;
      colIndices[0] = 0;
      ll.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ll = " << std::endl;
      //ll.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    {
      N_LAS_Matrix & lr = bW_exact->block(1,1);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "lr.epetraObj().Graph() = " << std::endl;
      //lr.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
      int row;
      int length = 1;
      double coeffs[1];
      int colIndices[1];
      row = 0;
      coeffs[0] = -beta;
      colIndices[0] = 0;
      lr.putRow(row,length,coeffs,colIndices);
      row = 1;
      colIndices[0] = 1;
      lr.putRow(row,length,coeffs,colIndices);
      row = 2;
      colIndices[0] = 2;
      lr.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "lr = " << std::endl;
      //lr.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    bW_exact->fillComplete();
    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "bW_exact->epetraObj().Graph() = " << std::endl;
    //bW_exact->epetraObj().Graph().Print(std::cout);
    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "bW_exact = " << std::endl;
    //bW_exact->printPetraObject();
    //std::cout << "-----------------------------------------------" << std::endl;

  }


	for ( int ti=0;ti<3; ++ti) {
 		 
		t = ti + 0.0001;	
		//std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "time is t = " << t << std::endl;
    //std::cout << "index is ti = " << ti << std::endl;
    //std::cout << "-----------------------------------------------" << std::endl;
    
		Epetra_Vector& xb = *xbar;
    Epetra_Vector& xbdot = *xbardot;
    Epetra_Vector& sb = *s;
    Epetra_Vector& sbdot = *sdot;
    Epetra_Vector& fb_exact = *fbar_exact;
    Epetra_Vector& fq_e = *fq_exact;
    //t = 0.0001; // 0.00025;
    //two blocks of three elements each.
    xb[0] = 2.0 ; // v_2
    xb[1] = 3.0; // i_vsrc
    xb[2] = 4.0; // v_1
    double v2 = xb[0]; // 2.0
    double i_vsrc = xb[1]; // 3.0
    double v1 = xb[2]; // 4.0

    xb[3] = 5.0; // z
    xb[4] = 6.0; 
    xb[5] = 7.0; 

    xbdot[0] = 8.0; // v_2 dot
    xbdot[1] = 9.0; // i_vsrc dot
    xbdot[2] = 10.0; // v_1 dot

    xbdot[3] = 11.0; // z dot
    xbdot[4] = 12.0; 
    xbdot[5] = 13.0; 

    sb[0] = 14.0; // s

    sbdot[0] = 17.0; // sdot

    alpha = 20;
    beta = 30;

    double C = 1.0e-7;
    double R = 1.0e6;
    double sinfreq = 1000;
    double pi = 4.0*atan(1.0);
    double vsrc = 12.0*sin(2*pi*sinfreq*t);
    // f:
    fq_e[0] = (v2-v1)/R;
    fq_e[1] = v1-vsrc;
    fq_e[2] = (v1-v2)/R+i_vsrc;
    // q:
    fq_e[3] = C*v2;
    fq_e[4] = 0.0;
    fq_e[5] = 0.0;

    fb_exact[0] = xbdot[3] + fq_e[0];
    fb_exact[1] = xbdot[4] + fq_e[1];
    fb_exact[2] = xbdot[5] + fq_e[2];
    fb_exact[3] = fq_e[3] - xb[3];
    fb_exact[4] = fq_e[4] - xb[4];
    fb_exact[5] = fq_e[5] - xb[5];


		/* Why does Todd do two evamodels???? */
  	// First we call the model to load the state vector s.
  	{
  	  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  	  inArgs.set_t(t);
  	  inArgs.set_x(xbar);
  	  inArgs.set_x_dot(xbardot);
  	  //inArgs.set_p(0,s);
  	  //inArgs.set_p(1,sdot);
  	  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
  	  
			/* Question: what is g?? */
			outArgs.set_g(0,s);
  	  //fbar->PutScalar(0.0);
  	  //outArgs.set_f(fbar);

  	  // evalModel
  	  model->evalModel(inArgs,outArgs);
  	}
  	// Verify model loaded s correctly:
  	// TODO
  	// Second we "calculate" sdot
  	// Third we call the model to load f
  	{
  	  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  	  inArgs.set_t(t);
  	  inArgs.set_x(xbar);
  	  inArgs.set_x_dot(xbardot);
  	  inArgs.set_p(0,s);
  	  inArgs.set_p(1,sdot);
  	  inArgs.set_alpha(alpha);
  	  inArgs.set_beta(beta);
  	  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
  	  fbar->PutScalar(0.0);
  	  outArgs.set_f(fbar);
  	  outArgs.set_W(W);
  	  outArgs.set_g(1,voltLimQ);
  	  outArgs.set_g(2,voltLimF);

  	  // evalModel
  	  model->evalModel(inArgs,outArgs);
  	}

  	//std::cout << "xbar = " << *xbar << std::endl;
  	//std::cout << "xbardot = " << *xbardot << std::endl;
  	//std::cout << "s = " << *s << std::endl;
  	//std::cout << "sdot = " << *sdot << std::endl;
  	//std::cout << "fbar = " << *fbar << std::endl;
  	//std::cout << "fq_exact = " << *fq_exact << std::endl;
  	//std::cout << "fbar_exact = " << *fbar_exact << std::endl;

  	// Verify model loaded F correctly
  	double tol = 1.0e-10;
  	TEST_FLOATING_EQUALITY( (*fbar)[0], (*fbar_exact)[0], tol );
  	TEST_FLOATING_EQUALITY( (*fbar)[1], (*fbar_exact)[1], tol );
  	TEST_FLOATING_EQUALITY( (*fbar)[2], (*fbar_exact)[2], tol );

  	TEST_FLOATING_EQUALITY( (*fbar)[3], (*fbar_exact)[3], tol );
  	TEST_FLOATING_EQUALITY( (*fbar)[4], (*fbar_exact)[4], tol );
  	TEST_FLOATING_EQUALITY( (*fbar)[5], (*fbar_exact)[5], tol );

  	// Verify model loaded W correctly
  	RCP<Epetra_CrsMatrix> W_matrix_exact = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_exact,true);
  	//std::cout << "W_matrix_exact = " << std::endl;
  	//W_matrix_exact->Print(std::cout);
  	//std::cout << "-----------------------------------------------" << std::endl;
  	RCP<Epetra_CrsMatrix> W_matrix = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W,true);
  	//std::cout << "W_matrix = " << std::endl;
  	//W_matrix->Print(std::cout);
  	//std::cout << "-----------------------------------------------" << std::endl;

  	for (int i=0 ; i<6 ; ++i) {
  	  int GlobalRow = i;
  	  if (i>2) {
  	    GlobalRow = (i-3)+10;
  	  }
  	  int Length = 6;
  	  int NumEntries;
  	  std::vector<double> Values; Values.resize(Length);
  	  std::vector<int> Indices ; Indices.resize(Length);
  	  W_matrix->ExtractGlobalRowCopy(GlobalRow, Length, NumEntries, &Values[0], &Indices[0]);

  	  int NumEntries_exact;
  	  std::vector<double> Values_exact; Values_exact.resize(Length);
  	  std::vector<int> Indices_exact ; Indices_exact.resize(Length);
  	  W_matrix_exact->ExtractGlobalRowCopy(GlobalRow, Length, NumEntries_exact, &Values_exact[0], &Indices_exact[0]);

  	  TEST_EQUALITY( NumEntries, NumEntries_exact );
  	  for (int j=0 ; j<NumEntries ; ++j) {
  	    TEST_FLOATING_EQUALITY( Values[j], Values_exact[j], tol );
  	  }
  	}
	}
}
#endif


#if 1 
TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, RCLoad_diff_x0_same_time ) {
  int iargs = 2;
  char *cargs[iargs];
  cargs[0] = "Xyce";
  cargs[1] = "modelEvaluatorTest.cir";
  // Use simple RC circuit netlist as input
  // Set up t,x,xdot,s,sdot,f,W
  RCP<N_ANP_ModelEvaluator> model = rcp(new N_ANP_ModelEvaluator());
  model->initialize(iargs,cargs);
  RCP<const Epetra_Map> x_map = model->get_x_map();
  RCP<const Epetra_Map> f_map = model->get_f_map();
  RCP<const Epetra_Map> s_map = model->get_p_map(0);
  RCP<const Epetra_Map> sdot_map = model->get_p_map(1);
  RCP<Epetra_Operator> W = model->create_W();
  RCP<Epetra_Operator> W_exact = model->create_W();
  
	double alpha,beta;
  double t;
  RCP<Epetra_Vector> xbar = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> xbardot = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> fbar = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> s = rcp(new Epetra_Vector(*s_map)); 
  RCP<Epetra_Vector> sdot = rcp(new Epetra_Vector(*sdot_map)); 
  RCP<Epetra_Vector> voltLimF = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> voltLimQ = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> fbar_exact = rcp(new Epetra_Vector(*f_map));
  RCP<Epetra_Vector> fq_exact = rcp(new Epetra_Vector(*f_map));
  {
    Epetra_Vector& xb = *xbar;
    Epetra_Vector& xbdot = *xbardot;
    Epetra_Vector& sb = *s;
    Epetra_Vector& sbdot = *sdot;
    Epetra_Vector& fb_exact = *fbar_exact;
    Epetra_Vector& fq_e = *fq_exact;
    t = 0.0001; // 0.00025;
    //two blocks of three elements each.
    xb[0] = 2.0; // v_2
    xb[1] = 3.0; // i_vsrc
    xb[2] = 4.0; // v_1
    double v2 = xb[0]; // 2.0
    double i_vsrc = xb[1]; // 3.0
    double v1 = xb[2]; // 4.0

    xb[3] = 5.0; // z
    xb[4] = 6.0; 
    xb[5] = 7.0; 

    xbdot[0] = 8.0; // v_2 dot
    xbdot[1] = 9.0; // i_vsrc dot
    xbdot[2] = 10.0; // v_1 dot

    xbdot[3] = 11.0; // z dot
    xbdot[4] = 12.0; 
    xbdot[5] = 13.0; 

    sb[0] = 14.0; // s

    sbdot[0] = 17.0; // sdot

    alpha = 20;
    beta = 30;

    double C = 1.0e-7;
    double R = 1.0e6;
    double sinfreq = 1000;
    double pi = 4.0*atan(1.0);
    double vsrc = 12.0*sin(2*pi*sinfreq*t);
    // f:
    fq_e[0] = (v2-v1)/R;
    fq_e[1] = v1-vsrc;
    fq_e[2] = (v1-v2)/R+i_vsrc;
    // q:
    fq_e[3] = C*v2;
    fq_e[4] = 0.0;
    fq_e[5] = 0.0;

    fb_exact[0] = xbdot[3] + fq_e[0];
    fb_exact[1] = xbdot[4] + fq_e[1];
    fb_exact[2] = xbdot[5] + fq_e[2];
    fb_exact[3] = fq_e[3] - xb[3];
    fb_exact[4] = fq_e[4] - xb[4];
    fb_exact[5] = fq_e[5] - xb[5];

    // df/dx:
    // df0/dv2, df0/di, df0/dv1 = [ 1/R, 0, -1/R ]
    // df1/dv2, df1/di, df1/dv1 = [ 0,   0,  1   ]
    // df2/dv2, df2/di, df2/dv1 = [-1/R, 1,  1/R ]

    // dq/dx:
    // dq0/dv2, dq0/di, dq0/dv1 = [ C, 0, 0 ]
    // dq1/dv2, dq1/di, dq1/dv1 = [ 0, 0, 0 ]
    // dq2/dv2, dq2/di, dq2/dv1 = [ 0, 0, 0 ]

    // W = [ beta*df/dx , alpha*I ]
    //     [ beta*dq/dx , -beta*I  ]
    std::string label = "N_LAS_BlockMatrix";
    RCP<N_LAS_BlockMatrix> bW_exact = Teuchos::get_extra_data<RCP<N_LAS_BlockMatrix> >(W_exact,label);

    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "bW_exact->epetraObj().Graph() = " << std::endl;
    //bW_exact->epetraObj().Graph().Print(std::cout);
    //std::cout << "-----------------------------------------------" << std::endl;
    bW_exact->put(0.0);
    {
      N_LAS_Matrix & ul = bW_exact->block(0,0);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ul.epetraObj().Graph() = " << std::endl;
      //ul.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
     
      /* Added by Prateek while playing with the code	*/
      N_LAS_Matrix & u2 = bW_exact->block(1,0);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "u2.epetraObj().Graph() = " << std::endl;
      //u2.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;

      int row;
      int length;
      double coeffs[3];
      int colIndices[3];
      row = 0;
      length=2;
      coeffs[0] = beta/R;
      coeffs[1] = -beta/R;
      colIndices[0] = 0;
      colIndices[1] = 2;
      ul.putRow(row,length,coeffs,colIndices);
      row = 1;
      length = 1;
      coeffs[0] = beta;
      colIndices[0] = 2;
      ul.putRow(row,length,coeffs,colIndices);
      row = 2;
      length=3;
      coeffs[0] = -beta/R;
      coeffs[1] = beta;
      coeffs[2] = beta/R;
      colIndices[0] = 0;
      colIndices[1] = 1;
      colIndices[2] = 2;
      ul.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ul = " << std::endl;
      //ul.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    {
      N_LAS_Matrix & ur = bW_exact->block(0,1);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ur.epetraObj().Graph() = " << std::endl;
      //ur.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
      int row;
      int length = 1;
      double coeffs[1];
      int colIndices[1];
      row = 0;
      coeffs[0] = alpha;
      colIndices[0] = 0;
      ur.putRow(row,length,coeffs,colIndices);
      row = 1;
      colIndices[0] = 1;
      ur.putRow(row,length,coeffs,colIndices);
      row = 2;
      colIndices[0] = 2;
      ur.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ur = " << std::endl;
      //ur.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    {
      N_LAS_Matrix & ll = bW_exact->block(1,0);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ll.epetraObj().Graph() = " << std::endl;
      //ll.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
      int row;
      int length;
      double coeffs[3];
      int colIndices[3];
      row = 0;
      length = 1;
      coeffs[0] = beta*C;
      colIndices[0] = 0;
      ll.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ll = " << std::endl;
      //ll.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    {
      N_LAS_Matrix & lr = bW_exact->block(1,1);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "lr.epetraObj().Graph() = " << std::endl;
      //lr.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
      int row;
      int length = 1;
      double coeffs[1];
      int colIndices[1];
      row = 0;
      coeffs[0] = -beta;
      colIndices[0] = 0;
      lr.putRow(row,length,coeffs,colIndices);
      row = 1;
      colIndices[0] = 1;
      lr.putRow(row,length,coeffs,colIndices);
      row = 2;
      colIndices[0] = 2;
      lr.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "lr = " << std::endl;
      //lr.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    bW_exact->fillComplete();
    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "bW_exact->epetraObj().Graph() = " << std::endl;
    //bW_exact->epetraObj().Graph().Print(std::cout);
    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "bW_exact = " << std::endl;
    //bW_exact->printPetraObject();
    //std::cout << "-----------------------------------------------" << std::endl;

  }


	for ( int xi= 0 ; xi < 3; ++xi) {
 		 
		Epetra_Vector& xb = *xbar;
    Epetra_Vector& xbdot = *xbardot;
    Epetra_Vector& sb = *s;
    Epetra_Vector& sbdot = *sdot;
    Epetra_Vector& fb_exact = *fbar_exact;
    Epetra_Vector& fq_e = *fq_exact;
    //t = 0.0001; // 0.00025;
    //two blocks of three elements each.
    
		// use the index xi and add some values to x0 to generate different x0 in
		// each iteration 
		xb[0] = sin(xi) + 2.0; // v_2
    xb[1] = cos(xi) + 3.0; // i_vsrc
    xb[2] = sin(xi)*cos(xi) + 4.0; // v_1
    double v2 = xb[0]; // 2.0
    double i_vsrc = xb[1]; // 3.0
    double v1 = xb[2]; // 4.0

    xb[3] = 5.0; // z
    xb[4] = 6.0; 
    xb[5] = 7.0; 

    xbdot[0] = 8.0; // v_2 dot
    xbdot[1] = 9.0; // i_vsrc dot
    xbdot[2] = 10.0; // v_1 dot

    xbdot[3] = 11.0; // z dot
    xbdot[4] = 12.0; 
    xbdot[5] = 13.0; 

    sb[0] = 14.0; // s

    sbdot[0] = 17.0; // sdot

    alpha = 20;
    beta = 30;

    double C = 1.0e-7;
    double R = 1.0e6;
    double sinfreq = 1000;
    double pi = 4.0*atan(1.0);
    double vsrc = 12.0*sin(2*pi*sinfreq*t);
    // f:
    fq_e[0] = (v2-v1)/R;
    fq_e[1] = v1-vsrc;
    fq_e[2] = (v1-v2)/R+i_vsrc;
    // q:
    fq_e[3] = C*v2;
    fq_e[4] = 0.0;
    fq_e[5] = 0.0;

    fb_exact[0] = xbdot[3] + fq_e[0];
    fb_exact[1] = xbdot[4] + fq_e[1];
    fb_exact[2] = xbdot[5] + fq_e[2];
    fb_exact[3] = fq_e[3] - xb[3];
    fb_exact[4] = fq_e[4] - xb[4];
    fb_exact[5] = fq_e[5] - xb[5];


		/* Why does Todd do two evamodels???? */
  	// First we call the model to load the state vector s.
  	{
  	  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  	  inArgs.set_t(t);
  	  inArgs.set_x(xbar);
  	  inArgs.set_x_dot(xbardot);
  	  //inArgs.set_p(0,s);
  	  //inArgs.set_p(1,sdot);
  	  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
  	  
			/* Question: what is g?? */
			outArgs.set_g(0,s);
  	  //fbar->PutScalar(0.0);
  	  //outArgs.set_f(fbar);

  	  // evalModel
  	  model->evalModel(inArgs,outArgs);
  	}
  	// Verify model loaded s correctly:
  	// TODO
  	// Second we "calculate" sdot
  	// Third we call the model to load f
  	{
  	  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  	  inArgs.set_t(t);
  	  inArgs.set_x(xbar);
  	  inArgs.set_x_dot(xbardot);
  	  inArgs.set_p(0,s);
  	  inArgs.set_p(1,sdot);
  	  inArgs.set_alpha(alpha);
  	  inArgs.set_beta(beta);
  	  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
  	  fbar->PutScalar(0.0);
  	  outArgs.set_f(fbar);
  	  outArgs.set_W(W);
  	  outArgs.set_g(1,voltLimQ);
  	  outArgs.set_g(2,voltLimF);

  	  // evalModel
  	  model->evalModel(inArgs,outArgs);
  	}

  	//std::cout << "xbar = " << *xbar << std::endl;
  	//std::cout << "xbardot = " << *xbardot << std::endl;
  	//std::cout << "s = " << *s << std::endl;
  	//std::cout << "sdot = " << *sdot << std::endl;
  	//std::cout << "fbar = " << *fbar << std::endl;
  	//std::cout << "fq_exact = " << *fq_exact << std::endl;
  	//std::cout << "fbar_exact = " << *fbar_exact << std::endl;

  	// Verify model loaded F correctly
  	double tol = 1.0e-10;
  	TEST_FLOATING_EQUALITY( (*fbar)[0], (*fbar_exact)[0], tol );
  	TEST_FLOATING_EQUALITY( (*fbar)[1], (*fbar_exact)[1], tol );
  	TEST_FLOATING_EQUALITY( (*fbar)[2], (*fbar_exact)[2], tol );

  	TEST_FLOATING_EQUALITY( (*fbar)[3], (*fbar_exact)[3], tol );
  	TEST_FLOATING_EQUALITY( (*fbar)[4], (*fbar_exact)[4], tol );
  	TEST_FLOATING_EQUALITY( (*fbar)[5], (*fbar_exact)[5], tol );

  	// Verify model loaded W correctly
  	RCP<Epetra_CrsMatrix> W_matrix_exact = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_exact,true);
  	//std::cout << "W_matrix_exact = " << std::endl;
  	//W_matrix_exact->Print(std::cout);
  	//std::cout << "-----------------------------------------------" << std::endl;
  	RCP<Epetra_CrsMatrix> W_matrix = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W,true);
  	//std::cout << "W_matrix = " << std::endl;
  	//W_matrix->Print(std::cout);
  	//std::cout << "-----------------------------------------------" << std::endl;

  	for (int i=0 ; i<6 ; ++i) {
  	  int GlobalRow = i;
  	  if (i>2) {
  	    GlobalRow = (i-3)+10;
  	  }
  	  int Length = 6;
  	  int NumEntries;
  	  std::vector<double> Values; Values.resize(Length);
  	  std::vector<int> Indices ; Indices.resize(Length);
  	  W_matrix->ExtractGlobalRowCopy(GlobalRow, Length, NumEntries, &Values[0], &Indices[0]);

  	  int NumEntries_exact;
  	  std::vector<double> Values_exact; Values_exact.resize(Length);
  	  std::vector<int> Indices_exact ; Indices_exact.resize(Length);
  	  W_matrix_exact->ExtractGlobalRowCopy(GlobalRow, Length, NumEntries_exact, &Values_exact[0], &Indices_exact[0]);

  	  TEST_EQUALITY( NumEntries, NumEntries_exact );
  	  for (int j=0 ; j<NumEntries ; ++j) {
  	    TEST_FLOATING_EQUALITY( Values[j], Values_exact[j], tol );
  	  }
  	}
	}
}
#endif

# if  1
TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, DiodeCapLoad ) {
  int iargs = 2;
  char *cargs[iargs];
  cargs[0] = "Xyce";
  cargs[1] = "forward_diode_model.cir";
  // Use simple Diode-Capacitor circuit netlist as input
  // Set up t,x,xdot,s,sdot,f,W
  RCP<N_ANP_ModelEvaluator> model = rcp(new N_ANP_ModelEvaluator());
  model->initialize(iargs,cargs);
  RCP<const Epetra_Map> x_map = model->get_x_map();
  RCP<const Epetra_Map> f_map = model->get_f_map();
  RCP<const Epetra_Map> s_map = model->get_p_map(0);
  RCP<const Epetra_Map> sdot_map = model->get_p_map(1);
  RCP<Epetra_Operator> W = model->create_W();
  RCP<Epetra_Operator> W_exact = model->create_W();
  double alpha,beta;
  double t;
  RCP<Epetra_Vector> xbar = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> xbardot = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> fbar = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> s = rcp(new Epetra_Vector(*s_map)); 
  RCP<Epetra_Vector> sdot = rcp(new Epetra_Vector(*sdot_map)); 
  RCP<Epetra_Vector> voltLimF = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> voltLimQ = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> fbar_exact = rcp(new Epetra_Vector(*f_map));
  RCP<Epetra_Vector> fq_exact = rcp(new Epetra_Vector(*f_map));
  {
    Epetra_Vector& xb = *xbar;
    Epetra_Vector& xbdot = *xbardot;
    Epetra_Vector& sb = *s;
    Epetra_Vector& sbdot = *sdot;
    Epetra_Vector& fb_exact = *fbar_exact;
    Epetra_Vector& fq_e = *fq_exact;
    t = 0.0001; // 0.00025;
    //two blocks of three elements each.
    xb[0] = 0.2; // v_2
    xb[1] = 0.3; // i_vsrc
    xb[2] = 0.4; // v_1
    double v2 = xb[0]; // 0.2
    double i_vsrc = xb[1]; // 0.3
    double v1 = xb[2]; // 0.4

    xb[3] = 0.5; // z
    xb[4] = 0.6; 
    xb[5] = 0.7; 

    xbdot[0] =0.8; // v_2 dot
    xbdot[1] =0.9; // i_vsrc dot
    xbdot[2] =0.10; // v_1 dot

    xbdot[3] =0.11; // z dot
    xbdot[4] =0.12; 
    xbdot[5] =0.13; 

    sb[0] = 0.14; // s

    sbdot[0] = 0.17; // sdot

    alpha = 20;
    beta = 30;
		
		// Capacitor Parms
    double C = 1.0e-12;
		// Diode Parms: Forward model of the diode 
    double IS = 1.0e-14; // Saturation current
    double N  = 1; // Mult Factor
    double Vth  = 0.025864186; // Threshold Factor
    double sinfreq = 1000;
    double pi = 4.0*atan(1.0);
    //double vsrc = 12.0*sin(2*pi*sinfreq*t);
    double vsrc = 0.2; 
		double R = 1e6;
    fq_e[0] = -IS*(exp((v1 - v2)/(N*Vth))-1);
    fq_e[1] =  v1-vsrc;
    fq_e[2] =  IS*(exp((v1 - v2)/(N*Vth))-1)+i_vsrc;
    // q:
    fq_e[3] = C*v2;
    fq_e[4] = 0.0;
    fq_e[5] = 0.0;

    fb_exact[0] = xbdot[3] + fq_e[0];
    fb_exact[1] = xbdot[4] + fq_e[1];
    fb_exact[2] = xbdot[5] + fq_e[2];
    fb_exact[3] = fq_e[3] - xb[3];
    fb_exact[4] = fq_e[4] - xb[4];
    fb_exact[5] = fq_e[5] - xb[5];

    // df/dx:
    // df0/dv2, df0/di, df0/dv1 = [ IS * exp ((V1 -V2)/N*Vth)*1/N*Vth,   0,    -IS*exp(v1-v2/N*Vth)]
    // df1/dv2, df1/di, df1/dv1 = [ 0,                                   0,                    1   ]
    // df2/dv2, df2/di, df2/dv1 = [ IS * exp ((V1 -V2)/N*Vth)*1/N*Vth,   1,    -IS*exp(v1-v2/N*Vth)]
    // dq/dx:
    // dq0/dv2, dq0/di, dq0/dv1 = [ C, 0, 0 ]
    // dq1/dv2, dq1/di, dq1/dv1 = [ 0, 0, 0 ]
    // dq2/dv2, dq2/di, dq2/dv1 = [ 0, 0, 0 ]

    // W = [ beta*df/dx , alpha*I ]
    //     [ beta*dq/dx , -beta*I  ]
    std::string label = "N_LAS_BlockMatrix";
    RCP<N_LAS_BlockMatrix> bW_exact = Teuchos::get_extra_data<RCP<N_LAS_BlockMatrix> >(W_exact,label);

    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "bW_exact->epetraObj().Graph() = " << std::endl;
    //bW_exact->epetraObj().Graph().Print(std::cout);
    //std::cout << "-----------------------------------------------" << std::endl;

    bW_exact->put(0.0);
    {
      N_LAS_Matrix & ul = bW_exact->block(0,0);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ul.epetraObj().Graph() = " << std::endl;
      //ul.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
     
      /* Added by Prateek while playing with the code	*/
      //N_LAS_Matrix & u2 = bW_exact->block(1,0);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "u2.epetraObj().Graph() = " << std::endl;
      //u2.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;

      int row;
      int length;
      double coeffs[3];
      int colIndices[3];
      row = 0;
      length=2;
      coeffs[0] = beta*IS*exp((v1 - v2)/N/Vth)*1/N/Vth;
      coeffs[1] = -beta*IS*exp((v1 - v2)/N/Vth)*1/N/Vth;
      colIndices[0] = 0;
      colIndices[1] = 2;
      ul.putRow(row,length,coeffs,colIndices);
      row = 1;
      length = 1;
      coeffs[0] = beta;
      colIndices[0] = 2;
      ul.putRow(row,length,coeffs,colIndices);
      row = 2;
      length=3;
      coeffs[0] = -beta*IS*exp((v1 - v2)/N/Vth)*1/N/Vth;
      coeffs[1] = beta;
      coeffs[2] = beta*IS*exp((v1 - v2)/N/Vth)*1/N/Vth;
      colIndices[0] = 0;
      colIndices[1] = 1;
      colIndices[2] = 2;
      ul.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ul = " << std::endl;
      //ul.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    {
      N_LAS_Matrix & ur = bW_exact->block(0,1);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ur.epetraObj().Graph() = " << std::endl;
      //ur.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
      int row;
      int length = 1;
      double coeffs[1];
      int colIndices[1];
      row = 0;
      coeffs[0] = alpha;
      colIndices[0] = 0;
      ur.putRow(row,length,coeffs,colIndices);
      row = 1;
      colIndices[0] = 1;
      ur.putRow(row,length,coeffs,colIndices);
      row = 2;
      colIndices[0] = 2;
      ur.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ur = " << std::endl;
      //ur.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    {
      N_LAS_Matrix & ll = bW_exact->block(1,0);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ll.epetraObj().Graph() = " << std::endl;
      //ll.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
      int row;
      int length;
      double coeffs[3];
      int colIndices[3];
      row = 0;
      length = 1;
      coeffs[0] = beta*C;
      colIndices[0] = 0;
      ll.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ll = " << std::endl;
      //ll.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    {
      N_LAS_Matrix & lr = bW_exact->block(1,1);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "lr.epetraObj().Graph() = " << std::endl;
      //lr.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
      int row;
      int length = 1;
      double coeffs[1];
      int colIndices[1];
      row = 0;
      coeffs[0] = -beta;
      colIndices[0] = 0;
      lr.putRow(row,length,coeffs,colIndices);
      row = 1;
      colIndices[0] = 1;
      lr.putRow(row,length,coeffs,colIndices);
      row = 2;
      colIndices[0] = 2;
      lr.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "lr = " << std::endl;
      //lr.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    bW_exact->fillComplete();
    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "bW_exact->epetraObj().Graph() = " << std::endl;
    //bW_exact->epetraObj().Graph().Print(std::cout);
    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "bW_exact = " << std::endl;
    //bW_exact->printPetraObject();
    //std::cout << "-----------------------------------------------" << std::endl;

  }



  /* Why does Todd do two evamodels???? */
  // First we call the model to load the state vector s.
  {
    EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
    inArgs.set_t(t);
    inArgs.set_x(xbar);
    inArgs.set_x_dot(xbardot);
    //inArgs.set_p(0,s);
    //inArgs.set_p(1,sdot);
    EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
    
		/* Question: what is g?? */
		outArgs.set_g(0,s);
    //fbar->PutScalar(0.0);
    //outArgs.set_f(fbar);

    // evalModel
    model->evalModel(inArgs,outArgs);
  }
  // Verify model loaded s correctly:
  // TODO
  // Second we "calculate" sdot
  // Third we call the model to load f
  {
    EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
    inArgs.set_t(t);
    inArgs.set_x(xbar);
    inArgs.set_x_dot(xbardot);
    inArgs.set_p(0,s);
    inArgs.set_p(1,sdot);
    inArgs.set_alpha(alpha);
    inArgs.set_beta(beta);
    EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
    fbar->PutScalar(0.0);
    outArgs.set_f(fbar);
    outArgs.set_W(W);
    outArgs.set_g(1,voltLimQ);
    outArgs.set_g(2,voltLimF);

    // evalModel
    model->evalModel(inArgs,outArgs);
  }

  //std::cout << "xbar = " << *xbar << std::endl;
  //std::cout << "xbardot = " << *xbardot << std::endl;
  //std::cout << "s = " << *s << std::endl;
  //std::cout << "sdot = " << *sdot << std::endl;
  //std::cout << "fbar = " << *fbar << std::endl;
  //std::cout << "fq_exact = " << *fq_exact << std::endl;
  //std::cout << "fbar_exact = " << *fbar_exact << std::endl;

  // Verify model loaded F correctly
  double tol = 1.0e-6;
  TEST_FLOATING_EQUALITY( (*fbar)[0], (*fbar_exact)[0], tol );
  TEST_FLOATING_EQUALITY( (*fbar)[1], (*fbar_exact)[1], tol );
  TEST_FLOATING_EQUALITY( (*fbar)[2], (*fbar_exact)[2], tol );

  TEST_FLOATING_EQUALITY( (*fbar)[3], (*fbar_exact)[3], tol );
  TEST_FLOATING_EQUALITY( (*fbar)[4], (*fbar_exact)[4], tol );
  TEST_FLOATING_EQUALITY( (*fbar)[5], (*fbar_exact)[5], tol );

  // Verify model loaded W correctly
  RCP<Epetra_CrsMatrix> W_matrix_exact = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_exact,true);
  //std::cout << "W_matrix_exact = " << std::endl;
  //W_matrix_exact->Print(std::cout);
  //std::cout << "-----------------------------------------------" << std::endl;
  RCP<Epetra_CrsMatrix> W_matrix = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W,true);
  //std::cout << "W_matrix = " << std::endl;
  //W_matrix->Print(std::cout);
  //std::cout << "-----------------------------------------------" << std::endl;

  for (int i=0 ; i<6 ; ++i) {
    int GlobalRow = i;
    if (i>2) {
      GlobalRow = (i-3)+10;
    }
    int Length = 6;
    int NumEntries;
    std::vector<double> Values; Values.resize(Length);
    std::vector<int> Indices ; Indices.resize(Length);
    W_matrix->ExtractGlobalRowCopy(GlobalRow, Length, NumEntries, &Values[0], &Indices[0]);

    int NumEntries_exact;
    std::vector<double> Values_exact; Values_exact.resize(Length);
    std::vector<int> Indices_exact ; Indices_exact.resize(Length);
    W_matrix_exact->ExtractGlobalRowCopy(GlobalRow, Length, NumEntries_exact, &Values_exact[0], &Indices_exact[0]);

    TEST_EQUALITY( NumEntries, NumEntries_exact );
    for (int j=0 ; j<NumEntries ; ++j) {
      TEST_FLOATING_EQUALITY( Values[j], Values_exact[j], tol );
    }
  }
}
#endif

# if 0 
TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, DiodeCapLoad_diff_times ) {
  int iargs = 2;
  char *cargs[iargs];
  cargs[0] = "Xyce";
  cargs[1] = "forward_diode_model.cir";
  // Use simple Diode-Capacitor circuit netlist as input
  // Set up t,x,xdot,s,sdot,f,W
  RCP<N_ANP_ModelEvaluator> model = rcp(new N_ANP_ModelEvaluator());
  model->initialize(iargs,cargs);
  RCP<const Epetra_Map> x_map = model->get_x_map();
  RCP<const Epetra_Map> f_map = model->get_f_map();
  RCP<const Epetra_Map> s_map = model->get_p_map(0);
  RCP<const Epetra_Map> sdot_map = model->get_p_map(1);
  RCP<Epetra_Operator> W = model->create_W();
  RCP<Epetra_Operator> W_exact = model->create_W();
  double alpha,beta;
  double t;
  RCP<Epetra_Vector> xbar = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> xbardot = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> fbar = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> s = rcp(new Epetra_Vector(*s_map)); 
  RCP<Epetra_Vector> sdot = rcp(new Epetra_Vector(*sdot_map)); 
  RCP<Epetra_Vector> voltLimF = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> voltLimQ = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> fbar_exact = rcp(new Epetra_Vector(*f_map));
  RCP<Epetra_Vector> fq_exact = rcp(new Epetra_Vector(*f_map));
  {
    Epetra_Vector& xb = *xbar;
    Epetra_Vector& xbdot = *xbardot;
    Epetra_Vector& sb = *s;
    Epetra_Vector& sbdot = *sdot;
    Epetra_Vector& fb_exact = *fbar_exact;
    Epetra_Vector& fq_e = *fq_exact;
    t = 0.0001; // 0.00025;
    //two blocks of three elements each.
    xb[0] = 0.2; // v_2
    xb[1] = 0.3; // i_vsrc
    xb[2] = 0.4; // v_1
    double v2 = xb[0]; // 0.2
    double i_vsrc = xb[1]; // 0.3
    double v1 = xb[2]; // 0.4

    xb[3] = 0.5; // z
    xb[4] = 0.6; 
    xb[5] = 0.7; 

    xbdot[0] =0.8; // v_2 dot
    xbdot[1] =0.9; // i_vsrc dot
    xbdot[2] =0.10; // v_1 dot

    xbdot[3] =0.11; // z dot
    xbdot[4] =0.12; 
    xbdot[5] =0.13; 

    sb[0] = 0.14; // s

    sbdot[0] = 0.17; // sdot

    alpha = 20;
    beta = 30;
		
		// Capacitor Parms
    double C = 1.0e-12;
		// Diode Parms: Forward model of the diode 
    double IS = 1.0e-14; // Saturation current
    double N  = 1; // Mult Factor
    double Vth  = 0.025864186; //
    double sinfreq = 1000;
    double pi = 4.0*atan(1.0);
    //double vsrc = 0.2*sin(2*pi*sinfreq*t);
    double vsrc = 0.2; 
		double R = 1e6;
    fq_e[0] = -IS*(exp((v1 - v2)/(N*Vth))-1);
    fq_e[1] =  v1-vsrc;
    fq_e[2] =  IS*(exp((v1 - v2)/(N*Vth))-1)+i_vsrc;
    // q:
    fq_e[3] = C*v2;
    fq_e[4] = 0.0;
    fq_e[5] = 0.0;

    fb_exact[0] = xbdot[3] + fq_e[0];
    fb_exact[1] = xbdot[4] + fq_e[1];
    fb_exact[2] = xbdot[5] + fq_e[2];
    fb_exact[3] = fq_e[3] - xb[3];
    fb_exact[4] = fq_e[4] - xb[4];
    fb_exact[5] = fq_e[5] - xb[5];

    // df/dx:
    // df0/dv2, df0/di, df0/dv1 = [ IS * exp ((V1 -V2)/N*Vth)*1/N*Vth,   0,    -IS*exp(v1-v2/N*Vth)]
    // df1/dv2, df1/di, df1/dv1 = [ 0,                                   0,                    1   ]
    // df2/dv2, df2/di, df2/dv1 = [ IS * exp ((V1 -V2)/N*Vth)*1/N*Vth,   1,    -IS*exp(v1-v2/N*Vth)]
    // dq/dx:
    // dq0/dv2, dq0/di, dq0/dv1 = [ C, 0, 0 ]
    // dq1/dv2, dq1/di, dq1/dv1 = [ 0, 0, 0 ]
    // dq2/dv2, dq2/di, dq2/dv1 = [ 0, 0, 0 ]

    // W = [ beta*df/dx , alpha*I ]
    //     [ beta*dq/dx , -beta*I  ]
    std::string label = "N_LAS_BlockMatrix";
    RCP<N_LAS_BlockMatrix> bW_exact = Teuchos::get_extra_data<RCP<N_LAS_BlockMatrix> >(W_exact,label);

    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "bW_exact->epetraObj().Graph() = " << std::endl;
    //bW_exact->epetraObj().Graph().Print(std::cout);
    //std::cout << "-----------------------------------------------" << std::endl;

    bW_exact->put(0.0);
    {
      N_LAS_Matrix & ul = bW_exact->block(0,0);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ul.epetraObj().Graph() = " << std::endl;
      //ul.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
     
      /* Added by Prateek while playing with the code	*/
      //N_LAS_Matrix & u2 = bW_exact->block(1,0);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "u2.epetraObj().Graph() = " << std::endl;
      //u2.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;

      int row;
      int length;
      double coeffs[3];
      int colIndices[3];
      row = 0;
      length=2;
      coeffs[0] = beta*IS*exp((v1 - v2)/N/Vth)*1/N/Vth;
      coeffs[1] = -beta*IS*exp((v1 - v2)/N/Vth)*1/N/Vth;
      colIndices[0] = 0;
      colIndices[1] = 2;
      ul.putRow(row,length,coeffs,colIndices);
      row = 1;
      length = 1;
      coeffs[0] = beta;
      colIndices[0] = 2;
      ul.putRow(row,length,coeffs,colIndices);
      row = 2;
      length=3;
      coeffs[0] = -beta*IS*exp((v1 - v2)/N/Vth)*1/N/Vth;
      coeffs[1] = beta;
      coeffs[2] = beta*IS*exp((v1 - v2)/N/Vth)*1/N/Vth;
      colIndices[0] = 0;
      colIndices[1] = 1;
      colIndices[2] = 2;
      ul.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ul = " << std::endl;
      //ul.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    {
      N_LAS_Matrix & ur = bW_exact->block(0,1);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ur.epetraObj().Graph() = " << std::endl;
      //ur.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
      int row;
      int length = 1;
      double coeffs[1];
      int colIndices[1];
      row = 0;
      coeffs[0] = alpha;
      colIndices[0] = 0;
      ur.putRow(row,length,coeffs,colIndices);
      row = 1;
      colIndices[0] = 1;
      ur.putRow(row,length,coeffs,colIndices);
      row = 2;
      colIndices[0] = 2;
      ur.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ur = " << std::endl;
      //ur.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    {
      N_LAS_Matrix & ll = bW_exact->block(1,0);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ll.epetraObj().Graph() = " << std::endl;
      //ll.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
      int row;
      int length;
      double coeffs[3];
      int colIndices[3];
      row = 0;
      length = 1;
      coeffs[0] = beta*C;
      colIndices[0] = 0;
      ll.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "ll = " << std::endl;
      //ll.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    {
      N_LAS_Matrix & lr = bW_exact->block(1,1);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "lr.epetraObj().Graph() = " << std::endl;
      //lr.epetraObj().Graph().Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
      int row;
      int length = 1;
      double coeffs[1];
      int colIndices[1];
      row = 0;
      coeffs[0] = -beta;
      colIndices[0] = 0;
      lr.putRow(row,length,coeffs,colIndices);
      row = 1;
      colIndices[0] = 1;
      lr.putRow(row,length,coeffs,colIndices);
      row = 2;
      colIndices[0] = 2;
      lr.putRow(row,length,coeffs,colIndices);
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "lr = " << std::endl;
      //lr.printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    bW_exact->fillComplete();
    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "bW_exact->epetraObj().Graph() = " << std::endl;
    //bW_exact->epetraObj().Graph().Print(std::cout);
    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "bW_exact = " << std::endl;
    //bW_exact->printPetraObject();
    //std::cout << "-----------------------------------------------" << std::endl;

  }



  /* Why does Todd do two evamodels???? */
  // First we call the model to load the state vector s.
  {
    EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
    inArgs.set_t(t);
    inArgs.set_x(xbar);
    inArgs.set_x_dot(xbardot);
    //inArgs.set_p(0,s);
    //inArgs.set_p(1,sdot);
    EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
    
		/* Question: what is g?? */
		outArgs.set_g(0,s);
    //fbar->PutScalar(0.0);
    //outArgs.set_f(fbar);

    // evalModel
    model->evalModel(inArgs,outArgs);
  }
  // Verify model loaded s correctly:
  // TODO
  // Second we "calculate" sdot
  // Third we call the model to load f
  {
    EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
    inArgs.set_t(t);
    inArgs.set_x(xbar);
    inArgs.set_x_dot(xbardot);
    inArgs.set_p(0,s);
    inArgs.set_p(1,sdot);
    inArgs.set_alpha(alpha);
    inArgs.set_beta(beta);
    EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
    fbar->PutScalar(0.0);
    outArgs.set_f(fbar);
    outArgs.set_W(W);
    outArgs.set_g(1,voltLimQ);
    outArgs.set_g(2,voltLimF);

    // evalModel
    model->evalModel(inArgs,outArgs);
  }

  std::cout << "xbar = " << *xbar << std::endl;
  std::cout << "xbardot = " << *xbardot << std::endl;
  std::cout << "s = " << *s << std::endl;
  std::cout << "sdot = " << *sdot << std::endl;
  std::cout << "fbar = " << *fbar << std::endl;
  std::cout << "fq_exact = " << *fq_exact << std::endl;
  std::cout << "fbar_exact = " << *fbar_exact << std::endl;

  // Verify model loaded F correctly
  double tol = 1.0e-6;
  TEST_FLOATING_EQUALITY( (*fbar)[0], (*fbar_exact)[0], tol );
  TEST_FLOATING_EQUALITY( (*fbar)[1], (*fbar_exact)[1], tol );
  TEST_FLOATING_EQUALITY( (*fbar)[2], (*fbar_exact)[2], tol );

  TEST_FLOATING_EQUALITY( (*fbar)[3], (*fbar_exact)[3], tol );
  TEST_FLOATING_EQUALITY( (*fbar)[4], (*fbar_exact)[4], tol );
  TEST_FLOATING_EQUALITY( (*fbar)[5], (*fbar_exact)[5], tol );

  // Verify model loaded W correctly
  RCP<Epetra_CrsMatrix> W_matrix_exact = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_exact,true);
  //std::cout << "W_matrix_exact = " << std::endl;
  //W_matrix_exact->Print(std::cout);
  //std::cout << "-----------------------------------------------" << std::endl;
  RCP<Epetra_CrsMatrix> W_matrix = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W,true);
  //std::cout << "W_matrix = " << std::endl;
  //W_matrix->Print(std::cout);
  //std::cout << "-----------------------------------------------" << std::endl;

  for (int i=0 ; i<6 ; ++i) {
    int GlobalRow = i;
    if (i>2) {
      GlobalRow = (i-3)+10;
    }
    int Length = 6;
    int NumEntries;
    std::vector<double> Values; Values.resize(Length);
    std::vector<int> Indices ; Indices.resize(Length);
    W_matrix->ExtractGlobalRowCopy(GlobalRow, Length, NumEntries, &Values[0], &Indices[0]);

    int NumEntries_exact;
    std::vector<double> Values_exact; Values_exact.resize(Length);
    std::vector<int> Indices_exact ; Indices_exact.resize(Length);
    W_matrix_exact->ExtractGlobalRowCopy(GlobalRow, Length, NumEntries_exact, &Values_exact[0], &Indices_exact[0]);

    TEST_EQUALITY( NumEntries, NumEntries_exact );
    for (int j=0 ; j<NumEntries ; ++j) {
      TEST_FLOATING_EQUALITY( Values[j], Values_exact[j], tol );
    }
  }
}
#endif

TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, RC_W_graph ) {
  // Graph for RC problem's W matrix should be as follows:
  // [ x, 0, x ] , [ x, 0, x ]
  // [ 0, x, x ] , [ 0, x, x ]
  // [ x, x, x ] , [ x, x, x ]
  // -------------------------
  // [ x, 0, x ] , [ x, 0, x ]
  // [ 0, x, x ] , [ 0, x, x ]
  // [ x, x, x ] , [ x, x, x ]
  int iargs = 2;
  char *cargs[iargs];
  cargs[0] = "Xyce";
  cargs[1] = "modelEvaluatorTest.cir";
  int numEntries = 6;
  RCP<N_ANP_ModelEvaluator> model = rcp(new N_ANP_ModelEvaluator());
  model->initialize(iargs,cargs);
  RCP<const Epetra_CrsGraph> graph;
  {
    RCP<Epetra_Operator> W = model->create_W();
    std::string label = "N_LAS_BlockMatrix";
    RCP<N_LAS_BlockMatrix> bA = Teuchos::get_extra_data<RCP<N_LAS_BlockMatrix> >(W,label);
    graph = Teuchos::rcpFromRef(bA->epetraObj().Graph());
  }
  int lenOfIndices = numEntries;
  int numIndices;
  std::vector<int> indices;  indices.resize(lenOfIndices);
  for (int globalRow=0 ; globalRow<=12 ; globalRow+=10) {
    graph->ExtractGlobalRowCopy(globalRow,lenOfIndices,numIndices,&indices[0]);
    TEST_EQUALITY_CONST( numIndices, 4 );
    TEST_EQUALITY_CONST( indices[0], 0 );
    TEST_EQUALITY_CONST( indices[1], 2 );
    TEST_EQUALITY_CONST( indices[2], 10 );
    TEST_EQUALITY_CONST( indices[3], 12 );
  }
  for (int globalRow=1 ; globalRow<=12 ; globalRow+=10) {
    graph->ExtractGlobalRowCopy(globalRow,lenOfIndices,numIndices,&indices[0]);
    TEST_EQUALITY_CONST( numIndices, 4 );
    TEST_EQUALITY_CONST( indices[0], 1 );
    TEST_EQUALITY_CONST( indices[1], 2 );
    TEST_EQUALITY_CONST( indices[2], 11 );
    TEST_EQUALITY_CONST( indices[3], 12 );
  }
  for (int globalRow=2 ; globalRow<=12 ; globalRow+=10) {
    graph->ExtractGlobalRowCopy(globalRow,lenOfIndices,numIndices,&indices[0]);
    TEST_EQUALITY_CONST( numIndices, 6 );
    TEST_EQUALITY_CONST( indices[0], 0 );
    TEST_EQUALITY_CONST( indices[1], 1 );
    TEST_EQUALITY_CONST( indices[2], 2 );
    TEST_EQUALITY_CONST( indices[3], 10 );
    TEST_EQUALITY_CONST( indices[4], 11 );
    TEST_EQUALITY_CONST( indices[5], 12 );
  }
}

TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, linearSolve ) {
  int iargs = 2;
  char *cargs[iargs];
  cargs[0] = "Xyce";
  cargs[1] = "modelEvaluatorTest.cir";
  int numEntries = 6;
  RCP<N_ANP_ModelEvaluator> model = rcp(new N_ANP_ModelEvaluator());
  model->initialize(iargs,cargs);
  RCP<const Epetra_Map> x_map = model->get_x_map();
  RCP<N_LAS_BlockMatrix> bA;
  {
    RCP<Epetra_Operator> W = model->create_W();
    std::string label = "N_LAS_BlockMatrix";
    bA = Teuchos::get_extra_data<RCP<N_LAS_BlockMatrix> >(W,label);
  }
  Epetra_CrsMatrix& A = bA->epetraObj();
  RCP<Epetra_Vector> x = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> b = rcp(new Epetra_Vector(*x_map)); 
  {
    // Fill bA with identity
    bA->put(0.0);
    int row;
    int length = 1;
    double coeffs[1];
    coeffs[0] = 1.0;
    int colIndices[1];
    for (row=0 ; row<numEntries ; ++row) {
      int globalRow = bA->epetraObj().GRID(row);
      colIndices[0] = globalRow;
      bA->putRow(globalRow,length,coeffs,colIndices);
    }
    bA->fillComplete();
    //bA->printPetraObject();
  }
  { 
    // Fill b with 1..numEntries
    for (int row=0 ; row<numEntries; ++row) {
      (*b)[row] = row;
    }
  }
  {
    // Solve linear system
    Epetra_LinearProblem problem(&A,&*x,&*b);
    Amesos_Klu solver(problem);
    solver.SetUseTranspose(false);
    solver.SymbolicFactorization(); // Repeat if nonzeros in A change
    solver.NumericFactorization(); // Repeat if values in A change
    solver.Solve(); // Repeat if values in A or b change
  }
  // Verify x == b
  for (int i=0 ; i<numEntries ; ++i) {
    TEST_EQUALITY( (*x)[i], (*b)[i] );
  }
}

#if 1
// Here we perform a nonlinear solve on the circuits. The flow of nonlinear
// solve is as below:
// 		Step 1: initialize the circuit with x0, x_i = x0
// 		Step 2: Call ModelEvaluator to compute the W and f.
// 		Step 3: Perform a linear system solve by calling some(?) Trilinos package
// 		to obtain Deltax 
// 		Step 4: update x, x_{i+1} = x_i + Deltax

TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, RCLoad_Newton_Raphson ) {
  int iargs = 2;
  char *cargs[iargs];
  cargs[0] = "Xyce";
  cargs[1] = "modelEvaluatorTest.cir";
  // Use simple RC circuit netlist as input
  // Set up t,x,xdot,s,sdot,f,W
  RCP<N_ANP_ModelEvaluator> model = rcp(new N_ANP_ModelEvaluator());
  model->initialize(iargs,cargs);
  RCP<const Epetra_Map> x_map = model->get_x_map();
  RCP<const Epetra_Map> f_map = model->get_f_map();
  RCP<const Epetra_Map> s_map = model->get_p_map(0);
  RCP<const Epetra_Map> sdot_map = model->get_p_map(1);
  RCP<const Epetra_Map> small_x_map = model->get_small_x_map();
  RCP<Epetra_Operator> W = model->create_W();
  
  double alpha,beta;
  double t;
  RCP<Epetra_Vector> xbar = rcp(new Epetra_Vector(*x_map)); 
  RCP<N_LAS_BlockVector> nlas_block_xbar = convertEpetraToNLASBlockVectorView(xbar,*small_x_map);
  RCP<Epetra_Vector> xbar_block_0 = Teuchos::rcpFromRef(*(nlas_block_xbar->block(0).epetraVector())); // This a view!
  RCP<Epetra_Vector> xbardot = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> fbar = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> s = rcp(new Epetra_Vector(*s_map)); 
  RCP<Epetra_Vector> sdot = rcp(new Epetra_Vector(*sdot_map)); 
  RCP<Epetra_Vector> voltLimF = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> voltLimQ = rcp(new Epetra_Vector(*f_map)); 
  
  Epetra_Vector& xb = *xbar;
  Epetra_Vector& xbdot = *xbardot;
  Epetra_Vector& sb = *s;
  Epetra_Vector& sbdot = *sdot;
  t = 0.0001; // 0.00025;
  //two blocks of three elements each.
  xb[0] = 2.0; // v_2
  xb[1] = 3.0; // i_vsrc
  xb[2] = 4.0; // v_1
  double v2 = xb[0]; // 2.0
  double i_vsrc = xb[1]; // 3.0
  double v1 = xb[2]; // 4.0

  xb[3] = 5.0; // z
  xb[4] = 6.0; 
  xb[5] = 7.0; 

  xbdot[0] = 0.0; // v_2 dot
  xbdot[1] = 0.0; // i_vsrc dot
  xbdot[2] = 0.0; // v_1 dot

  xbdot[3] = 0.0; // z dot
  xbdot[4] = 0.0; 
  xbdot[5] = 0.0; 

  sb[0] = 14.0; // s

  sbdot[0] = 0.0; // sdot

  alpha = 1;
  beta = 1;

  double C = 1.0e-7;
  double R = 1.0e6;
  double sinfreq = 1000;
  double pi = 4.0*atan(1.0);
  double vsrc = 12.0*sin(2*pi*sinfreq*t);

  //std::cout << "-----------------------------------------------" << std::endl;
  //std::cout << "vsrc = " << vsrc << std::endl;
  //std::cout << "-----------------------------------------------" << std::endl;
  // W = [ beta*df/dx , alpha*I ]
  //     [ beta*dq/dx , -beta*I  ]

  	
	// Now we will start the newton loop to perform the nonlinear solve
	int converged = 0;
	double newtonTol = 1e-6;
	int iter = 0;
	while (!converged) {
		
  	// First we call the model to load the state vector s.
  	{
  	  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  	  inArgs.set_t(t);
  	  inArgs.set_x(xbar);
  	  inArgs.set_x_dot(xbardot);
  	  //inArgs.set_p(0,s);
  	  //inArgs.set_p(1,sdot);
  	  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
  	  
			/* Question: what is g?? */
			outArgs.set_g(0,s);
  	  //fbar->PutScalar(0.0);
  	  //outArgs.set_f(fbar);

  	  // evalModel
  	  model->evalModel(inArgs,outArgs);
  	}
  	// Verify model loaded s correctly:
  	// TODO
  	// Second we "calculate" sdot
  	// Third we call the model to load f
  	{
  	  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  	  inArgs.set_t(t);
  	  inArgs.set_x(xbar);
  	  inArgs.set_x_dot(xbardot);
  	  inArgs.set_p(0,s);
  	  inArgs.set_p(1,sdot);
  	  inArgs.set_alpha(alpha);
  	  inArgs.set_beta(beta);
  	  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
  	  fbar->PutScalar(0.0);
  	  outArgs.set_f(fbar);
  	  outArgs.set_W(W);
  	  outArgs.set_g(1,voltLimQ);
  	  outArgs.set_g(2,voltLimF);
  	  // evalModel
  	  model->evalModel(inArgs,outArgs);
				
  	}

		
  	//std::cout << "xbar = " << *xbar << std::endl;
  	//std::cout << "xbardot = " << *xbardot << std::endl;
  	//std::cout << "s = " << *s << std::endl;
  	//std::cout << "sdot = " << *sdot << std::endl;
  	//std::cout << "fbar = " << *fbar << std::endl;
		
		// Now pull out the block vector from F(xbar, xbar_dot, t)  and block
		// matrix W which corresponds to f(x) and df_dx 
			
		// Pull f(x) 
		RCP<N_LAS_BlockVector> xyce_f = convertEpetraToNLASBlockVectorView(fbar,*small_x_map); 
		N_LAS_Vector & f = xyce_f->block(0);	
  	RCP<Epetra_Vector> ef = Teuchos::rcpFromRef(*f.epetraVector());

		//std::cout << "-----------------------------------------------" << std::endl;
  	//std::cout << "small f ="<< std::endl;
		//f.printPetraObject(); 
		
		double fNorm = *f.lpNorm(1); 
  	//std::cout << "norm(f) =" << fNorm  << std::endl;
  	//std::cout << "-----------------------------------------------" << std::endl;

	
		//std::cout << "-----------------------------------------------" << std::endl;
  	//std::cout << "small f converterd to Epetra vector ="<< std::endl;
		//ef->Print(std::cout);	
  	//std::cout << "-----------------------------------------------" << std::endl;
		
		if ( fNorm <= newtonTol) {
			converged = 1;
			//std::cout << "-----------------------------------------------" << std::endl;
			//std::cout << "***********************************************" << std::endl;
  		//std::cout << " Newton Raphson for RC load Converged in "<< iter << " iterations." <<  std::endl;
  		//std::cout << "Solution: xbar = " << *xbar << std::endl;
			//std::cout << "***********************************************" << std::endl;
  		//std::cout << "-----------------------------------------------" << std::endl;
		}		


  	RCP<Epetra_CrsMatrix> W_matrix = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W,true);
  	//std::cout << "W_matrix = " << std::endl;
  	//W_matrix->Print(std::cout);
  	//std::cout << "-----------------------------------------------" << std::endl;
  	
		std::string label = "N_LAS_BlockMatrix";
  	RCP<N_LAS_BlockMatrix> bW = Teuchos::get_extra_data<RCP<N_LAS_BlockMatrix> >(W,label);

  	N_LAS_Matrix & ul = bW->block(0,0);
  	//std::cout << "-----------------------------------------------" << std::endl;
  	//std::cout << "ul.epetraObj()=" << std::endl;
  	//ul.printPetraObject();
  	////ul.epetraObj().Graph().Print(std::cout);
  	//std::cout << "-----------------------------------------------" << std::endl;

		//Linear solve begins
  	RCP<Epetra_Vector> Deltax = rcp(new Epetra_Vector(*small_x_map)); 
  	Epetra_CrsMatrix& A = ul.epetraObj();
		Epetra_LinearProblem problem(&A,&*Deltax,&*ef);
    Amesos_Klu solver(problem);
    solver.SetUseTranspose(false);
    solver.SymbolicFactorization(); // Repeat if nonzeros in A change
    solver.NumericFactorization(); // Repeat if values in A change
    solver.Solve(); // Repeat if values in A or b change

		//std::cout << "-----------------------------------------------" << std::endl;
  	//std::cout << "Delta x="<< std::endl;
		//Deltax->Print(std::cout);	
  	//std::cout << "-----------------------------------------------" << std::endl;

		// update xbar  [ xbar = 1.0*xbar -1.0*Deltax ]
    xbar_block_0->Update(-1.0,*Deltax,1.0);
		iter++;
    // Example of Update with two vectors
    // xbar = 1.0*xbar + alpha*vecA + beta*vecB
    // xbar->update(alpha,*vecA,beta,*vecB,1.0);
    // Example of scaling a vector:
    // xbar = alpha*xbar
    // xbar->Scale(alpha);
	}
	
  // Verify that we computed the solution correctly
	// In steady state the the i_vsrc = 0, v1 = vsrc, v2 = vsrc  
  double tol = 1.0e-6;
  TEST_FLOATING_EQUALITY( (*xbar)[0], vsrc, tol ); // v2
  TEST_FLOATING_EQUALITY( (*xbar)[1], 0.0, tol );  //i_vsrc
  TEST_FLOATING_EQUALITY( (*xbar)[2], vsrc, tol ); //v1
 
}
#endif

#if 1 
// Here we perform a nonlinear solve on the circuits. The flow of nonlinear
// solve is as below:
// 		Step 1: initialize the circuit with x0, x_i = x0
// 		Step 2: Call ModelEvaluator to compute the W and f.
// 		Step 3: Perform a linear system solve by calling some(?) Trilinos package
// 		to obtain Deltax 
// 		Step 4: update x, x_{i+1} = x_i + Deltax
// Nonlinear solve of Diode-Capacitor circuit. 
TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, DiodeCapLoad_Newton_Raphson ) {
  int iargs = 2;
  char *cargs[iargs];
  cargs[0] = "Xyce";
  cargs[1] = "forward_diode_model.cir";
  // Use simple Diode-Capacitor netlist as input
  // Set up t,x,xdot,s,sdot,f,W
  RCP<N_ANP_ModelEvaluator> model = rcp(new N_ANP_ModelEvaluator());
  model->initialize(iargs,cargs);
  RCP<const Epetra_Map> x_map = model->get_x_map();
  RCP<const Epetra_Map> f_map = model->get_f_map();
  RCP<const Epetra_Map> s_map = model->get_p_map(0);
  RCP<const Epetra_Map> sdot_map = model->get_p_map(1);
  RCP<const Epetra_Map> small_x_map = model->get_small_x_map();
  RCP<Epetra_Operator> W = model->create_W();
  
  double alpha,beta;
  double t;
  RCP<Epetra_Vector> xbar = rcp(new Epetra_Vector(*x_map)); 
  RCP<N_LAS_BlockVector> nlas_block_xbar = convertEpetraToNLASBlockVectorView(xbar,*small_x_map);
  RCP<Epetra_Vector> xbar_block_0 = Teuchos::rcpFromRef(*(nlas_block_xbar->block(0).epetraVector())); // This a view!
  RCP<Epetra_Vector> xbardot = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> fbar = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> s = rcp(new Epetra_Vector(*s_map)); 
  RCP<Epetra_Vector> sdot = rcp(new Epetra_Vector(*sdot_map)); 
  RCP<Epetra_Vector> voltLimF = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> voltLimQ = rcp(new Epetra_Vector(*f_map)); 
  
  Epetra_Vector& xb = *xbar;
  Epetra_Vector& xbdot = *xbardot;
  Epetra_Vector& sb = *s;
  Epetra_Vector& sbdot = *sdot;
  t = 0.0001; // 0.00025;
  //two blocks of three elements each.

  //two blocks of three elements each.
  xb[0] = 0.2; // v_2
  xb[1] = 0.3; // i_vsrc
  xb[2] = 0.4; // v_1
  double v2 = xb[0]; // 0.2
  double i_vsrc = xb[1]; // 0.3
  double v1 = xb[2]; // 0.4

  xb[3] = 0.5; // z
  xb[4] = 0.6; 
  xb[5] = 0.7; 

  xbdot[0] =0; // v_2 dot
  xbdot[1] =0; // i_vsrc dot
  xbdot[2] =0; // v_1 dot

  xbdot[3] =0; // z dot
  xbdot[4] =0; 
  xbdot[5] =0; 

  sb[0] = 0.14; // s

  sbdot[0] = 0.17; // sdot

  alpha = 1;
  beta = 1;

	// Capacitor Parms
  double C = 1.0e-12;
	// Diode Parms: Forward model of the diode 
  double IS = 1.0e-14; // Saturation current
  double N  = 1; // Mult Factor
  double Vth  = 0.025864186; // Threshold Factor
  double sinfreq = 1000;
  double pi = 4.0*atan(1.0);
  //double vsrc = 12.0*sin(2*pi*sinfreq*t);
  double vsrc = 0.2; 
	double R = 1e6;

	// Now we will start the newton loop to perform the nonlinear solve
	int converged = 0;
	double newtonTol = 1e-18;
	int iter = 0;
	while (!converged) {
		
  	// First we call the model to load the state vector s.
  	{
  	  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  	  inArgs.set_t(t);
  	  inArgs.set_x(xbar);
  	  inArgs.set_x_dot(xbardot);
  	  //inArgs.set_p(0,s);
  	  //inArgs.set_p(1,sdot);
  	  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
  	  
			/* Question: what is g?? */
			outArgs.set_g(0,s);
  	  //fbar->PutScalar(0.0);
  	  //outArgs.set_f(fbar);

  	  // evalModel
  	  model->evalModel(inArgs,outArgs);
  	}
  	// Verify model loaded s correctly:
  	// TODO
  	// Second we "calculate" sdot
  	// Third we call the model to load f
  	{
  	  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  	  inArgs.set_t(t);
  	  inArgs.set_x(xbar);
  	  inArgs.set_x_dot(xbardot);
  	  inArgs.set_p(0,s);
  	  inArgs.set_p(1,sdot);
  	  inArgs.set_alpha(alpha);
  	  inArgs.set_beta(beta);
  	  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
  	  fbar->PutScalar(0.0);
  	  outArgs.set_f(fbar);
  	  outArgs.set_W(W);
  	  outArgs.set_g(1,voltLimQ);
  	  outArgs.set_g(2,voltLimF);
  	  // evalModel
  	  model->evalModel(inArgs,outArgs);
				
  	}

		
  	//std::cout << "xbar = " << *xbar << std::endl;
  	//std::cout << "xbardot = " << *xbardot << std::endl;
  	//std::cout << "s = " << *s << std::endl;
  	//std::cout << "sdot = " << *sdot << std::endl;
  	//std::cout << "fbar = " << *fbar << std::endl;
		
		// Now pull out the block vector from F(xbar, xbar_dot, t)  and block
		// matrix W which corresponds to f(x) and df_dx 
			
		// Pull f(x) 
		RCP<N_LAS_BlockVector> xyce_f = convertEpetraToNLASBlockVectorView(fbar,*small_x_map); 
		N_LAS_Vector & f = xyce_f->block(0);	
  	RCP<Epetra_Vector> ef = Teuchos::rcpFromRef(*f.epetraVector());

		//std::cout << "-----------------------------------------------" << std::endl;
  	//std::cout << "small f ="<< std::endl;
		//f.printPetraObject(); 
		
		double fNorm = *f.lpNorm(1); 
  	//std::cout << "norm(f) =" << fNorm  << std::endl;
  	//std::cout << "-----------------------------------------------" << std::endl;

	
		//std::cout << "-----------------------------------------------" << std::endl;
  	//std::cout << "small f converterd to Epetra vector ="<< std::endl;
		//ef->Print(std::cout);	
  	//std::cout << "-----------------------------------------------" << std::endl;
		
		if ( fNorm <= newtonTol) {
			converged = 1;
			//std::cout << "-----------------------------------------------" << std::endl;
			//std::cout << "***********************************************" << std::endl;
			//std::cout << " Newton Raphson for Diode Capacitor load Converged in "<<
			//iter << " iterations." <<  std::endl;
  		//std::cout << "Solution: xbar = " << *xbar << std::endl;
			//std::cout << "***********************************************" << std::endl;
  		//std::cout << "-----------------------------------------------" << std::endl;
		}		


  	RCP<Epetra_CrsMatrix> W_matrix = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W,true);
  	//std::cout << "W_matrix = " << std::endl;
  	//W_matrix->Print(std::cout);
  	//std::cout << "-----------------------------------------------" << std::endl;
  	
		std::string label = "N_LAS_BlockMatrix";
  	RCP<N_LAS_BlockMatrix> bW = Teuchos::get_extra_data<RCP<N_LAS_BlockMatrix> >(W,label);

  	N_LAS_Matrix & ul = bW->block(0,0);
  	//std::cout << "-----------------------------------------------" << std::endl;
  	//std::cout << "ul.epetraObj()=" << std::endl;
  	//ul.printPetraObject();
  	////ul.epetraObj().Graph().Print(std::cout);
  	//std::cout << "-----------------------------------------------" << std::endl;

		//Linear solve begins
  	RCP<Epetra_Vector> Deltax = rcp(new Epetra_Vector(*small_x_map)); 
  	Epetra_CrsMatrix& A = ul.epetraObj();
		Epetra_LinearProblem problem(&A,&*Deltax,&*ef);
    Amesos_Klu solver(problem);
    solver.SetUseTranspose(false);
    solver.SymbolicFactorization(); // Repeat if nonzeros in A change
    solver.NumericFactorization(); // Repeat if values in A change
    solver.Solve(); // Repeat if values in A or b change

		//std::cout << "-----------------------------------------------" << std::endl;
  	//std::cout << "Delta x="<< std::endl;
		//Deltax->Print(std::cout);	
  	//std::cout << "-----------------------------------------------" << std::endl;

		// update xbar [xbar = 1.0*xbar -1.0*Deltax]
		xbar_block_0->Update(-1.0,*Deltax,1.0);
		iter++;
	}
	
  // Verify that we computed the solution correctly
	// In steady state the the i_vsrc = 0, v1 = vsrc, v2 = vsrc  
  double tol = 1.0e-6;
  TEST_FLOATING_EQUALITY( (*xbar)[0], vsrc, tol ); // v2
  TEST_FLOATING_EQUALITY( (*xbar)[1], 0.0, tol );  //i_vsrc
  TEST_FLOATING_EQUALITY( (*xbar)[2], vsrc, tol ); //v1
 
}
#endif

#if 1
// Now we perform a transient analysis of simple circuits. As a first
// step we will use the Backward-Euler (BE) integration routine to
// perform the transient analysis. The standard form of the BE is as
// below: 
// 1) Start with the standard circuit DAE form:
//    q_dot(x) + f(x) + b(t) = 0; where except time everything is a
//    vector.
// 2) Discretize the DAE:
//    [q_dot(x_new)  - q_dot(x_old)]/Delta_t + f(x_new) + b(t_new) = 0; 
// 3) solve for the x_new using the standard Newton-Raphson Algorithm   

TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, RCLoad_Backward_Euler) {
  int iargs = 2;
  char *cargs[iargs];
  cargs[0] = "Xyce";
  cargs[1] = "modelEvaluatorTest.cir";

	// vars related with analytical solution of the BE for RC circuit
	int getInitICFlag =  0;
	double V2_n; //V2  at nth step 
	double Isrc_n; // I_src at nth step
	double V1_n; // V1 at nth step
	double omega ; // 1/Det(BE_LHS)
	double V2_next; 
	double Isrc_next; 

  // Use simple RC circuit netlist as input
  // Set up t,x,xdot,s,sdot,f,W
  RCP<N_ANP_ModelEvaluator> model = rcp(new N_ANP_ModelEvaluator());
  model->initialize(iargs,cargs);
  RCP<const Epetra_Map> x_map = model->get_x_map();
  RCP<const Epetra_Map> f_map = model->get_f_map();
  RCP<const Epetra_Map> s_map = model->get_p_map(0);
  RCP<const Epetra_Map> sdot_map = model->get_p_map(1);
  RCP<const Epetra_Map> small_x_map = model->get_small_x_map();
  RCP<Epetra_Operator> W = model->create_W();
  
  double alpha,beta;
  double t;

  RCP<Epetra_Vector> xbar = rcp(new Epetra_Vector(*x_map)); 
  RCP<N_LAS_BlockVector> nlas_block_xbar = convertEpetraToNLASBlockVectorView(xbar,*small_x_map);
  RCP<Epetra_Vector> xbar_block_0 = Teuchos::rcpFromRef(*(nlas_block_xbar->block(0).epetraVector())); // This a view!
  RCP<Epetra_Vector> xbardot = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> fbar = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> eq_old = rcp(new Epetra_Vector(*small_x_map));
  RCP<Epetra_Vector> s = rcp(new Epetra_Vector(*s_map)); 
  RCP<Epetra_Vector> sdot = rcp(new Epetra_Vector(*sdot_map)); 
  RCP<Epetra_Vector> voltLimF = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> voltLimQ = rcp(new Epetra_Vector(*f_map)); 
 	
	// This will be used to store the values at old time points
  RCP<Epetra_Vector> xbar_old = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> fbar_old = rcp(new Epetra_Vector(*f_map)); 
	
	// Make a temporary vector to compute the rhs of BE 
  RCP<Epetra_Vector> e_rhs_BE = rcp(new Epetra_Vector(*small_x_map));

  Epetra_Vector& xb = *xbar;
  Epetra_Vector& xbdot = *xbardot;
  Epetra_Vector& sb = *s;
  Epetra_Vector& sbdot = *sdot;

  ////two blocks of three elements each.
  xb[0] = 4.0; // v_2
  xb[1] = 4.0e-6; // i_vsrc
  xb[2] = 8.0; // v_1
  //two blocks of three elements each.
  //xb[0] = 0.0; // v_2
  //xb[1] = 0.0e-6; // i_vsrc
  //xb[2] = 0.0; // v_1
  double v2 = xb[0]; // 2.0
  double i_vsrc = xb[1]; // 3.0
  double v1 = xb[2]; // 4.0
	// Note: These are required to be set to zero, becuase we want to
	// extract the f_func and q_func from F.
  xb[3] = 0.0; // z
  xb[4] = 0.0; 
  xb[5] = 0.0; 

  xbdot[0] = 0.0; // v_2 dot
  xbdot[1] = 0.0; // i_vsrc dot
  xbdot[2] = 0.0; // v_1 dot

  xbdot[3] = 0.0; // z dot
  xbdot[4] = 0.0; 
  xbdot[5] = 0.0; 

  sb[0] = 14.0; // s

  sbdot[0] = 0.0; // sdot


  double C = 1.0e-7;
  double R = 1.0e6;
  double sinfreq = 1000;
  double pi = 4.0*atan(1.0);
	
	// Define constants for the time integration
	double t_final = 1;
	double t_step = t_final/10;
	double t_init = 0*1/sinfreq;
	double t_new = t_step + t_init;
 

  // W = [ beta*df/dx , alpha*I ]
  //     [ beta*dq/dx , -beta*I  ]
	// As we need df/dx and dq/dx we need to set beta=1 
	alpha = 1;
  beta = 1;

	// Start the time integration
	while (t_new <= t_final ) {
  	
		// Print IC:
  	//std::cout << t_init <<"\t" ;
  	//for (int i=0; i< 3; ++i) {
    //	std::cout << (*xbar)[i] << "\t";
 	 	//}
  	//std::cout<<std::endl;

		double t_old = t_new - t_step;
  	//std::cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << std::endl;
  	//std::cout << "NEW TIME POINT" << std::endl;
  	//std::cout << "t_new = " << t_new << std::endl;
  	//std::cout << "t_old= " << t_old << std::endl;
		
		double vsrc = 12.0*sin(2*pi*sinfreq*t_new);
  	//std::cout << "vsrc = " << vsrc << std::endl;
  	//std::cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << std::endl;
		
		// Function evaluations at old time point:
  	//----------------------------------------------------
		// First we call the model to load the state vector s.
  	{
			// update xbar_old [ xbar_old = 0.0*xbar_old + 1.0*xbar ]
			xbar_old->Update(1.0,*xbar,0.0);	
  		//std::cout << "xbar_old = " << *xbar_old << std::endl;
			EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  	  inArgs.set_t(t_old);
  	  inArgs.set_x(xbar_old);
  	  inArgs.set_x_dot(xbardot);
  	  //inArgs.set_p(0,s);
  	  //inArgs.set_p(1,sdot);
  	  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
			/* Question: what is g?? */
			outArgs.set_g(0,s);
  	  //fbar->PutScalar(0.0);
  	  //outArgs.set_f(fbar);
  	  // evalModel
  	  model->evalModel(inArgs,outArgs);
  	}
  	// Verify model loaded s correctly:
  	// TODO
  	// Second we "calculate" sdot
  	// Third we call the model to load f
  	{
  	  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  	  inArgs.set_t(t_old);
  	  inArgs.set_x(xbar_old);
  	  inArgs.set_x_dot(xbardot);
  	  inArgs.set_p(0,s);
  	  inArgs.set_p(1,sdot);
  	  inArgs.set_alpha(alpha);
  	  inArgs.set_beta(beta);
  	  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
  	  fbar->PutScalar(0.0);
  	  outArgs.set_f(fbar);
  	  outArgs.set_W(W);
  	  outArgs.set_g(1,voltLimQ);
  	  outArgs.set_g(2,voltLimF);
  	  // evalModel
  	  model->evalModel(inArgs,outArgs);
  	}
		
		// Pull q(x_old) i.e. at old time point 
		RCP<N_LAS_BlockVector> xyce_q_old = convertEpetraToNLASBlockVectorView(fbar,*small_x_map); 
		N_LAS_Vector & q_old = xyce_q_old->block(1);	
    eq_old->Update(1.0,*q_old.epetraVector(),0.0);
  	//RCP<Epetra_Vector> eq_old = Teuchos::rcpFromRef(*q_old.epetraVector());
		
		//std::cout << "-----------------------------------------------" << std::endl;
  	//std::cout << "q(x_old) vector ="<< std::endl;
		//eq_old->Print(std::cout);	
  	//std::cout << "-----------------------------------------------" << std::endl;

		// update e_rhs_BE  [ e_rhs_BE = 0.0*e_rhs_BE +1.0*eq_old ]
  	e_rhs_BE->Update(1.0,*eq_old,0.0);
		
		// Now e_rhs_BE is q(x_old)
		//std::cout << "-----------------------------------------------" << std::endl;
  	//std::cout << "e_rhs_BE when is equal to q(x_old) "<< std::endl;
		//e_rhs_BE->Print(std::cout);	
  	//std::cout << "-----------------------------------------------" << std::endl;
		
		// Now we will start the NEWTON LOOP to perform the nonlinear solve
		int converged = 0;
		double newtonTol = 1e-10;
		int iter = 0;

		while (!converged) {

      //std::cout << "----------------------------------------------" << std::endl;
      //std::cout << "BEGINNING OF NEWTON LOOP" << std::endl;
      //std::cout << "----------------------------------------------" << std::endl;
			// Function evaluations at new time point:
  		// First we call the model to load the state vector s.
  		{
  		  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  		  inArgs.set_t(t_new);
  		  inArgs.set_x(xbar);
  		  inArgs.set_x_dot(xbardot);
  		  //inArgs.set_p(0,s);
  		  //inArgs.set_p(1,sdot);
  		  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
				/* Question: what is g?? */
				outArgs.set_g(0,s);
  		  //fbar->PutScalar(0.0);
  		  //outArgs.set_f(fbar);
  		  // evalModel
  		  model->evalModel(inArgs,outArgs);
  		}
  		// Verify model loaded s correctly:
  		// TODO
  		// Second we "calculate" sdot
  		// Third we call the model to load f
  		{
  		  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  		  inArgs.set_t(t_new);
  		  inArgs.set_x(xbar);
  		  inArgs.set_x_dot(xbardot);
  		  inArgs.set_p(0,s);
  		  inArgs.set_p(1,sdot);
  		  inArgs.set_alpha(alpha);
  		  inArgs.set_beta(beta);
  		  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
  		  fbar->PutScalar(0.0);
  		  outArgs.set_f(fbar);
  		  outArgs.set_W(W);
  		  outArgs.set_g(1,voltLimQ);
  		  outArgs.set_g(2,voltLimF);
  		  // evalModel
  		  model->evalModel(inArgs,outArgs);
					
  		}
  		//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "xbar = " << *xbar << std::endl;
  		//std::cout << "xbardot = " << *xbardot << std::endl;
  		//std::cout << "s = " << *s << std::endl;
  		//std::cout << "sdot = " << *sdot << std::endl;
  		//std::cout << "fbar = " << *fbar << std::endl;
  		//std::cout << "-----------------------------------------------" << std::endl;
			
			// Now pull out the block vector from F(xbar, xbar_dot, t)  and block
			// matrix W which corresponds to f(x) and df_dx 
				
			// Pull f(x_new) 
			RCP<N_LAS_BlockVector> xyce_f_new = convertEpetraToNLASBlockVectorView(fbar,*small_x_map); 
			N_LAS_Vector & f_new = xyce_f_new->block(0);	
  		RCP<Epetra_Vector> ef_new = Teuchos::rcpFromRef(*f_new.epetraVector());

			// Pull q(x_new) i.e. at new time point 
			RCP<N_LAS_BlockVector> xyce_q_new= convertEpetraToNLASBlockVectorView(fbar,*small_x_map); 
			N_LAS_Vector & q_new = xyce_q_new->block(1);	
  		RCP<Epetra_Vector> eq_new = Teuchos::rcpFromRef(*q_new.epetraVector());
  		
			std::string label = "N_LAS_BlockMatrix";
  		RCP<N_LAS_BlockMatrix> bW = Teuchos::get_extra_data<RCP<N_LAS_BlockMatrix> >(W,label);

  		N_LAS_Matrix & ul = bW->block(0,0);
  		//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "ul.epetraObj()=" << std::endl;
  		//ul.printPetraObject();
  		////ul.epetraObj().Graph().Print(std::cout);
  		//std::cout << "-----------------------------------------------" << std::endl;
  		
			N_LAS_Matrix & ll = bW->block(1,0);
  		//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "ll.epetraObj()=" << std::endl;
  		//ll.printPetraObject();
  		////ll.epetraObj().Graph().Print(std::cout);
  		//std::cout << "-----------------------------------------------" << std::endl;


			// Assemble the rhs for BE method: (q(x_new) - q(x_old))/Delta_t
			// + f(x_new) = 0
			//=========================================================================
			// update e_rhs_BE [ e_rhs_BE = - 1.0*e_rhs_BE + 1.0*eq_new ]
      
      // update e_rhs_BE  [ e_rhs_BE = 0.0*e_rhs_BE +1.0*eq_old ]
      e_rhs_BE->Update(1.0,*eq_old,0.0);
      // Now e_rhs_BE is q(x_old)
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "e_rhs_BE when is equal to q(x_old) "<< std::endl;
      //e_rhs_BE->Print(std::cout);	
      //std::cout << "-----------------------------------------------" << std::endl;
			
  	  e_rhs_BE->Update(1.0,*eq_new,-1.0);
  	  e_rhs_BE->Scale(1/t_step);
			//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "diff of q(x_old) - q(x(new) ="<< std::endl;
			//e_rhs_BE->Print(std::cout); 
			//std::cout << "-----------------------------------------------" << std::endl;
			// finally update e_rhs_BE [ e_rhs_BE = 1.0*e_rhs_BE + 1.0*f_new ]
  	  e_rhs_BE->Update(1.0,*ef_new,1.0);
			//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "[q(x_old) - q(x(new)]/Delta_t + f(x_new) ="<< std::endl;
			//e_rhs_BE->Print(std::cout); 
			//std::cout << "-----------------------------------------------" << std::endl;
			
			double e_rhs_BE_Norm;
      e_rhs_BE->Norm1(&e_rhs_BE_Norm);
  		//std::cout << "norm(e_rhs_BE) =" << e_rhs_BE_Norm  << std::endl;
  		//std::cout << "-----------------------------------------------" << std::endl;
      {
				//std::cout << "-----------------------------------------------" << std::endl;
				//std::cout << "***********************************************" << std::endl;
  			//std::cout << "Newton Raphson for RC load iteration = "<< iter <<  std::endl;
  			//std::cout << "xbar = " << *xbar << std::endl;
				//std::cout << "***********************************************" << std::endl;
  			//std::cout << "-----------------------------------------------" << std::endl;
      }
			if ( e_rhs_BE_Norm <= newtonTol) {
				//std::cout << "-----------------------------------------------" << std::endl;
				//std::cout << "***********************************************" << std::endl;
  			//std::cout << "Newton Raphson for RC load Converged in "<< iter << " iterations." <<  std::endl;
  			//std::cout << "Solution: xbar = " << *xbar << std::endl;
				//std::cout << "***********************************************" << std::endl;
  			//std::cout << "-----------------------------------------------" << std::endl;
				converged = 1;
				break;
			}		
      if (iter>10) {
				std::cout << "-----------------------------------------------" << std::endl;
				std::cout << "***********************************************" << std::endl;
  			std::cout << "Newton Raphson for RC load did not converge in max iterations = "<< iter << std::endl;
  			std::cout << "Solution: xbar = " << *xbar << std::endl;
				std::cout << "***********************************************" << std::endl;
  			std::cout << "-----------------------------------------------" << std::endl;
        break;
      }

			// Assemble the d_rhs/d_x for BE method:
			// [dq(x_new)/dx_new]/Delta_t + df_x_new/dx_new
  		Epetra_CrsMatrix& drhs_dx = ll.epetraObj();
			//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "dq(x_new)/dx ="<< std::endl;
  		//drhs_dx.Print(std::cout);
			//std::cout << "-----------------------------------------------" << std::endl;
			// Scale the dq(x_new)/dx by 1/Delta_t
  		drhs_dx.Scale(1/t_step);
			//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "[dq(x_new)/dx]*1/Delta_t ="<< std::endl;
  		//drhs_dx.Print(std::cout);
			//std::cout << "-----------------------------------------------" << std::endl;
			// Now add the df(x_new)/dx 
  		Epetra_CrsMatrix& df_dx = ul.epetraObj();
			EpetraExt::MatrixMatrix::Add(df_dx, 0, 1.0, drhs_dx, 1.0);  // 0 = no transpose
			//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "[dq(x_new)/dx]*1/Delta_t + df(x_new)/dx = "<< std::endl;
  		//drhs_dx.Print(std::cout);
			//std::cout << "-----------------------------------------------" << std::endl;

			//Linear solve begins
  		RCP<Epetra_Vector> Deltax = rcp(new Epetra_Vector(*small_x_map)); 
			Epetra_LinearProblem problem(&drhs_dx,&*Deltax,&*e_rhs_BE);
  	  Amesos_Klu solver(problem);
  	  solver.SetUseTranspose(false);
  	  solver.SymbolicFactorization(); // Repeat if nonzeros in A change
  	  solver.NumericFactorization(); // Repeat if values in A change
  	  solver.Solve(); // Repeat if values in A or b change

			//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "Delta x="<< std::endl;
			//Deltax->Print(std::cout);	
  		//std::cout << "-----------------------------------------------" << std::endl;
      double e_dx_BE_Norm;
      Deltax->Norm1(&e_dx_BE_Norm);
  		//std::cout << "norm(e_dx_BE) =" << e_dx_BE_Norm  << std::endl;

			// update xbar  [ xbar = 1.0*xbar -1.0*Deltax ]
  	  xbar_block_0->Update(-1.0,*Deltax,1.0);
			iter++;
      
			// std::cout << "----------------------------------------------" << std::endl; 
			// std::cout << "END OF NEWTON LOOP, returning to top" << std::endl; 
			// std::cout << "----------------------------------------------" << std::endl;
		}
		
		// Here we compute the analytical solution for the first time step
		// and compare it with the one obtained from ModelEvaluator
		// interface.
		// The analytical solution is as below:
		//  omega = 1/ (C/Delta_t + 1/R)
		// | V2_{n+1} | = omega* [ C*V2_{n}/Delta_t + Vsrc_{n+1}/R]
		// | Isrc_{n+1} | = omega* ([ C*V2_{n}/Delta_t + Vsrc_{n+1}/R]*1/R - 
		// [C/Delta_t + 1/R]*Vsrc_{n+1}/R)
	
  double tol = 1.0e-8;
	
	{ 
		omega =  1/ (C/t_step + 1/R);  // 1/Det(BE_LHS)
		// during the first iteration get IC from the user.
		if(!getInitICFlag) {
			//std::cout << "-----------------------------------------------" << std::endl; 
			//std::cout << "Compare the solution with analytical solution" << std::endl;
			V2_n = (*xbar_old)[0]; 
			Isrc_n = (*xbar_old)[1]; 
			V1_n = (*xbar_old)[2]; 
			getInitICFlag = 1;
		}
		else {
			V2_n = V2_next; 
		}
		V2_next = omega* (C*V2_n/t_step + vsrc/R);
		Isrc_next = omega* (( C*V2_n/t_step + vsrc/R)*1/R - (C/t_step + 1/R)*vsrc/R);
 		//std::cout << t_new <<"\t" ;
		//std::cout << (*xbar)[0] <<"\t"<< V2_next; //V2 
		//std::cout << "\t" << (*xbar)[1] <<"\t"<< Isrc_next;  //Iscr
		//std::cout << "\t" << (*xbar)[2] <<"\t"<< vsrc<<std::endl ; //V1
		TEST_FLOATING_EQUALITY( (*xbar)[0], V2_next, tol ); // v2
  	TEST_FLOATING_EQUALITY( (*xbar)[1], Isrc_next, tol ); //I_src
	}
	t_new = t_step + t_new;
	//std::cout << "-------Next ti--------------------" << std::endl;
	}
}

#endif

# if 1
TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, DiodeCapLoad_Backward_Euler) {
  int iargs = 2;
  char *cargs[iargs];
  cargs[0] = "Xyce";
  cargs[1] = "forward_diode_model.cir";
  // Use simple Diode-Capacitor circuit netlist as input
  // Set up t,x,xdot,s,sdot,f,W
  RCP<N_ANP_ModelEvaluator> model = rcp(new N_ANP_ModelEvaluator());
  model->initialize(iargs,cargs);
  RCP<const Epetra_Map> x_map = model->get_x_map();
  RCP<const Epetra_Map> f_map = model->get_f_map();
  RCP<const Epetra_Map> s_map = model->get_p_map(0);
  RCP<const Epetra_Map> sdot_map = model->get_p_map(1);
  RCP<const Epetra_Map> small_x_map = model->get_small_x_map();
  RCP<Epetra_Operator> W = model->create_W();
  
  double alpha,beta;
  double t;

  RCP<Epetra_Vector> xbar = rcp(new Epetra_Vector(*x_map)); // stores x and z 
  RCP<N_LAS_BlockVector> nlas_block_xbar = convertEpetraToNLASBlockVectorView(xbar,*small_x_map);
  RCP<Epetra_Vector> xbar_block_0 = Teuchos::rcpFromRef(*(nlas_block_xbar->block(0).epetraVector())); // This a view!
  RCP<Epetra_Vector> xbardot = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> fbar = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> eq_old = rcp(new Epetra_Vector(*small_x_map));
  RCP<Epetra_Vector> s = rcp(new Epetra_Vector(*s_map)); 
  RCP<Epetra_Vector> sdot = rcp(new Epetra_Vector(*sdot_map)); 
  RCP<Epetra_Vector> voltLimF = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> voltLimQ = rcp(new Epetra_Vector(*f_map)); 
 	
	// This will be used to store the values at old time points
  RCP<Epetra_Vector> xbar_old = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> fbar_old = rcp(new Epetra_Vector(*f_map)); 
	
	// Make a temporary vector to compute the rhs of BE 
  RCP<Epetra_Vector> e_rhs_BE = rcp(new Epetra_Vector(*small_x_map));

  Epetra_Vector& xb = *xbar;
  Epetra_Vector& xbdot = *xbardot;
  Epetra_Vector& sb = *s;
  Epetra_Vector& sbdot = *sdot;

  //two blocks of three elements each.
  xb[0] = 0; // v_2
  xb[1] = 0; // i_vsrc
  xb[2] = 0.2; // v_1
  double v2 = xb[0]; // 0.2 
  double i_vsrc = xb[1]; // 4.0e-6
  double v1 = xb[2]; // 0.4

	// Note: These are required to be set to zero, becuase we want to
	// extract the f_func and q_func from F.
	// F(xbar,x,t)|_xbar=0 => F(0,x,t) = [f_func; q_func];
  xb[3] = 0.0; // z
  xb[4] = 0.0; 
  xb[5] = 0.0; 

  xbdot[0] = 0.0; // v_2 dot
  xbdot[1] = 0.0; // i_vsrc dot
  xbdot[2] = 0.0; // v_1 dot

  xbdot[3] = 0.0; // z dot
  xbdot[4] = 0.0; 
  xbdot[5] = 0.0; 

  sb[0] = 14.0; // s

  sbdot[0] = 0.0; // sdot

	// Capacitor Parms
  double C = 1.0e-12;
	// Diode Parms: Forward model of the diode 
  double IS = 1.0e-14; // Saturation current
  double N  = 1; // Mult Factor
  double Vth  = 0.025864186; // Threshold Factor
  double pi = 4.0*atan(1.0);
	double vsrc = 0.2; 
	double R = 1e6;

	// Define constants for the time integration
	double t_final = 10e-6;
	double t_step = t_final/10;
	double t_init = 0;
	double t_new = t_step + t_init;

	// The v2 is computed from Xyce and is stored in a STL map and used
	// for comparison with what ModelEvaluator computes:
	// TODO:FIXME: In a sense this is hardcoded, but C++ does not
	// provide accurate values while computing exp(.).
	std::map< int, double > diode_sim_data ;  
	
	diode_sim_data[1]	= 2.27874483e-05;
  diode_sim_data[2]	= 4.55548464e-05;
  diode_sim_data[3]	= 6.83022208e-05;
  diode_sim_data[4]	= 9.10296066e-05;
  diode_sim_data[5]	= 1.13737039e-04;
  diode_sim_data[6]	= 1.36424553e-04;
  diode_sim_data[7]	= 1.59092183e-04;
  diode_sim_data[8]	= 1.81739964e-04;
  diode_sim_data[9]	= 2.04367931e-04;
  diode_sim_data[10]= 2.26976119e-04;
	
  // W = [ beta*df/dx , alpha*I ]
  //     [ beta*dq/dx , -beta*I  ]
	// As we need df/dx and dq/dx we need to set beta=1 
	alpha = 1;
  beta = 1;
	// Start the time integration
	while (t_new <= t_final ) {
		//// Print IC:
  	//std::cout << t_init <<"\t" ;
  	//for (int i=0; i< 3; ++i) {
    //	std::cout << (*xbar)[i] << "\t";
 	 	//}
  	//std::cout<<std::endl;

		double t_old = t_new - t_step;
  	//std::cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << std::endl;
  	//std::cout << "NEW TIME POINT" << std::endl;
  	//std::cout << "t_new = " << t_new << std::endl;
  	//std::cout << "t_old= " << t_old << std::endl;
		
		// Function evaluations at old time point:
  	//----------------------------------------------------
		// First we call the model to load the state vector s.
  	{
			// update xbar_old [ xbar_old = 0.0*xbar_old + 1.0*xbar ]
			xbar_old->Update(1.0,*xbar,0.0);	
  		//std::cout << "xbar_old = " << *xbar_old << std::endl;
			EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  	  inArgs.set_t(t_old);
  	  inArgs.set_x(xbar_old);
  	  inArgs.set_x_dot(xbardot);
  	  //inArgs.set_p(0,s);
  	  //inArgs.set_p(1,sdot);
  	  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
			/* Question: what is g?? */
			outArgs.set_g(0,s);
  	  //fbar->PutScalar(0.0);
  	  //outArgs.set_f(fbar);
  	  // evalModel
  	  model->evalModel(inArgs,outArgs);
  	}
  	// Verify model loaded s correctly:
  	// TODO
  	// Second we "calculate" sdot
  	// Third we call the model to load f
  	{
  	  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  	  inArgs.set_t(t_old);
  	  inArgs.set_x(xbar_old);
  	  inArgs.set_x_dot(xbardot);
  	  inArgs.set_p(0,s);
  	  inArgs.set_p(1,sdot);
  	  inArgs.set_alpha(alpha);
  	  inArgs.set_beta(beta);
  	  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
  	  fbar->PutScalar(0.0);
  	  outArgs.set_f(fbar);
  	  outArgs.set_W(W);
  	  outArgs.set_g(1,voltLimQ);
  	  outArgs.set_g(2,voltLimF);
  	  // evalModel
  	  model->evalModel(inArgs,outArgs);
  	}
		
		// Pull q(x_old) i.e. at old time point 
		RCP<N_LAS_BlockVector> xyce_q_old = convertEpetraToNLASBlockVectorView(fbar,*small_x_map); 
		N_LAS_Vector & q_old = xyce_q_old->block(1);	
    eq_old->Update(1.0,*q_old.epetraVector(),0.0);
  	//RCP<Epetra_Vector> eq_old = Teuchos::rcpFromRef(*q_old.epetraVector());
		
		//std::cout << "-----------------------------------------------" << std::endl;
  	//std::cout << "q(x_old) vector ="<< std::endl;
		//eq_old->Print(std::cout);	
  	//std::cout << "-----------------------------------------------" << std::endl;

		// update e_rhs_BE  [ e_rhs_BE = 0.0*e_rhs_BE +1.0*eq_old ]
  	e_rhs_BE->Update(1.0,*eq_old,0.0);
		
		// Now e_rhs_BE is q(x_old)
		//std::cout << "-----------------------------------------------" << std::endl;
  	//std::cout << "e_rhs_BE when is equal to q(x_old) "<< std::endl;
		//e_rhs_BE->Print(std::cout);	
  	//std::cout << "-----------------------------------------------" << std::endl;
		
		// Now we will start the NEWTON LOOP to perform the nonlinear solve
		int converged = 0;
		double newtonTol = 1e-18;
		int iter = 0;

		while (!converged) {

      //std::cout << "----------------------------------------------" << std::endl;
      //std::cout << "BEGINNING OF NEWTON LOOP" << std::endl;
      //std::cout << "----------------------------------------------" << std::endl;
			// Function evaluations at new time point:
  		// First we call the model to load the state vector s.
  		{
  		  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  		  inArgs.set_t(t_new);
  		  inArgs.set_x(xbar);
  		  inArgs.set_x_dot(xbardot);
  		  //inArgs.set_p(0,s);
  		  //inArgs.set_p(1,sdot);
  		  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
				/* Question: what is g?? */
				outArgs.set_g(0,s);
  		  //fbar->PutScalar(0.0);
  		  //outArgs.set_f(fbar);
  		  // evalModel
  		  model->evalModel(inArgs,outArgs);
  		}
  		// Verify model loaded s correctly:
  		// TODO
  		// Second we "calculate" sdot
  		// Third we call the model to load f
  		{
  		  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  		  inArgs.set_t(t_new);
  		  inArgs.set_x(xbar);
  		  inArgs.set_x_dot(xbardot);
  		  inArgs.set_p(0,s);
  		  inArgs.set_p(1,sdot);
  		  inArgs.set_alpha(alpha);
  		  inArgs.set_beta(beta);
  		  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
  		  fbar->PutScalar(0.0);
  		  outArgs.set_f(fbar);
  		  outArgs.set_W(W);
  		  outArgs.set_g(1,voltLimQ);
  		  outArgs.set_g(2,voltLimF);
  		  // evalModel
  		  model->evalModel(inArgs,outArgs);
					
  		}
  		//std::cout << "xbar = " << *xbar << std::endl;
  		//std::cout << "xbardot = " << *xbardot << std::endl;
  		//std::cout << "s = " << *s << std::endl;
  		//std::cout << "sdot = " << *sdot << std::endl;
  		//std::cout << "fbar = " << *fbar << std::endl;
			
			// Now pull out the block vector from F(xbar, xbar_dot, t)  and block
			// matrix W which corresponds to f(x) and df_dx 
				
			// Pull f(x_new) 
			RCP<N_LAS_BlockVector> xyce_f_new = convertEpetraToNLASBlockVectorView(fbar,*small_x_map); 
			N_LAS_Vector & f_new = xyce_f_new->block(0);	
  		RCP<Epetra_Vector> ef_new = Teuchos::rcpFromRef(*f_new.epetraVector());

			// Pull q(x_new) i.e. at new time point 
			RCP<N_LAS_BlockVector> xyce_q_new= convertEpetraToNLASBlockVectorView(fbar,*small_x_map); 
			N_LAS_Vector & q_new = xyce_q_new->block(1);	
  		RCP<Epetra_Vector> eq_new = Teuchos::rcpFromRef(*q_new.epetraVector());
  		
			std::string label = "N_LAS_BlockMatrix";
  		RCP<N_LAS_BlockMatrix> bW = Teuchos::get_extra_data<RCP<N_LAS_BlockMatrix> >(W,label);

  		N_LAS_Matrix & ul = bW->block(0,0);
  		//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "ul.epetraObj()=" << std::endl;
  		//ul.printPetraObject();
  		////ul.epetraObj().Graph().Print(std::cout);
  		//std::cout << "-----------------------------------------------" << std::endl;
  		
			N_LAS_Matrix & ll = bW->block(1,0);
  		//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "ll.epetraObj()=" << std::endl;
  		//ll.printPetraObject();
  		////ll.epetraObj().Graph().Print(std::cout);
  		//std::cout << "-----------------------------------------------" << std::endl;

			// Assemble the rhs for BE method: (q(x_new) - q(x_old))/Delta_t
			// + f(x_new) = 0
			//=========================================================================
			// update e_rhs_BE [ e_rhs_BE = - 1.0*e_rhs_BE + 1.0*eq_new ]
      
      // update e_rhs_BE  [ e_rhs_BE = 0.0*e_rhs_BE +1.0*eq_old ]
      e_rhs_BE->Update(1.0,*eq_old,0.0);
      // Now e_rhs_BE is q(x_old)
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "e_rhs_BE when is equal to q(x_old) "<< std::endl;
      //e_rhs_BE->Print(std::cout);	
      //std::cout << "-----------------------------------------------" << std::endl;
			
  	  e_rhs_BE->Update(1.0,*eq_new,-1.0);
  	  e_rhs_BE->Scale(1/t_step);
			//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "diff of q(x_old) - q(x(new) ="<< std::endl;
			//e_rhs_BE->Print(std::cout); 
			//std::cout << "-----------------------------------------------" << std::endl;
			// finally update e_rhs_BE [ e_rhs_BE = 1.0*e_rhs_BE + 1.0*f_new ]
  	  e_rhs_BE->Update(1.0,*ef_new,1.0);
			//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "[q(x_old) - q(x(new)]/Delta_t + f(x_new) ="<< std::endl;
			//e_rhs_BE->Print(std::cout); 
			//std::cout << "-----------------------------------------------" << std::endl;
			
			double e_rhs_BE_Norm;
      e_rhs_BE->Norm1(&e_rhs_BE_Norm);
  		//std::cout << "norm(e_rhs_BE) =" << e_rhs_BE_Norm  << std::endl;
  		//std::cout << "-----------------------------------------------" << std::endl;
      {
				//std::cout << "-----------------------------------------------" << std::endl;
				//std::cout << "***********************************************" << std::endl;
  			//std::cout << "Newton Raphson for DiodeCap load iteration = "<< iter <<  std::endl;
  			//std::cout << "xbar = " << *xbar << std::endl;
				//std::cout << "***********************************************" << std::endl;
  			//std::cout << "-----------------------------------------------" << std::endl;
      }
			if ( e_rhs_BE_Norm <= newtonTol) {
				//std::cout << "-----------------------------------------------" << std::endl;
				//std::cout << "***********************************************" << std::endl;
  			//std::cout << "Newton Raphson for DiodeCap load Converged in "<< \ 
				//  iter << " iterations." <<  std::endl;
  			//std::cout << "Solution: xbar = " << *xbar << std::endl;
				//std::cout << "***********************************************" << std::endl;
  			//std::cout << "-----------------------------------------------" << std::endl;
				converged = 1;
				break;
			}		
      if (iter>100) {
				//std::cout << "-----------------------------------------------" << std::endl;
				//std::cout << "***********************************************" << std::endl;
  			std::cout << "Newton Raphson for DiodeCap load did not converge in \
					max iterations = "<< iter << std::endl;
  			//std::cout << "Solution: xbar = " << *xbar << std::endl;
				//std::cout << "***********************************************" << std::endl;
  			//std::cout << "-----------------------------------------------" << std::endl;
        break;
      }

			// Assemble the d_rhs/d_x for BE method:
			// [dq(x_new)/dx_new]/Delta_t + df_x_new/dx_new
  		Epetra_CrsMatrix& drhs_dx = ll.epetraObj();
			//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "dq(x_new)/dx ="<< std::endl;
  		//drhs_dx.Print(std::cout);
			//std::cout << "-----------------------------------------------" << std::endl;
			// Scale the dq(x_new)/dx by 1/Delta_t
  		drhs_dx.Scale(1/t_step);
			//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "[dq(x_new)/dx]*1/Delta_t ="<< std::endl;
  		//drhs_dx.Print(std::cout);
			//std::cout << "-----------------------------------------------" << std::endl;
			// Now add the df(x_new)/dx 
  		Epetra_CrsMatrix& df_dx = ul.epetraObj();
			EpetraExt::MatrixMatrix::Add(df_dx, 0, 1.0, drhs_dx, 1.0);  // 0 = no transpose
			//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "[dq(x_new)/dx]*1/Delta_t + df(x_new)/dx = "<< std::endl;
  		//drhs_dx.Print(std::cout);
			//std::cout << "-----------------------------------------------" << std::endl;

			//Linear solve begins
  		RCP<Epetra_Vector> Deltax = rcp(new Epetra_Vector(*small_x_map)); 
			Epetra_LinearProblem problem(&drhs_dx,&*Deltax,&*e_rhs_BE);
  	  Amesos_Klu solver(problem);
  	  solver.SetUseTranspose(false);
  	  solver.SymbolicFactorization(); // Repeat if nonzeros in A change
  	  solver.NumericFactorization(); // Repeat if values in A change
  	  solver.Solve(); // Repeat if values in A or b change

			//std::cout << "-----------------------------------------------" << std::endl;
  		//std::cout << "Delta x="<< std::endl;
			//Deltax->Print(std::cout);	
  		//std::cout << "-----------------------------------------------" << std::endl;
      double e_dx_BE_Norm;
      Deltax->Norm1(&e_dx_BE_Norm);
  		//std::cout << "norm(e_dx_BE) =" << e_dx_BE_Norm  << std::endl;
			// update xbar  [ xbar = 1.0*xbar -1.0*Deltax ]
  	  xbar_block_0->Update(-1.0,*Deltax,1.0);
			iter++;
			// std::cout << "----------------------------------------------" << std::endl; 
			// std::cout << "END OF NEWTON LOOP, returning to top" << std::endl; 
			// std::cout << "----------------------------------------------" << std::endl;
		}
		
		// Here we compute the analytical solution for the first time step
		// and compare it with the one obtained from ModelEvaluator
		// interface.
	
  	double tol = 1.0e-6;
		{ 
  		//std::cout << "t_new = " << t_new << std::endl;
			//std::cout << (*xbar)[0] <<"\t"<< diode_sim_data[ (int) (t_new/t_step)] << std::endl; //V2 
			TEST_FLOATING_EQUALITY( (*xbar)[0], diode_sim_data[(int )(t_new/t_step)], tol ); // v2
		}
		t_new = t_step + t_new;
	}
}
#endif

# if 1
TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, BlockMatrixSolve ) {
  int iargs = 2;
  char *cargs[iargs];
  cargs[0] = "Xyce";
  cargs[1] = "modelEvaluatorTest.cir";
  int numEntries = 6;
  // R and C values of the modelEvaluatorTest.cir to calculate the
  // A_exact for testing 
  double C = 1.0e-7;
  double R = 1.0e6;


  // Use simple RC circuit netlist as input
  // Set up t,x,xdot,s,sdot,f,W
  RCP<N_ANP_ModelEvaluator> model = rcp(new N_ANP_ModelEvaluator());
  model->initialize(iargs,cargs);
  RCP<const Epetra_Map> x_map = model->get_x_map();
  RCP<const Epetra_Map> f_map = model->get_f_map();
  RCP<const Epetra_Map> s_map = model->get_p_map(0);
  RCP<const Epetra_Map> sdot_map = model->get_p_map(1);
  RCP<const Epetra_Map> small_x_map = model->get_small_x_map();
  RCP<Epetra_Operator> W = model->create_W();
  RCP<Epetra_Operator> W_exact = model->create_W();
  
  double alpha,beta;
  double t;
  RCP<Epetra_Vector> xbar = rcp(new Epetra_Vector(*x_map)); 
  RCP<N_LAS_BlockVector> nlas_block_xbar = convertEpetraToNLASBlockVectorView(xbar,*small_x_map);
  RCP<Epetra_Vector> xbar_block_0 = Teuchos::rcpFromRef(*(nlas_block_xbar->block(0).epetraVector())); // This a view!
  RCP<Epetra_Vector> xbardot = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> fbar = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> s = rcp(new Epetra_Vector(*s_map)); 
  RCP<Epetra_Vector> sdot = rcp(new Epetra_Vector(*sdot_map)); 
  RCP<Epetra_Vector> voltLimF = rcp(new Epetra_Vector(*f_map)); 
  RCP<Epetra_Vector> voltLimQ = rcp(new Epetra_Vector(*f_map)); 
  
  Epetra_Vector& xb = *xbar;
  Epetra_Vector& xbdot = *xbardot;
  Epetra_Vector& sb = *s;
  Epetra_Vector& sbdot = *sdot;
  t = 0.0000; // 0.00025;
  //two blocks of three elements each.
  xb[0] = 2.0; // v_2
  xb[1] = 3.0; // i_vsrc
  xb[2] = 4.0; // v_1
  double v2 = xb[0]; // 2.0
  double i_vsrc = xb[1]; // 3.0
  double v1 = xb[2]; // 4.0

  xb[3] = 5.0; // z
  xb[4] = 6.0; 
  xb[5] = 7.0; 

  xbdot[0] = 8.0; // v_2 dot
  xbdot[1] = 9.0; // i_vsrc dot
  xbdot[2] = 10.0; // v_1 dot

  xbdot[3] = 11.0; // z dot
  xbdot[4] = 12.0; 
  xbdot[5] = 13.0; 

  sb[0] = 14.0; // s

  sbdot[0] = 0.0; // sdot

  alpha = 1;
  beta = 1;

{
  	  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  	  inArgs.set_t(t);
  	  inArgs.set_x(xbar);
  	  inArgs.set_x_dot(xbardot);
  	  //inArgs.set_p(0,s);
  	  //inArgs.set_p(1,sdot);
  	  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
  	  
			/* Question: what is g?? */
			outArgs.set_g(0,s);
  	  //fbar->PutScalar(0.0);
  	  //outArgs.set_f(fbar);

  	  // evalModel
  	  model->evalModel(inArgs,outArgs);
  	}
  	// Verify model loaded s correctly:
  	// TODO
  	// Second we "calculate" sdot
  	// Third we call the model to load f
  	{
  	  EpetraExt::ModelEvaluator::InArgs inArgs = model->createInArgs();
  	  inArgs.set_t(t);
  	  inArgs.set_x(xbar);
  	  inArgs.set_x_dot(xbardot);
  	  inArgs.set_p(0,s);
  	  inArgs.set_p(1,sdot);
  	  inArgs.set_alpha(alpha);
  	  inArgs.set_beta(beta);
  	  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();
  	  fbar->PutScalar(0.0);
  	  outArgs.set_f(fbar);
  	  outArgs.set_W(W);
  	  outArgs.set_g(1,voltLimQ);
  	  outArgs.set_g(2,voltLimF);
  	  // evalModel
  	  model->evalModel(inArgs,outArgs);
      //std::cout << "-----------------------------------------------" << std::endl;
      //RCP<Epetra_CrsMatrix> W_matrix = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W,true);
      //std::cout << "W_matrix = " << std::endl;
      //W_matrix->Print(std::cout);
      //std::cout << "-----------------------------------------------" << std::endl;
  	}
  	
  RCP<Epetra_Vector> x = rcp(new Epetra_Vector(*x_map)); 
  RCP<Epetra_Vector> b = rcp(new Epetra_Vector(*x_map)); 
  { 
    // Fill b with 1..numEntries
    for (int row=0 ; row<numEntries; ++row) {
      (*b)[row] = row;
    }
  }
  // Fill F and W with evalModel
  RCP<Epetra_CrsMatrix> A = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W,true);
  {
    // Solve for bogus b
    Epetra_LinearProblem problem(&*A,&*x,&*b);
    Amesos_Klu solver(problem);
    solver.SetUseTranspose(false);
    solver.SymbolicFactorization(); // Repeat if nonzeros in A change
    solver.NumericFactorization(); // Repeat if values in A change
    solver.Solve(); // Repeat if values in A or b change
  } 
  
  //----Check that the system can be solved
  // Now compute the actual solution manually by composing block matrix and
  // solving it.
  // Ideally we should have following structure of W from the ModelEvaluator
  // df/dx:
  // df0/dv2, df0/di, df0/dv1 = [ 1/R, 0, -1/R ]
  // df1/dv2, df1/di, df1/dv1 = [ 0,   0,  1   ]
  // df2/dv2, df2/di, df2/dv1 = [-1/R, 1,  1/R ]

  // dq/dx:
  // dq0/dv2, dq0/di, dq0/dv1 = [ C, 0, 0 ]
  // dq1/dv2, dq1/di, dq1/dv1 = [ 0, 0, 0 ]
  // dq2/dv2, dq2/di, dq2/dv1 = [ 0, 0, 0 ]

  // W_exact = [ beta*df/dx , alpha*I ]
  //           [ beta*dq/dx , -beta*I  ]
  std::string label = "N_LAS_BlockMatrix";
  RCP<N_LAS_BlockMatrix> bW_exact = Teuchos::get_extra_data<RCP<N_LAS_BlockMatrix> >(W_exact,label);
  //std::cout << "-----------------------------------------------" << std::endl;
  //std::cout << "bW_exact->epetraObj().Graph() = " << std::endl;
  //bW_exact->epetraObj().Graph().Print(std::cout);
  //std::cout << "-----------------------------------------------" << std::endl;
  bW_exact->put(0.0);
  {
    N_LAS_Matrix & ul = bW_exact->block(0,0);
    int row;
    int length;
    double coeffs[3];
    int colIndices[3];
    row = 0;
    length=2;
    coeffs[0] = beta/R;
    coeffs[1] = -beta/R;
    colIndices[0] = 0;
    colIndices[1] = 2;
    ul.putRow(row,length,coeffs,colIndices);
    row = 1;
    length = 1;
    coeffs[0] = beta;
    colIndices[0] = 2;
    ul.putRow(row,length,coeffs,colIndices);
    row = 2;
    length=3;
    coeffs[0] = -beta/R;
    coeffs[1] = beta;
    coeffs[2] = beta/R;
    colIndices[0] = 0;
    colIndices[1] = 1;
    colIndices[2] = 2;
    ul.putRow(row,length,coeffs,colIndices);
  }
  {
    N_LAS_Matrix & ur = bW_exact->block(0,1);
    int row;
    int length = 1;
    double coeffs[1];
    int colIndices[1];
    row = 0;
    coeffs[0] = alpha;
    colIndices[0] = 0;
    ur.putRow(row,length,coeffs,colIndices);
    row = 1;
    colIndices[0] = 1;
    ur.putRow(row,length,coeffs,colIndices);
    row = 2;
    colIndices[0] = 2;
    ur.putRow(row,length,coeffs,colIndices);
  }
  {
    N_LAS_Matrix & ll = bW_exact->block(1,0);
    int row;
    int length;
    double coeffs[3];
    int colIndices[3];
    row = 0;
    length = 1;
    coeffs[0] = beta*C;
    colIndices[0] = 0;
    ll.putRow(row,length,coeffs,colIndices);
  }
  {
    N_LAS_Matrix & lr = bW_exact->block(1,1);
    int row;
    int length = 1;
    double coeffs[1];
    int colIndices[1];
    row = 0;
    coeffs[0] = -beta;
    colIndices[0] = 0;
    lr.putRow(row,length,coeffs,colIndices);
    row = 1;
    colIndices[0] = 1;
    lr.putRow(row,length,coeffs,colIndices);
    row = 2;
    colIndices[0] = 2;
    lr.putRow(row,length,coeffs,colIndices);
  }
  bW_exact->fillComplete();
  //std::cout << "-----------------------------------------------" << std::endl;
  //std::cout << "bW_exact->epetraObj().Graph() = " << std::endl;
  //bW_exact->epetraObj().Graph().Print(std::cout);
  //std::cout << "-----------------------------------------------" << std::endl;
  //std::cout << "bW_exact = " << std::endl;
  //bW_exact->printPetraObject();
  //std::cout << "-----------------------------------------------" << std::endl;
  RCP<Epetra_CrsMatrix> A_exact = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_exact,true);
  //std::cout << "A_exact = " << std::endl;
  //A_exact->Print(std::cout);
  //std::cout << "-----------------------------------------------" << std::endl;
  
  RCP<Epetra_Vector> x_exact = rcp(new Epetra_Vector(*x_map)); 
  
  // Solve for bogus b
  {
    Epetra_LinearProblem problem_exact(&*A_exact,&*x_exact,&*b);
    Amesos_Klu solver(problem_exact);
    solver.SetUseTranspose(false);
    solver.SymbolicFactorization(); // Repeat if nonzeros in A change
    solver.NumericFactorization(); // Repeat if values in A change
    solver.Solve(); // Repeat if values in A or b change
  }

  double tol = 1.0e-6;
  for( int i=0; i<numEntries; ++i) {
    //std::cout<< " i = " << i << "\t";
    //std::cout<< " x(i) = " << (*x)[i]<< "\t";
    //std::cout<< " x_exact = " << (*x_exact)[i]<<endl;
		TEST_FLOATING_EQUALITY( (*x)[i],(*x_exact)[i], tol ); 
	}

}

#endif

#ifdef Xyce_TRILINOS_DEV_RYTHMOS
TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, rythmos_BE ) {
  
  int iargs = 2;
  char *cargs[iargs];
  cargs[0] = "Xyce";
  // Use simple RC circuit netlist as input
  cargs[1] = "modelEvaluatorTest.cir";
  
 // vars related with analytical solution of the BE for RC circuit
	int getInitICFlag =  0;
	double V2_n; //V2  at nth step 
	double Isrc_n; // I_src at nth step
	double V1_n; // V1 at nth step
	double omega ; // 1/Det(BE_LHS)
	double V2_next; 
	double Isrc_next; 
  double C = 1.0e-7;
  double R = 1.0e6;
  double sinfreq = 1000;
  double pi = 4.0*atan(1.0);
  
  RCP<N_ANP_ModelEvaluator> xyceModel = rcp(new N_ANP_ModelEvaluator());
  xyceModel->initialize(iargs,cargs);
  
  RCP<N_ANP_ModelEvaluator_Stateless> emodel = 
			N_ANP_modelEvaluator_Stateless(xyceModel);

  // Pull in Stratimikos to create linear solver
  Stratimikos::DefaultLinearSolverBuilder lowsfCreator;
  {
    RCP<Teuchos::ParameterList> sPL = Teuchos::parameterList();
    //sPL->set("Linear Solver Type","Amesos");
    //sPL->set("Preconditioner Type","None");
    lowsfCreator.setParameterList(sPL);
  }

  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > 
    W_factory = lowsfCreator.createLinearSolveStrategy("");
  TEST_ASSERT( Teuchos::nonnull(W_factory) );
  
  // Convert EpetraExt::ModelEvaluator to Thyra::ModelEvaluator
  RCP<Thyra::ModelEvaluator<double> > tmodel = 
    rcp(new Thyra::EpetraModelEvaluator(emodel,W_factory));
  
  // Create nonlinear solver
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver = 
    Rythmos::timeStepNonlinearSolver<double>();
  
  // Create Rythmos stepper
  RCP<Rythmos::BackwardEulerStepper<double> > stepper = 
    Rythmos::backwardEulerStepper<double>(tmodel,nlSolver);
  //stepper->setVerbLevel(Teuchos::VERB_EXTREME);
  
  // Set initial condition on stepper
  // Convert Epetra_Vector to Thyra::VectorBase
  // Set up initial conditions for Stepper:
  RCP<Epetra_Vector> xbar = rcp(new Epetra_Vector(*(emodel->get_x_map()))); 
  RCP<Epetra_Vector> xbardot = rcp(new Epetra_Vector(*(emodel->get_x_map()))); 
  
  // The initial conditions should be consistent in the following sense
  // explained below:
  // Note: var_dot == d(var)/dt
  // The inplicit-DAE is F(xbar_dot, xbar, t) = 0 
  // The DAE which Xyce uses is q_dot(x,t) + f(x,t) = 0;
  // So in implicit form x_bar = [x ; z]; 
  //                                  ^
  //      _                       _   z is the new variable such that q(x) - z = 0;
  // F = | F(1,1) = z_dot + f(x,t) |
  //     | F(2,1) = q(x,t) - z     |
  //      --                      -- 
  // Now, in the code we set initial condition on x_bar by setting xb 
  // therefore xb(3:5) should be such that it q(xb(0:2),t) = xb(3:5)
  // i.e. it satisfies F(2,1) and similarly for F(1,1) 
  {
    Epetra_Vector& xb = *xbar;
    xb[0] = 4.0; // v_2
    xb[1] = 4.0e-6; // i_vsrc
    xb[2] = 8.0; // v_1
    xb[3] = 4.0e-07;
    xb[4] = 0.0;
    xb[5] = 0.0;
    Epetra_Vector& xbdot = *xbardot;
    xbdot[0] = 0.0;
    xbdot[1] = 0.0;
    xbdot[2] = 0.0;
    xbdot[3] = 4.0e-6;
    xbdot[4] = -8;
    xbdot[5] = -8.0e-6;
    // the initial condition below actually work; they are not consistent
    //xbdot[3] = -4.0;
    //xbdot[4] = -4.0e-6;
    //xbdot[5] = -8.0;
  }
  // NOT REQUIRED NOW
  //RCP<Epetra_Vector> sbar = rcp(new Epetra_Vector(*(emodel->get_p_map(0)))); 
  //RCP<Epetra_Vector> sbardot = rcp(new Epetra_Vector(*(emodel->get_p_map(1)))); 
  //sbar->PutScalar(0.0);
  //sbardot->PutScalar(0.0);
  
  RCP<Thyra::VectorBase<double> > tx_init = // View!
    Thyra::create_Vector( xbar, tmodel->get_x_space() );
  RCP<Thyra::VectorBase<double> > txdot_init =  // View!
    Thyra::create_Vector( xbardot, tmodel->get_x_space() );

  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<double> model_ic = tmodel->createInArgs();
  model_ic.set_t(0.0);
  model_ic.set_x(tx_init);
  model_ic.set_x_dot(txdot_init);

  // WE DON'T NEED THEM ANYMORE 
  //model_ic.set_p(0,ts_init); // Constant for RC problem
  //model_ic.set_p(1,tsdot_init); // Constant for RC problem
  stepper->setInitialCondition(model_ic);
  // Take steps
  int N = 10;
  double dt = 0.1;
  double t_final = N*dt;
  for (int i=0 ; i<N ; ++i) {

    double dt_taken = stepper->takeStep(dt,Rythmos::STEP_TYPE_FIXED);
    TEST_ASSERT( dt_taken == dt );
    
    // Pull intermediate values out
    RCP<const Thyra::VectorBase<double> > x_intermediate = 
      Thyra::createMember(tmodel->get_x_space());
    x_intermediate = Rythmos::get_x(*stepper,dt +i*dt);
    // Convert Thrya::VectorBase to Epetra_Vector
    RCP<const Epetra_Vector> ex_intermediate = 
      Thyra::get_Epetra_Vector(*(emodel->get_x_map()),x_intermediate);
    //std::cout << "Intermediate Solution at final time t = " << dt + i*dt <<std::endl;
    //ex_intermediate->Print(std::cout);
    
    // Compare to exact solution
		
    // Here we compute the analytical solution for the BE 
		// and compare it with the one obtained from Rythmos 
		// The analytical solution is as below:
		//  omega = 1/ (C/Delta_t + 1/R)
		// | V2_{n+1} | = omega* [ C*V2_{n}/Delta_t + Vsrc_{n+1}/R]
		// | Isrc_{n+1} | = omega* ([ C*V2_{n}/Delta_t + Vsrc_{n+1}/R]*1/R - 
		// [C/Delta_t + 1/R]*Vsrc_{n+1}/R)
    double tol = 1.0e-6;
	  { 
		  double vsrc = 12.0*sin(2*pi*sinfreq*(i+1)*dt);
		  omega =  1/ (C/dt + 1/R);  // 1/Det(BE_LHS)
		  // during the first iteration get IC from the user.
		  if(!getInitICFlag) {
		  	std::cout << "-----------------------------------------------" << std::endl; 
		  	std::cout << "Compare the solution with analytical solution" << std::endl;
		  	V2_n = (*xbar)[0]; 
		  	Isrc_n = (*xbar)[1]; 
		  	V1_n = (*xbar)[2]; 
		  	getInitICFlag = 1;
		  }
		  else {
		  	V2_n = V2_next; 
		  }
		  V2_next = omega* (C*V2_n/dt + vsrc/R);
		  Isrc_next = omega* (( C*V2_n/dt + vsrc/R)*1/R - (C/dt + 1/R)*vsrc/R);
 		  //std::cout <<"t = "<< (i+1)*dt <<"\t" ;
		  //std::cout << "V2_next = " <<"\t"<< V2_next; //V2 
		  //std::cout << "\t" << "Isrc_next = " <<"\t"<< Isrc_next;  //Iscr
		  //std::cout << "\t" << "vsrc = " <<"\t"<< vsrc<<std::endl ; //V1
		  TEST_FLOATING_EQUALITY( (*ex_intermediate)[0], V2_next, tol ); // v2
  	  TEST_FLOATING_EQUALITY( (*ex_intermediate)[1], Isrc_next, tol ); //I_src
  	  //TEST_FLOATING_EQUALITY( (*ex_intermediate)[2], vsrc, tol ); //I_src
	  }
  }
  
  // Pull final answer out
  RCP<const Thyra::VectorBase<double> > x_final = Thyra::createMember(tmodel->get_x_space());
  x_final = Rythmos::get_x(*stepper,t_final);
  
  // Convert Thrya::VectorBase to Epetra_Vector
  RCP<const Epetra_Vector> ex_final = Thyra::get_Epetra_Vector(*(emodel->get_x_map()),x_final);
  //std::cout << "Final Solution at final time = " << std::endl;
  //ex_final->Print(std::cout);
  //TEST_ASSERT(false);
}
#endif // Xyce_TRILINOS_DEV_RYTHMOS


#ifdef Xyce_TRILINOS_DEV_RYTHMOS
TEUCHOS_UNIT_TEST( N_ANP_ModelEvaluator, rythmos_IRK_RCLoad ) {
  int iargs = 2;
  char *cargs[iargs];
  cargs[0] = "Xyce";
  cargs[1] = "modelEvaluatorTest.cir";
  // Use simple RC circuit netlist as input
  RCP<N_ANP_ModelEvaluator> xyceModel = rcp(new N_ANP_ModelEvaluator());
  xyceModel->initialize(iargs,cargs);
  	
  RCP<N_ANP_ModelEvaluator_Stateless> emodel = 
			N_ANP_modelEvaluator_Stateless(xyceModel);
	
  // Pull in Stratimikos to create linear solver
  Stratimikos::DefaultLinearSolverBuilder lowsfCreator;
  {
    RCP<Teuchos::ParameterList> sPL = Teuchos::parameterList();
    //sPL->set("Linear Solver Type","Amesos");
    //sPL->set("Preconditioner Type","None");
    lowsfCreator.setParameterList(sPL);
  }

  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > W_factory; 
  {
    //W_factory = lowsfCreator.createLinearSolveStrategy("");
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    RCP<Teuchos::ParameterList> stratPL = Teuchos::parameterList();
    stratPL->set("Linear Solver Type","AztecOO");
    stratPL->set("Preconditioner Type","None");
    linearSolverBuilder.setParameterList(stratPL);
    W_factory = Thyra::createLinearSolveStrategy(linearSolverBuilder);
  }
  TEST_ASSERT( Teuchos::nonnull(W_factory) );
  
  // Convert EpetraExt::ModelEvaluator to Thyra::ModelEvaluator
  RCP<Thyra::ModelEvaluator<double> > tmodel = 
    rcp(new Thyra::EpetraModelEvaluator(emodel,W_factory));
  
  // Create nonlinear solver
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver = 
    Rythmos::timeStepNonlinearSolver<double>();
  
  // Create Rythmos stepper
  //RCP<Rythmos::BackwardEulerStepper<double> > stepper = 
  //  Rythmos::backwardEulerStepper<double>(tmodel,nlSolver);
  //stepper->setVerbLevel(Teuchos::VERB_EXTREME);
  
  //RCP<Rythmos::RKButcherTableauBase<double> > rkbt = Rythmos::createRKBT<double>("Backward Euler");
  RCP<Rythmos::RKButcherTableauBase<double> > rkbt = Rythmos::createRKBT<double>("Implicit 1 Stage 2nd order Gauss");
  //RCP<Rythmos::RKButcherTableauBase<double> > rkbt = Rythmos::createRKBT<double>("Implicit 1 Stage 1st order Radau left");
  //RCP<Rythmos::RKButcherTableauBase<double> > rkbt = Rythmos::createRKBT<double>("Implicit 1 Stage 1st order Radau right");
  //RCP<Rythmos::RKButcherTableauBase<double> > rkbt = Rythmos::createRKBT<double>("Diagonal IRK 2 Stage 3rd order");  // not yet
  //RCP<Rythmos::RKButcherTableauBase<double> > rkbt = Rythmos::createRKBT<double>("Implicit 2 Stage 4th order Gauss"); // not yet
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory; 
  {
    //irk_W_factory = lowsfCreator.createLinearSolveStrategy("");
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    RCP<Teuchos::ParameterList> stratPL = Teuchos::parameterList();
    stratPL->set("Linear Solver Type","AztecOO");
    stratPL->set("Preconditioner Type","None");
    linearSolverBuilder.setParameterList(stratPL);
    irk_W_factory = Thyra::createLinearSolveStrategy(linearSolverBuilder);
  }
  RCP<Rythmos::ImplicitRKStepper<double> > stepper = 
    Rythmos::implicitRKStepper<double>(tmodel,nlSolver,irk_W_factory,rkbt);
  
  // Set initial condition on stepper
  //   Convert Epetra_Vector to Thyra::VectorBase
  // Set up initial conditions for Stepper:
  RCP<Epetra_Vector> xbar = rcp(new Epetra_Vector(*(emodel->get_x_map()))); 
  RCP<Epetra_Vector> xbardot = rcp(new Epetra_Vector(*(emodel->get_x_map()))); 
  
  // The initial conditions should be consistent in the following sense
  // explained below:
  // Note: var_dot == d(var)/dt
  // The inplicit-DAE is F(xbar_dot, xbar, t) = 0 
  // The DAE which Xyce uses is q_dot(x,t) + f(x,t) = 0;
  // So in implicit form x_bar = [x ; z]; 
  //                                  ^
  //      _                       _   z is the new variable such that q(x) - z = 0;
  // F = | F(1,1) = z_dot + f(x,t) |
  //     | F(2,1) = q(x,t) - z     |
  //      --                      -- 
  // Now, in the code we set initial condition on x_bar by setting xb 
  // therefore xb(3:5) should be such that it q(xb(0:2),t) = xb(3:5)
  // i.e. it satisfies F(2,1) and similarly for F(1,1) 
  
  {
    Epetra_Vector& xb = *xbar;
    xb[0] = 4.0; // v_2
    xb[1] = 4.0e-6; // i_vsrc
    xb[2] = 8.0; // v_1
    xb[3] = 4.0e-07;
    xb[4] = 0.0;
    xb[5] = 0.0;
    Epetra_Vector& xbdot = *xbardot;
    xbdot[0] = 0.0;
    xbdot[1] = 0.0;
    xbdot[2] = 0.0;    
    xbdot[3] = 4.0e-6;
    xbdot[4] = -8;
    xbdot[5] = -8.0e-6;
    // the initial condition below actually work; they are not consistent
    //xbdot[3] = -4.0;
    //xbdot[4] = -4.0e-6;
    //xbdot[5] = -8.0;
  }
  // NOT REQUIRED NOW
  //RCP<Epetra_Vector> sbar = rcp(new Epetra_Vector(*(emodel->get_p_map(0)))); 
  //RCP<Epetra_Vector> sbardot = rcp(new Epetra_Vector(*(emodel->get_p_map(1)))); 
  //sbar->PutScalar(0.0);
  //sbardot->PutScalar(0.0);
  
  RCP<Thyra::VectorBase<double> > tx_init = // View!
    Thyra::create_Vector( xbar, tmodel->get_x_space() );
  RCP<Thyra::VectorBase<double> > txdot_init =  // View!
    Thyra::create_Vector( xbardot, tmodel->get_x_space() );

  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<double> model_ic = tmodel->createInArgs();
  model_ic.set_t(0.0);
  model_ic.set_x(tx_init);
  model_ic.set_x_dot(txdot_init);

  // WE DON'T NEED THEM ANYMORE 
  //model_ic.set_p(0,ts_init); // Constant for RC problem
  //model_ic.set_p(1,tsdot_init); // Constant for RC problem
  stepper->setInitialCondition(model_ic);
  // Take steps
  int N = 20;
  double dt = 0.05;
  double t_final = N*dt;
  for (int i=0 ; i<N ; ++i) {
    double dt_taken = stepper->takeStep(dt,Rythmos::STEP_TYPE_FIXED);
    TEST_ASSERT( dt_taken == dt );
    // Pull intermediate values out
    RCP<const Thyra::VectorBase<double> > x_intermediate = 
      Thyra::createMember(tmodel->get_x_space());
    x_intermediate = Rythmos::get_x(*stepper,dt +i*dt);
    // Convert Thrya::VectorBase to Epetra_Vector
    RCP<const Epetra_Vector> ex_intermediate = 
      Thyra::get_Epetra_Vector(*(emodel->get_x_map()),x_intermediate);
    //std::cout << "Intermediate Solution at final time t = " << dt + i*dt <<std::endl;
    //ex_intermediate->Print(std::cout);
  }
  
  // Pull final answer out
  RCP<const Thyra::VectorBase<double> > x_final = Thyra::createMember(tmodel->get_x_space());
  x_final = Rythmos::get_x(*stepper,t_final);
  
  // Convert Thrya::VectorBase to Epetra_Vector
  RCP<const Epetra_Vector> ex_final = Thyra::get_Epetra_Vector(*(emodel->get_x_map()),x_final);
  
  // Compare to exact solution
  //std::cout << "Final Solution at final time = " << std::endl;
  
  //ex_final->Print(std::cout);
  TEST_ASSERT(false);
}
#endif // Xyce_TRILINOS_DEV_RYTHMOS

#endif // Xyce_TRILINOS_DEV
