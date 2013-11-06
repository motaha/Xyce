///////////////////////////////////////////////////////////////////////////////
// Experimental MOS1 continuation work.
// Eric Keiter 9233
// 8/29/03.
//
// The purpose of this program is to test out the application of 2 homotopy
// parameters (alpha1 and alpha2) to the level 1 MOSFET equations.  This
// general idea comes from a preprint by Jaijeet Roychowdhury (and 
// verbal discussion).  In his paper, he calls his homotopy parameters
// lambda1 and lambda2.  Unfortunately, there is another parameter called
// lambda in the MOSFET level 1 model, so I'm calling these parameters
// alpha1 and alpha2 instead.
//
// In a homotopy simulation, alpha1 is varried from 0 to 1, and 
// effects the shape of the output 
// characteristic.  The output characteristic is the behavior of
// Ids with respect to Vds, for a fixed Vgs.   There should be two regions
// in an output characteristic curve.  The linear region, for Vds < Vgst,
// and the saturation region, for Vds >= Vgst.  (Vgst = Vgs - Von).
//
// Normally (alpha1=1) , the linear region is very steep, and the
// saturation region is horizontal on the plot.  This makes the MOSFET
// highly nonlinear.  alpha1 is used so that when alpha1=0, this
// nonlinearity is mitigated.  At alpha1=0, the device still saturates, but
// at a much higher value of Vds, and the linear region slope is much more
// gentle.   The approach I used in this code was to scale the apparent
// Vds downward, so that from the perspective of the drain current
// calculation, it would appear that Vds was much smaller than it was.
// This appears to have the correct effect.
//
// The other parameter, alpha2, is also varried from 0 to 1, and affects the
// gain of the device.  At alpha2=0, the device has no gain, and acts as
// though Vgs is constant.  At alpha2=1, the device has the gain of a 
// regular MOSFET, and Ids (drain-source current) is linearly dependent
// upon Vgs.  The curve showing Ids as a function of Vgs is known as the
// "transfer characteristic",  and for alpha2=0, is just a horizontal line.
// In Jaijeets paper, he appears to choose the value of Ids that
// corresponds to a Vgs of 3.0, for the alpha2=1 case.
//
// Since both alpha1 and alpha2 are applied to the input variables (vds and
// vgs), rather than directly to the drain current equations, it is hoped
// that these adjustments  can also be applied for other, more complex
// MOSFET models, which probably have much more complicated expressions for
// calculating cdrain, etc.
// 
// Ironically, it turns out that for a real homotopy simulation, which
// actually works, it is important to vary alpha2 from 0 to 1 first, and
// then vary alpha1 from 0 to 1 second.
//


#include <iostream>
#include <cstdio>
#include <math.h>

// the default for VON in Xyce/Spice is 2.0.  Jaijeet's paper seems to use
// 0.0.
#define VON 0.0

using namespace std;

bool idrain (double vds, double vgs,  // inputs
   double & cdrain1, double & gm1, double & gds1, double & gmbs1) // outputs
{
  double mode = 1.0;
  if (vds >= 0)
  {
    mode = 1.0;
  }
  else
  {
    mode = -1.0;
  }
  double cdrain, gm, gds, gmbs, betap;

  double dtype = 1.0; // I think this is n-type.

  double Von=VON; // this is typical.  vbi = vt0 = 2.0, if room temp.

  double vgst=vgs-Von;

  // these are all just spice default for constants.  Assuming room temp.
  double oxideThickness =  1.0e-7;
  double oxideCapFactor = 3.9 * 8.854214871e-12/oxideThickness;
  double surfaceMobility= 600.0;
  double transconductance = surfaceMobility * oxideCapFactor * 1e-4;
  double w = 1.0e-4;  // defw
  double l = 1.0e-4;  // defl
  double Beta = transconductance * w/l;

  double lambda = 0.0; // channel length modulation (not continuation)

  double arg = 0.0;  // check this.  What is this?


  if (vgst <= 0) 
  {
    //
    //      cutoff region
    //
    cdrain=0;
    gm=0;
    gds=0;
    gmbs=0;
  } 
  else
  {
    //
    //     saturation region
    //
    betap=Beta*(1.0 + lambda*(vds*mode));
    if (vgst <= (vds*mode))
    {
      cdrain=betap*vgst*vgst*0.5;
      gm=betap*vgst;
      gds=lambda*Beta*vgst*vgst*0.5;
      gmbs=gm*arg;
    } 
    else 
    {
      //
      //     linear region
      //
      cdrain=betap*(vds*mode)* (vgst-0.5*(vds*mode));

      gm=betap*(vds*mode);

      gds=betap*(vgst-(vds*mode))+ 
	lambda*Beta* (vds*mode)* (vgst-0.5*(vds*mode));

      gmbs=gm*arg;
    }
  }

  cdrain1 = cdrain;
  gm1 = gm;
  gds1 = gds;
  gmbs1 = gmbs;

  return true;
}

int main (int arg, char cargs[])
{
  double cdrain, cdrain_test;
  double gm, gm_test;
  double gds, gds_test;
  double gmbs, gmbs_test;
  double vds, vds_test;
  double vgs, vgs_test;

  //if (arg < 2)
  //{
    //printf("Sorry, not enough args.  Need at least 2\n");
    //exit(0);
  //}

  double alpha1 = 0.0; // continaution alpha.  (nonlinearity parameter)
  double alpha2 = 1.0; // continuation alpha.  (gain parameter)

  fprintf(stdout,"Input alpha1 (nonlinearity param): ");
  fscanf(stdin, "%lf", &alpha1);

  fprintf(stdout,"Input alpha2 (gain param): ");
  fscanf(stdin, "%lf", &alpha2);

  // This number, min, is used in the scaling of Vds.  The point is to not
  // let the "effective" Vds ever be zero.  At alpha1 = 0, 
  double min = 0.3;

  int num_vds = 31;
  int num_vgs = 31;

  double vgs_homotopy_const = VON + 3.0; // to match jaijeet.

  double vgs_min =  -1.0 + VON;  // for Jaijeet, VON=0.0
  double vgs_max =  +3.0 + VON;
  double vds_min =  +0.0;
  double vds_max = +10.0;

  double delta_vds = (vds_max-vds_min)/(static_cast<double>(num_vds-1));
  double delta_vgs = (vgs_max-vgs_min)/(static_cast<double>(num_vgs-1));

  FILE *fp1 = fopen("idrain.dat","w");

  fprintf(fp1, " TITLE = \"MOS1 Continuation: alpha1=%f alpha2=%f\",\n", 
    alpha1, alpha2);
  fprintf(fp1,"\tVARIABLES = \"Vds  \",\"Vgs \",\n");
  fprintf(fp1,"\t    \"cdrain \",\n");
  fprintf(fp1,"\t    \"gm     \",\n");
  fprintf(fp1,"\t    \"gds    \",\n");
  fprintf(fp1,"\t    \"gmbs   \",\n");

  fprintf(fp1,"\tZONE I=%d J=%d F=POINT\n", num_vds, num_vgs);

  for (int ivds=0;ivds<num_vds;++ivds)
  {

    vds = (static_cast<double>(ivds))*delta_vds + vds_min;
    vds_test = vds * (alpha1*(1.0-min) + min);

    for (int ivgs=0;ivgs<num_vgs;++ivgs)
    {
      vgs = (static_cast<double>(ivgs))*delta_vgs + vgs_min;
      vgs_test = (alpha2)*vgs + (1.0-alpha2)*vgs_homotopy_const;

      idrain(vds, vgs, cdrain, gm, gds, gmbs);

      idrain(vds_test, vgs_test, cdrain_test,  gm_test, gds_test, gmbs_test);

      fprintf(fp1,"%12.4e ", vds);
      fprintf(fp1,"%12.4e %12.4e %12.4e %12.4e %12.4e", 
	  vgs, cdrain_test, gm_test, gds_test, gmbs_test);
      fprintf(fp1,"\n");

    } // end  of for loop.
  }

  fclose(fp1);
}

