======================================================================================
======================================================================================

  ---------------------------
  Verilog-A definition of PSP
  ---------------------------


  (c) Copyright 2009, All Rights Reserved, NXP Semiconductors


  Version: PSP 103.1 (including JUNCAP2 200.3.3), May 2009

  This version of PSP is contained in SiMKit 3.3

======================================================================================
======================================================================================

 Authors: G.D.J. Smit, A.J. Scholten, and D.B.M. Klaassen (NXP Semiconductors Research)
          R. van Langevelde (Philips Research)
          G. Gildenblat, X. Li, and W. Wu (Arizona State University)



The most recent version of the model code, the documentation, and contact information
can be found on:

         http://pspmodel.asu.edu/

======================================================================================
======================================================================================

This package consists of several files:

     - readme.txt                     This file

     - psp103.va                      Main file for PSP model
     - psp103_nqs.va                  Main file for PSP with NQS-effects
     - juncap200.va                   Main file for JUNCAP2 stand-alone model

     - Common103_macrodefs.include    Common macro definitions
     - PSP103_macrodefs.include       Macro definitions for PSP
     - PSP103_module.include          Actual model code for intrinsic MOS model
     - PSP103_SPCalculation.include   Surface potential and related calculations
     - PSP103_binning.include         Geometry scaling equation for binning
     - PSP103_binpars.include         Parameterlist for global PSP binning model
     - PSP103_nqs_macrodefs.include   Macro definitions for PSP-NQS
     - PSP103_InitNQS.include         PSP-NQS initialization code
     - PSP103_ChargesNQS.include      Calculation of NQS-charge contributions
     - JUNCAP200_macrodefs.include    Macro definitions for JUNCAP2 model
     - JUNCAP200_parlist.include      JUNCAP2 parameter list
     - JUNCAP200_varlist1.include     JUNCAP2 variable declarations
     - JUNCAP200_varlist2.include     JUNCAP2 variable declarations
     - JUNCAP200_InitModel.include    JUNCAP2 model initialization code

======================================================================================
======================================================================================

Usage
-----

Depending which model one wants to use, one should compile one of the three .va-files
(psp103.va, psp103_nqs.va, and juncap200.va). The module names are "PSP103VA" and "PSPNQS103VA"
(for QS and NQS, respectively), and "JUNCAP200" for the JUNCAP2-model.


======================================================================================
======================================================================================

Release notes vA-code of PSP 103.1 (May 2009)
---------------------------------------------

Changes:

- Added external sheet resistance RSHD for drain diffusion (used when SWJUNASYM=1)
- Extended NUD-model to allow for retrograde profiles (GFACNUD > 1)
- Added noise source labeling (vA-code only)
- Added value of gate resistance to OP-output
- Bugfix and minor implementation change in NUD-model
- Minor bug fix in conditional for SP-calculation of overlap areas.

PSP 103.1 is backwards compatible with the previous version, PSP 103.0.

Remark
------
PSP 103.1 verilog-A code is LRM 2.2-compliant. Yet, the PSP-team has identified some PSP
compilation issues (mostly related to the ddx()-operator) in a number of vA-compilers. The
EDA-vendors in question have been contacted. At this moment, they have already solved most of
the underlying limitations/bugs in the compilers and the fixes are available in the most recent
or near-future updates of their products.
Details will be provided on request by PSP-team.

Comments and code-adaptations related to these compilation issues are still present in this
release of the PSP-code, but will be removed in future updates.

======================================================================================
======================================================================================

The authors want to thank Laurent Lemaitre and Colin McAndrew (Freescale)
for their help with ADMS and the implementation of the model code. Geoffrey
Coram (Analog Devices) is acknowledged for input concerning the Verilog-A
implementation of the model.
