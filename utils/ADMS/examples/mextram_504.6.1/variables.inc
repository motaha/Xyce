// Declaration of variables 

// Model constants

real An, Bn;

// Temperature scaling variables

real Tk, Trk, tN, tNa, Tki, Tamb;
real Vt, Vtr, VtINV, VtrINV, VdtINV;
real Vdt;

real UdeT, VDE_T, UdcT, VDC_T;
real CJE_T, CJC_T, XP_T;
real CJCscale, CJCscaleINV;

real RE_T, RBV_T, RBC_T, RCC_T, RCV_T;
real BF_T, BRI_T;

real IS_T, IK_T, IBF_T, IBR_T, VEF_T, VER_T;

real TAUE_T, TAUB_T, TEPI_T, TAUR_T;
real BnT, DEG_T, Tk300;

`ifdef SELFHEATING
real RTH_Tamb;
`endif

`ifdef SUBSTRATE
real UdsT, VDS_T, CJS_T, ISS_T, IKS_T;
`endif

// MULT - scaling variables

real invMULT;
real IS_TM, IK_TM, IBF_TM, IBR_TM, IHC_M;
real CJE_TM, CJC_TM;

real RE_TM, RBC_TM, RBV_TM, RCC_TM, RCV_TM, SCRCV_M;
real KF_M, KFN_M;

`ifdef SELFHEATING
real RTH_Tamb_M, CTH_M;
`endif

`ifdef SUBSTRATE
real ISS_TM, IKS_TM, CJS_TM;
`endif


// Epilayer model variables

real K0, Kw, pW, Ec, Ic1c2;
real Vqs_th, Vqs, Iqs;
real alpha, vyi, yi, xi_w, xi_w1;
real gp0, gp02, p0star, eVb2c2star;
real B1, B2, Vxi0, Vch, Icap, pav;

// Effective emitter and collector junction bias variables    

real Vfe, Vje, Vte;
real Vjunc, bjc, Vfc, Vjc, fI, Vcv, Vtc;

// Transfer current variables

real If0, f1, f2, n0, nB;
real q0I, q1I, qBI, Ir, If, In;

// Base and substrate current(s) variables    

real Xext1;
real Ib1, Ib1_s, Ib2, Ib3;
real Ibf0, Iex, Isub;
real g1, g2, pWex, nBex;
real Xg1, XnBex, XIMex, XIMsub, Vex, VBex, Fex, XIex;

`ifdef SUBSTRATE
real XIsub, Isf;
`endif

// Distributed base effects variables

real q0Q, q1Q, qBQ, Rb2, Ib1b2;
real dVteVb2e1, dVteVje, dVjeVb2e1;
real dQteVb2e1, dQbeVb2e1, dQeVb2e1;
real dn0Vb2e1;

// Weak-avalanche current variables

real dEdx0, xd, Weff, Wd, Eav, E0, Em, SHw, Efi, Ew;
real lambda, Gem, Gmax, Iavl;
real Icap_IHC;

`ifdef SELFHEATING
real Vb2c2star, power;
`endif

// Charges and capacitances variables

real Qte, Vje_s, Qte_s;
real Qtc;
real Qb0, Qbe, Qbc, Qb1b2;
real Qbe_qs, Qbc_qs;
real Vjcex, Vtexv, Qtex, XVjcex, XVtexv, XQtex;

`ifdef SUBSTRATE
real Vfs, Vjs, Qts;
`endif

real Qe0, Qe;
real Qepi0, Qepi, Xg2, XpWex, XQex;
real Qex;
real CBEO_M, CBCO_M;

// Biases and exponential terms variables

real Vb2c1, Vb2c2, Vb2e1, Vb1e1, Vb1b2, Vb1c1, Vc1c2;
real Vsc1, Vee1, Vbb1, Vbc1, Vcc1, Vbe, Vbc;
real eVb2c2, eVb2e1, eVb1e1, eVb1b2, eVb1c1, eVbc1, eVsc1;
real eVb1c1VDC, eVb2c2VDC, eVbc1VDC, eVb2c1VDC;

// Help variables

real dxa, sqr_arg;
real eps2, x2;
real alpha1, vdif, Ic1c2_Iqs, gp0_help;
real EmEav_Em, Vb2e1Vfe, termE, termC;
real Vex_bias;
real eps_VDC, a_VDE, a_VDC;
real xtmp;

real tmpExp, tmpV;


`ifdef  SUBSTRATE
real a_VDS;
`endif

// Noise variables 

real common; 
real powerREC, powerRBC, powerRCC, powerRBV; 
real powerCCS;  
real powerFBCS, powerFBC1fB1, exponentFBC1fB2, powerFBC1fB2; 
real powerEBSCS, powerEBSC1f; 
real powerRBCS, powerRBC1f; 
real powerExCS, powerExCSMOD, powerExC1f, powerExC1fMOD; 

`ifdef SUBSTRATE
real powerSubsCS_B1S, powerSubsCS_BS; 
`endif

real twoqIavl, powerII;

`ifdef __VAMS_COMPACT_MODELING__
   real w_n_1;
`else
   electrical w_n_1;
`endif


//
// Operation point (output) variables
//

// The external currents and current gain
`OPP(OP_ic,        A, External DC collector current)
`OPP(OP_ib,        A, External DC base current)
`OPP(OP_betadc,     , External DC current gain Ic/Ib)

// The internal biases
`OPP(OP_vb2e1,     V, Internal base-emitter bias)
`OPP(OP_vb2c2,     V, Internal base-collector bias)
`OPP(OP_vb2c1,     V, Internal base-collector bias including epilayer)
`OPP(OP_vb1c1,     V, External base-collector bias without contact resistances)
`OPP(OP_ve1e,      V, Bias over emitter resistance)

// The actual currents
`OPP(OP_in,        A,   Main current)
`OPP(OP_ic1c2,     A,   Epilayer current)
`OPP(OP_ib1b2,     A,   Pinched-base current)
`OPP(OP_ib1,       A,   Ideal forward base current)
`OPP(OP_isb1,      A,   Ideal side-wall base current)
`OPP(OP_ib2,       A,   Non-ideal forward base current)
`OPP(OP_ib3,       A,   Non-ideal reverse base current)
`OPP(OP_iavl,      A,   Avalanche current)
`OPP(OP_iex,       A,   Extrinsic reverse base current)

`OPP(OP_xiex,      A,   Extrinsic reverse base current)
`OPP(OP_isub,      A,   Substrate current)
`OPP(OP_xisub,     A,   Substrate current)
`OPP(OP_isf,       A,   Substrate failure current)
`OPP(OP_ire,       A,   Current through emitter resistance)
`OPP(OP_irbc,      A,   Current through constant base resistance)
`OPP(OP_ircc,      A,   Current through constant collector resistance)

//The actual charges
`OPP(OP_qe,        C,   Emitter charge or emitter neutral charge)
`OPP(OP_qte,       C,   Base-emitter depletion charge)
`OPP(OP_qste,      C,   Sidewall base-emitter depletion charge)
`OPP(OP_qbe,       C,   Base-emitter diffusion charge)
`OPP(OP_qbc,       C,   Base_collector diffusion charge)
`OPP(OP_qtc,       C,   Base-collector depletion charge)
`OPP(OP_qepi,      C,   Epilayer diffusion charge)
`OPP(OP_qb1b2,     C,   AC current crowding charge)
`OPP(OP_qtex,      C,   Extrinsic base-collector depletion charge)
`OPP(OP_xqtex,     C,   Extrinsic base-collector depletion charge)
`OPP(OP_qex,       C,   Extrinsic base-collector diffusion charge)
`OPP(OP_xqex,      C,   Extrinsic base-collector diffusion charge)
`OPP(OP_qts,       C,   Collector-substrate depletion charge)

//Small signal equivalent circuit conductances and resistances
`OPP(OP_gx,        S,   Forward transconductance)
`OPP(OP_gy,        S,   Reverse transconductance)
`OPP(OP_gz,        S,   Reverse transconductance)
`OPP(OP_gspi,      S,   Conductance sidewall b-e junction)
`OPP(OP_gpix,      S,   Conductance floor b-e junction)
`OPP(OP_gpiy,      S,   Early effect on recombination base current)
`OPP(OP_gpiz,      S,   Early effect on recombination base current)
`OPP(OP_gmux,      S,   Early effect on avalanche current limiting)
`OPP(OP_gmuy,      S,   Conductance of avalanche current)
`OPP(OP_gmuz,      S,   Conductance of avalanche current)
`OPP(OP_gmuex,     S,   Conductance of extrinsic b-c junction)
`OPP(OP_xgmuex,    S,   Conductance of extrinsic b-c junction)
`OPP(OP_grcvy,     S,   Conductance of epilayer current)
`OPP(OP_grcvz,     S,   Conductance of epilayer current)
`OPP(OP_rbv,       Ohm, Base resistance)
`OPP(OP_grbvx,     S,   Early effect on base resistance)
`OPP(OP_grbvy,     S,   Early effect on base resistance)
`OPP(OP_grbvz,     S,   Early effect on base resistance)
`OPP(OP_re,        Ohm, Emitter resistance)
`OPP(OP_rbc,       Ohm, Constant base resistance)
`OPP(OP_rcc,       Ohm, Constant collector resistance)
`OPP(OP_gs,        S,   Conductance parasistic PNP transistor)
`OPP(OP_xgs,       S,   Conductance parasistic PNP transistor)
`OPP(OP_gsf,       S,   Conductance substrate failure current)

//Small signal equivalent circuit capacitances
`OPP(OP_csbe,      F,   Capacitance sidewall b-e junction)
`OPP(OP_cbex,      F,   Capacitance floor b-e junction)
`OPP(OP_cbey,      F,   Early effect on b-e diffusion charge)
`OPP(OP_cbez,      F,   Early effect on b-e diffusion charge)
`OPP(OP_cbcx,      F,   Early effect on b-c diffusion charge)
`OPP(OP_cbcy,      F,   Capacitance floor b-c junction)
`OPP(OP_cbcz,      F,   Capacitance floor b-c junction)
`OPP(OP_cbcex,     F,   Capacitance extrinsic b-c junction)
`OPP(OP_xcbcex,    F,   Capacitance extrinsic b-c junction)
`OPP(OP_cb1b2,     F,   Capacitance AC current crowding)
`OPP(OP_cb1b2x,    F,   Cross-capacitance AC current crowding)
`OPP(OP_cb1b2y,    F,   Cross-capacitance AC current crowding)
`OPP(OP_cb1b2z,    F,   Cross-capacitance AC current crowding)
`OPP(OP_cts,       F,   Capacitance s-c junction)

//Approximate small signal equivalent circuit
`OPP(OP_gm,        S,transconductance)
`OPP(OP_beta,       ,   Current amplification)
`OPP(OP_gout,      S,   Output conductance)
`OPP(OP_gmu,       S,   Feedback transconductance)
`OPP(OP_rb,        Ohm, Base resistance)
`OPP(OP_cbe,       C,   Base-emitter capacitance)
`OPP(OP_cbc,       C,   Base-collector capacitance)

//quantities to describe internal state of the model
`OPP(OP_ft,         ,   Good approximation for cut-off frequency)
`OPP(OP_iqs,       A,   Current at onset of quasi-saturation)
`OPP(OP_xiwepi,    m,   Thickness of injection layer)
`OPP(OP_vb2c2star, V, Physical value of internal base-collector bias)

//self-heating
`OPP(OP_pdiss,     W, Dissipation)
`OPP(OP_tk,        K, Actual temperature)

//help variables
real dydx, dydz, gpi;
real gammax, gammay, gammaz, gbfx, gbfy, gbfz, alpha_ft;
real rx, ry, rz, rb1b2, rex, xrex, taut;
real ddxInVb2, ddxIb1Vb2, ddxIavlVb2, ddxIc1c2Vb2, ddxIb1b2Vb2;
real ddxQbeVb2, ddxQcciVb2, ddxQb1b2Vb1;
