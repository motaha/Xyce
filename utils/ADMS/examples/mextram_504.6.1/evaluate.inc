// Evaluate model equations 

begin  // Currents and sharges
// Nodal biases

   Vb2c1 = TYPE * V(b2, c1);
   Vb2c2 = TYPE * V(b2, c2);
   Vb2e1 = TYPE * V(b2, e1);
   Vb1e1 = TYPE * V(b1, e1);
   Vb1b2 = TYPE * V(b1, b2);
   Vb1c1 = TYPE * V(b1, c1);
   Vbc1  = TYPE * V(b,  c1);
`ifdef SUBSTRATE
   Vsc1  = TYPE * V(s,  c1);
`endif
   Vc1c2 = TYPE * V(c1, c2);
   Vee1  = TYPE * V(e,  e1);
   Vbb1  = TYPE * V(b,  b1);
   Vcc1  = TYPE * V(c,  c1);
   Vbe   = TYPE * V(b,  e);
   Vbc   = TYPE * V(b,  c);

// Exponential bias terms

   `expLin(eVb2c2, Vb2c2 * VtINV)
   `expLin(eVb2e1, Vb2e1 * VtINV)
   `expLin(eVb1e1, Vb1e1 * VtINV)
   `expLin(eVb1c1, Vb1c1 * VtINV)
   `expLin(eVb1b2, Vb1b2 * VtINV)
   `expLin(eVbc1,  Vbc1  * VtINV)
`ifdef SUBSTRATE
   `expLin(eVsc1,  Vsc1  * VtINV)
`endif

   `expLin(eVbc1VDC, (Vbc1  - VDC_T) * VtINV)
   `expLin(eVb1c1VDC,(Vb1c1 - VDC_T) * VtINV)
   `expLin(eVb2c2VDC,(Vb2c2 - VDC_T) * VtINV)
   `expLin(eVb2c1VDC,(Vb2c1 - VDC_T) * VtINV)

// Governing equations

   // Epilayer model

   K0 = sqrt(1.0 + 4.0 * eVb2c2VDC);
   Kw = sqrt(1.0 + 4.0 * eVb2c1VDC);
   pW = 2.0 *  eVb2c1VDC / (1.0 + Kw);
   if (pW < `TEN_M40) pW = 0;
   Ec = Vt * (K0 - Kw - ln((K0 + 1.0) / (Kw + 1.0)) );
   Ic1c2 =  (Ec + Vc1c2) / RCV_TM;

   if (Ic1c2 > 0.0) begin

     `linLog(tmpV,Vb2c1,100.0);
     Vqs_th = VDC_T + 2.0 * Vt *
              ln(0.5 * Ic1c2 * RCV_TM * VtINV + 1.0) - tmpV;
     eps_VDC = 0.2 * VDC_T;
     `max_hyp0(Vqs, Vqs_th, eps_VDC);  
     Iqs = Vqs * (Vqs + IHC_M * SCRCV_M) / (SCRCV_M * (Vqs + IHC_M * RCV_TM));

     Ic1c2_Iqs = Ic1c2 / Iqs;
     `max_logexp(alpha1, Ic1c2_Iqs, 1.0, AXI);
     alpha = alpha1 / (1.0 + AXI * ln(1.0 + exp(-1.0 / AXI)));
     vyi = Vqs / (IHC_M * SCRCV_M);
     yi = (1.0 + sqrt(1.0 + 4.0 * alpha * vyi * (1.0 + vyi))) /
          (2.0 * alpha * (1.0 + vyi));

     xi_w = 1.0 - yi / (1.0 + pW * yi);
     gp0 = 0.5 * Ic1c2 * RCV_TM * xi_w * VtINV;

     gp0_help = 2.0 * gp0 + pW * (pW + gp0 + 1.0);
     gp02 = 0.5 * (gp0 - 1.0);
     sqr_arg = gp02 * gp02 + gp0_help;
     if (gp0 >= 1.0)
        p0star =  gp02 + sqrt(sqr_arg);
     else
        p0star = gp0_help / (sqrt(sqr_arg) - gp02);
     if (p0star < `TEN_M40) p0star = 0.0;


     eVb2c2star = p0star * (p0star + 1.0) * exp(VDC_T * VtINV);
     B1 = 0.5 * SCRCV_M * (Ic1c2 - IHC_M); 
     B2 = SCRCV_M * RCV_TM * IHC_M * Ic1c2;
     Vxi0 = B1 + sqrt(B1 * B1 + B2);
     Vch = VDC_T * (0.1 + 2.0 * Ic1c2 / (Ic1c2 + Iqs));
     Icap = IHC_M * Ic1c2 / (IHC_M + Ic1c2);
     Icap_IHC = IHC_M / (IHC_M + Ic1c2);

   end else begin

     p0star = 2.0 * eVb2c2VDC / (1.0 + K0);
     eVb2c2star = eVb2c2;
     if ((abs(Vc1c2) < 1.0e-5 * Vt) ||
         (abs(Ec) < `TEN_M40 * Vt * (K0 + Kw)))
        begin 
          pav = 0.5 * (p0star + pW);
          xi_w = pav / (pav + 1.0);
        end

     else
        begin 
          xi_w = Ec / (Ec + Vb2c2 - Vb2c1);
        end

     Vxi0 = Vc1c2;
     Vch = 0.1 * VDC_T;
     Icap = Ic1c2;
     Icap_IHC = 1.0 - Icap / IHC_M;

   end

   // Effective emitter junction capacitance bias

   Vfe = VDE_T * (1.0 - pow(`AJE , -1.0 / PE));
   a_VDE = 0.1 * VDE_T;
   `min_logexp(Vje, Vb2e1, Vfe, a_VDE);
   Vte = VDE_T / (1.0 - PE) * (1.0 - pow(1.0 - Vje / VDE_T, 1.0 - PE)) +
        `AJE * (Vb2e1 - Vje);

   // Effective collector junction capacitance bias

   Vjunc = Vb2c1 + Vxi0;
   bjc = (`AJC - XP_T) / (1.0 - XP_T);
   Vfc = VDC_T * (1.0 - pow(bjc, -1.0 / PC));
   `min_logexp(Vjc, Vjunc, Vfc, Vch);
   fI = pow(Icap_IHC, MC);
   Vcv = VDC_T / (1.0 - PC) * (1.0 - fI * pow(1.0 - Vjc / VDC_T, 1.0 - PC)) +
         fI * bjc * (Vjunc - Vjc);
   Vtc = (1.0 - XP_T) * Vcv + XP_T * Vb2c1;

   // Transfer current

   If0 = 4.0 * IS_TM / IK_TM;
   f1 =  If0 * eVb2e1;
   n0 =  f1 / (1.0 + sqrt(1.0 + f1));
   f2 =  If0 * eVb2c2star;
   nB =  f2 / (1.0 + sqrt(1.0 + f2));

   if (DEG == 0.0)
        q0I = 1.0 + Vte / VER_T + Vtc / VEF_T; 
   else
      begin
        termE = (Vte / VER_T + 1.0) * DEG_T * VtINV;
        termC = -Vtc / VEF_T * DEG_T * VtINV;
        q0I = (exp(termE) - exp(termC)) /
              (exp(DEG_T * VtINV) - 1.0);
      end

   `max_hyp0(q1I, q0I, 0.1);
   qBI = q1I * (1.0 + 0.5 * (n0 + nB));

   Ir = IS_TM *  eVb2c2star;
   If = IS_TM * eVb2e1;
   In = (If - Ir) / qBI;

   // Base and substrate current(s)

   Ibf0 = IS_TM / BF_T;
   if (XREC == 0.0)
      Ib1 = (1.0 - XIBI) * Ibf0 * (eVb2e1 - 1.0);
   else
      Ib1 = (1.0 - XIBI) * Ibf0 * ((1.0 - XREC) * (eVb2e1 - 1.0) +
            XREC * (eVb2e1 + eVb2c2star - 2.0) * (1.0 + Vtc / VEF_T));

   Ib1_s = XIBI * Ibf0 * (eVb1e1 - 1.0);
   `expLin(tmpExp,Vb2e1 * VtINV / MLF)
   Ib2 = IBF_TM * (tmpExp - 1.0) + GMIN * Vb2e1;
   `expLin(tmpExp,0.5 * Vb1c1 * VtINV)
   Ib3 = IBR_TM * (eVb1c1 - 1.0) /
         (tmpExp + exp(0.5 * VLR * VtINV)) +
         GMIN  * Vb1c1;

   // Iex, Isub (XIex, XIsub)

   g1 = If0 * eVb1c1;
   g2 = 4.0 * eVb1c1VDC;
   nBex = g1 / (1.0 + sqrt(1.0 + g1));
   pWex = g2 / (1.0 + sqrt(1.0 + g2));
   Iex = (1.0 / BRI_T) * (0.5 * IK_TM * nBex - IS_TM);

`ifdef SUBSTRATE
   Isub = 2.0 * ISS_TM * (eVb1c1 - 1.0) /
          (1.0 + sqrt(1.0 + 4.0 * (IS_TM / IKS_TM) * eVb1c1));
   Isf =  ISS_TM * (eVsc1 - 1.0);
`endif

   XIex =0.0;

`ifdef SUBSTRATE
   XIsub = 0.0;
`endif

   if (EXMOD == 1) begin

      Iex =   Iex *  Xext1;

`ifdef SUBSTRATE
      Isub  = Isub * Xext1;
`endif

      Xg1 = If0 * eVbc1;
      XnBex = Xg1 / (1.0 + sqrt(1.0 + Xg1));
      XIMex = XEXT * (0.5 * IK_TM * XnBex - IS_TM) / BRI_T;

`ifdef SUBSTRATE
      XIMsub = XEXT * 2.0 * ISS_TM * (eVbc1 - 1.0) /
               (1.0 + sqrt(1.0 + 4.0 * IS_T / IKS_T * eVbc1));
      Vex_bias = XEXT * (IS_TM / BRI_T + ISS_TM) * RCC_TM;
`else
      XIMsub = 0.0;
      Vex_bias = XEXT * (IS_TM / BRI_T) * RCC_TM;
`endif

      Vex = Vt * (2.0 - ln( Vex_bias * VtINV));
      vdif = Vbc1 - Vex;
      `max_hyp0(VBex, vdif, 0.11);

      Fex = VBex /(Vex_bias + (XIMex + XIMsub) * RCC_TM + VBex);
      XIex = Fex * XIMex;

`ifdef SUBSTRATE
      XIsub = Fex * XIMsub;
`endif

   end

   // Variable base resistance

   q0Q = 1.0 + Vte / VER_T + Vtc / VEF_T;
   `max_hyp0(q1Q, q0Q, 0.1);
   qBQ = q1Q * (1.0 + 0.5 * (n0 + nB));

   Rb2 = 3.0 * RBV_TM / qBQ;
   Ib1b2 =  (2.0 * Vt * (eVb1b2 - 1.0) + Vb1b2) / Rb2;

   // Weak-avalanche current

   Iavl = 0.0;
   Gem  = 0.0;
   if ((Ic1c2 > 0.0) && (Vb2c1 < VDC_T)) begin

      dEdx0 = 2.0 * VAVL / (WAVL * WAVL);
      sqr_arg = (VDC_T - Vb2c1) / Icap_IHC;
      xd = sqrt(2.0 * sqr_arg / dEdx0);
      if (EXAVL == 0.0)
         Weff = WAVL;
      else
         begin
           xi_w1 = 1.0 - 0.5 * xi_w;
           Weff = WAVL * xi_w1 * xi_w1;
         end
      Wd = xd * Weff / sqrt(xd * xd + Weff * Weff);
      Eav = (VDC_T - Vb2c1) / Wd;
      E0 = Eav + 0.5 * Wd * dEdx0 * Icap_IHC;

      if (EXAVL == 0)
         Em = E0;
      else
         begin
           SHw = 1.0 + 2.0 * SFH * (1.0 + 2.0 * xi_w);
           Efi = (1.0 + SFH) / (1.0 + 2.0 * SFH);
           Ew = Eav - 0.5 * Wd * dEdx0 * (Efi - Ic1c2 / (IHC_M * SHw));
           sqr_arg = (Ew - E0) * (Ew - E0) + 0.1 * Eav * Eav * Icap / IHC_M;
           Em = 0.5 * (Ew + E0 + sqrt(sqr_arg));
         end

      EmEav_Em = (Em - Eav) / Em;
      if (abs(EmEav_Em) > `TEN_M07)
         begin
           lambda = 0.5 * Wd / EmEav_Em;
           Gem = An / BnT * Em * lambda *
                (exp(-BnT / Em) - exp(-BnT / Em * (1.0 + Weff / lambda)) );
         end
      else
         Gem = An * Weff * exp(-BnT / Em);

      Gmax = Vt / (Ic1c2 * (RBC_TM + Rb2)) +  qBI / BF_T +
             RE_TM / (RBC_TM + Rb2);
      Iavl = Ic1c2 * Gem  / (Gem +Gem / Gmax + 1.0);
   end

`ifdef SELFHEATING
   // Power dissipation

   if (eVb2c2star > 0.0)
      Vb2c2star = Vt * ln(eVb2c2star);
   else
      Vb2c2star = Vb2c2;

   power =  In * (Vb2e1 - Vb2c2star) +
            Ic1c2 * (Vb2c2star - Vb2c1) -
            Iavl * Vb2c2star +
            Vee1 * Vee1 / RE_TM +
            Vcc1 * Vcc1 / RCC_TM +
            Vbb1 * Vbb1 / RBC_TM +
            Ib1b2 * Vb1b2 +
            (Ib1 + Ib2) * Vb2e1 +
            Ib1_s * Vb1e1 +
`ifdef SUBSTRATE
            (Iex + Isub + Ib3) * Vb1c1 +
            (XIex + XIsub) * Vbc1 -
            (XIsub + Isub - Isf) * Vsc1;
`else
            (Iex + Ib3) * Vb1c1 + XIex * Vbc1;
`endif

`endif


   // Charges

   Qte = (1.0 - XCJE) * CJE_TM * Vte;
   `min_logexp(Vje_s, Vb1e1, Vfe, a_VDE);
   Qte_s = XCJE * CJE_TM * (VDE_T / (1.0 - PE) *
           (1.0 - pow(1.0 - Vje_s / VDE_T, 1.0 - PE)) +
            `AJE * (Vb1e1 - Vje_s));

   Qtc = XCJC * CJC_TM * Vtc;
   Qb0 = TAUB_T * IK_TM;
   Qbe_qs = 0.5 * Qb0 * n0 * q1Q;
   Qbc_qs = 0.5 * Qb0 * nB * q1Q;

   a_VDC = 0.1 * VDC_T;
   `min_logexp(Vjcex, Vb1c1, Vfc, a_VDC);
   Vtexv = VDC_T / (1.0 - PC) * (1.0 - pow(1.0 - Vjcex / VDC_T, 1.0 - PC)) +
           bjc * (Vb1c1 - Vjcex);
   Qtex = CJC_TM * ((1.0 - XP_T) * Vtexv + XP_T * Vb1c1) *
          (1.0 - XCJC) * (1.0 - XEXT);

   `min_logexp(XVjcex, Vbc1, Vfc, a_VDC);
   XVtexv = VDC_T / (1.0 - PC) * (1.0 - pow(1.0 - XVjcex / VDC_T, 1.0 - PC)) +
            bjc * (Vbc1 - XVjcex);
   XQtex = CJC_TM * ((1.0 - XP_T) * XVtexv + XP_T * Vbc1) *
           (1.0 - XCJC) * XEXT;

`ifdef SUBSTRATE
   a_VDS = 0.1 * VDS_T;
   Vfs = VDS_T * (1.0 - pow(`AJS , -1.0 / PS));
   `min_logexp(Vjs, Vsc1, Vfs, a_VDS);
   Qts = CJS_TM * (VDS_T / (1.0 - PS) *
         (1.0 - pow(1.0 - Vjs / VDS_T, 1.0 - PS)) +  `AJS * (Vsc1 - Vjs));
`endif

   Qe0 = TAUE_T * IK_TM * pow(IS_TM / IK_TM, 1.0 / MTAU);
   `expLin(tmpExp,Vb2e1 / (MTAU * Vt))
   Qe = Qe0 * (tmpExp - 1.0);

   Qepi0 = 4.0 * TEPI_T * Vt / RCV_TM;
   Qepi = 0.5 * Qepi0 * xi_w * (p0star + pW + 2.0);

   Qex = TAUR_T * 0.5 * (Qb0 * nBex + Qepi0 * pWex) / (TAUB_T + TEPI_T);
   XQex = 0.0;

   if (EXMOD == 1) begin

      Qex = Qex * (1.0 - XEXT);
      Xg2 = 4.0 * eVbc1VDC;
      XpWex = Xg2 / (1.0 + sqrt(1.0 + Xg2));
      XQex = 0.5 * Fex * XEXT * TAUR_T *
             (Qb0 * XnBex + Qepi0 * XpWex) / (TAUB_T + TEPI_T);

   end

   Qb1b2 = 0.0;
   if (EXPHI == 1)
      begin
        dVteVje = pow(1.0 - Vje / VDE_T, -PE) - `AJE;
        Vb2e1Vfe = (Vb2e1 - Vfe) / a_VDE;
        if (Vb2e1Vfe < 0.0)
           dVjeVb2e1 = 1.0 / (1.0 + exp(Vb2e1Vfe));
        else
           dVjeVb2e1 = exp(- Vb2e1Vfe) / (1.0 + exp(- Vb2e1Vfe));

        dVteVb2e1 = dVteVje * dVjeVb2e1 + `AJE;
        dQteVb2e1 = (1.0 - XCJE) * CJE_TM * dVteVb2e1;

        dn0Vb2e1 = If0 * eVb2e1 * VtINV * (0.5 / sqrt(1.0 + f1));
        dQbeVb2e1 = 0.5 * Qb0 * q1Q * dn0Vb2e1;

        dQeVb2e1 = (Qe + Qe0) / (MTAU * Vt);

        Qb1b2 = 0.2 * Vb1b2 * (dQteVb2e1 + dQbeVb2e1 + dQeVb2e1);

        Qbc = Qbe_qs / 3.0 + Qbc_qs;
        Qbe = 2.0 * Qbe_qs / 3.0;
      end
    else
      begin
        Qbe = Qbe_qs;
        Qbc = Qbc_qs;
      end


// Add branch current contributions

   // Static currents
   I(c1, c2) <+ TYPE * Ic1c2;
   I(c2, e1) <+ TYPE * In;
   I(b1, e1) <+ TYPE * Ib1_s;
   I(b2, e1) <+ TYPE * (Ib1 + Ib2);
   I(b1, c1) <+ TYPE * (Ib3 + Iex);
   I(b,  c1) <+ TYPE * XIex;
`ifdef SUBSTRATE
   I(b1, s)  <+ TYPE * Isub;
   I(b,  s)  <+ TYPE * XIsub;
   I(s,  c1) <+ TYPE * Isf;
`endif
   I(b1, b2) <+ TYPE * Ib1b2;
   I(b2, c2) <+ TYPE * (-1.0 * Iavl);
   I(e,  e1) <+ TYPE * Vee1 / RE_TM;
   I(b,  b1) <+ TYPE * Vbb1 / RBC_TM;
   I(c,  c1) <+ TYPE * Vcc1 / RCC_TM;
 
`ifdef SELFHEATING
   // Electrical equivalent for the thermal network 
   I(dt) <+ V(dt) / RTH_Tamb_M;
   I(dt) <+ ddt(CTH_M * V(dt));
   I(dt) <+ -1.0 *  power;
`endif


   // Dynamic currents 
   I(b2, e1) <+ ddt(TYPE * (Qte + Qbe + Qe));
   I(b1, e1) <+ ddt(TYPE * (Qte_s));
   I(b2, c2) <+ ddt(TYPE * (Qtc + Qbc + Qepi));
   I(b1, c1) <+ ddt(TYPE * (Qtex + Qex));
   I(b,  c1) <+ ddt(TYPE * (XQtex + XQex));
`ifdef SUBSTRATE
   I(s,  c1) <+ ddt(TYPE * Qts);
`endif
   I(b1, b2) <+ ddt(TYPE * Qb1b2);
   I(b,   e) <+ ddt(TYPE * CBEO_M * Vbe);
   I(b,   c) <+ ddt(TYPE * CBCO_M * Vbc);

   end  // Currents and charges 


// Noise model
begin // Noise model
   // Thermal noise 
   common = 4.0 * `KB * Tk;
   powerREC = common / RE_TM;  // Emitter resistance
   powerRBC = common / RBC_TM; // Base resistance
   powerRCC = common / RCC_TM; // Collector resistance
   powerRBV = common / Rb2 * (4.0 * eVb1b2 + 5.0) / 3.0; // Variable base resistance
   I(e, e1)   <+ white_noise(powerREC);
   I(b, b1)   <+ white_noise(powerRBC);
   I(c, c1)   <+ white_noise(powerRCC);
   I(b1, b2)  <+ white_noise(powerRBV);

   // Collector current shot noise 
powerCCS = 2.0 * `QQ * (If + Ir) / qBI; 
`ifdef __VAMS_COMPACT_MODELING__
   w_n_1 = white_noise(1.0);
   I(c2, e1)  <+ powerCCS * w_n_1;
`else
   I(w_n_1)   <+ V(w_n_1);
   I(w_n_1)   <+ white_noise(1.0);
   I(c2, e1)  <+ powerCCS * V(w_n_1);
`endif

   // Excess noise due to the avalanche multiplication
   twoqIavl = KAVL * 2.0 * `QQ * Iavl;
   powerII  = twoqIavl * (Gem + 1.0); // Intrinsic avalanche multiplication noise
  // Adding intrinsic avalanche noise and amplified main current noise
`ifdef __VAMS_COMPACT_MODELING__
   I(b2, c2) <+ white_noise(powerII) + twoqIavl * w_n_1;
`else
   I(b2, c2) <+ white_noise(powerII) + twoqIavl * V(w_n_1);
`endif

   // Forward base current shot noise and 1/f noise 
   powerFBCS = 2.0 * `QQ * (abs(Ib1) + abs(Ib2));
   powerFBC1fB1 = (1.0 - XIBI) * pow((abs(Ib1) / (1 - XIBI)), AF) * KF_M; 
   exponentFBC1fB2 = (2.0 * (MLF - 1.0)) + (AF * (2.0 - MLF)); 
   powerFBC1fB2 = KFN_M * pow(abs(Ib2), exponentFBC1fB2); 
   I(b2, e1)  <+ white_noise(powerFBCS);
   I(b2, e1)  <+ flicker_noise(powerFBC1fB1, 1);
   I(b2, e1)  <+ flicker_noise(powerFBC1fB2, 1);

  // Emitter-base sidewall current shot and 1/f noise 
   powerEBSCS = 2.0 * `QQ * abs(Ib1_s);
   if (XIBI == 0)
       powerEBSC1f = 0.0;
   else
       powerEBSC1f = KF_M * XIBI * pow((abs(Ib1_s / XIBI)), AF);
   I(e1, b1)  <+ white_noise(powerEBSCS);
   I(e1, b1)  <+ flicker_noise(powerEBSC1f, 1);

  // Reverse base current shot noise and 1/f noise
   powerRBCS = 2.0 * `QQ * abs(Ib3);
   powerRBC1f = KF_M * pow(abs(Ib3), AF);
   I(b1, c1)  <+ white_noise(powerRBCS);
   I(b1, c1)  <+ flicker_noise(powerRBC1f, 1);

   // Extrinsic current shot noise and 1/f noise
   powerExCS = 2.0 * `QQ * abs(Iex);
   powerExC1f = KF_M * (1 - (EXMOD * XEXT)) *
                pow((abs(Iex) / (1 - (EXMOD * XEXT))), AF);
   powerExCSMOD = 2.0 * `QQ * abs(XIex) * EXMOD;
   if (XEXT == 0.0)
       powerExC1fMOD = 0.0;
   else
       powerExC1fMOD = KF_M * EXMOD * XEXT * pow((abs(XIex) / XEXT), AF);
   I(b1, c1)  <+ white_noise(powerExCS);
   I(b1, c1)  <+ flicker_noise(powerExC1f, 1);
   I(b, c1)   <+ white_noise(powerExCSMOD);
   I(b, c1)   <+ flicker_noise(powerExC1fMOD, 1);

`ifdef SUBSTRATE
   // Substrate current shot noise (between nodes B1 and S, resp. B and S)
   powerSubsCS_B1S = 2.0 * `QQ * abs(Isub);
   powerSubsCS_BS  = 2.0 * `QQ * abs(XIsub);
   I(b1, s)   <+ white_noise(powerSubsCS_B1S);
   I(b, s)    <+ white_noise(powerSubsCS_BS);
`endif

end // Noise model

// Evaluate the operating point (outout) variables
begin 

`ifdef __VAMS_COMPACT_MODELING__

// The external currents and the current gain
OP_ic     = I(<c>);        // External DC collector current
OP_ib     = I(<b>);        // External DC base Current
OP_betadc = OP_ic / OP_ib;  // External DC Current gain

// The internal voltage differences
OP_vb2e1 = Vb2e1;  // Internal base-emiter bias
OP_vb2c2 = Vb2c2;  // Internal base-emiter bias 
OP_vb2c1 = Vb2c1;  // Internal base-collector bias including epilayer
OP_vb1c1 = Vb1c1;  // External base-collector bias without contact resistances
OP_ve1e  = - Vee1; // Bias over emiter resistance

// The branch currents
OP_in    = In;            // Main current
OP_ic1c2 = Ic1c2;         // Epilayer current
OP_ib1b2 = Ib1b2;         // Pinched-base current
OP_ib1   = Ib1;           // Ideal forward base current
OP_isb1  = Ib1_s;         // Ideal side-wall base current
OP_ib2   = Ib2;           // Non-ideal forward base current
OP_ib3   = Ib3;           // Non-ideal reverse base current
OP_iavl  = Iavl;          // Avalanche current
OP_iex   = Iex;           // Extrinsic reverse base current
OP_xiex  = XIex;          // Extrinsic reverse base current
OP_isub  = Isub;          // Substrate current
OP_xisub = XIsub;         // Substrate current
OP_isf   = Isf;           // Substrate failure current
OP_ire   = Vee1 / RE_TM;  // Current through emiter resistance
OP_irbc  = Vbb1 / RBC_TM; // Current through constant base resistance
OP_ircc  = Vcc1 / RCC_TM; // Current through constant collector resistance

// The branch charges
OP_qe    = Qe;    // Emitter charge or emitter neutral charge
OP_qte   = Qte;   // Base-emiter depletion charge
OP_qste  = Qte_s; // Sidewall base-emiter depletion charge
OP_qbe   = Qbe;   // Base-emiter diffusion charge
OP_qbc   = Qbc;   // Base-collector diffusion charge
OP_qtc   = Qtc;   // Base-colector depletion charge
OP_qepi  = Qepi;  // Epilayer diffusion charge 
OP_qb1b2 = Qb1b2; // AC current crowding charge
OP_qtex  = Qtex;  // Extrinsic base-collector depletion charge
OP_xqtex = XQtex; // Extrinsic base-collector depletion charge
OP_qex   = Qex;   // Extrinsic base-collector diffusion charge
OP_xqex  = XQex;  // Extrinsic base-collector diffusion charge
OP_qts   = Qts;   // Collector substrate depletion charge

// Small signal equivalent circuit conductances and resistances
ddxInVb2 = ddx(In, V(b2));
OP_gx     = ddxInVb2            - ddx(In, V(e1));      // Forward transconductance
OP_gy     = ddxInVb2            - ddx(In, V(c2));      // Reverse transconductance
OP_gz     = ddxInVb2            - ddx(In, V(c1));      // Reverse transconductance
OP_gspi   = ddx(Ib1_s, V(b1))   - ddx(Ib1_s, V(e1));   // Conductance sidewal b-e junction
OP_gpix   = ddx(Ib1+Ib2, V(b2)) - ddx(Ib1+Ib2, V(e1)); // Conductance floor b-e junction
ddxIb1Vb2 = ddx(Ib1, V(b2));
OP_gpiy   = ddxIb1Vb2     - ddx(Ib1, V(c2));  // Early effect on recombination base current
OP_gpiz   = ddxIb1Vb2     - ddx(Ib1, V(c1));  // Early effect on recombination base current
ddxIavlVb2 = - ddx(Iavl, V(b2));
OP_gmux   = ddxIavlVb2   - ddx(-Iavl, V(e1)); // Early effect on avalanche current limitting
OP_gmuy   = ddxIavlVb2   - ddx(-Iavl, V(c2)); // Conductance of avalanche current
OP_gmuz   = ddxIavlVb2   - ddx(-Iavl, V(c1)); // Conductance of avalanche current
OP_gmuex  = ddx(Iex+Ib3, V(b1)) - ddx(Iex+Ib3, V(c1)); // Conductance extrinsic b-c current
OP_xgmuex = ddx(XIex, V(b))     - ddx(XIex, V(c1));    // Conductance extrinsic b-c current
ddxIc1c2Vb2 = ddx(Ic1c2, V(b2));
OP_grcvy  = ddxIc1c2Vb2  - ddx(Ic1c2, V(c2)); // Conductance of epilayer current
OP_grcvz  = ddxIc1c2Vb2  - ddx(Ic1c2, V(c1)); // Conductance of epilayer current
OP_rbv    = 1.0 / (ddx(Ib1b2, V(b1)) - ddx(Ib1b2, V(b2)));  // Base resistance
ddxIb1b2Vb2 = ddx(Ib1b2, V(b2));
OP_grbvx  = ddxIb1b2Vb2  - ddx(Ib1b2, V(e1)); // Early effect on base resistance
OP_grbvy  = ddxIb1b2Vb2  - ddx(Ib1b2, V(c2)); // Early effect on base resistance
OP_grbvz  = ddxIb1b2Vb2  - ddx(Ib1b2, V(c1)); // Early effect on base resistance
OP_re     = RE_TM;  // Emiter resistance
OP_rbc    = RBC_TM; // Constant base resistance
OP_rcc    = RCC_TM; // Constant collector resistance
OP_gs     = ddx(Isub, V(b1)) - ddx(Isub, V(c1));  // Conductance parasitic PNP transitor
OP_xgs    = ddx(XIsub, V(b)) - ddx(XIsub, V(c1)); // Conductance parasitic PNP transistor
OP_gsf    = ddx(Isf, V(s))   - ddx(Isf, V(c1));   // Conductance substrate failure current

// Small signal equivalent circuit capacitances
OP_csbe   = ddx(Qte_s, V(b1)) - ddx(Qte_s, V(e1)); // Capacitance sidewall b-e junction
OP_cbex   = ddx(Qte + Qbe + Qe, V(b2)) -           // Capacitance floor b-e junction
            ddx(Qte + Qbe + Qe, V(e1));
ddxQbeVb2 = ddx(Qbe, V(b2));
OP_cbey   = ddxQbeVb2        - ddx(Qbe, V(c2)); // Early effect on b-e diffusion junction
OP_cbez   = ddxQbeVb2        - ddx(Qbe, V(c1)); // Early effect on b-e diffusion junction
OP_cbcx   = ddx(Qbc, V(b2))  - ddx(Qbc, V(e1)); // Early effect on b-c diffusion junction
ddxQcciVb2 = ddx(Qtc + Qbc + Qepi, V(b2));
OP_cbcy   = ddxQcciVb2 - ddx(Qtc + Qbc + Qepi, V(c2)); // Capacitance floor b-c junction
OP_cbcz   = ddxQcciVb2 - ddx(Qtc + Qbc + Qepi, V(c1)); // Capacitance floor b-c junction
OP_cbcex  = ddx(Qtex + Qex,V(b1)) -   // Capacitance extrinsic b-c junction 
            ddx(Qtex + Qex,V(c1)); 
OP_xcbcex = ddx(XQtex + XQex, V(b)) - // Capacitance extrinsic b-c junction
            ddx(XQtex + XQex, V(c1));
ddxQb1b2Vb1 = ddx(Qb1b2, V(b1));
OP_cb1b2  = ddxQb1b2Vb1 - ddx(Qb1b2, V(b2));  // Capacitance AC current crowding
OP_cb1b2x = ddxQb1b2Vb1 - ddx(Qb1b2, V(e1));  // Cross-capacitance AC current crowding
OP_cb1b2y = ddxQb1b2Vb1 - ddx(Qb1b2, V(c2));  // Cross-capacitance AC current crowding
OP_cb1b2z = ddxQb1b2Vb1 - ddx(Qb1b2, V(c1));  // Cross-capacitance AC current crowding
OP_cts    = ddx(Qts, V(s)) - ddx(Qts, V(c1)); // Capacitance s-c junction

// Approximate small signal equivalent circuit
dydx = (OP_gx - OP_gmux)            / (OP_grcvy + OP_gmuy - OP_gy);
dydz = (OP_gz - OP_grcvz - OP_gmuz) / (OP_grcvy + OP_gmuy - OP_gy);
gpi  = OP_gspi + OP_gpix + OP_gmux + OP_gpiz + OP_gmuz +    
       (OP_gpiy + OP_gmuy) * (dydx + dydz);
OP_gm   = (OP_grcvy * (OP_gx - OP_gmux +                   // Transconductance
           OP_gz - OP_gmuz) - OP_grcvz * 
          (OP_gy - OP_gmuy)) / (OP_grcvy + OP_gmuy - OP_gy);
OP_beta = OP_gm / gpi;                                     // Current amplification
OP_gout = ((OP_gy - OP_gmuy) * OP_grcvz -                  // Output conductance
          (OP_gz - OP_gmuz) * OP_grcvy) / 
          (OP_grcvy + OP_gmuy - OP_gy);
OP_gmu  = OP_gpiz + OP_gmuz + (OP_gpiy + OP_gmuy) * dydz + // Feedback transconductance
          OP_gmuex + OP_xgmuex;
OP_rb   = RBC_TM + OP_rbv;                                 // Base resistance
OP_cbe  = OP_cbex + OP_csbe + OP_cbcx +                    // Base-emitter capacitance
          (OP_cbey + OP_cbcy) * dydx + CBEO_M;
OP_cbc  = (OP_cbey + OP_cbcy) * dydz + OP_cbcz +           // Base-collector capacitance
           OP_cbcex + OP_xcbcex + CBCO_M;

// Quantities to describe internal state of the model
gammax = (OP_gpix + OP_gmux - OP_grbvx) * OP_rbv;
gammay = (OP_gpiy + OP_gmuy - OP_grbvy) * OP_rbv;
gammaz = (OP_gpiz + OP_gmuz - OP_grbvz) * OP_rbv;
gbfx   =  OP_gpix + OP_gspi * (1.0 + gammax);
gbfy   =  OP_gpiy + OP_gspi * gammay;
gbfz   =  OP_gpiz + OP_gspi * gammaz;
alpha_ft  = (1.0 + (OP_grcvy * dydx * RCC_TM) +
            (OP_gx + gbfx + (OP_gy + gbfy) * dydx) * RE_TM)/
            (1.0 - (OP_grcvz + OP_grcvy * dydz) * RCC_TM -
            (OP_gz + gbfz + (OP_gy + gbfy) * dydz) * RE_TM);
rx = pow((OP_grcvy * dydx + alpha_ft * (OP_grcvz + OP_grcvy * dydz)), -1);
rz = alpha_ft * rx;
ry = (1.0 - OP_grcvz * rz) / OP_grcvy;
rb1b2 = gammax * rx + gammay * ry + gammaz * rz;
rex = rz + rb1b2;
xrex = rex + RBC_TM * ((gbfx + OP_gmux) * rx + (gbfy + OP_gmuy) * ry +
       (gbfz + OP_gmuz) * rz);
taut = OP_csbe * (rx + rb1b2) + (OP_cbex + OP_cbcx) * rx + (OP_cbey + OP_cbcy) *
          ry + (OP_cbez + OP_cbcz) * rz + OP_cbcex * rex + OP_xcbcex * xrex +
          (CBEO_M + CBCO_M) * (xrex - RCC_TM);

OP_ft = 1.0 / (2.0 * `PI * taut); // Good approximation for cut-off frequency
OP_iqs = Iqs;                     // Current at onset of quasi-saturation
OP_xiwepi = xi_w;                 // Thickness of injection layer
OP_vb2c2star = Vb2c2star;         // Physical value of internal base-collector bias

//self-heating
OP_pdiss = power;                 // Dispipatiom
OP_tk = Tk;                       // Actual temperature

`endif
end
