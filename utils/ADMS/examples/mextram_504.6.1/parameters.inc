// Mextram parameters

`PAI(Model level, )
     LEVEL = 504      from [504:505);
`PAR(Reference temperature, Celsius)
     TREF  = 25.0     from [-273.0:inf);
`PAR(Difference between the local and global ambient temperatures, Celsius)
     DTA   = 0.0;
`PAI(Flag for extended modeling of the reverse current gain, )
     EXMOD = 1        from [0:1];
`PAI(Flag for the distributed high-frequency effects in transient, )
     EXPHI = 1        from [0:1];
`PAI(Flag for extended modeling of avalanche currents, )
     EXAVL = 0        from [0:1];

`PAR(Collector-emitter saturation current, A)
     IS = 22.0a       from (0.0:inf);
`PAR(Collector-emitter high injection knee current, A)
     IK = 0.1         from [1.0p:inf);
`PAR(Reverse Early voltage, V)
     VER   = 2.5      from [0.01:inf);
`PAR(Forward Early voltage, V)
     VEF   = 44.0     from [0.01:inf);
`PAR(Ideal forward current gain, )
     BF    = 215.0    from [0.1m:inf);
`PAR(Saturation current of the non-ideal forward base current, A)
     IBF   = 2.7f     from [0.0:inf);
`PAR(Non-ideality factor of the non-ideal forward base current, )
     MLF   = 2.0      from [0.1:inf);
`PAR(Part of ideal base current that belongs to the sidewall, )
     XIBI  = 0.0      from [0.0:1.0];
`PAR(Ideal reverse current gain, )
     BRI   = 7.0      from [1.0e-4:inf);
`PAR(Saturation current of the non-ideal reverse base current, A)
     IBR   = 1.0f     from [0.0:inf);
`PAR(Cross-over voltage of the non-ideal reverse base current, V)
     VLR   = 0.2;
`PAR(Part of currents and charges that belong to extrinsic region, )
     XEXT  = 0.63     from [0.0:1.0];

`PAR(Epilayer thickness used in weak-avalanche model, m)
     WAVL  = 1.1u     from [1.0n:inf);
`PAR(Voltage determining curvature of avalanche current, V)
     VAVL  = 3.0      from [0.01:inf);
`PAR(Current spreading factor of avalanche model when EXAVL=1, )
     SFH   = 0.3      from [0.0:inf);

`PAR(Emitter resistance, Ohm)
     RE    = 5.0      from [1.0u:inf);
`PAR(Constant part of the base resistance, Ohm)
     RBC   = 23.0     from [1.0u:inf);
`PAR(Zero-bias value of the variable part of the base resistance, Ohm)
     RBV   = 18.0     from [1.0u:inf);
`PAR(Constant part of the collector resistance, Ohm)
     RCC   = 12.0     from [1.0u:inf);
`PAR(Resistance of the un-modulated epilayer, Ohm)
     RCV   = 150.0    from [1.0u:inf);
`PAR(Space charge resistance of the epilayer, Ohm)
     SCRCV = 1250.0   from [1.0u:inf);
`PAR(Critical current for velocity saturation in the epilayer, A)
     IHC   = 4.0m     from [1.0p:inf);
`PAR(Smoothness parameter for the onset of quasi-saturation, )
     AXI   = 0.3      from [0.02:inf);

`PAR(Zero-bias emitter-base depletion capacitance, F)
     CJE   = 73.0f    from [0.0:inf);
`PAR(Emitter-base diffusion voltage, V)
     VDE   = 0.95     from [0.05:inf);
`PAR(Emitter-base grading coefficient, )
     PE    = 0.4      from [0.01:0.99);
`PAR(Sidewall fraction of the emitter-base depletion capacitance, )
     XCJE  = 0.4 from [0.0:1.0];
`PAR(Emitter-base overlap capacitance, F)
     CBEO  = 0.0      from [0.0:inf);

`PAR(Zero-bias collector-base depletion capacitance, F)
     CJC   = 78.0f    from [0.0:inf);
`PAR(Collector-base diffusion voltage, V)
     VDC   = 0.68     from [0.05:inf);
`PAR(Collector-base grading coefficient, V)
     PC    = 0.5      from [0.01:0.99);
`PAR(Constant part of Cjc, )
     XP    = 0.35     from [0.0:0.99);
`PAR(Coefficient for current modulation of CB depletion capacitance, )
     MC    = 0.5      from [0.0:1.0);
`PAR(Fraction of CB depletion capacitance under the emitter, )
     XCJC  = 32.0m    from [0.0:1.0];
`PAR(Collector-base overlap capacitance, F)
     CBCO  = 0.0      from [0.0:inf);

`PAR(Non-ideality factor of the emitter stored charge, )
     MTAU  = 1.0      from [0.1:inf);
`PAR(Minimum transit time of stored emitter charge, s)
     TAUE  = 2.0p     from [0.0:inf);
`PAR(Transit time of stored base sharge, s)
     TAUB  = 4.2p     from (0.0:inf);
`PAR(Transit time of stored epilayer charge, s)
     TEPI  = 41.0p    from (0.0:inf);
`PAR(Transit time of reverse extrinsic stored base charge, s)
     TAUR  = 520.0p   from [0.0:inf);

`PAR(Bandgap difference over the base, eV)
     DEG   = 0.0;
`PAR(Pre-factor of the recombination part of Ib1, )
     XREC  = 0.0      from [0.0:inf);

`PAR(Temperature coefficient of the zero-bias base charge, )
     AQBO  = 0.3;
`PAR(Temperature coefficient of the resistivity of the emitter, )
     AE    = 0.0;
`PAR(Temperature coefficient of the resistivity of the base, )
     AB    = 1.0;
`PAR(Temperature coefficient of the resistivity of the epilayer, )
     AEPI  = 2.5;
`PAR(Temperature coefficient of the resistivity of the extrinsic base, )
     AEX   = 0.62;
`PAR(Temperature coefficient of the resistivity of the buried layer, )
     AC    = 2.0;
`PAR(Fine tuning of temperature dependence of C-E saturation current, )
     DAIS  = 0.0;
`PAR(Band-gap voltage difference of the forward current gain, V)
     DVGBF = 50.0m;
`PAR(Band-gap voltage difference of the reverse current gain, V)
     DVGBR = 45.0m;
`PAR(Band-gap voltage of the base, V)
     VGB   = 1.17    from [0.1:inf);
`PAR(Band-gap voltage of the collector, V)
     VGC   = 1.18    from [0.1:inf);
`PAR(Band-gap voltage recombination emitter-base junction, V)
     VGJ   = 1.15    from [0.1:inf);
`PAR(Band-gap voltage difference of emitter stored charge, V)
     DVGTE = 0.05;


`PAR(Exponent of the Flicker-noise, )
     AF    = 2.0     from [0.01:inf);
`PAR(Flicker-noise coefficient of the ideal base current, )
     KF    = 20.0p   from [0.0:inf);
`PAR(Flicker-noise coefficient of the non-ideal base current, )
     KFN   = 20.0p   from [0.0:inf);
`PAI(Switch for white noise contribution due to avalanche, )
     KAVL  = 0       from [0:1];


`ifdef SUBSTRATE
`PAR(Base-substrate saturation current, A)
     ISS   = 48.0a   from [0.0:inf);
`PAR(Base-substrate high injection knee current, A)
     IKS   = 250.0u  from [1.0p:inf);
`PAR(Zero-bias collector-substrate depletion capacitance, F)
     CJS   = 315.0f  from [0:inf);
`PAR(Collector-substrate diffusion voltage, V)
     VDS   = 0.62    from (0.05:inf);
`PAR(Collector-substrate grading coefficient, )
     PS    = 0.34    from (0.01:0.99);
`PAR(Band-gap voltage of the substrate, V)
     VGS   = 1.20    from [0.1:inf);
`PAR(Substrate temperature coefficient, )
     AS    = 1.58;
`endif

`ifdef SELFHEATING
`PAR(Thermal resistance, Celsius/W)
     RTH   = 300.0   from (0.0:inf);
`PAR(Thermal capacitance, J/Celsius)
     CTH   = 3.0n    from [0.0:inf);
`PAR(Temperature coefficient of the thermal resistance, )
     ATH   = 0.0;
`endif

`PAR(Multiplication factor, )
MULT       = 1.0     from (0.0:inf);


// Non-standard model parameters
// (introduced for users' convenience)

`PAI(Flag for NPN (1) or PNP (-1) transistor type, )
TYPE       = 1       from [-1:1] exclude 0;
`PAR(Minimum conductance, S)
GMIN       = 1.0e-12  from (0:1e-10];

