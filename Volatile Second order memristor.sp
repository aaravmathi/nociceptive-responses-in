*Second order memristor in Volatile Nociceptor regime
.SUBCKT MEM_Second_order TE BE XSV G R T_CF T_bulk

* THERMAL PARAMETERS  
.PARAM Cp=5e-11 Kth=4.0e-5 T_amb=300 kB=8.617e-5
.PARAM I0=0.02  g0=0.5e-9

* --- GEOMETRY ---
.PARAM g_min=0.2e-9 g_max=3.0e-9
.PARAM r_min=1.0e-9 r_max=10.0e-9

* --- DYNAMICS ---
.PARAM Relax_Rate=2.0e4

* --- FUNCTIONS ---
.func S(x) = 0.5*(1+tanh(x))
.func AbsS(x) = x*tanh(10*x)
.func T_Safe(T) = limit(T, 200, 2500)

.func Vp_eff(T) = Vp0 - gammaV*(T - T_amb)
.func Vn_eff(T) = Vn0 - gammaV*(T - T_amb)

* Gates
.func GateP(V,T) = S((V - Vp_eff(T))/Vs)
.func GateN(V,T) = S((-V - Vn_eff(T))/Vs)

* Masks
.func Mprog(V) = S((AbsS(V) - Vprog)/Vs)
.func Midle(V) = S((Vidle - AbsS(V))/Vs)

.func Rate_SET(T)   = exp(-Ea_set   /(kB*T_Safe(T)))
.func Rate_RESET(T) = exp(-Ea_reset /(kB*T_Safe(T)))

* Current Equation
.func I_Phys(V,g,r) = I0*((r/r_min)**2)*exp(-(g-g_min)/g0)*sinh(1.2*V)
.func Wg(g) = ((g-g_min)*(g_max-g))/((g_max-g_min)**2)
.func Wr(r) = ((r-r_min)*(r_max-r))/((r_max-r_min)**2)

* --- STATES ---
C_TCF T_CF 0 1
R_TCF T_CF 0 1G
C_TB  T_bulk 0 1
R_TB  T_bulk 0 1G
Cg G 0 1
R_G G 0 1G
Cr R 0 1
R_R R 0 1G
Cx XSV 0 1
R_XSV XSV 0 1G
.ic V(T_CF)={T_amb} V(T_bulk)={T_amb} V(G)=2.95e-9 V(R)=1.0e-9

* --- EQUATIONS ---
Gm TE BE VALUE = { limit(I_Phys(V(TE,BE), V(G), V(R)), -200m, 200m) }
G_TCF 0 T_CF VALUE = { (1/Cp) * ( AbsS(V(TE,BE)*I(Gm)) - Kth*(V(T_CF)-V(T_bulk)) ) }
G_TB  0 T_bulk VALUE = { 4e-3*(V(T_CF)-V(T_bulk)) - 8e-4*(V(T_bulk)-T_amb) }

* Gap Dynamics (Volatile)
* +V (GateP) -> Rate_SET (Closes Gap/Increases G) -> Sensitization
.func dG_raw(V,T,G) = Mprog(V) * Wg(G) * ( - GateP(V,T) * Rate_SET(T) * Ap_g * (G - g_min) + GateN(V,T) * Rate_RESET(T) * An_g * (g_max - G) ) + Midle(V) * Relax_Rate * (g_max - G)
Ggap 0 G VALUE = { limit(dG_raw(V(TE,BE), V(T_CF), V(G)), -500, 500) }

* Radius Dynamics (Persistent)
.func dR_raw(V,T,R) = Mprog(V) * Wr(R) * ( - GateP(V,T) * Rate_RESET(T) * An_r * (R - r_min) + GateN(V,T) * Rate_SET(T) * Ap_r * (r_max - R) )
Gr 0 R VALUE = { limit(dR_raw(V(TE,BE), V(T_CF), V(R)), -10, 10) }

Gx 0 XSV VALUE = { limit( 1 - (V(G)/g_max), 0, 1) }
.ENDS MEM_Second_order
