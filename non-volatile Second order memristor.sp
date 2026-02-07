*Second order memristor(non-Volatile)
.SUBCKT MEM_Second_order TE BE XSV G R T_CF T_bulk

* --- PARAMETERS ---
.PARAM Vp0=0.4 Vn0=0.4
.PARAM Vp_min=0.20 Vn_min=0.20
.PARAM gammaV=0.0030

.PARAM T_amb=300 kB=8.617e-5 T_limit=2500
.PARAM Cp=2e-14         
.PARAM CpB=2e-12         
.PARAM Kth=1.2e-7          
.PARAM Kamb=2e-3        
.PARAM g_min=0.2e-9 g_max=3.0e-9 r_min=1.0e-9 r_max=10.0e-9
.PARAM Ea_set=0.20 Ea_reset=0.18
.PARAM Ap_g=9.0e5  An_g=8.0e5  Ap_r=9.0e5  An_r=8.0e5
.PARAM I0=2e-4 g0=0.5e-9
.PARAM Vs=0.02

* --- FUNCTIONS ---
.func T_Safe(T) = limit(T, 200, T_limit)
.func Rate_SET(T)   = exp(-Ea_set   /(kB*T_Safe(T)))
.func Rate_RESET(T) = exp(-Ea_reset /(kB*T_Safe(T)))
.func Vp_eff(T) = max(Vp_min, Vp0 - gammaV*(T - T_amb))
.func Vn_eff(T) = max(Vn_min, Vn0 - gammaV*(T - T_amb))
.func GateP_T(V,T) = 1/(1+exp(-(V - Vp_eff(T))/Vs))
.func GateN_T(V,T) = 1/(1+exp(-((-V) - Vn_eff(T))/Vs))  ; activates when V is sufficiently negative
.func I_Phys(V,g,r) = I0*((r/r_min)**2)*exp(-(g-g_min)/g0)*sinh(1.5*V)
.func Wg(g) = ((g-g_min)*(g_max-g))/((g_max-g_min)**2)
.func Wr(r) = ((r-r_min)*(r_max-r))/((r_max-r_min)**2)
.func ReadMask(V) = if(abs(V) <= 0.25, 0, 1)
* --- STATE STORAGE 
C_TCF  T_CF   0 {Cp}
R_TCF  T_CF   0 1G
C_TB   T_bulk 0 {CpB}
R_TB   T_bulk 0 1G
Cg G 0 1
R_G G 0 1G
Cr R 0 1
R_R R 0 1G
Cx XSV 0 1
R_XSV XSV 0 1G
.ic V(T_CF)={T_amb} V(T_bulk)={T_amb} V(G)=1.6e-9 V(R)=5.0e-9
Gm TE BE VALUE = { limit(I_Phys(V(TE,BE), V(G), V(R)), -8e-3, 8e-3) }
G_TCF 0 T_CF VALUE = {
+ ( ReadMask(V(TE,BE))*abs(V(TE,BE)*I(Gm))
+   - Kth*(V(T_CF)-V(T_bulk)) )
+}
G_TB  0 T_bulk VALUE = {
+ ( Kth*(V(T_CF)-V(T_bulk))
+   - Kamb*(V(T_bulk)-T_amb) )
+}
Ggap 0 G VALUE = {
+ ReadMask(V(TE,BE)) * Wg(V(G)) *
+ ( GateP_T(V(TE,BE),V(T_CF)) * Rate_RESET(V(T_CF)) * An_g * (g_max - V(G))
+ - GateN_T(V(TE,BE),V(T_CF)) * Rate_SET(V(T_CF))   * Ap_g * (V(G) - g_min) )
+}
Gr 0 R VALUE = {
+ ReadMask(V(TE,BE)) * Wr(V(R)) *
+ ( - GateP_T(V(TE,BE),V(T_CF)) * Rate_RESET(V(T_CF)) * An_r * (V(R) - r_min)
+ + GateN_T(V(TE,BE),V(T_CF)) * Rate_SET(V(T_CF))    * Ap_r * (r_max - V(R)) )
+}
Gx 0 XSV VALUE = { limit( 1 - (V(G)/g_max), 0, 1) }
.ENDS MEM_Second_order