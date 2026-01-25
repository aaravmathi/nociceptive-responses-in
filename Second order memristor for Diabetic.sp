*Second order memristor for Diabetic Foot Ulcer Detection

.SUBCKT MEM_Second_order TE BE XSV G R T_CF T_bulk
* Model Parameters *
.PARAM a1=0.2 a2=0.3 b=0.083 Vp=0.4 Vn=0.4
.PARAM Ap=10000 An=10000 xp=0.25 xn=0.5 alphap=1 alphan=2
.PARAM xo=0.2 eta=1
.PARAM Ea=0.15 kB=8.617e-5 T_amb=300
.PARAM Cp=0.15e-6 k_th=1e-3 delta=1e-4
.PARAM H=1
.PARAM L=6e-9 P=2.5e-9 r0=2.2e-9 rm=0.6e-9 alpha_r=2e4 beta_r=1e4
.PARAM gm=0.6e-9

* State Variables *
Cx XSV 0 1
.ic V(XSV) = xo
Cg G 0 1
.ic V(G) = 0.2e-9
Cr R 0 1
.ic V(R) = 5.5e-9
.ic V(T_CF) = T_amb
.ic V(T_bulk) = T_amb

* Multiplicative Functions *
.func wp(V) = xp/(1-xp) - V/(1-xp) + 1
.func wn(V) = V/(1-xn)
.func kT(V) = kB * V / 300  

* Device Threshold Function G(V) *
.func G(V) = IF(V <= Vp, IF(V >= -Vn, 0, -An*(exp(-V)-exp(Vn))), Ap*(exp(V)-exp(Vp)))

* First Order State Variable Motion Function F(V1, V2) *
.func F(V1, V2) = IF(eta*V1 >= 0, IF(V2 >= xp, exp(-alphap*(V2-xp))*wp(V2), 1), IF(V2 <= (1-xn), exp(alphan*(V2+xn-1))*wn(V2), 1))

* IV Relationship - Hyperbolic Sine Model *
.func IVRel(V1, V2) = IF(V1 >= 0, a1*V2*sinh(b*V1), a2*V2*sinh(b*V1))

* State Equation for x *
Gx 0 XSV VALUE = {eta * F(V(TE,BE), V(XSV)) * G(V(TE,BE))}

* Gap and Radius Dynamics *
Ggap 0 G VALUE = {
+ IF(V(TE,BE) > 0,
+   V(G) - alpha_r * exp(-Ea/(kB*V(T_CF))) * (V(G) - gm),  ; SET: Close gap
+   V(G) + beta_r * exp(-Ea/(kB*V(T_CF))) * (L - P - V(G)) ; RESET: Open gap
+)
+}

Gr 0 R VALUE = {
+ IF(V(TE,BE) < 0,
+   min(r0, V(R) + beta_r * exp(-Ea/(kB*max(V(T_CF), 250))) * log(1 + (r0 - V(R)))),
+   IF(V(TE,BE) > 0,
+      max(rm, V(R) - alpha_r * exp(-Ea/(kB*max(V(T_CF), 250))) * log(1 + (V(R) - rm))),
+      V(R)
+   )
+)
+}

* CF Temperature Dynamics *
G_TCF 0 T_CF VALUE = {
+  (1/Cp) * (V(TE,BE) * I(Gm) - k_th * (V(T_CF) - V(T_bulk)))
+}
CTCF T_CF 0 1e-4
RTCF T_CF 0 1G
.ic V(T_CF) = 300

* Bulk Temperature Dynamics *
G_Tbulk 0 T_bulk VALUE = {
+  (1/Cp) * (H * V(TE,BE) * I(Gm) - delta * (V(T_bulk) - 300))
+}
CTbulk T_bulk 0 1e-4
RTbulk T_bulk 0 1G
.ic V(T_bulk) = 300

* Memristor Current Response *
Gm TE BE VALUE = {IVRel(V(TE,BE), V(XSV))}

* High-value resistors for floating nodes *
R1 XSV 0 1G
R2 G 0 1G
R3 R 0 1G
R4 T_CF 0 1G
R5 T_bulk 0 1G

.ENDS MEM_Second_order