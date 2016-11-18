nx[a_,b_,c_] := 1/2 * a*b*c * NIntegrate[1/(s+a^2)/Sqrt[(s + a^2)*(s + b^2)*(s + c^2)],{s,0,Infinity},PrecisionGoal->10]
ny[a_,b_,c_] := 1/2 * a*b*c * NIntegrate[1/(s+b^2)/Sqrt[(s + a^2)*(s + b^2)*(s + c^2)],{s,0,Infinity},PrecisionGoal->10]
nz[a_,b_,c_] := 1/2 * a*b*c * NIntegrate[1/(s+c^2)/Sqrt[(s + a^2)*(s + b^2)*(s + c^2)],{s,0,Infinity},PrecisionGoal->10]

Em[a_,b_,c_,Bm_,mu_] := -(mu-1)/12*Bm*(1/(1+(mu-1)*nx[a,b,c]) + 1/(1+(mu-1)*ny[a,b,c]))
(* we use constraint a*b*c=1; gamma=1; *)

(* constraining a>b>c *)
(* checked with http://planetcalc.com/149/ *)
(* https://www.researchgate.net/post/Briefly_how_do_you_calculate_the_surface_area_of_ellipsoids *)

AreaEl[c_,b_,a_] := (
        cosphi = c/a;
        phi = ArcCos[cosphi];
        sinphi = Sqrt[1 - cosphi^2];
        k = Sqrt[a^2*(b^2-c^2)/b^2/(a^2-c^2)];
        2*Pi*c^2 + 2*Pi*a*b/sinphi * (sinphi^2 * EllipticE[phi,k] + cosphi^2 * EllipticF[phi,k])
        )

AreaEl2[c_,b_,a_] := (
        x = Sqrt[1-c^2/a^2];
        k = Sqrt[1 - c^2/b^2]/x;
        2*Pi*(c^2 + a*b*x*EllipticE[x,k] + b*c^2/a/x*EllipticF[x,k])
        )

(* http://www-elsa.physik.uni-bonn.de/~dieckman/SurfaceEllipsoid/SurfEll.html *)
AreaEl3[c_,b_,a_] := (
        t = ArcCos[c/a];
        s = ArcCos[c/b];
        k = Sin[s]/Sin[t];
        2 Pi c^2 + 2 Pi a b/Sin[t] (Cos[t]^2 * EllipticF[t,k^2] + Sin[t]^2 * EllipticE[t,k^2])
        )


Etot[a_,b_,c_,Bm_,mu_] := Apply[AreaEl3,Sort[{a,b,c}]] + Em[a,b,c,Bm,mu]

EtotDim[a_?NumericQ,b_?NumericQ,Bm_,mu_] := (Print[{a,b}];(Apply[AreaEl,Sort[{a,b,1/a/b}]] + Em[a,b,1/a/b,Bm,mu])/4/Pi)


mu = 10;

Bm = {}
ac = {}
bc = {}
ab = {}
En = {}

Bmi = 1;
ai = 1.001;
bi = 1.002;

Do[
        Print[Bmi];
        result = FindMinimum[{EtotDim[a,b,Bmi,mu],0.1<a<=b<100.},{a,ai},{b,bi},MaxIterations->10];
        AppendTo[Bm,Bmi]
        AppendTo[En,result[[1]]];
        AppendTo[ac,a /. result[[2]]];
        AppendTo[bc,b /. result[[2]]];
        AppendTo[ab,a/b /. result[[2]]];
        ai = a /. result[[2]]; bi = b /. result[[2]],
                {Bmi,1,30,1}
]

fig = ListPlot[Transpose[{Bm,ab}],PlotRange->All]
Export["minimalelipsoid.pdf",fig]

(* Testing formula for surface area *)

k1[a_] := ArcSin[Sqrt[1 - a^2]];
m1[a_, b_] := (1 - b^2)/(1 - a^2);
ss[a_, b_] := (a*b)^(2/3)*(1 + 
    a/b/Sqrt[1 - a^2]*(Re[EllipticE[k1[a], m1[a, b]]]*(1/a^2 - 1) + 
       Re[EllipticF[k1[a], m1[a, b]]]))

R0 = 1.;
a = 3.;
b = 8.;
c = R0^3/a/b;

ss[c/a,c/b]*2*Pi
AreaEl3[c,a,b]

(* Testing of energy Formula *)
nz[a_, b_] := 
 nz[a, b] = 
  NIntegrate[
   1/2/((u + 1)^(3/2)*Sqrt[(1 + u^2*a^2)*(1 + u^2*b^2)]), {u, 0, 
    Infinity}, PrecisionGoal -> 10]
nx[a_, b_] := nz[1/a, b/a];
ny[a_, b_] := nz[a/b, 1/b];

enrot[a_, b_, Bm_, chi_] := 
 Which[a == b == 1, 
  2. - Bm*chi/6*(1/(1 + 4*Pi*chi/3) + 1/(1 + 4*Pi*chi/3)), True, 
  ss[a, b] - 
   Bm*chi/6*(1/(1 + 4*Pi*chi*nx[a, b]) + 1/(1 + 4*Pi*chi*ny[a, b]))]


R0 = 1.;
a = 1.;
b = 2.;
c = R0^3/a/b;

enrot[c/a,c/b,10,(10-1)/4/Pi]/2

EtotDim[a,b,10,10]

(* Testing demagnetization coeficients *)

R0 = 1.;
a = 1.;
b = 2.;
c = R0^3/a/b;

nz[c/a,c/b]

nz[a,b,c]

nx[0.2,0.1] + ny[0.2,0.1] + nz[0.2,0.1] 
