(* 
# Principal symbol of the three-fold interactions 

To run this in a notebook use  
```
Get[NotebookDirectory[] <> "prin_symb.wl"]
```


## Signature of the Minkowski metric 
*)
m = {-1, 1, 1, 1}; 
(* 


## Representation of the Clifford algebra 

*)
sigmas = {
   {{0, 1},
    {1, 0}},
   
   {{0, -I},
    {I, 0}},
   
   {{1, 0},
    {0, -1}}
   };
(* 

Gamma matrices with lower indices 
*)
Gammas = {ArrayFlatten[
     I {{0, IdentityMatrix[2]},
       {IdentityMatrix[2], 0}}
     ]}~Join~Table[ArrayFlatten[
     I {{0, sigmas[[k]]},
       {-sigmas[[k]], 0}}
     ], {k, 1, 3}];
(* 


## Incoming and outgoing directions 

All the directions are defined with lower indices

The eta and xi directions are given in Section 4.3 
*)
a[r_] := Sqrt[1 - r^2];
eta = {1, -a[r], r, 0};
xis = {
       {1, 1, 0, 0},
       {1, a[s], s, 0},
       {1, a[s], -s, 0}
   };
(* 

The directions eta in Section 5.2 
*)
ks = LinearSolve[Transpose[xis], eta];
etas = DiagonalMatrix[ks] . xis;
(* 


## Incoming principal symbols 

We ignore the factors on V that don't play a role in the computation.
In other words, we drop the b factors.

The omega 1-forms are given in Section 5.2
*)
omegas = {
       {0, 0, 1, 0},
       {s, 0, 1, 0},
       {-s, 0, 1, 0}
   };
(* 

As factors on V are ignored, Ij is treated as a vector in Delta = C^4.
I2 and I3 are I1 in the leading order 
*)
I1 = {I10, I11, I12, I13};
I2p = {I2p0, I2p1, I2p2, I2p3};
I3p = {I3p0, I3p1, I3p2, I3p3};
I2pp = {I2pp0, I2pp1, I2pp2, I2pp3};
I3pp = {I3pp0, I3pp1, I3pp2, I3pp3};
Is = {I1, I1 + s I2p + s^2 I2pp, I1 + s I3p + s^2 I3pp};
(* 

Formula for the incoming symbol hat phi as in Section 5.2.2 
*)
phis = Table[
   -1/2 Sum[Gammas[[a]] m[[a]] omegas[[j, a]], {a, 1, 4}] .
     Sum[Gammas[[b]] m[[b]] xis[[j, b]], {b, 1, 4}] . Is[[j]]
   , {j, 1, 3}];
(* 


## Twofold interactions

Implemented using the formulas in Section 5.2.1. 
*)
N2IdPrimitive[j_, k_] := Sum[
  2 omegas[[j, a]] m[[a]] I etas[[k, a]]
  , {a, 1, 4}]
N2GammaPrimitive[j_] := Sum[
  Gammas[[a]] . Gammas[[b]] m[[a]] m[[b]] I etas[[j, a]] omegas[[j, b]]
  , {a, 1, 4}, {b, 1, 4}]
N2[j_, k_] :=
 N2IdPrimitive[j, k] phis[[k]] + N2IdPrimitive[k, j] phis[[j]] +
  N2GammaPrimitive[j] . phis[[k]] + N2GammaPrimitive[k] . phis[[j]]

sigmaBox[j_, k_] := With[{e = etas[[j]] + etas[[k]]},
  Sum[e[[a]] m[[a]] e[[a]], {a, 1, 4}]]
phi2[j_, k_] := 1/sigmaBox[j, k] N2[j, k]
(* 


## Threefold interactions

Implemented using the formulas in Section 5.2.2.

Sum over S3 
*)
SumS3[s_] := Total[Map[Function[pi,
    s[Splice[pi, _]]],
   Permutations[Range[3]]]]
(* 

Threefold 1+1+1 interactions 
*)
N3111Primitive[j_, k_, l_] :=
 Sum[m[[a]] omegas[[j, a]] omegas[[k, a]] phis[[l]], {a, 1, 4}]
N3111 = SumS3[N3111Primitive];
(* 

Threefold wedge interactions 
*)
N3wedgePrimitive[j_, k_, l_] := With[{e = etas[[k]] + etas[[l]]},
  Sum[m[[a]] omegas[[j, a]] I e[[a]] phi2[k, l], {a, 1, 4}]
  ]
N3wedge = SumS3[N3wedgePrimitive];
(* 

Threefold bullet interactions 
*)
N3bulletPrimitive[j_, k_, l_] := Sum[
  1/2 m[[a]] m[[b]] I etas[[j, a]] omegas[[j, b]] *
   Gammas[[a]] . Gammas[[b]] . phi2[k, l]
  , {a, 1, 4}, {b, 1, 4}]
N3bullet = SumS3[N3bulletPrimitive];

N3 = N3111 + N3wedge + N3bullet;
(* 


## Verification of the formulas for N_123^D in Section 5.2.2 

Gamma matrices with upper indices 
*)
G0 = m[[1]] Gammas[[1]];
G1 = m[[2]] Gammas[[2]];
G2 = m[[3]] Gammas[[3]];

Print["1. Should give 0: ", Norm[Normal[Series[
   -4/(3 s^2 r) N3
    - (
     1/r (G0 - G1) . (G2 . I1 + G1 . (I2p - I3p))
      + I1 - 1/2 (G0 + G1) . G2 . (I2p - I3p)
     ), {s, 0, 0}, {r, 0, 0}]]]]

Print["2. Should give 0: ", Norm[Normal[Series[
   -4/(3 s^2 r) (G0 + G1) . (N3 - (N3 /. r -> 0))
    -
    (G0 + G1) . I1
   , {s, 0, 0}, {r, 0, 0}]]]]
