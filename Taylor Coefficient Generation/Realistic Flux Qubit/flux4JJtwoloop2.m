(* ::Package:: *)

flux4JJtwoloop[qbps1_,fluxvals_,LNnlevs_,L0nvals_]:=
Module[
{Jc,JJAs1,ICa1,ICb1,EJa1,EJb1,CJa1,CJb1,Cmat,JJmat0,gaugeflux,fluxvars,Lbasis,Lmat,islands,L1partition,LNpartition,L1nlevsJJ,out1,out0,fE1,fD1,Csh1,Lq1,LD1},
{Jc,JJAs1,{Csh1},{Lq1,LD1}}=qbps1;
ICa1=Jc*JJAs1[[1]]*10^-6;
ICb1=Jc*JJAs1[[2]]*10^-6;
EJa1=phi0*ICa1/(2*Pi)/h/10^9;
EJb1=phi0*ICb1/(2*Pi)/h/10^9;
CJa1=Sc*JJAs1[[1]]*1000;
CJb1=Sc*JJAs1[[2]]*1000;

Cmat=({
 {0, Csh1, CJb1, 0, 0},
 {Csh1, CJb1, 0, CJa1, CJa1},
 {CJb1, 0, 0, 0, 0},
 {0, CJa1, 0, 0, 0},
 {0, CJa1, 0, 0, 0}
});
JJmat0=({
 {0, 0, EJb1, 0, 0},
 {0, EJb1, 0, EJa1/2, EJa1/2},
 {EJb1, 0, 0, 0, 0},
 {0, EJa1/2, 0, 0, 0},
 {0, EJa1/2, 0, 0, 0}
});
gaugeflux=({
 {0, 0, 0, 0, 0},
 {0, 0, 0, fE1, fE1+fD1},
 {0, 0, 0, 0, 0},
 {0, fE1, 0, 0, 0},
 {0, fE1+fD1, 0, 0, 0}
});
fluxvars={fE1,fD1};
Lbasis=({
 {-1, 0, 0, 0, 1},
 {-1, 0, 0, 1, 0},
 {0, 0, 1, 0, 0}
}); (* matrix whose rows are phases across linear inductors in node basis *)
Lmat=({
 {LD1/2, 0, 0},
 {0, LD1/2, 0},
 {0, 0, Lq1}
});  (* matrix whose diagonal entries are linear inductances and off-diagonals are mutual inductances *)
islands={};

L1partition={{"Q1",{1,2,3,4,5}}}; (* partition of circuit into sub-circuits for hierarchical diagonalization *)
LNpartition={};
L1nlevsJJ={100};


out0=JJcircuitFL0[Cmat,Lmat,Lbasis,L1partition,islands,JJmat0,gaugeflux];
out1=JJcircuitFLN[out0,LNpartition,fluxvars,fluxvals,L0nvals,L1nlevsJJ,{{LNnlevs}},False,False];
out1[[1;;2]]
]
