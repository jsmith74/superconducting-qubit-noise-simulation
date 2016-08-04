(* Mathematica Package *)

BeginPackage["JJcircuit`"]

JJcircuitFL0::usage = "calculates Hamiltonian for JJ circuits"
JJcircuitFLN::usage = "diagonalizes full JJ circuit"
outinfo::usage="displays info about Hamiltonian construction from JJcircuitFL0"
islcutF::usage="removes states with island charges other than zero from full system"
Vdiptot::usage="gives the total effective Voltage decay matrix element from a given exc state assuming equal, independent charge noise at each physical circuit node"
Idiptot::usage="gives the total effective current decay matrix element from a given exc state assuming equal, independent flux noise in each inductor listed in Lbasis"
decFV::usage="lists the matrix elements between initstate and each other state, for each of the node voltage operators; returns pairs of frequency in GHz, and list of matrix elements"
decFI::usage="lists the matrix elements between initstate and each other state, for each of the inductor current operators specified by Lbasis; returns pairs of frequency in GHz, and list of matrix elements"


phi0 = 2.07 10^-15;
h = 6.626 10^-34;
e0 = 1.602 10^-19;	
JJfixed = 0.2;
Sc = 60 10^-3;
L0 = 10^-12;
C0 = 10^-15;
nu0 = 10^9;

Begin["`Private`"] 

Needs["JJfuncts`"]


(* solves individual level 1 systems (green blocks) composed of coupled physical (L0) variables *)
JJcircuitFL1[info0_,fluxvars_,fluxvals_,L0nvals_,L1nlevsJJ_,LNnlevs_,prediagin_]:=Module[
	{vecp,subvecp,Dsubvecp,qsubvecp,fsubvecp,qvec,qvecp,fvec,fvecp,chopvec,subsinfo,subdims,zerop,outHsubs,flag,subN,Rmat,subECs1,subELs1,preinfo,
	subH,subrepsF,nsubs,subvals,subvecs,subvals0,subvecs0,fluxreps,subHNO,Htot,HtotNO,varnames,oscinfo,prediag,nlevssub,L1partition,JJmat0,gaugeflux,Cinv,Linv},
	nlevssub=LNnlevs[[1]];
	
	{Htot,varnames,JJmat0,gaugeflux,Cinv,Linv,preinfo,oscinfo}=info0;
	{L1partition,Rmat,subECs1,subELs1}=preinfo;
	
	nsubs=Length[L1partition];
	
	If[Length[fluxvars]!=Length[fluxvals],Print["not all flux values specified"];Abort[]];
	fluxreps=Table[fluxvars[[j]]->fluxvals[[j]],{j,1,Length[fluxvals]}];
	If[Apply[And,Map[NumericQ,DeleteCases[Flatten[gaugeflux/.fluxreps],0]]],Null,Print["not all flux values specified"];Abort[]];

(* checks properties of L1partition *)
	Do[If[Length[L0nvals[[subN]]]<Length[L1partition[[subN,2]]],Print["dimensions not specified for all of subsystem ",subN," variables"];Abort[],Null];
	If[Length[L0nvals[[subN]]]>Length[L1partition[[subN,2]]],Print["too many dimensions specified for subsystem ",subN];Abort[],Null],{subN,1,nsubs}];

	fluxreps=Table[fluxvars[[j]]->fluxvals[[j]],{j,1,Length[fluxvals]}];
	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;
	
	subsinfo=Table[Table[varinfo[subvecp[[subN]][[i]],L1partition][[2]],{i,1,Length[subvecp[[subN]]]}],{subN,1,nsubs}]; (* lists which type of coordinate each one is: 1:oscillator, 2:JJ, 3: island *)
	subdims=Table[Table[Which[subsinfo[[subN,i]]==1,L0nvals[[subN]][[i]],subsinfo[[subN,i]]==2,2*L0nvals[[subN]][[i]]+1,subsinfo[[subN,i]]==3,2*L0nvals[[subN]][[i]]+1],{i,1,Length[subsinfo[[subN]]]}],{subN,1,nsubs}];(* resulting dimension of each subspace *)

	prediag=Which[Length[prediagin]==0,Table[prediagin,{subN,1,nsubs}],Length[prediagin]==nsubs,prediagin,Length[prediagin]!=nsubs,Null];
	If[prediag==Null,Print["Incorrect number of prediag values for number of subsystems"];Abort[]];
	flag=Table[And[prediag[[subN]],MemberQ[subsinfo[[subN]],2]],{subN,1,nsubs}]; (* flag==True means there is at least one JJ degree of freedom, so that initial diagonalization of only these will be carried out *)

	If[Apply[Or,flag],HtotNO=L0prediagH[JJmat0,gaugeflux,varnames,preinfo],Null]; (* if prediagonalization is to be used for any subsystem, calculate the appropriate Hamiltonians without oscillators *)

	zerop=Table[0,{k,1,nsubs}];
	subvals0=zerop;subvecs0=zerop;subvals=zerop;subvecs=zerop;
	
	For[subN=1;outHsubs=zerop;subrepsF=zerop,
		subN<=nsubs,
		subN++,
		subH=Chop[matrixconv[Simplify[Chop[Htot[[subN,subN]]]],varnames]/.fluxreps];
		(* for each subsystem, choose whether to prediagonalize or not *)
		If[flag[[subN]],
			subHNO=Chop[matrixconv[Simplify[Chop[HtotNO[[subN]]]],varnames]/.fluxreps];
			{subrepsF[[subN]],outHsubs[[subN]]}=subHFP[L1partition,subHNO,subH,subN,varnames,L0nvals,Htot,oscinfo,subsinfo[[subN]],subdims[[subN]],L1nlevsJJ[[subN]]];
			,
			{subrepsF[[subN]],outHsubs[[subN]]}=subHF[L1partition,subH,subN,varnames,L0nvals,Htot,oscinfo,subsinfo[[subN]],subdims[[subN]]];
		];
		
	];

(* diagonalize subsystem Hamiltonians, truncating at nlevssub *)
	Do[{subvals0[[subN]],subvecs0[[subN]]}=Eigensystem[outHsubs[[subN]],-nlevssub[[subN]]];{subvals[[subN]],subvecs[[subN]]}={Reverse[subvals0[[subN]]],Reverse[subvecs0[[subN]]]},{subN,1,nsubs}];
	{subvals,subvecs,subrepsF}
];


leviterateN[slev_,LNpartition_,LNnlevs_,Clev_,varnames_,nsys_,Htot_,info_,debug_]:=Module[
	{partitionsys,nlevssys,nlevssub,subopsH,HCsubreps3,zerop,sysvals3,sysvecs3,HCreps3,Hcouple,Hsystem,subpart,nsubs,IDmatssub,Hcmat,Hsubmats,sysvals0,sysvecs0,HCreps,sysvals,sysvecs},

	{sysvals,sysvecs,HCreps}=info;
	
	partitionsys=LNpartition[[slev-1,2]];
	
If[debug,Print[""];	
	Print["partitionsys: ",partitionsys];
	Print["nsys: ",nsys[[slev]]]];
	
	nlevssys=LNnlevs[[slev]];
	nlevssub=LNnlevs[[slev-1]];

	subopsH=subopsF[LNpartition,Clev,Htot,varnames,slev-1]; (* all of the required coupling operators indexed by L2 subsystem (green blocks) *)
	HCsubreps3=Table[Association[Table[subopsH[[subN,i]]->Chop[vecfunct[subopsH[[subN,i]],HCreps[[subN]][subopsH[[subN,i]]],sysvecs[[subN]],varnames,nlevssub[[subN]]]],{i,1,Length[subopsH[[subN]]]}]],{subN,1,nsys[[slev-1]]}];
	(* steps through all operators in coupling Hamiltonian, indexed by L2 subsystem (L1 systems, green blocks, typically qubits) and expresses each in the L2 subsystem eigenbasis subvecs *)

If[debug,Print["Keys[HCsubreps3]:",Keys[HCsubreps3]]];

	zerop=Table[0,{i,1,nsys[[slev]]}];
	sysvals3=zerop;sysvecs3=zerop;HCreps3=zerop;Hcouple=zerop;Hsystem=zerop;
	Hcouple=HCfunct[slev,LNpartition,Htot,varnames];

(* this loop constructs the coupling Hamiltonian between L2 subsystems (L1 systems, typically qubits), to diagonalize each L2 system (blue blocks) *)
	Do[
		subpart=partitionsys[[sysnum]]; (* set of subsystems that comprises system sysnum *)
		nsubs=Length[subpart]; (* number of subsystems that comprises system sysnum *)

If[debug,Print["slev= ",slev,", sysnum= ",sysnum," - nsubs= ",nsubs," subpart=",subpart]];

	(* if only a single subsystem, carry info forward from previous level *)
		If[nsubs==1,
			HCreps3[[sysnum]]=HCreps[[subpart[[1]]]];		
			{sysvals3[[sysnum]],sysvecs3[[sysnum]]}={sysvals[[subpart[[1]]]],sysvecs[[subpart[[1]]]]};

If[debug,Print["subpart[[1]]: ",subpart[[1]]]];
If[debug,Print["Keys[HCreps3[[sysnum]]]: ",Keys[HCreps3[[sysnum]]]]];	
			,
			IDmatssub=Table[IdentityMatrix[nlevssub[[subpart[[subN]]]],SparseArray],{subN,1,nsubs}]; (* identity matrices for L2 subsystems within sysnum^th L2 system (blue block), each corresponding to the chosen dimension of a green sub block *)

			(* goes through all operators in HCouple for syssnum^th L2 system (blue block), picks out L2 subsystem representation from HCsubreps, and takes appropriate Kronecker product with identity matrices for L2 subsystems within sysnum^th L2 system *)
			HCreps3[[sysnum]]=Association[Flatten[Table[Table[Keys[HCsubreps3[[subpart]]][[subN,i]]->Apply[KroneckerProduct,Table[If[j==subN,Values[HCsubreps3[[subpart]]][[subN,i]],IDmatssub[[j]]],{j,1,nsubs}]],{i,1,Length[Keys[HCsubreps3[[subpart]]][[subN]]]}],{subN,1,nsubs}]]];

If[debug,Print["Keys[HCreps3[[sysnum]]]: ",Keys[HCreps3[[sysnum]]]]];				

			Hcmat=Hcouple[[sysnum]]/.HCreps3[[sysnum]]; 

If[debug,Print["operators in Hcouple[[sysnum]]: ",opfind[Hcouple[[sysnum]],varnames]]];

			Hsubmats=Sum[Apply[KroneckerProduct,Table[If[i==subN,SparseArray[Band[{1,1}]->sysvals[[subpart[[i]]]]],IDmatssub[[i]]],{i,1,nsubs}]],{subN,1,nsubs}]; (* subsystem Hamiltonians in eigenbasis *)
			Hsystem=Hsubmats+Hcmat; (* Hamiltonian for sysnum^th L2 system (blue blocks) *)

			{sysvals0,sysvecs0}=Chop[Eigensystem[Hsystem,-nlevssys[[sysnum]]]];
			{sysvals3[[sysnum]],sysvecs3[[sysnum]]}={Reverse[sysvals0],Reverse[sysvecs0]}	
		];
		,
		{sysnum,1,nsys[[slev]]}
	];
	Chop[{sysvals3,sysvecs3,HCreps3}]
]

JJcircuitFL0[Cmat_,Lmat_,Lbasis_,L1partitionin_,islands_,JJmat0_,gaugeflux_]:=Module[
	{Rmat,L1partition,varnames,vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp,Htot,subH0oscs,subp0oscs,subx0oscs,wls,subJJpots,
	nsubs,subH0oscslist,ECmat,ELmat,ECoscmat,x0oscs,p0oscs,subECs0,subELs0,tempb,subECs1,subELs1,temp2b,temp3,subECs,subELs,JJpot,Cinv,Linv},
	
	Rmat=convertL[Lbasis,islands,gaugeflux];
	L1partition=convpart[Rmat,L1partitionin];
	nsubs=Length[L1partition];

	varnames=makevars[Lbasis,islands,L1partition];(* varnames={vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}*)
	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;

	Cinv=Inverse[convertCmat[Cmat]];
	Linv=Inverse[Lmat];
	ECmat=e0^2/2/(h*nu0)*FullSimplify[Rmat.Cinv.Transpose[Rmat]]/C0; (* in GHz, C values in fF *)
	ELmat=(phi0/(2*Pi))^2/(h*nu0)*FullSimplify[Linv]/L0; (* in GHz, L values in pH *)
	ECoscmat=Take[ECmat[[1;;Length[ELmat],1;;Length[ELmat]]]];
	
	wls=Table[N[Sqrt[8*ECoscmat[[i,i]]*ELmat[[i,i]]]],{i,1,Length[Lmat]}];
	x0oscs=Table[(ECoscmat[[i,i]]/(2*ELmat[[i,i]]))^(1/4)/Pi,{i,1,Length[Lmat]}]; (* plasma oscillator flux in units of phi/phi_0 - nonangular. phase is then 2pi*x*)
	p0oscs=Table[1/(2*Pi*x0oscs[[i]]),{i,1,Length[Lmat]}];(* plasma oscillator charge in units of 2e - nonangular. *)
	subp0oscs=Table[PadRight[p0oscs,Length[Cmat]][[L1partition[[i,2]]]],{i,1,Length[L1partition]}];(* padded with zeroes for non-oscillator variables *)
	subx0oscs=Table[PadRight[x0oscs,Length[Cmat]][[L1partition[[i,2]]]],{i,1,Length[L1partition]}]; (* padded with zeroes for non-oscillator variables *)

	subECs0=Table[ECmat[[L1partition[[k,2]],L1partition[[l,2]]]],{k,1,Length[L1partition]},{l,1,Length[L1partition]}];
	subELs0=Table[PadRight[ELmat,{Length[ECmat],Length[ECmat]}][[L1partition[[k,2]],L1partition[[l,2]]]],{k,1,Length[L1partition]},{l,1,Length[L1partition]}];

	tempb=Table[Table[varinfo[vecp[[L1partition[[pnum,2,i]]]],L1partition][[2]]==1,{i,1,Length[L1partition[[pnum,2]]]}],{pnum,1,Length[L1partition]}]; (* true or false for each variable, if it is or is not an oscillator *)
	subECs1=Table[If[pnum==pnump,ReplacePart[subECs0[[pnum,pnum]],Flatten[Table[If[tempb[[pnum,i]]==True,{i,i}->0,{}],{i,1,Length[tempb[[pnum]]]}],1]],subECs0[[pnum,pnump]]],{pnum,1,Length[L1partition]},{pnump,1,Length[L1partition]}];
	subELs1=Table[If[pnum==pnump,ReplacePart[subELs0[[pnum,pnum]],Flatten[Table[If[tempb[[pnum,i]]==True,{i,i}->0,{}],{i,1,Length[tempb[[pnum]]]}],1]],subELs0[[pnum,pnump]]],{pnum,1,Length[L1partition]},{pnump,1,Length[L1partition]}];

	subH0oscslist=Table[Table[If[tempb[[pnum,i]],ToExpression[StringJoin["H0",ToString[vecp[[L1partition[[pnum,2,i]]]]]]],0],{i,1,Length[L1partition[[pnum,2]]]}],{pnum,1,Length[L1partition]}];
	subH0oscs=DiagonalMatrix[Table[Total[subH0oscslist[[pnum]]],{pnum,1,Length[L1partition]}]];
	
	temp2b=Table[Table[varinfo[vecp[[L1partition[[pnum,2,i]]]],L1partition][[2]]!=2,{i,1,Length[L1partition[[pnum,2]]]}],{pnum,1,Length[L1partition]}];
	temp3=Table[Table[If[Or[temp2b[[i,j]],temp2b[[i,k]]],0,1],{j,1,Length[temp2b[[i]]]},{k,1,Length[temp2b[[i]]]}],{i,1,Length[L1partition]}];

	(* zero out oscillator and island variables *)
	subECs=Simplify[Chop[Table[4*qvecp[[L1partition[[i,2]]]].subECs1[[i,j]].qvecp[[L1partition[[j,2]]]],{i,1,Length[L1partition]},{j,1,Length[L1partition]}]]]; (* charge in cooper pairs, and EC multiplies this directly. No factor of 1/2 because EC already includes this *)
	subELs=Simplify[Chop[Table[(4*Pi^2)/2*fvecp[[L1partition[[i,2]]]].subELs1[[i,j]].fvecp[[L1partition[[j,2]]]],{i,1,Length[L1partition]},{j,1,Length[L1partition]}]]]; (* 4Pi^2 is because flux is defined in units of phi0, but EL multiplies (angular phase)^2; 1/2 because EL does not include this *)
	JJpot=JJHfunct[JJmat0,gaugeflux,Rmat,varnames,L1partition,True];
	subJJpots=partitionJJmat[L1partition,JJpot,chopvec,Global`IDmat,gaugeflux];
	Htot=Chop[Simplify[subECs+subELs+subJJpots+subH0oscs]];
	
	{Htot,varnames,JJmat0,gaugeflux,Cinv,Linv,{L1partition,Rmat,subECs1,subELs1},{subp0oscs,subx0oscs,wls,subH0oscslist,subH0oscs}}
]

(* displays info about Hamiltonian construction from JJcircuitFL0 *)  
outinfo[info0_]:=Module[
	{Htot,varnames,L1partition,subp0oscs,subx0oscs,wls,subH0oscslist,subH0oscs,gaugeflux,Rmat,JJmat0,subECs1,subELs1,vecp,subvecp,qvec,qvecp,
		fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp,Cinv,Linv},
	
	{Htot,varnames,JJmat0,gaugeflux,Cinv,Linv,{L1partition,Rmat,subECs1,subELs1},{subp0oscs,subx0oscs,wls,subH0oscslist,subH0oscs}}=info0;
    {vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;
    Print["Quantum variable definitions:"];
    Print[Table[{fvecp[[i]]->FullSimplify[(Rmat.fvec)[[i]]],qvecp[[i]]->FullSimplify[(Transpose[Inverse[Rmat]].qvec)[[i]]]},{i,1,Length[fvecp]}]];
]  


(* removes states with island charges other than zero from full system *)
islcutF[info0_,info1_,L1partition_]:=Module[
	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp,subN,islchg,islchgs,inds,indsout,outvals,outvecs,
		HCreps,subsinfo,nsubs,varnames,Htot,JJmat0,gaugeflux,Cinv,Linv,Rmat,subECs1,subELs1,subp0oscs,subx0oscs,wls,subH0oscslist,subH0oscs},
	
	{outvals,outvecs,HCreps}=info1;
	{Htot,varnames,JJmat0,gaugeflux,Cinv,Linv,{L1partition,Rmat,subECs1,subELs1},{subp0oscs,subx0oscs,wls,subH0oscslist,subH0oscs}}=info0;
	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;

	nsubs=Length[L1partition];
	subsinfo=Table[Table[varinfo[subvecp[[subN]][[i]],L1partition][[2]],{i,1,Length[subvecp[[subN]]]}],{subN,1,nsubs}]; (* lists which type of coordinate each one is: 1:oscillator, 2:JJ, 3: island *)

	islchg=Table[0,{i,1,nsubs}];
	islchgs=Table[0,{i,1,nsubs}];
	inds=Table[0,{i,1,nsubs}];

	Do[
		islchg[[subN]]=Pick[qsubvecp[[subN]],Table[subsinfo[[subN,i]]==3,{i,1,Length[subsinfo[[subN]]]}]] ;(* finds the island charge variables for each subsystem *)
		If[Length[islchg[[subN]]]>0,
			islchgs[[subN]]=Chop[Table[Table[Conjugate[outvecs[[j]]].HCreps[islchg[[subN,i]]].outvecs[[j]],{j,1,Length[outvals]}],{i,1,Length[islchg[[subN]]]}]];
			inds[[subN]]=DeleteCases[Table[If[FreeQ[Table[islchgs[[subN,i,j]]==0,{i,1,Length[islchgs[[subN]]]}],False],j,Null],{j,1,Length[islchgs[[subN,1]]]}],Null],(* list of state indices with zero charge on all islands *)
			inds[[subN]]=Table[i,{i,1,Length[outvecs[[subN]]]}](* all states chosen if no islands in subsystem *)
		];
		,{subN,1,nsubs}
	];
	indsout=Apply[Intersection,inds];
	Chop[{outvals[[indsout]],outvecs[[indsout]]}]
]

(* lists the matrix elements between states i and j of all of the node charge operators *)
VM2[info0_,info1_,i_,j_]:=Module[
 	{vals,vecs,Htot,varnames,JJmat0,gaugeflux,L1partition,subECs1,subELs1,subp0oscs,subx0oscs,wls,subH0oscslist,subH0oscs,
 		vecp,subvecp,qvec,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp,reps,Rmat,qvecp,Cinv,Linv},

  	{Htot,varnames,JJmat0,gaugeflux,Cinv,Linv,{L1partition,Rmat,subECs1,subELs1},{subp0oscs,subx0oscs,wls,subH0oscslist,subH0oscs}}=info0;
  	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;
  	{vals,vecs,reps}=info1;

  	{vals[[i]]-vals[[j]],Table[Abs[2*e0/C0*Conjugate[vecs[[i]]].((Cinv.Transpose[Rmat].qvecp)[[k]] /. reps).vecs[[j]]],{k,1,Length[qvecp]}]}
]

(* lists the matrix elements between states i and j of all of the inductor flux operators specified by Lbasis *)
IM2[info0_,info1_,i_,j_]:=Module[
 	{vals,vecs,Htot,varnames,JJmat0,gaugeflux,L1partition,subECs1,subELs1,subp0oscs,subx0oscs,wls,subH0oscslist,subH0oscs,
 		vecp,subvecp,qvec,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp,reps,Rmat,qvecp,Cinv,Linv},

  	{Htot,varnames,JJmat0,gaugeflux,Cinv,Linv,{L1partition,Rmat,subECs1,subELs1},{subp0oscs,subx0oscs,wls,subH0oscslist,subH0oscs}}=info0;
  	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;
  	{vals,vecs,reps}=info1;

	If[Dimensions[Linv]!={1,1},
  		{vals[[i]]-vals[[j]],Table[Abs[phi0/L0*Conjugate[vecs[[i]]].((Linv.fvecp[[1;;Length[wls]]])[[k]]/.reps).vecs[[j]]],{k,1,Length[wls]}]},
  		{vals[[i]]-vals[[j]],{Abs[phi0/L0*Conjugate[vecs[[i]]].((Linv[[1,1]]*fvecp[[1]])/.reps).vecs[[j]]]}}
	]
]

(* gives the total effective Voltage decay matrix element from a given exc state assuming equal, independent charge noise at each physical circuit node *)
Vdiptot[info0_,info1_,excstate_]:=Sqrt[Sum[Norm[VM2[info0,info1,excstate,k][[2]]]^2,{k,excstate-1,1,-1}]]

(* gives the total effective current decay matrix element from a given exc state assuming equal, independent flux noise in each inductor listed in Lbasis *)
Idiptot[info0_,info1_,excstate_]:=Sqrt[Sum[Norm[IM2[info0,info1,excstate,k][[2]]]^2,{k,excstate-1,1,-1}]]


(* lists the matrix elements between initstate and each other state, for each of the node voltage operators; returns pairs of frequency in GHz, and list of matrix elements *)
decFV[info0_,info1_,initstate_]:=Module[
 	{vals,vecs,Htot,varnames,JJmat0,gaugeflux,L1partition,subECs1,subELs1,subp0oscs,subx0oscs,wls,subH0oscslist,subH0oscs,
 		vecp,subvecp,qvec,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp,reps,Rmat,qvecp,Cinv,Linv},

  	{Htot,varnames,JJmat0,gaugeflux,Cinv,Linv,{L1partition,Rmat,subECs1,subELs1},{subp0oscs,subx0oscs,wls,subH0oscslist,subH0oscs}}=info0;
  	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;
  	{vals,vecs,reps}=info1;

  	Chop[DeleteCases[Table[If[initstate!=j,{vals[[initstate]]-vals[[j]],Table[Abs[2*e0/C0*Conjugate[vecs[[initstate]]].((Cinv.Transpose[Rmat].qvecp)[[k]] /. reps).vecs[[j]]],{k,1,Length[qvecp]}]},Null],{j,1,Length[vals]}],Null]]
]

(* lists the matrix elements between initstate and each other state, for each of the inductor current operators specified by Lbasis; returns pairs of frequency in GHz, and list of matrix elements *)
decFI[info0_,info1_,initstate_]:=Module[
 	{vals,vecs,Htot,varnames,JJmat0,gaugeflux,L1partition,subECs1,subELs1,subp0oscs,subx0oscs,wls,subH0oscslist,subH0oscs,
 		vecp,subvecp,qvec,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp,reps,Rmat,qvecp,Cinv,Linv},

  	{Htot,varnames,JJmat0,gaugeflux,Cinv,Linv,{L1partition,Rmat,subECs1,subELs1},{subp0oscs,subx0oscs,wls,subH0oscslist,subH0oscs}}=info0;
  	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;
  	{vals,vecs,reps}=info1;

  	Chop[DeleteCases[Table[If[initstate!=j,{vals[[initstate]]-vals[[j]],Table[Abs[phi0/L0*Conjugate[vecs[[initstate]]].((Linv.fvecp[[1;;Length[wls]]])[[k]]/.reps).vecs[[j]]],{k,1,Length[wls]}]},Null],{j,1,Length[vals]}],Null]]
]

JJcircuitFLN[info0_,LNpartition_,fluxvars_,fluxvals_,L0nvals_,L1nlevsJJ_,LNnlevs_,prediag_,debug_]:=Module[
	{outN,Htot,varnames,nsys,slev,levsys,Clev,oscinfo,preinfo,JJmat0,gaugeflux,info1,out,Cinv,Linv},
	
	{Htot,varnames,JJmat0,gaugeflux,Cinv,Linv,preinfo,oscinfo}=info0;

	info1=JJcircuitFL1[info0,fluxvars,fluxvals,L0nvals,L1nlevsJJ,LNnlevs,prediag];
	
	If[Length[LNpartition]!=0,
		Clev=findClev[LNpartition,varnames,Htot];(*association which gives the maximum level at which each operator in Hcoupletot is nonzero*)
		levsys=Join[Table[Flatten[Transpose[LNpartition][[2,k]]],{k,1,Length[LNpartition]}],{{1}}];(*list of systems by level*)
		nsys=Table[Length[levsys[[i]]],{i,1,Length[levsys]}];(*list of number of systems at each level*)
		outN=Join[{info1},Table[0,{i,1,Length[LNpartition]}]];
		Do[
 			outN[[slev]]=leviterateN[slev,LNpartition,LNnlevs,Clev,varnames,nsys,Htot,outN[[slev-1]],debug],
 			{slev,2,Length[LNpartition]+1}
 		];
 		out=Extract[outN[[Length[outN]]],{{1,1},{2,1},{3,1}}]; (* take first system (whole system at top level) at Nth level, and elements 1 and 2 (eigenvalues and eigenvectors) *)
 		,
 		out=Extract[info1,{{1,1},{2,1},{3,1}}];
	];
 	Chop[out]
]


End[] 

EndPackage[]