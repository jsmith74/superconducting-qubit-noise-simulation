

BeginPackage["JJfuncts`"]

convertJJmat::usage=""
convertCmat::usage=""
checkpartition::usage=""
checkLinfo::usage=""
checkJJinfo::usage=""
checkJJdisp::usage=""
findothers::usage=""
convertL::usage=""
JJHfunct::usage=""
excludeF::usage=""
includeF::usage=""
includetwoF::usage=""
chopvecF::usage=""
zeroiseJJ::usage=""
partitionJJmat::usage=""
matrixconv::usage=""
varinfo::usage=""
oscpopF::usage=""
oscxopF::usage=""
qmat::usage=""
fmat::usage=""
qfmatS::usage=""
qDmat::usage=""
qDmatS::usage=""
suboscHmatS::usage=""
makevars::usage=""
convpart::usage=""
optype::usage=""
vecfunct::usage=""
checkptot::usage=""
subpartF::usage=""
findsysN::usage=""
findClev::usage=""
subopsF::usage=""
opfind::usage=""
HCfunct::usage=""
L0prediagH::usage=""
subHFP::usage=""
subHF::usage=""

Begin["`Private`"] (* Begin Private Context *) 

Needs["expmat`"]

checkpartition[partition_,chopvec_]:=
	If[
	Sort[Flatten[Transpose[partition][[2]]]]==Table[i,{i,1,Length[chopvec]}],
	Null,
	Print["Invalid partition (checkpartition)"];Abort[]
] (* validates that the specified circuit partition contains all variables *)

checkLinfo[Lbasis_,Lmat_]:=
	If[
	And[FreeQ[Table[Or[And[Count[Lbasis[[i]],1]==1,Count[Lbasis[[i]],-1]==1,Count[Lbasis[[i]],0]==Length[Lbasis[[i]]]-2],And[Count[Lbasis[[i]],1]==1,Count[Lbasis[[i]],0]==Length[Lbasis[[i]]]-1]],{i,1,Length[Lbasis]}],False],Lmat==Transpose[Lmat]],
	Null,
	Print["Invalid partition (checkLinfo)"];Abort[]
] (* validates that Lmat and Lbasis are correctly specified *)

checkJJinfo[JJmat0_,gaugeflux_]:=
	If[
	And[JJmat0==Transpose[JJmat0],gaugeflux==Transpose[gaugeflux],Length[JJmat0]==Length[gaugeflux],Length[JJmat0]==Length[JJmat0[[1]]],Length[gaugeflux]==Length[gaugeflux[[1]]]],
	Null,
	Print["Invalid partition (checkJJinfo)"];Abort[]
] (* validates that JJmat0 and gaugeflux are correctly specified *)

checkJJdisp[JJpot_,chopvec_]:=Module[{list,inds,nonoscinds}, (* checks that only oscillator variables have nonunity charge displacements in chosen representation *)
	list=Table[Position[JJpot,chopvec[[k]]],{k,1,Length[chopvec]}];
	inds=Table[Table[ReplacePart[list[[k,i]],Length[list[[k,i]]]->list[[k,i,Length[list[[k,i]]]]]+1],{i,1,Length[list[[k]]]}],{k,1,Length[chopvec]}];
	nonoscinds=Drop[inds,Length[chopvec]];
	If[FreeQ[Flatten[Table[DeleteDuplicates[Abs[Table[Chop[Extract[JJpot,nonoscinds[[k,i]]]],{i,1,Length[nonoscinds[[k]]]}]]],{k,1,Length[chopvec]-Length[chopvec]}]],False],Null,Print["nonunity charge displacements for JJs or islands"];Abort[];];
]

(* validates specification of hierarchical partition *)
checkptot[LNpartition_, nlevstot_] := Module[
 	{temp},
  	temp=Transpose[LNpartition];
  	If[
   		And[Length[temp[[1]]]==Length[temp[[2]]],
    	Length[nlevstot]==Length[temp[[2]]],
    	Apply[And,Table[Length[nlevstot[[i]]]==Length[temp[[2,i]]],{i,1,Length[temp[[2]]]}]],
    	Apply[And,Table[Length[temp[[2, i]]]==Length[Flatten[temp[[2,i+1]]]],{i,1,Length[temp[[2]]]-1}]]],
   		Null,
   		Print["LNpartition invalid"];
   		Abort[]
   	];
 ]

convertJJmat[JJmat0_,gaugeflux_]:=Module[{temp,temp2,temp3},
	temp=UpperTriangularize[JJmat0,1]+Transpose[UpperTriangularize[JJmat0,1]];
	temp2=Total[temp];
	temp3=DiagonalMatrix[Diagonal[JJmat0]+temp2]-temp;
	temp3*Exp[2*Pi*I*gaugeflux]
];
				
convertCmat[JJmat0_]:=Module[{temp,temp2},
	temp=UpperTriangularize[JJmat0,1]+Transpose[UpperTriangularize[JJmat0,1]];
	temp2=Total[temp];
	DiagonalMatrix[Diagonal[JJmat0]+temp2]-temp
];

findothers[Lbasis_,temp1_]:=IdentityMatrix[Length[Lbasis[[1]]]][[Flatten[Position[Table[If[Total[Abs[Transpose[temp1][[i]]]]==0,1,0],{i,1,Length[temp1[[1]]]}],1]]]];  (* finds any coordinates not yet included *)

convertL[Lbasis_,islands_,gaugeflux_]:=Module[{temp1,qvec,fvec,LbasisR,numLs,leftqs,Cvec,lincomb,islq,islvars,Rislmat,Rmat,zvec,inds,newrow}, 
	(* generates basis rotation matrix acting on fvec which transforms to coordinates including inductor branch fluxes and island charges *)
	fvec=Table[ToExpression[StringJoin["f",ToString[i]]],{i,1,Length[Lbasis[[1]]]}];
	qvec=Table[ToExpression[StringJoin["q",ToString[i]]],{i,1,Length[Lbasis[[1]]]}];
	temp1=Lbasis;
(* adds a row to temp1 corresponding to the average phase across any inductor that is not connected to any other inductor *)
	Do[If[And[{Count[Lbasis[[i]],1],Count[Lbasis[[i]],-1]}=={1,1},Drop[Transpose[temp1][[Position[Lbasis[[i]],1][[1,1]]]],{i}]===Table[0,{p,1,Length[temp1]-1}],Drop[Transpose[temp1][[Position[Lbasis[[i]],-1][[1,1]]]],{i}]===Table[0,{p,1,Length[temp1]-1}]],temp1=Join[temp1,{Abs[Lbasis[[i]]]/2}],Null],{i,1,Length[Lbasis]}]; 
 
 (* adds rows to Rmat for any nodes that are completely absent in Lbasis *)
	temp1=Join[temp1,findothers[Lbasis,Lbasis]];

(* if still one row missing, pick a single node which makes Rmat the correct Rank *)
	If[Length[temp1]==Length[temp1[[1]]]-1,
		zvec=Table[0,{i,1,Length[temp1[[1]]]}];
		inds=Flatten[Position[Table[MatrixRank[Join[temp1,{ReplacePart[zvec,i->1]}]],{i,1,Length[temp1[[1]]]}],Length[temp1[[1]]]]];
		newrow=Pick[inds,Table[Total[gaugeflux][[inds[[i]]]]===zvec[[i]],{i,1,Length[inds]}]][[1]];
		LbasisR=Join[temp1,{ReplacePart[zvec,newrow->1]}],
		LbasisR=temp1
	];
(* if not a square matrix, abort *)
	If[Length[LbasisR]!=Length[LbasisR[[1]]],Print["problem with convertL"];Abort[],Null];

(* if there are islands in the circuit, rotate non-oscillator sub-basis to include island charge variables *)
	If[Length[islands]!=0,
		numLs=Length[Lbasis]; (* number of oscillator variables *)
		leftqs=Drop[Inverse[Transpose[LbasisR]].qvec,numLs]; (* non-oscillator variables *)
		Cvec=Table[ToExpression[StringJoin["C",ToString[i]]],{i,1,Length[leftqs]}];
		lincomb=leftqs.Cvec;
		islq=islands.qvec;
		islvars=Table[PadLeft[Cvec/.Table[D[lincomb->islq[[k]],qvec[[i]]],{i,1,Length[qvec]}],Length[qvec]],{k,1,Length[islq]}]; (* solve for transformation which results in island charge variables *)
		Rislmat=ReplacePart[IdentityMatrix[Length[qvec]],Table[i+numLs->islvars[[i]],{i,1,Length[islq]}]]; (* construct square matrix which transforms previous variables to new *)
		Rmat=Inverse[Transpose[Rislmat]].LbasisR,
		Rmat=LbasisR
	];
	Rmat
]

JJHfunct[JJmat0_,gaugeflux_,Rmat_,varnames_,partition_,flag_]:=Module[{eat4,fvec,qvec,fvecp,qvecp,rep,chgdisp,eat2,eat3,JJpot,JJmat,eat1,
list,inds,temp2,temp,reps,elist,elist2,eat5,dlist,dlist2,chopvec,Dsubvecp,qsubvecp,fsubvecp,vecp,subvecp},
	checkJJinfo[JJmat0,gaugeflux];
	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;

	rep=Table[fvec[[i]]->(Inverse[Rmat].fvecp)[[i]],{i,1,Length[JJmat0]}];
	JJmat=Chop[JJmat0*Exp[2*Pi*I*gaugeflux]];

	chgdisp[in_]:=Module[{temp3},
		temp3=Table[D[in,fvecp[[i]]],{i,1,Length[JJmat0]}];
		Table[If[temp3[[i]]==0,1,Exp[I*temp3[[i]]*fvecp[[i]]]],{i,1,Length[JJmat0]}]
	];

	eat1=Table[If[i==j,(chgdisp[fvec[[i]]/.rep]),chgdisp[fvec[[i]]/.rep]*chgdisp[-fvec[[j]]/.rep]],{i,1,Length[JJmat]},{j,i,Length[JJmat]}];
	list=Position[eat1,E];
	inds=Table[ReplacePart[list[[i]],Length[list[[i]]]->list[[i,Length[list[[i]]]]]+1],{i,1,Length[list]}];
	temp=Table[Chop[Extract[eat1,inds[[i]]]],{i,1,Length[inds]}];
	temp2=Table[Chop[Cases[Table[D[temp[[i]],fvecp[[k]]],{i,1,Length[temp]}],Except[_Integer]]],{k,1,Length[fvecp]}]; (* reps: extracting exponents of displacement operators to replace them with symbolic objects chopvec[disp] *)
	reps=DeleteDuplicates[Chop[Flatten[Table[Table[Exp[temp2[[k,i]]*fvecp[[k]]]->If[Or[flag==True,varinfo[vecp[[k]],partition][[2]]==2],chopvec[[k]][Chop[-I*temp2[[k,i]]]],1],{i,1,Length[temp2[[k]]]}],{k,1,Length[fvecp]}]]]];
	eat2=eat1/.reps;
	eat3=Chop[Table[If[i==j,(JJmat[[i,j+i-1]]*Product[eat2[[i,j,k]],{k,1,Length[Rmat]}])/2,JJmat[[i,j+i-1]]/2*Product[eat2[[i,j,k]],{k,1,Length[Rmat]}]],{i,1,Length[JJmat]},{j,1,Length[JJmat0]-i+1}]];
	eat4=Chop[Total[Flatten[eat3]]]; (* one half of Hamiltonian *)

	elist=Position[eat4,E]; (* this part builds the complex conjugate of eat4 *)
	elist2=Table[ReplacePart[elist[[i]],Length[elist[[i]]]->elist[[i,Length[elist[[i]]]]]+1],{i,1,Length[elist]}];
	eat5=eat4/.Table[Extract[eat4,elist2[[i]]]->ComplexExpand[Conjugate[Extract[eat4,elist2[[i]]]]],{i,1,Length[elist2]}];
	dlist=Flatten[Table[Position[eat4,chopvec[[i]]],{i,1,Length[chopvec]}],1];
	dlist2=Table[ReplacePart[dlist[[i]],Length[dlist[[i]]]->dlist[[i,Length[dlist[[i]]]]]+1],{i,1,Length[dlist]}];
	eat5=ReplacePart[eat5,Table[dlist2[[i]]->-Extract[eat4,dlist2[[i]]],{i,1,Length[dlist2]}]];
	JJpot=Chop[Expand[-eat4-eat5]];
	checkJJdisp[JJpot,chopvec];
	JJpot
]

excludeF[JJpot_,chopvec_,exclude_]:=Module[{termpos,allterms,exterms}, (* excludes terms from Josephson potential involving the charge operators in the list exclude *)
	termpos=Quiet[Table[Transpose[Position[JJpot,chopvec[[i]]]][[1]],{i,1,Length[chopvec]}]];
	allterms=Table[i,{i,1,Length[JJpot]}];
	exterms=Complement[allterms,Flatten[termpos[[exclude]]]];
	JJpot[[exterms]]
]

includeF[JJpot_,chopvec_,include_]:=Module[{interms,termpos,allchops,exclude,exterms,outterms}, (* extracts terms from Josephson potential involving ONLY the charge operators in the list include *)
	allchops=Table[i,{i,1,Length[chopvec]}];
	termpos=Quiet[Table[Transpose[Position[JJpot,chopvec[[i]]]][[1]],{i,1,Length[chopvec]}]];
	interms=Flatten[termpos[[include]]]; (* terms in expression which have one of the specified operators *)
	exclude=Complement[allchops,include]; (* operators to be excluded *)
	exterms=Flatten[termpos[[exclude]]]; (* terms in expression which have one of the non-specified operators *)
	outterms=Complement[interms,exterms];
	JJpot[[outterms]]
]

includetwoF[JJpot_,chopvec_,include1_,include2_]:=Module[{interms1,termpos,interms2,outterms}, (* extracts terms from Josephson potential involving ONLY the charge operators in the list include *)
	termpos=Quiet[Table[Transpose[Position[JJpot,chopvec[[i]]]][[1]],{i,1,Length[chopvec]}]];
	interms1=Flatten[termpos[[include1]]]; (* terms in expression which have one of the specified operators in set 1 *)
	interms2=Flatten[termpos[[include2]]]; (* terms in expression which have one of the specified operators in set 2 *)
	outterms=Intersection[interms1,interms2];
	JJpot[[outterms]]
]

chopvecF[Lbasis_,islands_,partition_]:=Module[{temp,temp2,inds,vecp,subvecp}, (* creates labels for rotated coordinates according to subsystem label and variable type *)
	temp=Table[StringJoin[Apply[Which,Flatten[Table[{MemberQ[partition[[i,2]],inval],partition[[i,1]]},{i,1,Length[partition]}]]],Which[inval<=Length[Lbasis],"o",0<inval-Length[Lbasis]<=Length[islands],"i",inval>Length[Lbasis]+Length[islands],"j"]],{inval,1,Length[Lbasis[[1]]]}];
	temp2=Pick[Transpose[Tally[temp]][[1]],Table[Transpose[Tally[temp]][[2]][[i]]>1,{i,1,Length[Transpose[Tally[temp]][[2]]]}]];
	inds=Table[Flatten[Position[temp,temp2[[k]]]],{k,1,Length[temp2]}];
	vecp=ToExpression[ReplacePart[temp,Flatten[Table[inds[[i,j]]->StringJoin[temp[[inds[[i,j]]]],ToString[j]],{i,1,Length[inds]},{j,1,Length[inds[[i]]]}]]]];
	subvecp=Table[vecp[[partition[[i,2]]]],{i,1,Length[partition]}];
	{vecp,subvecp}
]

zeroiseJJ[JJpot_,chopvec_,gaugeflux_]:=Module[{list,inds,displist,nodisp,fluxes,fzero,fluxes0,JJpotin}, (* adds back energy offset to make minimum Josephson potential zero *)
	JJpotin=ReplacePart[JJpot,Table[{i,1}->Abs[JJpot[[i,1]]],{i,1,Length[JJpot]}]];
	list=Table[Position[JJpotin,chopvec[[k]]],{k,1,Length[chopvec]}];
	inds=Table[Table[ReplacePart[list[[k,i]],Length[list[[k,i]]]->list[[k,i,Length[list[[k,i]]]]]+1],{i,1,Length[list[[k]]]}],{k,1,Length[chopvec]}];
	displist=Table[Table[Chop[Extract[JJpotin,inds[[k,i]]]],{i,1,Length[inds[[k]]]}],{k,1,Length[chopvec]}];
	nodisp=Flatten[Table[If[Length[displist[[i]]]>0,{chopvec[[i]][displist[[i,j]]]->1,chopvec[[i]][-displist[[i,j]]]->1},Null],{i,1,Length[chopvec]},{j,1,Length[displist[[i]]]}]];
	fluxes0=Cases[Flatten[UpperTriangularize[gaugeflux]],Except[0]];
	fluxes=Table[If[Length[fluxes0[[i]]]==2,Extract[fluxes0[[i]],{2}],fluxes0[[i]]],{i,1,Length[fluxes0]}];
	fzero=Join[nodisp,Table[fluxes[[i]]->0,{i,1,Length[fluxes]}]];
	JJpotin /. fzero
]

partitionJJmat[partition_,JJpot_,chopvec_,IDmat_,gaugeflux_]:=Module[{offsets,subJJpots}, (* splits up Josephson potential Hamiltonian using specified partition*)
	checkpartition[partition,chopvec];
	offsets=Table[zeroiseJJ[If[k==l,includeF[JJpot,chopvec,partition[[k,2]]],includetwoF[JJpot,chopvec,partition[[k,2]],partition[[l,2]]]],chopvec,gaugeflux]*IDmat,{k,1,Length[partition]},{l,1,Length[partition]}];
	subJJpots=Table[If[k==l,includeF[JJpot,chopvec,partition[[k,2]]]+offsets[[k,l]],includetwoF[JJpot,chopvec,partition[[k,2]],partition[[l,2]]]+offsets[[k,l]]],{k,1,Length[partition]},{l,1,Length[partition]}];
	subJJpots
]

matrixconv[Hintemp_,varnames_]:=Module[
{oplist,vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp,temp2,list2,list3,listsq,factreps2,sqreps,Dfactreps,termlist,kvals,maxnum,listDfact3,
	DfactrepsF,Hin,dummy},
	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;

	Hin=Expand[Hintemp+dummy];

	oplist=Flatten[Join[Dsubvecp,qsubvecp,fsubvecp]]; (* list of all operators *)
	temp2=Sort[Flatten[Table[Position[Hin,oplist[[i]]],{i,1,Length[oplist]}],1]]; (* positions of each operator in Hin *)

	list2=Pick[temp2,Table[Length[temp2[[i]]]==2,{i,1,Length[temp2]}]]; (* list of terms where operator appears at second level *)
	list3=Pick[temp2,Table[Length[temp2[[i]]]==3,{i,1,Length[temp2]}]];(* list of terms where operator appears at third level *)
	listsq=Pick[list3,Table[list3[[i,3]]==1,{i,1,Length[list3]}]]; (* terms with an operator squared *)
	termlist=Table[list3[[i,1]],{i,1,Length[list3]}]; (* list of terms where each operator appears *)

	factreps2=Table[Apply[Times,Extract[Hin,{list2[[i]],list2[[i+1]]}]]->Apply[Dot,Extract[Hin,{list2[[i]],list2[[i+1]]}]],{i,1,Length[list2],2}];

	sqreps=Table[Extract[Hin,Drop[listsq[[i]],-1]]->Dot[Extract[Hin,listsq[[i]]],Extract[Hin,listsq[[i]]]],{i,1,Length[listsq]}];

	maxnum=5; (* maximum number of operator factors to be considered *)
	kvals=Table[k,{k,2,maxnum}];
	Dfactreps=Table[Null,{i,1,Length[kvals]}];

	Do[
		listDfact3=Pick[list3,Table[And[list3[[i,3]]==0,Count[termlist,termlist[[i]]]==kvals[[k]]],{i,1,Length[list3]}]];
		Dfactreps[[k]]=If[Length[listDfact3]>0,Table[Apply[Times,Table[Extract[Hin,Drop[listDfact3[[i+j]],-1]],{j,0,kvals[[k]]-1}]]->Apply[Dot,Table[Extract[Hin,Drop[listDfact3[[i+j]],-1]],{j,0,kvals[[k]]-1}]],{i,1,Length[listDfact3],kvals[[k]]}],Null],{k,1,Length[kvals]-1}
	];

	DfactrepsF=DeleteCases[Flatten[Dfactreps],Null];

	(Hin/.Join[factreps2,sqreps,DfactrepsF])/.{dummy->0}
]

varinfo[invar_,partition_]:=Module[{whichsub,temp,vartypes,vartype,temp2,varnum},
	whichsub=Position[Table[Length[StringPosition[ToString[invar],partition[[i,1]]]]>0,{i,1,Length[partition]}],True][[1,1]]; (* extracts which subsystem a given operator is part of *)
	temp=StringTrim[ToString[invar],partition[[whichsub,1]]];
	vartypes={"o","j","i"};
	vartype=Position[Table[MemberQ[Characters[temp],vartypes[[i]]],{i,1,Length[vartypes]}],True][[1,1]];
	temp2=StringTrim[temp,vartypes[[vartype]]];
	varnum=If[Length[temp2]>0,temp2,0];
	{whichsub,vartype,varnum} (* returns number of subsystem, type of variable - 1:oscillator, 2:junction, 3:island, and variable number where 0 indicates only one of that type  *)
]

oscpopF[numosc_]:=Module[{oscvals},
	oscvals=Table[Sqrt[i/2],{i,1,numosc-1}];
	SparseArray[{Band[{1,2}]->I*oscvals,Band[{2,1}]->-I*oscvals},{numosc,numosc}]
] (* generates the oscillator momentum operator of dimension numosc *)

oscxopF[numosc_]:=Module[{oscvals},
	oscvals=Table[Sqrt[i/2],{i,1,numosc-1}];
	SparseArray[{Band[{1,2}]->oscvals,Band[{2,1}]->oscvals},{numosc,numosc}]
] (* generates the oscillator momentum operator of dimension numosc *)

qmat[vartype_,dim_,subp0osc_]:=Which[vartype==1,oscpopF[dim]*subp0osc,MemberQ[{2,3},vartype],SparseArray[{Band[{1,1}]->Table[i,{i,-dim,dim}]},{2*dim+1,2*dim+1}]]
 (* generates the charge operators for a generic variable type; flag is true if including oscillator variables *)
 
fmat[vartype_,dim_,subx0osc_]:=Which[vartype==1,oscxopF[dim]*subx0osc,MemberQ[{2,3},vartype],SparseArray[{Band[{1,1}]->0},{2*dim+1,2*dim+1}]]
 (* generates the flux operators for a generic variable type, where set to zero unless oscillator variable; flag is true if including oscillator variables *)
 
qfmatS[i_,subqmat_,subIDmats_]:=If[Length[subIDmats]>1,Apply[KroneckerProduct,DeleteCases[Table[If[j==i,subqmat,subIDmats[[j]]],{j,1,Length[subIDmats]}],{}]],subqmat];(* builds full subsystem operator around output of qmat or fmat *)

qDmat[vartype_,dim_,subx0osc_,disp_]:=Which[vartype==1,Chop[expmat[subx0osc*disp,dim]],MemberQ[{2,3},vartype],Which[disp==-1,SparseArray[{Band[{1,2}]->1},{dim,dim}],disp==1,SparseArray[{Band[{2,1}]->1},{dim,dim}]]]
(* generates the charge displacement operators for a generic variable type; flag is true if including oscillator variables*)
(* i steps through charge variable list for a subsystem, l is displacement value for each charge variable, and j is internal index used to construct Kronecker product for each displacement operator *)

qDmatS[i_,subsinfo_,subdims_,subsx0_,disp_,subIDmats_]:=If[Length[subIDmats]>1,Chop[Apply[KroneckerProduct,DeleteCases[Table[If[i==j,qDmat[subsinfo[[i]],subdims[[i]],subsx0[[i]],disp],subIDmats[[j]]],{j,1,Length[subsinfo]}],{}]]],qDmat[subsinfo[[1]],subdims[[1]],subsx0[[1]],disp]](* builds full subsystem operator around output of qDmat ; flag is true if including oscillator variables *)

suboscHmatS[j_,subH0numlist_,subIDmats_,subdims_,wls_]:=
	If[Length[subIDmats]>1,Apply[KroneckerProduct,Table[If[i!=j,subIDmats[[i]],SparseArray[Band[{1,1}]->Table[wls[[subH0numlist[[i]]]]*(v+1/2),{v,0,subdims[[i]]-1}],{subdims[[i]],subdims[[i]]}]],{i,1,Length[subdims]}]],SparseArray[Band[{1,1}]->Table[wls[[subH0numlist[[1]]]]*(v+1/2),{v,0,subdims[[1]]-1}],{subdims[[1]],subdims[[1]]}]] (* builds full subsystem oscillator operator  *)

makevars[Lbasis_,islands_,partition_]:=Module[{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp},
	{vecp,subvecp}=chopvecF[Lbasis,islands,partition]; (* labels for new coordinates *)
	qvec=Table[ToExpression[StringJoin["q",ToString[i]]],{i,1,Length[Lbasis[[1]]]}];
	qvecp=Table[ToExpression[StringJoin[ToString[vecp[[i]]],"q"]],{i,1,Length[vecp]}];
	fvec=Table[ToExpression[StringJoin["f",ToString[i]]],{i,1,Length[Lbasis[[1]]]}];
	fvecp=Table[ToExpression[StringJoin[ToString[vecp[[i]]],"f"]],{i,1,Length[vecp]}];
	chopvec=Table[ToExpression[StringJoin[ToString[vecp[[i]]],"D"]],{i,1,Length[vecp]}];(* charge displacement operators for new coordinates *)
	Dsubvecp=Table[Table[ToExpression[StringJoin[ToString[subvecp[[j,i]]],"D"]],{i,1,Length[subvecp[[j]]]}],{j,1,Length[partition]}];
	qsubvecp=Table[Table[ToExpression[StringJoin[ToString[subvecp[[j,i]]],"q"]],{i,1,Length[subvecp[[j]]]}],{j,1,Length[partition]}];
	fsubvecp=Table[Table[ToExpression[StringJoin[ToString[subvecp[[j,i]]],"f"]],{i,1,Length[subvecp[[j]]]}],{j,1,Length[partition]}];
	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}
]

convpart[Rmat_,inpartition_]:=Module[
	{zvec,pvec,temppart,varnums,numocc,twonum,twoocc,reparray,pairlist,newpairlist,truepos,outpartition},
	zvec=Table[0,{i,1,Length[Rmat]}];
	pvec=Table[i,{i,1,Length[Rmat]}];
	temppart=Table[{inpartition[[l,1]],Pick[pvec,Table[Sum[Abs[Transpose[Rmat][[inpartition[[i,2]]]]][[j]],{j,1,Length[inpartition[[i,2]]]}][[k]]>0,{i,1,Length[inpartition]},{k,1,Length[Rmat]}][[l]]]},{l,1,Length[inpartition]}];
	varnums=Flatten[Transpose[temppart][[2]]];
	numocc=Table[Count[varnums,pvec[[i]]],{i,1,Length[pvec]}];
(* this part accounts for when a primed variable is shared across subsystems in temppart; common variable is removed from larger subsystem *)
	If[Apply[Times,numocc]!=1,
		If[Count[numocc,3]>0,Print[];Abort[],Null];
		twonum=Flatten[Position[numocc,2]];
		twoocc=Table[Table[MemberQ[Transpose[temppart][[2]][[i]],twonum[[j]]],{i,1,Length[Transpose[temppart][[2]]]}],{j,1,Length[twonum]}];
		If[MemberQ[Table[Count[twoocc[[j]],True]==2,{j,1,Length[twoocc]}],False],Print["problem in partition processing"];Abort[],Null];
		reparray=Table[0,{i,1,Length[twoocc]}];
		Do[
			pairlist=Reverse[SortBy[Pick[Transpose[temppart][[2]],twoocc[[j]]],Length]];
			newpairlist={DeleteCases[pairlist[[1]],Apply[Intersection,pairlist][[1]]],pairlist[[2]]};
			truepos=Flatten[Position[twoocc[[j]],True]];
			reparray[[j]]={{truepos[[1]],2}->newpairlist[[1]],{truepos[[2]],2}->newpairlist[[2]]},
		{j,1,Length[twoocc]}
		];
		outpartition=ReplacePart[temppart,Flatten[reparray]],
		outpartition=temppart;
	];
	outpartition
]

optype[tempe_,varnames_]:=Module[
	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp},

	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;
	Which[Length[Flatten[Table[Position[tempe,chopvec[[i]]],{i,1,Length[chopvec]}],1]]>0,"D",Length[Flatten[Table[Position[tempe,qvecp[[i]]],{i,1,Length[chopvec]}],1]]>0,"Q",Length[Flatten[Table[Position[tempe,Flatten[fsubvecp][[i]]],{i,1,Length[chopvec]}],1]]>0,"F"]
]

vecfunct[inop_,inmat_,vecs_,varnames_,nlevs2_]:=Module[
	{temp,temp2,type,out},

	type=optype[inop,varnames];
	out=0;
	If[type=="D",out=Table[Conjugate[vecs[[i]]].inmat.vecs[[j]],{i,1,Length[vecs]},{j,1,Length[vecs]}],Null]; (* charge displacement operators are antisymmetric *)
	If[Or[type=="Q",type=="F"],
	temp=SparseArray[Flatten[Table[{i,j}->Conjugate[vecs[[i]]].inmat.vecs[[j]],{i,1,Length[vecs]},{j,i,Length[vecs]}],1]];
	temp2=Conjugate[Transpose[UpperTriangularize[temp,1]]];
	out=temp+temp2,Null];(* charge and flux operators are hermitian *)
	If[type==Null,out=SparseArray[Band[{1,1}]->1,{nlevs2,nlevs2}]]; (* identity *)
	If[out==0,Print["vecfunct error"];Abort[],Null];
	out
]
 
 (* looks up all physical level subsystems contained in a given subsystem at hierarchy level slevin *)
 subpartF[partitiontot_, slevin_, spart_] := Module[
  	{inds, slev},
  	For[inds=partitiontot[[slevin-1,2,spart]];
  		slev=slevin-2,
  		slev>0,
  		slev--,
  		inds=Flatten[partitiontot[[slev, 2]][[inds]]]
  	];
  	inds
]

findsysN[LNpartition_,L1num_]:=Module[
	{LNnum, lev},(* for a given L1 system, finds its system number at higher levels *)

  	LNnum=L1num;
  	lev=Table[0,{k,1,Length[LNpartition]}];
  	Do[
   		lev[[i]]=Position[LNpartition[[i,2]],LNnum][[1,1]];
   		LNnum=lev[[i]],
   		{i,1,Length[LNpartition]}
   	];
  	Join[{L1num},lev]
  ]

findClev[LNpartition_,varnames_,Htot_]:=Module[
	{nsubsL1,subops,HC,temp,temp2,temp4,temp5,temp6,temp3,vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp,maxlevF,L1subsysF},

	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;

  	nsubsL1=Length[Flatten[LNpartition[[1,2]]]];
  	subops=Table[Flatten[Transpose[{Dsubvecp,qsubvecp,fsubvecp}][[subN]]],{subN,1,nsubsL1}];
  	HC=Htot*SparseArray[Band[{1,1}]->0,{Length[Htot],Length[Htot]},1]; 
  	(* coupling Hamiltonian between L1 systems (green blocks) *)
  	temp=Table[Flatten[Position[Table[Position[HC[[subN,j]],subops[[subN,k]]]!={},{j,1,Length[HC]}],True]],
  		{subN,1,Length[HC]},{k,1,Length[subops[[subN]]]}];
  	temp2=DeleteCases[Flatten[Table[If[temp[[subN,k]]!= {},{{subops[[subN,k]],subN},temp[[subN,k]]},Null],
  		{subN,1,Length[HC]},{k,1,Length[subops[[subN]]]}],1],Null];
  	(* lists each operator in HC and which L1 subsystems it couples to *)
  	temp3=Table[{temp2[[i,1,1]],Table[{temp2[[i,1,2]],temp2[[i,2,k]]},{k,1,Length[temp2[[i,2]]]}]},{i,1,Length[temp2]}];
  	temp6=Table[0,{tnum,1,Length[temp3]}];
  	Do[
  		temp4=Transpose[Table[Transpose[{findsysN[LNpartition,temp3[[tnum,2,i,1]]],findsysN[LNpartition,temp3[[tnum,2,i,2]]]}],{i,1,Length[temp3[[tnum,2]]]}]];
   		temp5=Table[MemberQ[Table[temp4[[levt,i,1]]==temp4[[levt,i,2]],{i,1,Length[temp4[[levt]]]}],False],{levt,1,Length[LNpartition]+1}];
   		temp6[[tnum]]=temp3[[tnum,1]]->FirstPosition[temp5,False][[1]]-1,{tnum,1,Length[temp3]}
   	];
   	maxlevF=Association[temp6];(* association which takes as input an operator in HC, and gives out the maximum level at which it is present in inter-system coupling *)
   	L1subsysF=Association[Table[temp2[[i,1,1]]->findsysN[LNpartition,temp2[[i,1,2]]],{i,1,Length[temp2]}]];(* association which takes as input an operator in HC, and gives out a list of the systems it operates in at each level *)
  	{maxlevF,L1subsysF}
]

(* creates a list of coupling operators to be evaluated at level slev *)
subopsF[LNpartition_,{maxlevF_,L1subsysF_},Htot_,varnames_,slev_]:=Module[
	{inop,pos,vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp,subopsH,subopsHin,subopsH0,Hcoupletot},	
 
  	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;
  	Hcoupletot=Chop[Expand[Simplify[Sum[2*Htot[[i,j]],{i,1,Length[Htot]},{j,i+1,Length[Htot]}]]]]; (* total coupling Hamiltonian, factor of two since only summing upper half of matrix *)

	subopsH0=Table[Keys[L1subsysF][[Flatten[Position[Transpose[Values[L1subsysF]][[slev]],i]]]],{i,1,Length[Flatten[LNpartition[[slev,2]]]]}]; (* operators contained in Hcoupletot,ordered by L=slev subsystem *)
	subopsHin=DeleteCases[Table[If[maxlevF[subopsH0[[subN,i]]]>=slev,subopsH0[[subN,i]],Null],{subN,1,Length[subopsH0]},{i,1,Length[subopsH0[[subN]]]}],Null,2]; (* pick out the subset of these operators that are needed for levels>=slev  *)
  	subopsH=subopsHin;
  
  	Do[
  		If[Length[subopsHin[[subN]]]>0, 
  			inop=subopsHin[[subN, i]];
    		If[MemberQ[Flatten[Dsubvecp],inop], 
     			pos=Position[Hcoupletot, inop];
     			subopsH[[subN]]=DeleteDuplicates[DeleteCases[Join[subopsH[[subN]],Table[Extract[Hcoupletot,Drop[pos[[k]],-1]],{k,1,Length[pos]}]],inop]],
     			Null
     		],
     		Null
     	],
   		{subN,1,Length[subopsH]}, {i,1,Length[subopsHin[[subN]]]}
   	];
   	subopsH
]

(* lists operators contained in Hin, by L1 system *)
opfind[Hin_,varnames_]:=Module[
	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp,subops},
  	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;
  	subops=Table[Flatten[Transpose[{Dsubvecp, qsubvecp, fsubvecp}][[subN]]],{subN,1,Length[Dsubvecp]}];
  	DeleteCases[Table[If[Length[Position[Hin, subops[[subN,i]]]]>0,subops[[subN,i]]],{subN,1,7},{i,1,Length[subops[[subN]]]}],Null,2]
 ]

(* extracts coupling Hamiltonian from Htot, at a given level slev *)
HCfunct[slev_, LNpartition_, Htot_, varnames_] := 
 Module[{subparts,subparts2,out},
  	subparts=LNpartition[[slev-1,2]];
  	If[slev>2,
  		subparts2=Table[Table[subpartF[LNpartition,slev-1,subparts[[subN,i]]],{i,1,Length[subparts[[subN]]]}],{subN,1,Length[subparts]}]; 
   		out=Table[matrixconv[Chop[Expand[Simplify[Sum[2*Htot[[subparts2[[subN,i,iS]],subparts2[[subN,j,jS]]]],
   			{i,1,Length[subparts2[[subN]]]},
   				{j,i+1,Length[subparts2[[subN]]]},
   					{iS,1,Length[subparts2[[subN,i]]]},
   						{jS,1,Length[subparts2[[subN,j]]]}]]]],varnames],
   							{subN,1,Length[subparts]}],
   		Null
   	];
  	If[slev==2,
  		out=Table[matrixconv[Chop[Expand[Simplify[Sum[2*Htot[[subparts[[subN,i]],subparts[[subN,j]]]],
  			{i,1,Length[subparts[[subN]]]},{j,i+1,Length[subparts[[subN]]]}]]]],varnames],
  			{subN,1,Length[subparts]}],
  		Null
   	];
   	If[slev<2,Print["slev<2"];Abort[]];
  	out
  ]
  
L0prediagH[JJmat0_,gaugeflux_,varnames_,{L1partition_,Rmat_,subECs1_,subELs1_}]:=Module[
	{temp2b,temp3,subECsNO,subELsNO,JJpotNO,subJJpotsNO,vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp},
	
	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;
		
	temp2b=Table[Table[varinfo[vecp[[L1partition[[pnum,2,i]]]],L1partition][[2]]!=2,{i,1,Length[L1partition[[pnum,2]]]}],{pnum,1,Length[L1partition]}];
	temp3=Table[Table[If[Or[temp2b[[i,j]],temp2b[[i,k]]],0,1],{j,1,Length[temp2b[[i]]]},{k,1,Length[temp2b[[i]]]}],{i,1,Length[L1partition]}];

	(* zero out oscillator and island variables *)	
	subECsNO=Simplify[Chop[Table[4*qvecp[[L1partition[[i,2]]]].(subECs1[[i,i]]*temp3[[i]]).qvecp[[L1partition[[i,2]]]],{i,1,Length[L1partition]}]]];
	subELsNO=Simplify[Chop[Table[(4*Pi^2)/2*fvecp[[L1partition[[i,2]]]].(subELs1[[i,i]]*temp3[[i]]).fvecp[[L1partition[[i,2]]]],{i,1,Length[L1partition]}]]];
	JJpotNO=JJHfunct[JJmat0,gaugeflux,Rmat,varnames,L1partition,False];
	subJJpotsNO=partitionJJmat[L1partition,JJpotNO,chopvec,Global`IDmat,gaugeflux];
	Diagonal[Chop[Simplify[subECsNO+subELsNO+subJJpotsNO]]]
]

(* calculates subsystem Hamiltonian with prediagonalization *)
subHFP[L1partition_,subHNO_,subH_,subN_,varnames_,nvalsin_,Htot_,{subp0oscs_,subx0oscs_,wls_,subH0oscslist_,subH0oscs_},subsinfo_,subdims_,nlevsJJ_]:=Module[
	{subIDmats,totIDmat,subqmats,subfmats,displist,subH0list,subH0varlist,subH0list2,subH0numlist,Dreps,IDreps,qreps,freps,Hoscreps,subrepsF,
	vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp,nJJs,subdims0,JJinds,tempq,tempf,Dreps0,IDreps0,qreps0,freps0,subreps0,
	temp0,temp,temp2,flag0,subdims1,vals,vecs,vals0,vecs0,Dmat,Dtemp,Dtemp2,Dreps1,qreps1,freps1,subreps1},

	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;

	nJJs=Count[subsinfo,2];
	temp0=Table[And[subsinfo[[i]]!=2,MemberQ[subsinfo,2]],{i,1,Length[subsinfo]}];
	temp=Pick[subdims,temp0];
	temp2=Apply[Times,temp];(* dimension of non-JJ Hilbert space, set to 1 if no Josephson variables in subsystem *)

	flag0=Table[subsinfo[[i]]==2,{i,1,Length[subsinfo]}];
	subdims0=Table[2*Pick[nvalsin[[subN]],flag0][[i]]+1,{i,1,nJJs}] ;(* resulting dimension of each JJ subspace *)
	subdims1=Table[Which[subsinfo[[i]]==1,nvalsin[[subN]][[i]],subsinfo[[i]]==2,1,subsinfo[[i]]==3,2*nvalsin[[subN]][[i]]+1],{i,1,Length[subsinfo]}];

	JJinds=Flatten[Position[subsinfo,2]];
	subIDmats=Table[IdentityMatrix[subdims0[[i]],SparseArray],{i,1,nJJs}];
	totIDmat=IdentityMatrix[Apply[Times,subdims0],SparseArray];

	tempq=DeleteCases[Table[If[Length[Position[subHNO,qsubvecp[[subN,i]]]]>0,{qsubvecp[[subN,i]],qmat[2,nvalsin[[subN,i]],0]},Null],{i,1,Length[qsubvecp[[subN]]]}],Null];
	tempf=DeleteCases[Table[If[Length[Position[subHNO,fsubvecp[[subN,i]]]]>0,{qsubvecp[[subN,i]],fmat[2,nvalsin[[subN,i]],0]},Null],{i,1,Length[fsubvecp[[subN]]]}],Null];

	displist=Table[DeleteDuplicates[Extract[subHNO,Position[subHNO,Dsubvecp[[subN,JJinds[[i]]]]]/. 0->1]],{i,1,nJJs}];

	Dreps0=Association[Flatten[Chop[Table[If[Length[displist[[i]]]>0,Dmat=qDmatS[i,Table[2,{k,1,nJJs}],subdims0,subx0oscs[[subN]],displist[[i,1]],subIDmats];
		{Dsubvecp[[subN,JJinds[[i]]]][displist[[i,1]]]->Dmat,Dsubvecp[[subN,JJinds[[i]]]][displist[[i,2]]]->ConjugateTranspose[Dmat]},Null],{i,1,nJJs}]],1]];
	IDreps0=Association[{Global`IDmat->totIDmat}];
	qreps0=Association[If[Length[tempq]>0,Table[tempq[[i,1]]->qfmatS[i,tempq[[i,2]],subIDmats],{i,1,Length[tempq]}],{dummy->1}]];
	freps0=Association[If[Length[tempf]>0,Table[tempf[[i,1]]->qfmatS[i,tempf[[i,2]],subIDmats],{i,1,Length[tempf]}],{dummy->1}]];
	subreps0=Join[Dreps0,IDreps0,qreps0,freps0];

	{vals0,vecs0}=Eigensystem[Chop[subHNO/.subreps0],-nlevsJJ];
	vals=Chop[Reverse[vals0]];
	vecs=Chop[Reverse[vecs0]];
	Dtemp=Association[Table[Keys[Dreps0][[i]]->KroneckerProduct[vecfunct[Keys[Dreps0][[i]],Values[Dreps0][[i]],vecs,varnames,nlevsJJ],IdentityMatrix[temp2,SparseArray]],{i,1,Length[Dreps0],2}]];
	Dtemp2=Association[Table[Keys[Dreps0][[2*i]]->ConjugateTranspose[Values[Dtemp][[i]]],{i,1,Length[Dtemp]}]];
	Dreps1=Join[Dtemp,Dtemp2];
	qreps1=Association[If[Length[tempq]>0,Table[Keys[qreps0][[i]]->KroneckerProduct[vecfunct[Keys[qreps0][[i]],Values[qreps0][[i]],vecs,varnames,nlevsJJ],IdentityMatrix[temp2,SparseArray]],{i,1,Length[qreps0]}],{dummy->1}]];
	freps1=Association[If[Length[tempf]>0,Table[Keys[freps0][[i]]->KroneckerProduct[vecfunct[Keys[freps0][[i]],Values[freps0][[i]],vecs,varnames,nlevsJJ],IdentityMatrix[temp2,SparseArray]],{i,1,Length[freps0]}],{dummy->1}]];
	subreps1=KeyDrop[Join[Dreps1,qreps1,freps1],dummy];

(* this subsection takes JJ variables expressed in prediagonalized basis and writes full Hamiltonian including oscillators and islands  *)

	subIDmats=Table[IdentityMatrix[subdims1[[i]],SparseArray],{i,1,Length[subdims1]}];
	totIDmat=IdentityMatrix[nlevsJJ*Apply[Times,subdims1],SparseArray];
	subqmats=Table[If[subsinfo[[i]]!=2,qmat[subsinfo[[i]],nvalsin[[subN]][[i]],subp0oscs[[subN]][[i]]],Null],{i,1,Length[subsinfo]}];
	subfmats=Table[If[subsinfo[[i]]!=2,fmat[subsinfo[[i]],nvalsin[[subN]][[i]],subx0oscs[[subN]][[i]]],Null],{i,1,Length[subsinfo]}];

	displist=Table[If[subsinfo[[i]]!=2,DeleteDuplicates[Extract[Htot,Position[Htot,Dsubvecp[[subN,i]]]/. 0->1]],{}],{i,1,Length[subvecp[[subN]]]}];
	subH0list=PadRight[Diagonal[subH0oscs],Length[vecp]][[L1partition[[subN,2]]]];
	subH0varlist=Table[If[NumberQ[subH0oscslist[[subN,i]]],0,ToExpression[StringDrop[ToString[subH0oscslist[[subN,i]]],2]]],{i,1,Length[subH0oscslist[[subN]]]}];
	subH0list2=Table[If[NumberQ[subH0oscslist[[subN,i]]],0,subH0oscslist[[subN,i]]],{i,1,Length[subH0oscslist[[subN]]]}];
	subH0numlist=Table[If[NumberQ[subH0varlist[[i]]],0,Position[vecp,subH0varlist[[i]]][[1,1]]],{i,1,Length[subH0list]}];

	Dreps=Association[Flatten[Chop[Table[If[And[subsinfo[[i]]!=2,Length[displist[[i]]]>0],Dsubvecp[[subN,i]][displist[[i,l]]]->KroneckerProduct[IdentityMatrix[nlevsJJ,SparseArray],
		qDmatS[i,subsinfo,subdims1,subx0oscs[[subN]],displist[[i,l]],subIDmats]],Null],{i,1,Length[subvecp[[subN]]]},{l,1,Length[displist[[i]]]}]],1]];
	IDreps=Association[{Global`IDmat->totIDmat}];
	qreps=Association[DeleteCases[Table[If[subsinfo[[i]]!=2,qsubvecp[[subN,i]]->KroneckerProduct[IdentityMatrix[nlevsJJ,SparseArray],qfmatS[i,subqmats[[i]],subIDmats]],Null],{i,1,Length[subvecp[[subN]]]}],Null]];
	freps=Association[DeleteCases[Table[If[subsinfo[[i]]!=2,fsubvecp[[subN,i]]->KroneckerProduct[IdentityMatrix[nlevsJJ,SparseArray],qfmatS[i,subfmats[[i]],subIDmats]],Null],{i,1,Length[subvecp[[subN]]]}],Null]];
	Hoscreps=Association[DeleteCases[Flatten[Table[If[subH0numlist[[j]]==0,Null,subH0list2[[j]]->KroneckerProduct[IdentityMatrix[nlevsJJ,SparseArray],suboscHmatS[j,subH0numlist,subIDmats,subdims1,wls]]],{j,1,Length[subH0list]}]],Null]];
	
	subrepsF=Join[Dreps,IDreps,qreps,freps,Hoscreps,subreps1];
	{subrepsF,Chop[subH/.subrepsF]}
]

(* calculates subsystem Hamiltonian without prediagonalization *)
subHF[L1partition_,subH_,subN_,varnames_,nvalsin_,Htot_,{subp0oscs_,subx0oscs_,wls_,subH0oscslist_,subH0oscs_},subsinfo_,subdims_]:=Module[
	{subIDmats,totIDmat,subqmats,subfmats,displist,subH0list,subH0varlist,subH0list2,subH0numlist,Dreps,IDreps,qreps,freps,Hoscreps,subrepsF,
	vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp},

	{vecp,subvecp,qvec,qvecp,fvec,fvecp,chopvec,Dsubvecp,qsubvecp,fsubvecp}=varnames;

	subIDmats=Table[SparseArray[{Band[{1,1}]->1},{subdims[[i]],subdims[[i]]}],{i,1,Length[subsinfo]}];
	totIDmat=If[Length[subIDmats]>1,Apply[KroneckerProduct,Table[subIDmats[[i]],{i,1,Length[subIDmats]}]],subIDmats[[1]]];
	subqmats=Table[qmat[subsinfo[[i]],nvalsin[[subN,i]],subp0oscs[[subN]][[i]]],{i,1,Length[subsinfo]}];
	subfmats=Table[fmat[subsinfo[[i]],nvalsin[[subN,i]],subx0oscs[[subN]][[i]]],{i,1,Length[subdims]}];

	displist=Table[DeleteDuplicates[Extract[Htot,Position[Htot,Dsubvecp[[subN,i]]]/. 0->1]],{i,1,Length[subvecp[[subN]]]}];
	subH0list=PadRight[Diagonal[subH0oscs],Length[vecp]][[L1partition[[subN,2]]]];			
	subH0varlist=Table[If[NumberQ[subH0oscslist[[subN,i]]],0,ToExpression[StringDrop[ToString[subH0oscslist[[subN,i]]],2]]],{i,1,Length[subH0oscslist[[subN]]]}];
	subH0list2=Table[If[NumberQ[subH0oscslist[[subN,i]]],0,subH0oscslist[[subN,i]]],{i,1,Length[subH0oscslist[[subN]]]}];
	subH0numlist=Table[If[NumberQ[subH0varlist[[i]]],0,Position[vecp,subH0varlist[[i]]][[1,1]]],{i,1,Length[subH0list]}];
	
	Dreps=Association[Flatten[Chop[Table[If[Length[displist[[i]]]>0,Dsubvecp[[subN,i]][displist[[i,l]]]->qDmatS[i,subsinfo,subdims,subx0oscs[[subN]],displist[[i,l]],subIDmats],Null],{i,1,Length[subvecp[[subN]]]},{l,1,Length[displist[[i]]]}]],1]];
	IDreps=Association[{Global`IDmat->totIDmat}];
	qreps=Association[Table[qsubvecp[[subN,i]]->qfmatS[i,subqmats[[i]],subIDmats],{i,1,Length[subvecp[[subN]]]}]];
	freps=Association[Table[fsubvecp[[subN,i]]->qfmatS[i,subfmats[[i]],subIDmats],{i,1,Length[subvecp[[subN]]]}]];
	Hoscreps=Association[DeleteCases[Flatten[Table[If[subH0numlist[[j]]==0,Null,subH0list2[[j]]->suboscHmatS[j,subH0numlist,subIDmats,subdims,wls]],{j,1,Length[subH0list]}]],Null]];

	subrepsF=Join[Dreps,IDreps,qreps,freps,Hoscreps];
	{subrepsF,Chop[subH /.subrepsF]}
]

 

End[] (* End Private Context *)

EndPackage[]