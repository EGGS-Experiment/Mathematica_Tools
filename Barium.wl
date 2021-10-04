(* ::Package:: *)

(* ::Chapter::Closed:: *)
(*Begin*)


(*version: 0.9.1
last updated: Jun 7 2021*)


BeginPackage["Barium`"];


(*remove all conflicts*)
Unprotect["`*"];
ClearAll["`*"];


Barium::usage="Barium.m: ***";


(* ::Chapter::Closed:: *)
(*Declarations*)


(* ::Section::Closed:: *)
(*Config*)


ChangeAll::usage="Opens a Dialog Box to allow all variables to be changed";
ChangeVariable::usage="Changes the value of a variable";


(* ::Section::Closed:: *)
(*Variables*)


(* ::Subsection::Closed:: *)
(*Ion Values*)


mAtom::usage="Atom Mass";
mMol::usage="Molecule Mass";


(* ::Subsection::Closed:: *)
(*Trap Values*)


VRF::usage="Trap RF Voltage";
\[Omega]RF::usage="Trap RF Frequency";
U::usage="Trap Endcap Voltage";
r0::usage="Trap Radial Electrode Distance";
z0::usage="Trap Axial Electrode Distance";
\[Kappa]r::usage="Trap Radial Geometric Factor";
\[Kappa]z::usage="Trap Axial Geometric Factor";


qMol::usage="Trap q-factor for the molecule/qubit";
wsecMol::usage="Trap secular frequencies in the radial and axial direction";


(* ::Subsection::Closed:: *)
(*Normal Modes*)


wmvals::usage="Normal mode frequencies of the trap";
vmvals::usage="Normal mode eigenvectors of the trap";
\[Eta]::usage="Lamb-dicke parameters";


(* ::Subsection::Closed:: *)
(*Operators*)


Htrap::usage="Trap hamiltonian for individual modes";
Ut::usage="Time evolution operators for individual modes";
\[Sigma]p::usage="Qubit raising operator";
\[Sigma]m::usage="Qubit lowering operator";
\[Sigma]x::usage="Pauli x-operator";
ap::usage="Phonon creation operator";
am::usage="Phonon annihilation operator";
nh::usage="Phonon number operator";
Dni::usage="Nth order Debye-waller operator";


(* ::Subsection::Closed:: *)
(*Bases*)


Basis::usage="Holds the identities of all the states in the basis";
\[Psi]Basis::usage="Wavefunctions for the basis";


(* ::Section::Closed:: *)
(*Functions*)


(* ::Subsection::Closed:: *)
(*States*)


PureState::usage="Returns the desired qubit and phonon state";
GroundState::usage="Returns the motional ground state in the desired qubit state";
ThermalState::usage="Returns a thermal state at the given temperature and in the desired qubit state";
CoherentState::usage="Returns a coherent state with the given coherences in the desired qubit state";
(*BellState::usage="";*)
ToIntPic::usage="Takes a Hamiltonian into the specified interaction picture";
SolveSchrodinger::usage="Solves the interaction picture Schrodinger equation";


(* ::Subsection::Closed:: *)
(*Processing Data*)


GetPops::usage="Gets the population of each qubit state";
GetPhonons::usage="Gets the average number of phonons in each mode";
GetFidelity::usage="Gets the Bell-state fidelity";
GetPhononPops::usage="Gets the population of each phonon state";
GetPoints::usage="Gets a number of points from a function";


(* ::Subsection::Closed:: *)
(*Getting Values*)


ShowEigen::usage="Shows the normal mode data";
ShowTrap::usage="Shows the trap data";
ShowAll::usage="Shows all data";


(* ::Section::Closed:: *)
(*Hamiltonians*)


HCarrier::usage="Returns a carrier hamiltonian acting on a single qubit";
Hrsb::usage="Returns a red-detuned sideband hamiltonian acting on a single qubit";
Hbsb::usage="Returns a blue-detuned sideband hamiltonian acting on a single qubit";
Hbad::usage="Returns the bad terms acting on a single qubit";
HrsbN::usage="Returns the Nth order red-detuned sideband hamiltonian acting on a single qubit";
HbsbN::usage="Returns the Nth order blue-detuned sideband hamiltonian acting on a single qubit";
HAll::usage="Returns the carrier, red-detuned sideband, and blue-detuned sideband hamiltonians acting on a single qubit";


(* ::Chapter::Closed:: *)
(*Implementation*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*Setup*)


(* ::Subsection::Closed:: *)
(*Constants*)


(*fundamental*)
\[HBar]=1.0545718*^-34;qe=1.60217662*^-19;kb=1.38064852*^-23;ke2=2.306569*^-28;
(*defined*)
amu=1.66*^-27;wmhz=2\[Pi] 1.*^6;wkhz=2\[Pi] 1.*^3;debye=3.336*^-30;\[Mu]s=1.*^-6;mm=1.*^-3;


(* ::Subsection::Closed:: *)
(*Chain Configuration*)


mArr={1,1};
numQubits=2;
activeModes={{1,2},{}};
minPhonons=0;
maxPhonons=10;


(* ::Subsection::Closed:: *)
(*Config*)


numIons=0;
numModes=0;
numModesTOT=0;
numPhonons=0;
numStates=0;
numQubitStates=0;
numPhononStates=0;


(*setter function to update the simulation variables*)
UpdateConfig[]:=With[{},
(*allow variables to be changed*)
Unprotect[numIons,numModes,numModesTOT,numPhonons,numStates,numQubitStates,numPhononStates];
(*calculate number of ions/modes/states*)
numIons=Length[mArr];
numModes=Length/@activeModes;numModesTOT=Total[numModes];
numPhonons=maxPhonons-minPhonons+1;
numStates=2^numQubits*numPhonons^numModesTOT;
numQubitStates=2^numQubits;
numPhononStates=numPhonons^numModesTOT;
(*prevent variables from being changed by user*)
Protect[numIons,numModes,numModesTOT,numPhonons,numStates,numQubitStates,numPhononStates];]


UpdateConfig[];


(* ::Subsection::Closed:: *)
(*Bases*)


Basis={};
\[Psi]Basis={};
(*\[Psi]diff holds the symbols for the derivative of \[Psi]Basis and is only used in the schrodinger solver
this is just to make things a bit quicker for the schrodinger solver*)
\[Psi]diff={};


(*setter function to update the basis*)
CreateBasis[]:=With[{BasisHolder=Join[Array[a,numQubits],Array[n,numModesTOT]]},
(*allow variables to be changed*)
Unprotect[Basis,\[Psi]Basis,\[Psi]diff];
(*Basis tells us the identity of each state. For example, if we have 2 qubits/5 phonons/2 modes,
then the second state would be (0,0,0,1): both qubits in the ground state and one phonon in the second mode*)
Basis=Tuples[Join[ConstantArray[{0,1},numQubits],ConstantArray[Range[minPhonons,maxPhonons],numModesTOT]]];
(*define \[Psi]Basis in regular context to allow use*)
Block[{$Context="Barium`"},\[Psi]Basis=ToExpression[StringInsert["\[Psi][t]",ToString[#],2]&/@Range[numStates]];];
(*evaluate \[Psi]diff once and ahead of time to reduce schrodinger solver overhead*)
\[Psi]diff=ToExpression[StringInsert["\[Psi]'[t]",ToString[#],2]&/@Range[numStates]];
(*prevent variables from being changed by user*)
Protect[Basis,\[Psi]Basis,\[Psi]diff];]


CreateBasis[];


(* ::Section::Closed:: *)
(*Variables*)


(* ::Subsection::Closed:: *)
(*Ion Values*)


mMol=40*amu;
mAtom=40*amu;
\[Omega]0=80*wmhz;


(* ::Subsection::Closed:: *)
(*Trap Values*)


VRF=150;
\[Omega]RF=20*wmhz;
\[Kappa]r=1;
r0=0.5*mm;
U=20;
z0=2*mm;
\[Kappa]z=1;


(*calculate the q-parameter for the molecule*)
qMol:=(2*qe*VRF*\[Kappa]r)/(mMol*\[Omega]RF^2*r0^2);
(*calculates the q-parameter for the atom*)
qAtom:=(2*qe*VRF*\[Kappa]r)/(mAtom*\[Omega]RF^2*r0^2);
(*calculates the secular frequency for the molecule*)
wsecMol:=\[Omega]RF/2*Sqrt[{qMol^2/2-#,2.*#}]&[(4*qe*U*\[Kappa]z)/(z0^2*mMol*\[Omega]RF^2)];
(*calculates the secular frequency for atom*)
wsecAtom:=\[Omega]RF/2*Sqrt[{qAtom^2/2-#,2.*#}]&[(4*qe*U*\[Kappa]z)/(z0^2*mMol*\[Omega]RF^2)];


(* ::Subsection::Closed:: *)
(*Laser Values*)


(*store \[CapitalDelta]k as an array in case we need different splittings*)
\[CapitalDelta]k=ConstantArray[0,numQubits];
(*same for R0*)
R0=ConstantArray[0,numQubits];


(* ::Subsection::Closed:: *)
(*Normal Modes*)


(*the equilibrium separation between ions*)
l0:=(ke2/(mMol*wsecMol[[2]]^2))^(1/3);


wmvals={};
vmvals={};
\[Eta]={};


(*setter functions to calculate eigenmodes from trap parameters*)
UpdateEigen[]:=Module[{massArray=mArr/.{1->mMol,0->mAtom},wsecvalsMAT=Transpose[mArr/. {1->wsecMol,0->wsecAtom}],z0eqArr=Array[z0eq,numIons],L0EQlist,guessz0eq=N[Range[-#,#,2*#/(numIons-1)]&[(numIons)^(1/3)]],L0Mn,VMatX,VMatZ,eigenholder},
Unprotect[z0eq,wmvals,vmvals,\[Eta]];
(*calculate ion equilibrium distances*)
L0EQlist=z0eqArr[[#]]==(Sum[(z0eqArr[[#]]-z0eqArr[[k]])^(-2),{k,1,#-1}]-Sum[(z0eqArr[[#]]-z0eqArr[[k]])^(-2),{k,#+1,numIons}])&/@Range[numIons];
guessz0eq=Table[{z0eqArr[[i]],guessz0eq[[i]]},{i,numIons}];z0eqArr=z0eqArr/. Chop[FindRoot[L0EQlist,guessz0eq]];
(*get Subscript[l, 0] to the nth power*)
L0Mn[i_,j_,n_]:=If[i==j,0,Abs[z0eqArr[[i]]-z0eqArr[[j]]]^n];
(*create matrices for normal modes*)
VMatZ=Table[Sqrt[massArray[[i]]/massArray[[j]]]*wsecvalsMAT[[2,i]]^2*If[i==j,1+2*Sum[L0Mn[i,id,-3],{id,numIons}],-2*L0Mn[i,j,-3]],{i,numIons},{j,numIons}];
VMatX=Table[Sqrt[massArray[[i]]/massArray[[j]]]*wsecvalsMAT[[2,i]]^2*If[i==j,(wsecvalsMAT[[1,i]]/wsecvalsMAT[[2,i]])^2-Sum[L0Mn[i,id,-3],{id,numIons}],L0Mn[i,j,-3]],{i,numIons},{j,numIons}];
(*get eigendata*)
eigenholder={Eigensystem[VMatX],Eigensystem[VMatZ]};
wmvals=Sqrt[eigenholder[[;;,1]]];vmvals=eigenholder[[;;,2]];
(*calculate lamb-dicke*)
\[Eta]=Sqrt[\[HBar]/(2)]*Table[\[CapitalDelta]k[[i]]/Sqrt[massArray[[i]]*wmvals[[d,q]]],{d,2},{i,numIons},{q,numIons}];
Protect[wmvals,vmvals,\[Eta]];]


(* ::Subsection::Closed:: *)
(*Operators*)


(*holders are used to store the diagonal elements of the time-independent hamiltonians
since these values are reused often*)
Htrapholder:=Table[wmvals[[dir,activeModes[[dir,n]]]]*(Basis[[;;,numQubits+n+(dir-1)*numModes[[dir]]]]+1/2),{dir,2},{n,numModes[[dir]]}];
Htrap:=Map[SparseArray[{Band[{1,1}]->\[HBar]*#},{numStates,numStates}]&,Htrapholder,{2}];
Ut:=Map[SparseArray[{Band[{1,1}]->Exp[-I*t*#]},{numStates,numStates}]&,Htrapholder,{2}];


\[Sigma]p={};\[Sigma]m={};\[Sigma]x={};
ap={};am={};
nh={};
dp={};dm={};


\[Sigma]pholder[i_]:=With[{rstates=ConstantArray[ConstantArray[1,#1],#2]&[2^(numQubits-i),2^(i-1)],lstates=ConstantArray[ConstantArray[0,#1],#2-1]&[2^(numQubits-i),2^(i-1)]},
If[i>1,Riffle[rstates,lstates],rstates]//Flatten]


apholder[i_]:=With[{rstates=ConstantArray[ConstantArray[#,numPhonons^(numModesTOT-i)]&/@Sqrt[Range[minPhonons+1,maxPhonons]],numPhonons^(i-1)],lstates=ConstantArray[ConstantArray[0,#1],#2-1]&[numPhonons^(numModesTOT-i),numPhonons^(i-1)]},
Flatten[If[i>1,Riffle[rstates,lstates],rstates]]];


dpholder[i_]:=With[{rstates=ConstantArray[ConstantArray[#,numPhonons^(numModesTOT-i)]&/@Sqrt[1/Range[minPhonons+1,maxPhonons]],numPhonons^(i-1)],lstates=ConstantArray[ConstantArray[0,#1],#2-1]&[numPhonons^(numModesTOT-i),numPhonons^(i-1)]},
Flatten[If[i>1,Riffle[rstates,lstates],rstates]]];


Lnn[\[CapitalDelta]n_,ion_,mode_,dir_]:=SparseArray[{Band[{1,1}]->#},{numStates,numStates}]&[LaguerreL[#,\[CapitalDelta]n,(\[Eta][[dir,ion,mode]]*vmvals[[dir,mode,ion]])^2]&/@Basis[[;;,numQubits+mode+(dir-1)*numModes[[dir]]]]]

Dni[\[CapitalDelta]n_,ion_,mode_,dir_]:=Module[{modelist=Flatten[Position[activeModes[[dir]],#]&/@DeleteCases[activeModes[[dir]],mode]],L0diag,LMAT,prefactor},
L0diag=Times@@MapIndexed[LaguerreL[#1,0,(\[Eta][[dir,ion,modelist[[First[#2]]]]]*vmvals[[dir,modelist[[First[#2]]],ion]])^2]&,Transpose[Basis[[;;,(numQubits+modelist)]]],{2}];
LMAT=SparseArray[{Band[{1,1}]->L0diag},{numStates,numStates}] . Lnn[Abs[\[CapitalDelta]n],ion,mode,dir];
(I*If[\[CapitalDelta]n==0,1,\[Eta][[dir,ion,mode]]])^Abs[\[CapitalDelta]n]*Exp[-(Total[\[Eta][[dir,;;,mode]]^2]/2)]*(Dot@@#)&[If[\[CapitalDelta]n>0,Join[{LMAT},ConstantArray[dp[[dir,mode]],Abs[\[CapitalDelta]n]]],Join[ConstantArray[dm[[dir,mode]],Abs[\[CapitalDelta]n]],{LMAT}]]]];


(*setter function to define the qubit operators;
since most operators depend only on how the system is set up,
the setters for operators are only called when the config is changed (in ChangeAll)*)
CreateQubitOperators[]:=With[{\[Sigma]ptmp=SparseArray[{Band[{2^(numQubits-#)+1,1}]->\[Sigma]pholder[#]},{numQubitStates,numQubitStates}]&/@Range[numQubits]},
Unprotect[\[Sigma]p,\[Sigma]m,\[Sigma]x];
\[Sigma]p=KroneckerProduct[#,IdPhonon]&/@\[Sigma]ptmp;
\[Sigma]m=Transpose/@\[Sigma]p;
\[Sigma]x=\[Sigma]p+\[Sigma]m;
Protect[\[Sigma]p,\[Sigma]m,\[Sigma]x];]


(*setter function to define the mode operators
mode operators are a bit different: they are 2D arrays (instead of 1D arrays, like the qubit operators)
to separately hold the radial and axial operators - radial is the first row, and axial is the second*)
CreateModeOperators[]:=With[{aptmp=SparseArray[{Band[{numPhonons^(numModesTOT-#)+1,1}]->apholder[#]},{numPhononStates,numPhononStates}]&/@Range[numModesTOT],
nhtmp=SparseArray[{Band[{1,1}]->Basis[[;;,numQubits+#]]},{numStates,numStates}]&/@Range[numModesTOT]},
Unprotect[ap,am,nh];
ap={aptmp[[;;numModes[[1]]]],aptmp[[numModes[[1]]+1;;]]};
ap=Map[KroneckerProduct[IdQubit,#]&,ap,{2}];
am=Map[Transpose,ap,{2}];
nh={nhtmp[[;;numModes[[1]]]],nhtmp[[numModes[[1]]+1;;]]};
Protect[ap,am,nh];]


(*setter function to define the debye-waller operators
since the debye-waller operators act on the phonon modes, they are also 2D arrays like the regular mode operators*)
CreateDebyeWallerOperators[]:=With[{dptmp=SparseArray[{Band[{numPhonons^(numModesTOT-#)+1,1}]->dpholder[#]},{numPhononStates,numPhononStates}]&/@Range[numModesTOT]},
Unprotect[dp,dm];
dp={dptmp[[;;numModes[[1]]]],dptmp[[numModes[[1]]+1;;]]};
dp=Map[KroneckerProduct[IdQubit,#]&,dp,{2}];
dm=Map[Transpose,dp,{2}];
Protect[dp,dm];]


CreateOperators[]:=With[{},
Unprotect[IdPhonon,IdQubit];
IdPhonon=IdentityMatrix[numPhononStates,SparseArray];
IdQubit=IdentityMatrix[numQubitStates,SparseArray];
Protect[IdPhonon,IdQubit];
CreateQubitOperators[];
CreateModeOperators[];
CreateDebyeWallerOperators[];]


(* ::Section::Closed:: *)
(*Functions*)


(* ::Subsection::Closed:: *)
(*States*)


PureState[qubitstate_List,phononstate_List]:=With[{qubitInd=((2^Range[numQubits-1,0,-1]) . qubitstate)*numPhononStates,phononInd=(numPhonons^Range[numModesTOT-1,0,-1]) . phononstate},
(*use sparsearrays to stay compatible with everything else*)
SparseArray[{{qubitInd+phononInd+1}->1},numStates]]


GroundState[qubitstate_List]:=PureState[qubitstate,ConstantArray[0,numModesTOT]];


ThermalState[qubitstate_List,temp_]:=With[{startInd=((2^Range[numQubits-1,0,-1]) . qubitstate)*numPhononStates+1, (*first calculate index that corresponds to given qubit state*)
(*get phonon numbers from basis and multiply with mode frequencies, then sum and multiply by \[HBar]/Subscript[k, b] to get probability of state*)
(*square root since we have a probability*)
prob=Sqrt[Exp[-\[HBar]/(kb*temp)*(#1 . #2)]]&[Basis[[;;numPhononStates,numQubits+1;;]],Flatten[MapIndexed[wmvals[[First@#2,#1]]&,activeModes]]]},
(*normalize and return as sparsearray*)
SparseArray[Band[{startInd}]->Normalize[prob],numStates]]


CoherentState[qubitstate_List,\[Alpha]Arr_List]:=Module[{startInd=((2^Range[numQubits-1,0,-1]) . qubitstate)*numPhononStates+1,\[Alpha]Arr2=Flatten[\[Alpha]Arr],coherentfunc,prob},
(*used to convert \[Alpha] value into coefficient*)
coherentfunc[n_,\[Alpha]_]:=\[Alpha]^n/Sqrt[n!];
(*calculate coefficients for each state then multiply them all together*)
prob=Times@@Map[Thread[coherentfunc[Basis[[;;numPhononStates,numQubits+#]],\[Alpha]Arr2[[#]]]]&,Range[numModesTOT]];
(*normalize and return as sparsearray*)
SparseArray[Band[{startInd}]->Normalize[prob],numStates]]


BellState[ion1_Integer,ion2_Integer]:=1/Sqrt[2]*(SparseArray[{{1}->1,{2^(ion2-1)+2^(ion1-1)+1}->I},{numQubitStates}])


(*need to use simplify and specify t\[Element]Reals otherwise t is treated as potentially complex*)
ToIntPic[H_,U_]:=Simplify[ConjugateTranspose[U],t\[Element]Reals] . H . U;


(* ::Subsection::Closed:: *)
(*Diff Eqs*)


SolveSchrodinger[H_,initialstate_,t1_,t2_]:=With[{$Context="`Barium`"},(*need regular context to allow users to look at their results*)
(*dispatch basically makes replacement a lot faster*)
(*use stiffness switching since populations may oscillate a lot, and set method to automatic since mathematica knows best*)
Dispatch[First@NDSolve[#,\[Psi]Basis,{t,t1,t2},MaxSteps->\[Infinity],Method->{"StiffnessSwitching",Method->{Automatic,Automatic},"EquationSimplification"->"Solve"}]]&[Join[Thread[(\[Psi]Basis/.{t->t1})==initialstate],Thread[\[Psi]diff+I/\[HBar]*(#)==0]&[H . \[Psi]Basis]]]]


(* ::Subsection::Closed:: *)
(*Processing Data*)


(*take sum of squares within a qubit state to get population*)
GetPops[soln_]:=Total[Abs[ArrayReshape[\[Psi]Basis,{numQubitStates,numPhononStates}]]^2,{2}]/.soln;


(*apply the phonon number operator and take expectation value for each mode*)
GetPhonons[soln_,mode_,dir_]:=Re[Conjugate[\[Psi]Basis] . nh[[dir,mode]] . \[Psi]Basis]/.soln;


(*since we don't deal with mixed states, fidelity is just propability of desired state;
we effectively trace over phonon modes by calculating fidelity for each phonon state and summing them all up*)
GetFidelity[soln_,qubitstate_]:=Total[Abs[qubitstate . # ]^2]&[ArrayReshape[\[Psi]Basis,{numQubitStates,numPhononStates}]]/.soln;


(*get probability that a given phonon mode is populated*)
GetPhononPops[soln_]:=Total[Abs[#]^2,{2}]&[Transpose[ArrayReshape[\[Psi]Basis,{numQubitStates,numPhononStates}]]]/.soln;


GetPoints[func_,t1_,t2_,numPoints_:300]:=Module[{func2,tArr=Range[t1,t2,(t2-t1)/(numPoints)]},
(*defining a function that evaluates the interpolating functions speeds up the getting of points*)
func2[thold_]:=Evaluate[func/.{t->thold}];Thread[func2[tArr]]]


(* ::Subsection::Closed:: *)
(*Getting Values*)


ShowEigen[]:=With[{$Context="EGGS`Private`",eigenres=Table[{MatrixForm[vmvals[[d,j]]],wmvals[[d,j]]/wmhz},{d,2},{j,numIons}],wsecres=Transpose[{wsecMol,wsecAtom}]/wmhz},
Print[Labeled[TableForm[wsecres,TableHeadings->{{"Radial","Axial"},{"Molecule","Atom"}}],"Secular Frequencies (x2\[Pi] MHz)",Top],"  ",Labeled[TableForm[{TableForm[#,TableHeadings->{None,{"Eigenvector", "Mode Freq.(x2\[Pi] MHz)"}}]&/@eigenres},TableHeadings->{None,{"Radial","Axial"}}],"Normal Mode Data: "<>ToString[mArr/.{1->"Q",0->"A"}],Top]]];


ShowTrap[]:=With[{trapres={VRF,U,\[Omega]RF/wmhz,r0/mm,z0/mm,qMol},traplabels=If[ValueQ[U],{"\!\(\*SubscriptBox[\(V\), \(RF\)]\) (V)","\!\(\*SubscriptBox[\(V\), \(DC\)]\) (V)","\!\(\*SubscriptBox[\(\[Omega]\), \(RF\)]\) (x2\[Pi] MHz)","r0 (mm)","z0 (mm)","q (Mol)"},{"\!\(\*SubscriptBox[\(V\), \(RF\)]\) (V)","\!\(\*SubscriptBox[\(\[Omega]\), \(RF\)]\) (x2\[Pi] MHz)","r0 (mm)","q (Mol)"}],speciesres={mMol/amu,mAtom/amu,\[Omega]0/wmhz}},
Print[TableForm[{{TableForm[speciesres,TableHeadings->{{"mMol (amu)","mAtom (amu)","\[Omega]0 (x2\[Pi] MHz)"},None}]}},TableHeadings->{None,{"Ion Species"}}],"  ",TableForm[{{TableForm[trapres,TableHeadings->{traplabels,None}]}},TableHeadings->{None,{"Trap Parameters"}}]]];


ShowAll[]:=With[{},
ShowEigen[];
ShowTrap[];];


(* ::Section::Closed:: *)
(*Hamiltonians*)


(* ::Subsection::Closed:: *)
(*Carrier*)


HCarrier[\[CapitalOmega]_,\[Delta]\[Omega]_,\[CapitalDelta]\[Phi]_,ion_Integer,mode_Integer,dir_Integer]:=-\[HBar] \[CapitalOmega]/2*#1 . (#2*\[Sigma]p[[ion]]+#3*\[Sigma]m[[ion]])&[Dni[0,ion,mode,dir],Exp[I( \[CapitalDelta]k[[ion]]*R0[[ion]]-\[CapitalDelta]\[Phi]-\[Delta]\[Omega]*t)],Exp[-I( \[CapitalDelta]k[[ion]]*R0[[ion]]-\[CapitalDelta]\[Phi]-\[Delta]\[Omega]*t)]];


(* ::Subsection::Closed:: *)
(*RSB*)


Hrsb[\[CapitalOmega]_,\[Delta]\[Omega]_,\[CapitalDelta]\[Phi]_,ion_Integer,mode_Integer,dir_Integer]:=-\[HBar] \[CapitalOmega]/2*(#1 . \[Sigma]p[[ion]]+#2 . \[Sigma]m[[ion]])&[Exp[-I( \[CapitalDelta]k[[ion]]*R0[[ion]]- \[CapitalDelta]\[Phi] -\[Delta]\[Omega]*t)]*Dni[-1,ion,mode,dir],-Exp[I( \[CapitalDelta]k[[ion]]*R0[[ion]]- \[CapitalDelta]\[Phi] -\[Delta]\[Omega]*t)]*Dni[1,ion,mode,dir]];


(* ::Subsection::Closed:: *)
(*BSB*)


Hbsb[\[CapitalOmega]_,\[Delta]\[Omega]_,\[CapitalDelta]\[Phi]_,ion_,mode_Integer,dir_Integer]:=-\[HBar] \[CapitalOmega]/2*(#1 . \[Sigma]p[[ion]]+#2 . \[Sigma]m[[ion]])&[Exp[-I( \[CapitalDelta]k[[ion]]*R0[[ion]]- \[CapitalDelta]\[Phi] -\[Delta]\[Omega]*t)]*Dni[1,ion,mode,dir],-Exp[I( \[CapitalDelta]k[[ion]]*R0[[ion]]- \[CapitalDelta]\[Phi] -\[Delta]\[Omega]*t)]*Dni[-1,ion,mode,dir]];


(* ::Subsection::Closed:: *)
(*Bad Stuff*)


Hbad[\[CapitalOmega]_,ion_Integer]:=(-\[HBar]*\[CapitalOmega])/2*(Exp[I*\[Omega]0*t]*\[Sigma]p[[ion]]+Exp[-I*\[Omega]0*t]*\[Sigma]m[[ion]]);


(* ::Subsection::Closed:: *)
(*Higher-order Sidebands*)


HrsbN[\[CapitalOmega]_,\[Delta]\[Omega]_,\[CapitalDelta]\[Phi]_,ord_,ion_Integer,mode_Integer,dir_Integer]:=-\[HBar] \[CapitalOmega]/2*(#1 . \[Sigma]p[[ion]]+#2 . \[Sigma]m[[ion]])&[Exp[-I( \[CapitalDelta]k[[ion]]*R0[[ion]]- \[CapitalDelta]\[Phi] -\[Delta]\[Omega]*t)]*Dni[-ord,ion,mode,dir],(-1)^ord*Exp[I( \[CapitalDelta]k[[ion]]*R0[[ion]]- \[CapitalDelta]\[Phi] -\[Delta]\[Omega]*t)]*Dni[ord,ion,mode,dir]];


HbsbN[\[CapitalOmega]_,\[Delta]\[Omega]_,\[CapitalDelta]\[Phi]_,ord_,ion_,mode_Integer,dir_Integer]:=-\[HBar] \[CapitalOmega]/2*(#1 . \[Sigma]p[[ion]]+#2 . \[Sigma]m[[ion]])&[Exp[-I( \[CapitalDelta]k[[ion]]*R0[[ion]]- \[CapitalDelta]\[Phi] -\[Delta]\[Omega]*t)]*Dni[ord,ion,mode,dir],(-1)^ord*Exp[I( \[CapitalDelta]k[[ion]]*R0[[ion]]- \[CapitalDelta]\[Phi] -\[Delta]\[Omega]*t)]*Dni[-ord,ion,mode,dir]];


(* ::Subsection::Closed:: *)
(*QoL*)


HAll[\[CapitalOmega]_,\[Delta]\[Omega]_,\[CapitalDelta]\[Phi]_,ion_Integer,mode_Integer,dir_Integer]:=-\[HBar] \[CapitalOmega]/2*((\[Sigma]p[[ion]]*#1-\[Sigma]m[[ion]]*#2) . (Dni[1,ion,mode,dir]+Dni[-1,ion,mode,dir])+Dni[0,ion,mode,dir] . (\[Sigma]p[[ion]]*#2+\[Sigma]m[[ion]]*#1))&[Exp[-I( \[CapitalDelta]k[[ion]]*R0[[ion]]- \[CapitalDelta]\[Phi] -\[Delta]\[Omega]*t)],Exp[I( \[CapitalDelta]k[[ion]]*R0[[ion]]- \[CapitalDelta]\[Phi] -\[Delta]\[Omega]*t)]];


(* ::Section::Closed:: *)
(*Updating*)


(*triggers all the other setter functions; mostly for when user uses ChangeAll to change trap parameters*)
UpdateAll[]:=With[{},
UpdateConfig[];
CreateBasis[];
UpdateEigen[];
CreateOperators[];]


ChangeAll[]:=With[{$Context="`Barium`Private`",
(*create dynamic input boxes*)
chaingrid:=Column[{"Chain Config",Row[{"Chain: ",InputField[Dynamic[mArr],Expression]}],Row[{"Active Modes: ",InputField[Dynamic[activeModes],Expression]}],Row[{"Number of Qubits: ",InputField[Dynamic[numQubits],Number]}],Row[{"min Phonon #: ",InputField[Dynamic[minPhonons],Number]}],Row[{"max Phonon #: ",InputField[Dynamic[maxPhonons],Number]}]},":"],
iongrid:=Column[{"Ion Config",Row[{"Molecule Mass: ",InputField[Dynamic[mMol],Number]," amu"}],Row[{"Atom Mass: ",InputField[Dynamic[mAtom],Number]," amu"}],Row[{"\[Omega]0: ",InputField[Dynamic[\[Omega]0],Number]," x2\[Pi] MHz"}]},":"],
lasergrid:=Column[{"Laser Config",Row[{"\[CapitalDelta]k: ",InputField[Dynamic[\[CapitalDelta]k],Expression]," x2\[Pi] \!\(\*SuperscriptBox[\(\[Mu]m\), \(-1\)]\)"}],Row[{"\!\(\*SubscriptBox[\(R\), \(0\)]\): ",InputField[Dynamic[R0],Expression]," \[Mu]m"}]},":"],
trapgrid:=Column[{"Trap Config",Row[{"\!\(\*SubscriptBox[\(V\), \(RF\)]\): ",InputField[Dynamic[VRF],Number]," V"}],Row[{"\!\(\*SubscriptBox[\(\[Omega]\), \(RF\)]\): ",InputField[Dynamic[\[Omega]RF],Number]," x2\[Pi] MHz"}],Row[{"\!\(\*SubscriptBox[\(V\), \(EC\)]\): ",InputField[Dynamic[U],Number]," V"}],Row[{"\!\(\*SubscriptBox[\(r\), \(0\)]\): ",InputField[Dynamic[r0],Number]," mm"}],Row[{"\!\(\*SubscriptBox[\(z\), \(0\)]\): ",InputField[Dynamic[z0],Number]," mm"}],Row[{"\!\(\*SubscriptBox[\(\[Kappa]\), \(r\)]\): ",InputField[Dynamic[\[Kappa]r],Number]}],Row[{"\!\(\*SubscriptBox[\(\[Kappa]\), \(z\)]\): ",InputField[Dynamic[\[Kappa]z],Number]}]},":"]},
Unprotect[mArr,activeModes,numQubits,minPhonons,maxPhonons,mMol,mAtom,\[Omega]0,\[CapitalDelta]k,R0,VRF,\[Omega]RF,U,r0,z0,\[Kappa]r,\[Kappa]z];
(*edit variables for input*)
mMol/=amu;mAtom/=amu;\[Omega]0/=wmhz;\[CapitalDelta]k/=(\[Pi]*2*^6);\[Omega]RF/=wmhz;r0/=mm;z0/=mm;R0*=1*^6;
(*get input*)
DialogInput[{Row[{chaingrid,iongrid,lasergrid,trapgrid},Spacer[5]],Button["Proceed",DialogReturn[]]}];
(*edit variables back to original values*)
mMol*=amu;mAtom*=amu;\[Omega]0*=wmhz;\[CapitalDelta]k*=(\[Pi]*2*^6);\[Omega]RF*=wmhz;r0*=mm;z0*=mm;R0/=1*^6;
Protect[mArr,activeModes,numQubits,minPhonons,maxPhonons,mMol,mAtom,\[Omega]0,\[CapitalDelta]k,R0,VRF,\[Omega]RF,U,r0,z0,\[Kappa]r,\[Kappa]z];
UpdateAll[];]


(*prevents "var" argument in changevariable from being evaluated/converted to an rvalue*)
SetAttributes[ChangeVariable,HoldAll];


(*need to be in private context to change variables; also ensures that we are not simply creating a local variable in the scope of With[]*)
ChangeVariable[var_Symbol,val_]:=With[{$Context="`Barium`Private`",expr=Hold[Unevaluated[var]=val]},
Unprotect[Unevaluated[var]];
ReleaseHold[expr];
Protect[Unevaluated[var]];]


(* ::Section::Closed:: *)
(*Finish*)


Protect["`*"];
ChangeAll[];
End[];


(* ::Chapter::Closed:: *)
(*Finish*)


EndPackage[];
