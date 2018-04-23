(* ::Package:: *)

(**** INTRODUCTION ****)

(** This package serves as a simple Mathematica frontend to the Bertini numerical algebraic geometry software: **)
(**    Bertini(TM) version 1.5.1 **)
(**    Authors: D.J. Bates, J.D. Hauenstein, A.J. Sommese, C.W. Wampler **)
(**    Copyright (C) 2016 **)
(** (available at http://bertini.nd.edu) **)

(** This package is intended for use by those interested in using Bertini to analyze a potential landscape in e.g. string theory. It is capable of: **)
(**    producing input files that Bertini can interpret, **)
(**    reading output files that Bertini has created, **)
(**    performing simple analyses of the characteristics of the potential at its vacua. **)

(** There are doubtless bugs in this package. If you find one, I would appreciate you letting me know **)
(** This package was originally created by Jeremy M Wachter (github: fvnu ; web: www.jmwachter.com). It is provided without guarantee or warranty, and is released under the CC BY-SA 4.0 license (https://creativecommons.org/licenses/by-sa/4.0/) **)



BeginPackage["Eugenio`"]

EugenioVersion = "1.1.0"

(**** CONSTANTS ****)

OutputFiles={"raw_solutions","finite_solutions","real_finite_solutions","nonsingular_solutions","singular_solutions"};

(**** OPTIONS, USAGE, AND ERRORS FOR CERTAIN FUNCTIONS ****)

MakeInputFile::usage = "MakeInputFile[V] takes a potential function 'V' in N variables and produces a Bertini input file whose N functions are the gradient of 'V'. You may pass certain Bertini config options as options to this function."
Options[MakeInputFile]={FileName->"input",VariableGrouping->{},FinalTol->Power[10,-11],ImagThreshold->Power[10,-8],MPType->0,Precision->96,CoeffBound->1000,DegreeBound->5,RandomSeed->0,ScreenOut->0};
FileName::usage = "FileName an option for MakeInputFile, which is the name of the file produced. It defaults to 'input'."
VariableGrouping::usage = "VariableGrouping is an option for MakeInputFile, which instructs how the variables are to be grouped and if they are homogeneous or not. The proper format is: {N,{x_i..x_j},{x_j+1..x_k}..}, where N is 0 for nonhomogeneous groups and 1 is for homogeneous groups. Each list of variables belongs to a different group."
MakeInputFile::vgtype = "The first argument of VariableGrouping must be 0 (for nonhomogeneous variable groups) or 1 (for homogeneous)."

ReadOutputFile::usage = "ReadOutputFile[filename] reads a solutions file produced by Bertini. The file can have any name, but must be in one of the standard forms of: finite_solutions, nonsingular_solutions, raw_solutions, real_finite_solutions, or singular_solutions."
Options[ReadOutputFile]={AllReal->False};
AllReal::usage = "AllReal is an option for ReadOutputFile which can be False (default) or True. If set to True, it drops the imaginary parts of all vacua."

FindVacua::usage = "HOLD"
Options[FindVacua]={BertiniPath->"",SaveOutputFiles->{},RenameOutputFiles->""};
BertiniPath::usage = "BertiniPath is an option for FindVacua which tells Eugenio where the 'bertini.sh' file lives. It should be given without the trailing '/'."
SaveOutputFiles::usage = "SaveOutputFiles is an option for FindVacua which is a list of length 2. The first entry should be the location to which you wish to save the output files, and the second entry should be which of the output files you wish to save. This second entry should always use the list of possible output files known to Bertini (see the variable 'OutputFiles' for a full list), even if you are using the RenameOutputFiles option at the same time."
RenameOutputFiles::usage = "RenameOutputFiles is an option for FindVacua which appends to all of the Bertini output files (as given in the variable 'OutputFiles') the string given to it as an argument."
FindVacua::badsave = "The option 'SaveOutputFiles' is malformed. Check to ensure that it is a two-element list with the first element being the directory to save files to, and the second being a list of the Bertini output files you wish to save."
FindVacua::makefail = "FindVacua failed while trying to make an output file with input: `1`."
FindVacua::nobertini = "The command 'bertini' could not be executed. If you set 'BertiniPath', make sure it points to the folder containing the Bertini executable. If you did not set 'BertiniPath', make sure you have put the Bertini executable in your /bin/ folder."

FindEigenvalueList::usage = "FindEigenvalueList[{V,k},vacua,A] finds the generalized eigenvalues of the potential function 'V' with respect to the matrix 'k' at a list of minima 'vacua'. It is optional to provide a vector 'A', in which case the Hessian is constructed with the covariant derivative D_i = partial_i + A_i. You may also provide 'V' in place of '{V,k}', in which case the regular eigenvalues are found."
Options[FindEigenvalueList]={MatrixNorm->Identity};
MatrixNorm::usage = "MatrixNorm is an option for FindEigenvalueList, which maps the Hessian of V to another matrix of identical dimension which is then used as the LHS of the eigenvalue problem."

GetTadpole::usage = "GetTadpole[Q,n] returns a list of n integers which satisfy the condition that the sum of the squares of the integers is less than or equal to Q. Note that Q can be any real positive number."
GetTadpole::lown = "N must be an integer greater than 1";


(**** I/O FUNCTIONS ****)

(** Creates an input file for Bertini given a potential function 'V', and optionally a vector 'A' by which we form the covariant derivative D_i = partial_i + A_i **)
(** Alternately, accepts as 'V' a list of the polynomials representing the gradients and constraint equations which comprise the polynomial system to solve **)
(** The user may specify a subset of the full set of settable parameters of Bertini, which is to say, the ones I thought might be useful, or needed for some reason; they are named exactly as Bertini understands them, but in CamelCase **)
MakeInputFile[V_,opts:OptionsPattern[]]:=MakeInputFile[V,Table[0,Length[Variables@V]],Sequence[opts]];
MakeInputFile[V_,A_,OptionsPattern[]]:=Module[{fileName=OptionValue[FileName],pNames,X=Sort@Variables@V,P},
    If[FileExistsQ[fileName],
        Print["File '",fileName,"' already exists!"];Return[False,Module]];
    If[ListQ[V],
        P=V;X=Sort[Union@@Variables/@P],
        P=PotentialToGradient[V,A];X=Sort@Variables@V;
            If[MemberQ[Variables/@P,{}],
                Print["One or more of the gradients are constant. Exiting."];Return[False,Module]]];
    pNames=Table["f"<>ToString[i],{i,Length[P]}];
    s=CreateFile[fileName];
    WriteLine[s,"CONFIG"];
    If[OptionValue[FinalTol]>0,
        WriteLine[s,"FINALTOL: "<>ToString[FortranForm[OptionValue[FinalTol]]]<>";"],
        Print["FinalTol must be >0"]];
    If[OptionValue[ImagThreshold]>0,
        WriteLine[s,"IMAGTHRESHOLD: "<>ToString[FortranForm[OptionValue[ImagThreshold]]]<>";"],
        Print["ImagThreshold must be >0"]];
    If[MemberQ[{0,1,2},OptionValue[MPType]],
        WriteLine[s,"MPTYPE: "<>ToString[OptionValue[MPType]]<>";"];
        If[OptionValue[MPType]==1,
            If[And[IntegerQ[OptionValue[Precision]],64<=OptionValue[Precision]<=3328],
                WriteLine[s,"PRECISION: "<>ToString[OptionValue[Precision]]<>";"],
                Print["Precision must be an integer in [64,3328]"]]];
        If[OptionValue[MPType]==2,
            WriteLine[s,"COEFFBOUND: "<>ToString[OptionValue[CoeffBound]]<>";"];
            WriteLine[s,"DEGREEBOUND: "<>ToString[OptionValue[DegreeBound]]<>";"];
            If[Or[0>OptionValue[CoeffBound],0>OptionValue[DegreeBound]],
                Print["CoeffBound and DegreeBound must be >=0. Using default values"]]],
        Print["MPType must be 0, 1, or 2"]];
    If[OptionValue[RandomSeed]>=0,
        WriteLine[s,"RANDOMSEED: "<>ToString[FortranForm[OptionValue[RandomSeed]]]<>";"],
        Print["RandomSeed must be >=0"]];
    If[MemberQ[{0,1},OptionValue[ScreenOut]],
        WriteLine[s,"SCREENOUT: "<>ToString[OptionValue[ScreenOut]]<>";"],
        Print["ScreenOut must be 0 or 1"]];
    WriteLine[s,"END;"];
    WriteLine[s,"INPUT"];
    vg=OptionValue[VariableGrouping];
    If[vg=={},
        WriteLine[s,"variable_group "<>StringTake[#,{2,Length[#]-2}]&@ToString[X]<>";"],
        If[Or[vg[[1]]==0,vg[[1]]==1],
            Do[WriteLine[s,If[vg[[1]]==0,"variable_group ","hom_variable_group "]<>StringTake[#,{2,Length[#]-2}]&@ToString[vg[[i]]]<>";"],{i,2,Length[vg]}],
            Message[MakeInputFile::vgtype];Return[0,Module]]];
    WriteLine[s,"function "<>StringTake[#,{2,Length[#]-2}]&@ToString[pNames]<>";"];
    Do[WriteLine[s,pNames[[i]]<>" = "<>ToString[P[[i]],InputForm]<>";"],{i,Length[pNames]}];
    WriteLine[s,"END;"];
    Close[s];
    Return[True];
];


(** Reads the output file specified by 'fileName' and returns a list of the locations of the vacua reported by that output file **)
(** Because Bertini reports as real those points with all imaginary parts below IMAGTHRESHOLD, but does not write those imaginary parts as zero, you may set the option 'AllReal' to True in order to discard all imaginary components **)
(** Note that the output file must be in the current working directory, or must be specified such that OpenRead can find it **)
(** The output file can have any name, but it must be exactly in the format of the Bertini output files by name of: finite_solutions ; nonsingular_solutions ; raw_solutions ; real_finite_solutions ; singular_solutions **)
ReadOutputFile[fileName_,OptionsPattern[]]:=Module[{},
    If[FileExistsQ[fileName],
        s=OpenRead[fileName],
        Print["File '",fileName,"' doesn't exist!"];Return[False,Module]];
    results=ReadList[s,Number];
    Close[s];
    nVacua=results[[1]];
    If[nVacua==0,
        Print["No results in '",fileName,"'. Exiting."];Return[False,Module]];
    If[IntegerQ[results[[2]]],
        results=Delete[results,Table[{2+((Length[results]-1)/nVacua)*(i-1)},{i,nVacua}]]];
    nFields=(Length[results]-1)/(2*nVacua);
    listResults=ArrayReshape[Table[#[[2*k-1]]+I*#[[2*k]],{k,(Length[results]-1)/2}]&@results[[2;;Length[results]]],{nVacua,nFields}];
    If[OptionValue[AllReal],listResults=ImTrim[listResults,Infinity]];
    Return[listResults];
];

(** Removes the imaginary parts of each coordinate of a list of points LIST which are below CUTOFF **)
ImTrim[list_,cutoff_:Power[10,-8]]:=Map[If[Abs[Im[#]]<cutoff,Re[#]]&,list,{2}];

(** Creates an input file, runs Bertini on the input file, and renames the output file(s) according to the name of the input file **)
FindVacua[V_,opts:OptionsPattern[{FindVacua,MakeInputFile}]]:=FindVacua[V,Table[0,{Length[Variables@V]}],Sequence[opts]];
FindVacua[V_,A_,opts:OptionsPattern[{FindVacua,MakeInputFile}]]:=Module[{LocPrefix,Suffix,SaveFiles,SavePath},
    LocPrefix=OptionValue[BertiniPath];
    Suffix=OptionValue[RenameOutputFiles];
    If[LocPrefix!="/home/jmw/BertiniLinux64_v1.5.1/",Message[FindVacua::nobertini];Return[False,Module]];
    (* If the input file can't be created, warn and exit *)
    makeFlag=MakeInputFile[V,A,FilterRules[{opts},Options[MakeInputFile]]];
    If[Not[makeFlag],Message[FindVacua::makefail,{V,A}];Return[False,Module]];
    (* If the user doesn't specify anything for SaveOutputFiles, do the default behavior. Otherwise, get their choices. *)
    sof = OptionValue[SaveOutputFiles];
    If[sof=={},
        SavePath="";SaveFiles=OutputFiles,
        If[Length[sof]==2,
            SavePath=sof[[1]];SaveFiles=sof[[2]],
            Message[FindVacua::badsave];Return[False,Module]]];
    (* If the user provided any files names to SaveOutputFiles which aren't Bertini output files, warn them and exit *)
    If[Complement[SaveFiles,OutputFiles]!={},Message[FindVacua::badsave];Return[False,Module]];
    (* Who can ever remember if you're supposed to have a trailing / or not, anyhow? *)
    If[StringTake[LocPrefix,-1]!="/",LocPrefix=LocPrefix<>"/"];
    (* Actually run Bertini now, but warn and exit if it can't be found *)
    Quiet[Check[RunProcess[{LocPrefix<>"bertini",OptionValue[FileName]}],Message[FindVacua::nobertini];Return[False,Module],RunProcess::pnfd],RunProcess::pnfd];
    (* Rename all of the files if desired, changing the contents of SaveFiles as we go. *)
    If[Suffix!="",
        For[i=1,i<=Length[SaveFiles],i++,
            If[RunProcess[{"mv",SaveFiles[[i]],SaveFiles[[i]]<>Suffix}][[3]]!="",Message[FindVacua::badsave]];
            SaveFiles[[i]]=SaveFiles[[i]]<>Suffix]];
    If[SavePath!="",
        For[i=1,i<=Length[SaveFiles],i++,
            If[RunProcess[{"mv",SaveFiles[[i]],SavePath}][[3]]!="",Message[FindVacua::badsave]]]];
];

(**** ANALYSIS FUNCTIONS ****)

(** Given a potential function 'V', returns the gradient **)
(** If you provide an optional one-dimensional array 'A' of the correct size, instead performs the covariant derivative D_i = partial_i + A_i **)
PotentialToGradient[V_]:=PotentialToNthDerivative[V,1,Table[0,Length[Variables@V]]];
PotentialToGradient[V_,A_]:=PotentialToNthDerivative[V,1,A];

(** Given a potential function 'V', returns the Hessian **)
(** If you provide an optional one-dimensional array 'A' of the correct size, instead performs the covariant derivative D_i = partial_i + A_i **)
PotentialToHessian[V_]:=PotentialToNthDerivative[V,2,Table[0,Length[Variables@V]]];
PotentialToHessian[V_,A_]:=PotentialToNthDerivative[V,2,A];

(** Given a potential function 'V' and a vector 'A', returns the n-index covariant derivative array, where the covariant derivative is define as D_i = partial_i + A_i **)
PotentialToNthDerivative[V_,n_,A_]:=With[{X=Sort@Variables@V},Nest[Table[D[#,X[[i]]]+A[[i]]*#,{i,Length[X]}]&,V,n]];

(** Given a potential function 'V' and a list of vacuum points 'vacua', returns a list of the eigenvalues associated with each set of points **)
(** If you want to solve the generalized eigenvalue problem instead, you can provide the function with {V,k} instead of V, where 'k' is a matrix of the correct size **)
(** You may optionally specify the normalization of the Hessian matrix before its eigenvalues are found by providing the desired function via the 'MatrixNorm' option **)
FindEigenvalueList[Vk_,vacua_,opts:OptionsPattern[]]:=FindEigenvalueList[Vk,vacua,Table[0,Length[Variables@If[ListQ[Vk],Vk[[1]],Vk]]],FilterRules[{opts},Options[FindEigenvalueList]]];
FindEigenvalueList[Vk_,vacua_,A_,OptionsPattern[]]:=Module[{V,X,k,norm=OptionValue[MatrixNorm]},
    If[ListQ[Vk],
        V=Vk[[1]];X=Sort@Variables@V;k=Vk[[2]],
        V=Vk;X=Sort@Variables@V;k=IdentityMatrix[Length@X]];
    Eigenvalues[{norm@#,k}]&/@(PotentialToHessian[V,A]/.(Thread[X->#]&/@vacua))];


(**** GENESIS FUNCTIONS ****)

(** Given an (optional) upper limit 'Q' and an (optional) number of parameters 'n', returns a list of 'n' parameters which satisfy the tadpole condition: **)
(** \sum_i x_i^2 <= Q, i\in{1,n}, each x_i an integer **)
GetTadpole[Q_:100,n_:4]:=Module[{tp},
    If[n<2,Message[GetTadpole::lown];Return[False,Module]];
    While[True, 
        tp = Round[(((Sqrt[Q]+Sqrt[n])*Power[RandomReal[],1/n]/Norm[#])*#)&@RandomVariate[NormalDistribution[],n]];
        If[Norm[tp]<Sqrt[Q],Break[]]
    ];
    tp
];

(** Generates a list of for a toy model of supergravity of two complex scalar fields z and t (tau) representing the modulus and the axio-dilaton, respectively **)
(** The variables are thus z, t, their complex conjugates (zb and tb), and the auxillary real variables ax and ay (which prevent solutions with z=zb and t=tb) **)
(** Takes an argument 'Q' which is the upper limit of the tadpole constraint **)
ToyModelZT[pList_:GetTadpole[]]:=Module[{m,n,p,q,f1,f2,f3,f4,f5,f6},
    {m,n,p,q}=pList;
    f1 = (m-t*n)+(p+t*q)*zb+(1/2)*(q-t*p)*zb^2+(1/6)*(n+t*m)*zb^3;
    f2 = (m-tb*n)+(p+tb*q)*z+(1/2)*(q-tb*p)*z^2+(1/6)*(n+tb*m)*z^3;
    f3 = -3*(m-t*n)-(p+t*q)*(zb+2*z)-(q-t*p)*(zb+(1/2)*z)*z-(1/2)*(n+t*m)*z^2*zb;
    f4 = -3*(m-tb*n)-(p+tb*q)*(z+2*zb)-(q-tb*p)*(z+(1/2)*zb)*zb-(1/2)*(n+tb*m)*zb^2*z;
    f5 = 1-I*ax*(t-tb);
    f6 = 1-I*ay*(z-zb);
    {f1,f2,f3,f4,f5,f6}
];

EndPackage[]
