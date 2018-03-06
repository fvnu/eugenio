(* ::Package:: *)

(**** INTRODUCTION ****)

(** This package serves as a simple Mathematica frontend to the Bertini numerical algebraic geometry software (http://bertini.nd.edu) for those interested in using that program to analyze a potential landscape in e.g. string theory. It is capable of: **)
(**    producing input files that Bertini can interpret, **)
(**    reading output files that Bertini has created, **)
(**    performing simple analyses of the characteristics of the potential at its vacua. **)
(** There are doubtless bugs in the program. If you find one, I would appreciate you letting me know. But, I do not generally consider errors thrown by Mathematica due to incorrect user input to be bugs **)
(** This package was originally created by Jeremy M Wachter (www.jmwachter.com). It is provided without guarantee or warranty, and is released under the CC BY-SA 4.0 license (https://creativecommons.org/licenses/by-sa/4.0/) **)



BeginPackage["Eugenio`"]

EugenioVersion = "1.0.0"

(**** OPTIONS AND USAGE FOR CERTAIN FUNCTIONS ****)

MakeInputFile::usage = "MakeInputFile[V] takes a potential function 'V' in N variables and produces a Bertini input file whose N functions are the gradient of 'V'. You may pass certain Bertini config options as options to this function."
Options[MakeInputFile]={FileName->"input",FinalTol->Power[10,-11],ImagThreshold->Power[10,-8],MPType->0,Precision->96,CoeffBound->1000,DegreeBound->5,RandomSeed->0,ScreenOut->0};
FileName::usage = "FileName an option for MakeInputFile, which is the name of the file produced. It defaults to 'input'."

ReadOutputFile::usage = "ReadOutputFile[filename] reads a solutions file produced by Bertini. The file can have any name, but must be in one of the standard forms of: finite_solutions, nonsingular_solutions, raw_solutions, real_finite_solutions, or singular_solutions."
Options[ReadOutputFile]={AllReal->False};
AllReal::usage = "AllReal is an option for ReadOutputFile which can be False (default) or True. If set to True, it drops the imaginary parts of all vacua."

FindEigenvalueList::usage = "FindEigenvalueList[{V,k},vacua,A] finds the generalized eigenvalues of the potential function 'V' with respect to the matrix 'k' at a list of minima 'vacua'. It is optional to provide a vector 'A', in which case the Hessian is constructed with the covariant derivative D_i = partial_i + A_i. You may also provide 'V' in place of '{V,k}', in which case the regular eigenvalues are found."
Options[FindEigenvalueList]={MatrixNorm->Identity};
MatrixNorm::usage = "MatrixNorm is an option for FindEigenvalueList, which maps the Hessian of V to another matrix of identical dimension which is then used as the LHS of the eigenvalue problem."


(**** I/O FUNCTIONS ****)

(** Creates an input file for Bertini given a potential function 'V', and optionally a vector 'A' by which we form the covariant derivative D_i = partial_i + A_i **)
(** The user may specify a subset of the full set of settable parameters of Bertini, which is to say, the ones I thought might be useful, or needed for some reason; they are named exactly as Bertini understands them, but in CamelCase **)
MakeInputFile[V_,opts:OptionsPattern[]]:=MakeInputFile[V,Table[0,Length[Variables@V]],FilterRules[{opts},Options[MakeInputFile]]];
MakeInputFile[V_,A_,OptionsPattern[]]:=Module[{fileName=OptionValue[FileName],pNames,X=Sort@Variables@V,P},
    If[FileExistsQ[fileName],
        Print["File '",fileName,"' already exists!"];Return[0,Module]];
    P=PotentialToGradient[V,A];
    If[MemberQ[Variables/@P,{}],
        Print["One or more of the gradients are constant. Exiting."];Return[0,Module]];
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
    WriteLine[s,"variable_group "<>StringTake[#,{2,Length[#]-2}]&@ToString[X]<>";"];
    WriteLine[s,"function "<>StringTake[#,{2,Length[#]-2}]&@ToString[pNames]<>";"];
    Do[WriteLine[s,pNames[[i]]<>" = "<>ToString[P[[i]],InputForm]<>";"],{i,Length[pNames]}];
    WriteLine[s,"END;"];
    Close[s];
    Return[1];
];

(** Reads the output file specified by 'fileName' and returns a list of the locations of the vacua reported by that output file **)
(** Because Bertini reports as real those points with all imaginary parts below IMAGTHRESHOLD, but does not write those imaginary parts as zero, you may set the option 'AllReal' to True in order to discard all imaginary components **)
(** Note that the output file must be in the current working directory, or must be specified such that OpenRead can find it **)
(** The output file can have any name, but it must be exactly in the format of the Bertini output files by name of: finite_solutions ; nonsingular_solutions ; raw_solutions ; real_finite_solutions ; singular_solutions **)
ReadOutputFile[fileName_,OptionsPattern[]]:=Module[{},
    If[FileExistsQ[fileName],
        s=OpenRead[fileName],
        Print["File '",fileName,"' doesn't exist!"];Return[0,Module]];
    results=ReadList[s,Number];
    Close[s];
    nVacua=results[[1]];
    If[nVacua==0,
        Print["No results in '",fileName,"'. Exiting."];Return[0,Module]];
    If[IntegerQ[results[[2]]],
        results=Delete[results,Table[{2+((Length[results]-1)/nVacua)*(i-1)},{i,nVacua}]]];
    nFields=(Length[results]-1)/(2*nVacua);
    listResults=ArrayReshape[Table[#[[2*k-1]]+I*#[[2*k]],{k,(Length[results]-1)/2}]&@results[[2;;Length[results]]],{nVacua,nFields}];
    If[OptionValue[AllReal],listResults=ImTrim[listResults,Infinity]];
    Return[listResults];
];

(** Removes the imaginary parts of each coordinate of a list of points LIST which are below CUTOFF **)
ImTrim[list_,cutoff_]:=Map[If[Abs[Im[#]]<cutoff,Re[#]]&,list,{2}];



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

EndPackage[]
