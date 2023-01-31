(* ::Package:: *)

BeginPackage["ComputeFlow`"]

GenerateInitialCondition::usage = "blabala"

GenerateXGrid::usage = "blabla"

GenerateSolvingSymbols::usage = "blabla"

GenerateEquations::usage = "blabla"

getObservables::usage = "calculates all observables of a flow vector u"

flow::usage = "compute the flow"

Begin["`Private`"]

GenerateInitialCondition[initialConditionFunctionIn_, xGrid_, U_, tStart_
	] :=
	Module[{initialConditionFunction = initialConditionFunctionIn, initialConditionTable,
		 initialCondition},
		initialConditionTable = Map[initialConditionFunction, xGrid];
		initialCondition = Thread[Map[#[tStart]&, U] == initialConditionTable
			];
		initialCondition
	]

GenerateXGrid[sigmaMaxIn_, nIn_] :=
	Module[{sigmaMax = sigmaMaxIn, n = nIn, h, xGrid},
		h = sigmaMax / n;
		xGrid = Range[0, sigmaMax, h];
		xGrid
	]

GenerateSolvingSymbols[nIn_, uSymIn_] :=
	Module[{uSym = uSymIn, n = nIn},
		Table[Subscript[uSym, i], {i, 0, n}]
	]

GenerateEquations[uVector_, xGrid_, QIn_, SIn_, uSymb_, tSymb_] :=
	Module[{Q = QIn, S = SIn, uDerivative, qValues, heatFlux, source, n,
		 h},
		n = Length[xGrid];
		h = (xGrid[[-1]] - xGrid[[1]]) / n;
		uDerivative = ListCorrelate[{-1, 1} / h, Join[{-Subscript[uSymb, 1][
			tSymb]}, Map[#[tSymb]&, uVector], {2 * Subscript[uSymb, n - 1][tSymb]
			 - Subscript[uSymb, n - 2][tSymb]}]];
		qValues = Map[Q, uDerivative];
		heatFlux = ListCorrelate[{-1, 1} / h, qValues];
		source = Map[S, xGrid];
		Thread[D[Map[#[tSymb]&, uVector], tSymb] == source + heatFlux]
	]

argRootMin[vectorIn_] :=
	Module[{vector = vectorIn, root = 1},
		Do[
			If[vector[[i]] < 0 && vector[[i + 1]] > 0,
				root = i
			]
			,
			{i, 2, Length[vector] - 1}
		];
		root
	]

getMassSquare[uVectorIn_, dxIn_, rootIn_] :=
	Module[{uVector = uVectorIn, dx = dxIn, root = rootIn, result},
		If[root == 1,
			result = uVector[[2]] / dx
			,
			result = (uVector[[root + 1]] - uVector[[root - 1]]) / 2*dx
		];
		N[result]
	]

getMinimum[uVectorIn_, xGridIn_, rootIn_] :=
	Module[{uVector = uVectorIn, xGrid = xGridIn, root = rootIn, result,
		 lowerSum, upperSum},
		If[root == 1,
			result = xGrid[[1]]
			,
			lowerSum = Last @ FoldList[Plus, 0, uVector[[1 ;; root]]];
			upperSum = Last @ FoldList[Plus, 0, uVector[[1 ;; root + 1]]];
			If[(lowerSum + upperSum) / 2 < 0,
				result = (xGrid[[root]] + xGrid[[root + 1]]) / 2
				,
				result = xGrid[[1]]
			]
		];
		N[result]
	]

getPressure[uVectorIn_, dxIn_, rootIn_] :=
	Module[{uVector = uVectorIn, dx = dxIn, root = rootIn, result, lowerSum,
		 upperSum},
		If[root == 1,
				result = 0
				,
				lowerSum = Last @ FoldList[Plus, 0, uVector[[1 ;; root]]];
				upperSum = lowerSum + uVector[[root + 1]];
				result = -dx * (lowerSum + upperSum) / 2;
			];
		N[result]
	]

getObservables[uVectorIn_, xGridIn_] :=
	Module[{uVector = uVectorIn, xGrid = xGridIn, root = argRootMin[uVectorIn
		], dx = (xGridIn[[-1]] - xGridIn[[1]]) / Length[xGridIn], mass, sigmaMin
		},
		mass = getMassSquare[uVector, dx, root];
		sigmaMin = getMinimum[uVector, xGrid, root];
		pressure = getPressure[uVector, dx, root];
		<|"bosonMassSquare" -> mass, "bosonExpectationValue" -> sigmaMin, "pressure"
			 -> pressure|>
	]

End[]

EndPackage[]
