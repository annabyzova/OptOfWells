(* ::Package:: *)

BeginPackage["optnet`"];


makeTriangulation::usage = "makeTriangulation[f_, nSegments_Integer, area:{{x1_?NumericQ,y1_?NumericQ},{x2_?NumericQ,y2_?NumericQ}}] make Triangulation on area-Rectangle with lower left corner at {x1,y1} and upper right corner at {x2,y2} to n*n rectangles. Then add z=f(x,y)and get triangulation of surface";


Triangulation::usage="jnhb";


WellsOnTr::usage="zxcv";


putWells::usage ="zxcv";


graph::usage="jnhb";


area::usage="jnhb";


gridStepX::usage="jnhb";


gridStepY::usage="jnhb";


(* ::Code::Initialization::Bold:: *)
nodeFinder::usage="jnhb";


graphNodes::usage="jnhb";


distance::usage = "distance[tr_Triangulation][ pt1:{x1_?NumericQ, y1_?NumericQ}, pt2:{x2_?NumericQ, y2_?NumericQ}] calculate the distance between the points pt1 and pt2";


makeSpanTreeK::usage = "make spanning tree with n points and k fork points";


SpanTree::usage="njb";


paths::usage="njb";


distances::usage="njb";


pointslist::usage="njb";


adjmatr::usage="njb";


sppaths::usage="njb";


totalWeight::usage="njb";


sptreeWeight::usage="njb";


sptreeminimize::usage="hgfr";


drawPict::usage="make a 3d plot with surface, points of and paths. Red paths is optimal (with fork points) paths, blue isn't optimal ";


(*surface::usage ="jnhbg";*)


Begin["`Private`"];


makeTriangulation[f_, nSegments_Integer, area:{{x1_?NumericQ,y1_?NumericQ},{x2_?NumericQ,y2_?NumericQ}}]:=
Module[{g,grdist,graphNodes,node,xvert,yvert, edges,
n=nSegments+1},
graphNodes=Flatten@Array[node, {n,n}];
xvert=With[{i=Part[#,1], j=Part[#,2]},
x1+j*(x2-x1)/nSegments]&;
yvert=With[{i=Part[#,1], j=Part[#,2]},
x1+i*(x2-x1)/nSegments]&;
grdist=Function[{n1,n2},
N@Norm[{xvert[n1]-xvert[n2],
yvert[n1]-yvert[n2],
f[xvert[n1],yvert[n1]]-f[xvert[n2],yvert[n2]]}]];
edges = Join[
Flatten@Table[node[i,j]\[UndirectedEdge]node[i,j+1],{i,n},{j,1,n-1}],
Flatten@Table[node[i,j]\[UndirectedEdge]node[i+1,j],{i,1,n-1},{j,n}],
Flatten@Table[node[i,j]\[UndirectedEdge]node[i+1,j-1],{i,1,n-1},{j,2,n}]
];
g=Graph[edges,EdgeWeight->{e:UndirectedEdge[n1_,n2_]:>grdist[n1,n2]}];
Print[g];
(*Print[FindSpanningTree[g]];*)
(*HighlightGraph[g,FindSpanningTree[g],GraphHighlightStyle\[Rule]"Thick"];*)
(* Print[WeightedAdjacencyMatrix[g]//Normal//MatrixForm]; 
 Print@MatrixPlot[WeightedAdjacencyMatrix[g]]; *)
Triangulation[g, area, N@Abs[x2-x1]/nSegments,N@Abs[y2-y1]/nSegments, Function[{i,j},node[i,j]],node]
]

graph[tr_Triangulation]:=tr[[1]]
area[tr_Triangulation]:=tr[[2]]
gridStepX[tr_Triangulation]:=tr[[3]]
gridStepY[tr_Triangulation]:=tr[[4]]
nodeFinder[tr_Triangulation]:=tr[[5]]
graphNodes[tr_Triangulation]:=tr[[6]]



putWells[tr_Triangulation, points:{{_?NumericQ,_?NumericQ}..}]:= Module[{},

TrWells[tr, edges]
]


distance::ErrorArea="\:0412\:044b\:0445\:043e\:0434 \:0437\:0430 \:0433\:0440\:0430\:043d\:0438\:0446\:0443 \:043e\:0431\:043b\:0430\:0441\:0442\:0438 `1` \:043f\:0440\:0438 \:0432\:044b\:0447\:0438\:0441\:043b\:0435\:043d\:0438\:0438 \:0440\:0430\:0441\:0441\:0442\:043e\:044f\:043d\:0438\:044f \:043c\:0435\:0436\:0434\:0443 `2` \:0438 `3`.";

distance[tr_Triangulation][ pt1:{x1_?NumericQ, y1_?NumericQ}, pt2:{x2_?NumericQ, y2_?NumericQ}]:= 
distance[tr][pt1,pt2]=
With[{area=area[tr], nodeFinder=nodeFinder[tr]},
	If[Not@And[area[[1,1]]<= x1 <= area[[2,1]], area[[1,1]] <= x2 <= area[[2,1]], area[[1,2]] <= y1 <= area[[2,2]], area[[1,2]] <= y2 <= area[[2,2]]],
	Message[distance::ErrorArea, area, pt1, pt2]; Abort[]];
	(*cccc=cccc+1;*)
	With[{
	knear=Function[x, 1+Round[(x-area[[1,1]])/gridStepX[tr]]],
	lnear=Function[y, 1+Round[(y-area[[1,2]])/gridStepY[tr]]],
	node2point=With[{i=Part[#,1], j=Part[#,2]},{area[[1,1]]+j*gridStepX[tr], area[[1,2]] + i * gridStepY[tr]}]& 
	},
	{"path"-> Map[node2point, FindShortestPath[graph[tr], nodeFinder[lnear[y1], knear[x1]], nodeFinder[lnear[y2], knear[x2]]]],
	"distance"-> GraphDistance[graph[tr], nodeFinder[lnear[y1], knear[x1]], nodeFinder[lnear[y2], knear[x2]]]}
	]
]  




makeSpanTreeK[tr_Triangulation, points:{{_?NumericQ,_?NumericQ}..}, k: {{_?NumericQ,_?NumericQ}...}]:=
	Module[{edges,distances,g, paths, sppaths}, 
		With[{NewPoints = Join[points, k], nPoints = Length[points] + Length[k]},
			edges = Join@@Table[If[i!=j, UndirectedEdge[i, j]],{i, 1, nPoints}, {j, 1, nPoints}];
			edges = DeleteCases[edges, Null];
			distances = Map[With[{i=#[[1]], j=#[[2]]}, "distance"/.distance[tr][NewPoints[[i]], NewPoints[[j]]]]&, edges];
			paths = Map[With[{i=#[[1]], j=#[[2]]},"path"/.distance[tr][NewPoints[[i]], NewPoints[[j]]]]&, edges];
			g = Graph[edges, EdgeWeight->distances];
			h = FindSpanningTree[g];
			sppaths = Map[With[{i=#[[1]], j=#[[2]]},"path"/.distance[tr][NewPoints[[i]], NewPoints[[j]]]]&, EdgeList[h]];
			(*Total[AnnotationValue[{h,#},EdgeWeight]&/@EdgeList[h]];*)
			
			SpanTree[paths, distances, NewPoints, AdjacencyMatrix[h], sppaths, Total[AnnotationValue[{h,#}, EdgeWeight]&/@EdgeList[h]]]
		]
	]

paths[st_SpanTree]:=st[[1]]
distances[st_SpanTree]:=st[[2]]
pointslist[st_SpanTree]:=st[[3]]
adjmatr[st_SpanTree]:=st[[4]]
sppaths[st_SpanTree]:=st[[5]]
totalWeight[st_SpanTree]:=st[[6]]


sptreeWeight[tr_Triangulation,points:{{_?NumericQ,_?NumericQ}..}, k: {{_?NumericQ,_?NumericQ}...}]:=
With[{st=makeSpanTreeK[tr, points, k]},
totalWeight[st]
]


drawPict[st_SpanTree, optst_SpanTree, surface_, bounds:{{xl_?NumericQ, yl_?NumericQ},{xr_?NumericQ,yr_?NumericQ}}]:= 
	With[{paths=sppaths[st],len = Length[sppaths[st]], distances=distances[st], ptl=pointslist[st], 
		adjmatr = adjmatr[st], optpaths=sppaths[optst], optlen = Length[sppaths[optst]]},
		Print@Show[{
			Plot3D[surface[x,y],{x,xl,xr},{y,yl,yr}],
			Graphics3D[
			Map[{PointSize[Large],Black,Point[Join[#,{surface@@#}]]}&,ptl]],
			Table[{ ListLinePlot3D[ Map[Append[#,surface@@#]&,paths[[t]] ]]}, {t,1,len}],
			Table[{ ListLinePlot3D[ Map[Append[#,surface@@#]&,optpaths[[t]] ],  PlotStyle->{Red, Opacity[0.5]}]}, {t,1,optlen}]
		}];
]



sptreeminimize[tr_Triangulation, points_, borders:{{x0_, y0_}, {x1_, y1_}}, numCrossings_Integer/;numCrossings>0]:=
Module[{varPoints=Partition[Table[Unique[],{numCrossings*2}], 2],
	bounds
	},
	bounds=Join[
		Thread[LessEqual[x0, First/@varPoints, x1]],
		Thread[LessEqual[y0, Last/@varPoints, y1]]
		];
	Print["Bounds: ", bounds];
	NMinimize[{sptreeWeight[tr, points, varPoints],
	   		And@@bounds},
		Flatten@varPoints,
		Method->"SimulatedAnnealing",
		PrecisionGoal->2,
		AccuracyGoal->2,
		EvaluationMonitor:>Print[varPoints]]
];


(*Remove[surface];
surface[a_?NumericQ]:= Function[{x,y},
3Cos[a*x*y]+hill[5,5][x,y]
];*)


End[];


EndPackage[];
