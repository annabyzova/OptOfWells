(* ::Package:: *)

BeginPackage["optnet`"];


makeTriangulation::usage = "makeTriangulation[f_, nSegments_Integer, area:{{x1_?NumericQ,y1_?NumericQ},{x2_?NumericQ,y2_?NumericQ}}] make Triangulation on area-Rectangle with lower left corner at {x1,y1} and upper right corner at {x2,y2} to n*n rectangles. Then add z=f(x,y)and get triangulation of surface";


Triangulation::usage="Triangulation[g, area, N@Abs[x2-x1]/nSegments,N@Abs[y2-y1]/nSegments, Function[{i,j},node[i,j]],node]";


graph::usage="graph of triangulation";


area::usage="area of triangulation (bounds)";


gridStepX::usage="N@Abs[x2-x1]/nSegments grid step of triangulation for X";


gridStepY::usage="N@Abs[y2-y1]/nSegmentsgrid step of triangulation for Y";


vertex::usage="vertex of triangulation";


(* ::Code::Initialization::Bold:: *)
nodeFinder::usage="by (x,y) function get the vertex number of the graph tr";


graphNodes::usage="array of graph tr nodes ";


distance::usage = "distance[tr_Triangulation][ pt1:{x1_?NumericQ, y1_?NumericQ}, pt2:{x2_?NumericQ, y2_?NumericQ}] calculate the distance between the points pt1 and pt2";


makeSpanTreeK::usage = "make spanning tree with n points and k fork points";


SpanTree::usage="SpanTree[paths, distances, NewPoints, AdjacencyMatrix[h], sppaths, Total[AnnotationValue[{h,#}, EdgeWeight]&/@EdgeList[h]]]";


paths::usage="path from one vertex to other";


distances::usage="distance from one vertex to other";


pointslist::usage="the coordinates of wells";


adjmatr::usage="adjancy matrix of spanning tree";


sppaths::usage="paths in spanning tree";


totalWeight::usage="total weight of spanning tree";


sptreeWeight::usage="njb";


sptreeminimize::usage="minimize the weight of spanning tree";


drawPict::usage="make a 3d plot with surface, points of and paths. Red paths is optimal (with fork points) paths, blue isn't optimal ";


drawTriangl::usage="in process";


basin::usage="in process";


center::usage="in process";


theHill::usage="in process";



vertexPos::usage="in process";


hill::usage = 
  "\:0425\:043e\:043b\:043c \:0432 \:0437\:0430\:0434\:0430\:043d\:043d\:043e\:0439 \:0442\:043e\:0447\:043a\:0435. Height - \:0432\:044b\:0441\:043e\:0442\:0430, Width - \:043f\:0440\:0438\:043c\:0435\:0440\:043d\:044b\:0439 \:0440\:0430\:0434\:0438\:0443\:0441 \
\:043e\:0441\:043d\:043e\:0432\:0430\:043d\:0438\:044f.";


(*surface::usage ="jnhbg";*)


Begin["`Private`"];


makeTriangulation[f_, nSegments_Integer, area:{{x1_?NumericQ,y1_?NumericQ},{x2_?NumericQ,y2_?NumericQ}}]:=
Module[{g,grdist, graphNodes,node,xvert,yvert, edges,
n=nSegments+1, vertex, vertexPosFun},
graphNodes=Flatten@Array[node, {n,n}];
xvert=With[{i=Part[#,1], j=Part[#,2]},
x1+(j-1)*(x2-x1)/nSegments]&;
yvert=With[{i=Part[#,1], j=Part[#,2]},
x1+(i-1)*(x2-x1)/nSegments]&;

vertexPosFun=Function[n, {xvert[n], yvert[n]}];
vertex = Partition[Flatten@Table[ {x1+k*(x2-x1)/nSegments, y1+l*(y2-y1)/nSegments(*, f[x1+k*(x2-x1)/nSegments, y1+l*(y2-y1)/nSegments]*)}, {k, 1, n}, {l, 1, n}],2];

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
Print["This graph from makeTriangulation "];
Print[vertex];
Print[g];

(*Print[FindSpanningTree[g]];*)
(*HighlightGraph[g,FindSpanningTree[g],GraphHighlightStyle\[Rule]"Thick"];*)
(* Print[WeightedAdjacencyMatrix[g]//Normal//MatrixForm]; 
 Print@MatrixPlot[WeightedAdjacencyMatrix[g]]; *)
Triangulation[g, area, N@Abs[x2-x1]/nSegments, N@Abs[y2-y1]/nSegments, Function[{i,j},node[i,j]], node, vertexPosFun]
]

graph[tr_Triangulation]:=tr[[1]]
area[tr_Triangulation]:=tr[[2]]
gridStepX[tr_Triangulation]:=tr[[3]]
gridStepY[tr_Triangulation]:=tr[[4]]
nodeFinder[tr_Triangulation]:=tr[[5]]
graphNodes[tr_Triangulation]:=tr[[6]]
vertexPos[tr_Triangulation]:=tr[[7]]




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
	node2pointOld=With[{i=Part[#,1], j=Part[#,2]},{area[[1,1]]+j*gridStepX[tr], area[[1,2]] + i * gridStepY[tr]}]&,
	node2point=vertexPos[tr]
	},
	{"path"-> Map[node2point, FindShortestPath[graph[tr], nodeFinder[lnear[y1], knear[x1]], nodeFinder[lnear[y2], knear[x2]]]],
	"distance"-> GraphDistance[graph[tr], nodeFinder[lnear[y1], knear[x1]], nodeFinder[lnear[y2], knear[x2]]]}
	]
]  




drawTriangl[tr_Triangulation, surface_,  bounds:{{xl_?NumericQ, yl_?NumericQ},{xr_?NumericQ,yr_?NumericQ}}]:=
Module[{vpf = vertexPos[tr], g = graph[tr], edges, ptl},
	edges = Map[Line[Function[n, With[{xy=vpf[n]}, Append[xy, surface@@xy]]]/@Apply[List, #]]&, 
				EdgeList[g]];
	ptl = vpf/@VertexList[g];
	Print@Show[{
		Plot3D[surface[x,y],{x,xl,xr},{y,yl,yr}, 
				Mesh->None,
				PlotStyle->{LightYellow, Opacity[0.5]},
				ImageSize->Large],
		Graphics3D[{
			Map[{PointSize[Large],Cyan,Point[Join[#,{surface@@#}]]}&,ptl],
			edges
		}]
	}];
	{
	Graphics3D[{
		Map[{PointSize[Medium],Black,Point[Join[#,{surface@@#}]]}&,ptl],
		edges
	}]
	}
]




makeSpanTreeK[tr_Triangulation, points:{{_?NumericQ,_?NumericQ}..}, k: {{_?NumericQ,_?NumericQ}...}]:=
	Module[{edges,distances,g, h, paths, sppaths}, 
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
		adjmatr = adjmatr[st], optpaths=sppaths[optst], optlen = Length[sppaths[optst]], optptl=pointslist[optst]},
		Print["Red line is paths of optimal graph, Magenta points is fork-points"];
		{
			Plot3D[surface[x,y],{x,xl,xr},{y,yl,yr},
					PlotStyle->{Opacity[0.5]},
					Mesh->None,
					ImageSize->Large,
					ColorFunction->Function[{x,y,z}, Hue[(Sin^2)[x]]]],
			Graphics3D[{
			Map[{PointSize[Large],Cyan,Point[Join[#,{surface@@#}]]}&,ptl],
			Map[{PointSize[Medium],Magenta,Point[Join[#,{surface@@#}]]}&,optptl]}],
			Table[{ ListLinePlot3D[ Map[Append[#,surface@@#]&, paths[[t]]   ],  PlotStyle->{Blue, Dashed}]},{t,1,len}],
			Table[{ ListLinePlot3D[ Map[Append[#,surface@@#]&, optpaths[[t]]],  PlotStyle->{Red, Opacity[0.5], Thick}]}, {t,1,optlen}]
			}
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


Options[hill] = {"Height" -> 10, "Width" -> 0.2};

hill[center : {cx_?NumericQ, cy_?NumericQ}, OptionsPattern[] ] :=
   Function[{x, y}, 
   With[{r = Norm[{x - cx, y - cy}],
     h = OptionValue["Height"], w = OptionValue["Width"]},
    h*Exp[-(r)^2/(0.5*w)]]];
    


basin[center : {cx_?NumericQ, cy_?NumericQ}] :=
 With[{theHill = hill[center, "Width" -> 5]},
  Function[{x, y}, -theHill[x, y]]];


End[];


EndPackage[];
