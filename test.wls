#!/usr/bin/env wolframscript
hill[cx_?NumericQ, cy_?NumericQ] :=
  Function[{x, y}, 
   With[{r = Norm[{x - cx, y - cy}]}, 20*Sin[r]^2*Exp[-r + 1]]];
<< "snet.wl"

With[{surface = Function[{x, y},
    3 Cos[0.4*x*y] + hill[5, 5][x, y]
    ], x0 = 0, y0 = 0, x1 = 20, y1 = 20, gridSize = 50},
 
 Module[{points},
  points = RandomReal[{Max[x0, y0], Min[x1, y1]}, {20, 2}];
  tr = makeTriangulation[surface, gridSize, {{x0, y0}, {x1, y1}}];
  st = makeSpanTreeK[tr, points, 
    RandomReal[{Max[x0, y0], Min[x1, y1]}, {5, 2}]];
  Print["Триангулировано!"];
  
  bestK1 = sptreeminimize[tr, points, {{x0, y0}, {x1, y1}}];
  
  tempfile = FileNameJoin[{$TemporaryDirectory, "saved.wl"}];
  Save[tempfile, bestK1];
  bestK = Get[tempfile];
  (*bestK1={37.95096022172817`,{a\[Rule]8.78037488177198`,
  b\[Rule]13.158424837351651`,c\[Rule]8.18574239113796`,
  d\[Rule]14.067974598082426`}};*)
  originalst = makeSpanTreeK[tr, points, {}];
  st1 = makeSpanTreeK[tr, points, {{a, b}, {c, d}} /. bestK[[2]]];
  
  drawPict[originalst, st1,  surface, {{x0, y0}, {x1, y1}}]
  
  
  ]
 ]

