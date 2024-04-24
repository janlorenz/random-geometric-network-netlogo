to setup-compute-deg-from-N-r ; computes the deg from N and r before setup
  set deg (N - 1) * prob-length-less-r r ; (N - 1) * (1 / 2 * r ^ 4 - 8 / 3 * r ^ 3 + pi * r ^ 2)
  setup
end

; PROBABLY THIS IS THE INTERESTING CASE
to setup-compute-r-from-N-deg  ; computes the r from N and deg before setup
  set r (r-from-N-deg N deg)
  setup
end

to setup ; produces the network
  clear-all
  ask patches [set pcolor white]
  create-turtles N [ setxy random-xcor random-ycor ]
  ask turtles [
    create-links-with other turtles with [distance myself < r * (max-pxcor + 1)] [set hidden? hide-links]
  ] ;; We use (max-pxcor + 1) because this is the length of the square of the world
    ;; That way r * (max-pxcor + 1) matches distances of r if the world was a unit square as in
    ;; the standard random geometric graph
  update-plots
end

to-report r-from-N-deg [N_ deg_]
  ; This implements the Newton's method (https://en.wikipedia.org/wiki/Newton%27s_method) for finding the root of f
  ; which represent how to compute the radius r from N and the desired average degree
  let r_ sqrt (deg_ / (pi * N_)) ; start with this proxy for r
  let rn r_ - (f r_ N_ deg_) / (df r_) ; one step of Newton's method to improve toward a root
  while [rn - r_ > 10 ^ -10] [ ; repeat until "convergence" (distance to new value very small)
    set r_ rn
    set rn r_ - (f r_ N_ deg_) / (df r_)
  ]
  report rn
end


to-report prob-length-less-r [r_]
  ; referred to as Prob(l-r) in Info tab
  report (r_ ^ 4 / 2 - 8 / 3 * r_ ^ 3 + pi * r_ ^ 2)
end


to-report f [r_ N_ deg_]
  ; referred to as f in Info tab
  report prob-length-less-r r_ - (deg_ / (N_ - 1))
end

to-report df [x]
  ; The derivative of f to be used for the Newton method in r-from-N-deg
  report 2 * x ^ 3 - 8 * x ^ 2 + 2 * pi * x
end
@#$#@#$#@
GRAPHICS-WINDOW
215
10
734
530
-1
-1
12.463415
1
10
1
1
1
0
0
0
1
0
40
0
40
1
1
1
ticks
30.0

BUTTON
5
110
205
143
NIL
setup-compute-deg-from-N-r\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
5
10
205
43
N
N
0
2500
2500.0
50
1
NIL
HORIZONTAL

MONITOR
5
335
205
380
realized average degree
mean [count link-neighbors] of turtles
2
1
11

PLOT
5
380
205
530
degree distribution
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" "let max-degree max [count link-neighbors] of turtles\nplot-pen-reset  ;; erase what we plotted before\nset-plot-x-range 1 (max-degree + 1)  ;; + 1 to make room for the width of the last bar\nhistogram [count link-neighbors] of turtles"

SWITCH
20
265
147
298
hide-links
hide-links
1
1
-1000

TEXTBOX
40
305
145
323
Hide links for speed
9
0.0
1

INPUTBOX
5
160
205
220
deg
10.0
1
0
Number

BUTTON
5
220
205
251
NIL
setup-compute-r-from-N-deg\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
5
50
205
110
r
0.036247765075687326
1
0
Number

@#$#@#$#@
## WHAT IS IT?

Code to produces a [**random geometric network**](https://en.wikipedia.org/wiki/Random_geometric_graph) in the world. This is a network with **spatial structure**. 

**Definition:** The random geometric network for *N* agents is constructed by placing  nodes randomly in a unit square (with side length one). Each node is connected to all nodes with a radius *r* our its position. x and y coordinates are from a uniform distribution. Euclidean distance is used. 

The outcome is an undirected network with a spatial structure like a regular two-dimensional grid but with less artificial and more heterogeneous degree distribution. Such networks can be a good proxy for all real-world networks for which spatiality is relevant. Examples are the spreading of diseases or information in space through (static) social contacts.


**Problem:** In many  agent-based models, you typically want to analyze your society with a certain average number of links. In a random geometric network, the average number of links is a function of the number of agents *N* and the radius *r*. So, typically you want to make the number *N* as large as possible such that simulations still run at a reasonable speed. So, how to find the right radius *r* such that for a given *N* the average degree is a specific desired number? 

**Solution:** Implemented in this model. You can explore it by fixing *N* with the slider and writing the desired average degree into the variable *deg* and click the adjacent button. The *r* will be computed and put in the variable input box and the network is created. Check the average degree monitor and the link distribution. It also works the other way around. 



How does it work? Here is the "math":

The question was raised here https://math.stackexchange.com/questions/796746/average-degree-of-a-random-geometric-graph

From that we can deduce an equation for the average degree given N and radius r as

deg = (N - 1) * Prob(l < r)

where Prob(l < r) is the probability that a line between two random points in a unit square has a length less than r.

The factor (N-1) instead of N is because the average degree in an undirected graph is 2*L/N with L being the number of links. (In an undirected graph, we have to count each link twice because each link contributes to the degree of two nodes.) The total number of possible undirected links is L = choose(N,2)=N*(N-1)/2 (choose 2 out of N without replacement). All together the degree when all links are present is 2*(N*(N-1)/2)/N = N-1.

What remains is to find an estimate for Prob(l < r), the probability that a link's length in Euclidean space is lower than r. The answer is related to [*Square Line Picking*](https://mathworld.wolfram.com/SquareLinePicking.html). The resource 
provides an equation for the cumulative distribution function for the length of a line l between two random points in the unit square (uniform distribution on both coordinates).

Prob(l < r) = 1/2 * r^4 - 8/3 * r^3 + pi * r^2 

(This is only valid for r <= 1 which is reasonable in most cases, see https://mathworld.wolfram.com/SquareLinePicking.html to extend to r>1.)

So, Prob(l < r) stands for the probability that a line between two random points is shorter than r. 

So we have 

Prob(l < r) = 1/2 * l=r^4 - 8/3 * r^3 + pi * r^2     and 
Prob(l < r) = deg / (N-1)

If we want to compute r given deg and N we have to solve for r in 

(1/2 * r^4 - 8/3 * r^3 + pi * r^2) = deg / (N-1)

As it is a polynomial of degree 4 it is hard to give a closed form solution. So, we  solve using Newton's method to find the root of the function 

f(r) = (1/2 * l=r^4 - 8/3 * r^3 + pi * r^2) - deg / (N-1). 

(Newton's method requires that we also need to know the derivative of f, but this is simple in this case.)


## CREDITS AND REFERENCES

Developed by Jan Lorenz (http://janlo.de) in 2023
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.4.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
1
@#$#@#$#@
