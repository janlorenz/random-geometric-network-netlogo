# random-geometric-network-netlogo
Provides code to produce a random geometric network in the NetLogo world. In particular, it solves the problem of getting the right radius to realize a certain average degree for a given number of agents.

Programmed in NetLogo 6.4

![image](https://github.com/janlorenz/random-geometric-network-netlogo/assets/7503499/b8aa9ca7-f777-47e1-bb17-b29f6a22dd33)

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
