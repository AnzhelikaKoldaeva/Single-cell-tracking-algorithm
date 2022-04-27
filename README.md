# Single-cell-tracking-algorithm
The algorithm for tracking every single bacteria in a dense growing population

INPUT: movies with recorded growing populations of bacteria E.coli in microchannels. 
OUTPUT: spatial trajectories of every single bacteria with captured mother-daughter reationships 

The algorithm is based on the idea of constructing local structures using Delaunay triangulations. Since it is a dense population (with no or very small distances between the individuals), relative positions of each cell with respect to other cells are more informative than their absolute positions. Therefore, the construction of local structures is a more beneficial approach here. 

The algorithm encompasses a forward and a backward step:

• In the forward step, we associated a cell on a frame at time t with the one on the next frame at time t + ∆t by minimizing
the difference between their local structures.

• In the backward step, all bacteria from a frame at time t + ∆t that were not associated with any bacteria on t frame
were compared to them in the same way in order to determine their mothers.

See panel (a) in the figure below for the reference. As the result, we construct spatial trajectories for the cells from the initial populations along with all the offsprings (see panel (a) in the figure below). This king of trajectories are also called spatial geneology trees.

More details of the algorithm and its application can be found in the following publication 

\href{https://www.pnas.org/doi/abs/10.1073/pnas.2120821119 }{Koldaeva et al. “Population genetics in microchannels.” Proceedings of the National Academy of Sciences, 2022 }


![track_alg2](https://user-images.githubusercontent.com/44256309/165434754-14ee2353-7cd1-40a7-9930-364a9efdb038.jpg)
