% Computes the average degree of a node in a graph, defined as 2*num_edges 
% divided by the num_nodes (every edge is counted in degrees twice).
% Other routines used: numnodes.m, numedges.m
% cited from http://strategic.mit.edu/downloads.php?page=matlab_networks (Matlab Tools for Network Analysis ) 

function k=average_degree(adj)

k=2*numedges(adj)/numnodes(adj);