% counts the number of self-loops in the graph
% INPUT: adjacency matrix
% OUTPUT: interger, number of self-loops
% cited from http://strategic.mit.edu/downloads.php?page=matlab_networks (Matlab Tools for Network Analysis ) 

function sl=selfloops(adj)

sl=sum(diag(adj));