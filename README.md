# Who Should Get Vaccinated? Individualized Allocation of Vaccines Over SIR Network

Replications file for "Who Should Get Vaccinated? Individualized Allocation of Vaccines Over SIR Network" by Toru Kitagawa and Guanyi Wang

The files in this archive replicate the results reported in the "Simulation Exercises" section of the paper.

### Softwares

- MATLAB R2020b

### Analysis replication files

Network_Generation: Generate the random network

- Random_Graph_Generation.m: We totally generate nine random networks. We choose three number of nodes: N = 200, 500, 800 and three number of edges for each number of node: edges = 7N, 10N, Full.
- Small_Network_Generation.m: We totally generate nine random networks. We choose three number of nodes: N = 25, 30, 35 and three number of edges for each number of node: edges = 7N, 10N, Full.
- Multiple_Graph_Generation.m: For each pair number of nodes and number of edges, we generate 20 different random networks. There are totally 9 pairs of number of nodes and number of edges in large network setting. 
- Multiple_small_generation.m : For each pair number of nodes and number of edges, we generate 20 different random networks. There are totally 9 pairs of number of nodes and number of edges in small network setting. 

Greedy: Greedy Algorithm

- Greedy.m: The greedy algorithm for a single large network.
- Greedy_small: The greedy algorithm for a single small network.
- Greedy_multi: The greedy algorithm for multiple (in here we choose 20) large networks.
- Greedy_samll_multi: The greedy algorithm for multiple (in here we choose 20) small networks.
- obj.m : Objective Function

Global_opt: Brute force serach

- Global_search.m: The global optimal solution for a single small network.
- Global_search_multi.m: The global optimal solution for multiple (in here we choose 20) small networks.
- obj.m: Objective Function

Random_Search: Random allocation rule

- Random_search.m: The random allocation rules for a single large network.
- Random_search_multi.m:  The random allocation rules for multiple (in here we choose 20) large networks.
- obj.m: Objective Function

Allocation_nonet: Targeting without network information

- Allocation_without_network.m: The allocation rule without network information for a single large network.
- Allocation_without_network_multi.m: The allocation rule without network information for multiple (in here we choose 20) large networks.
- obj.m: Objective Function

Graph: Draw the figure

- figure_compare_greedy_nonet.m: Draw the figure for the comparison between greedy algorithm and allocation without network information.
- figure_compare_greedy_random.m: Draw the figure for the comparision between greedy algorithm and random allocation.
- figure_compare_random_nonet.m: Draw the figure for the comparision between random allocation and allocation without network information.
