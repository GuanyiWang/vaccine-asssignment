# Optimal Vaccine Allocation

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

- Greedy.m:
- Greedy_small:
- Greedy_multi:
- Greedy_samll_multi:
- obj.m : Objective Function

Global_opt: Brute force serach

- Global_search.m:
- Global_search_multi.m:
- obj.m: Objective Function

Random_Search: Random allocation rule

- Random_search.m:
- Random_search_multi.m:
- obj.m: Objective Function

Allocation_nonet: Targeting without network information

