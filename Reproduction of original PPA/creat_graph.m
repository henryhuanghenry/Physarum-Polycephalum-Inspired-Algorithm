%% Main file
% Execute with function findspin_PPA.m and lyameba_basic.m
%-------------------------------------Introduction-------------------------------------------------
%Author: Yusheng Huang (Created Date: 10.24 2020)
%Use this file to compare the original PPA with the MATLAB function in the shortest path problem.
%--------------------------------------------------------------------------------------------------
%--------------------Warning--------------------------------
% Every time you run this file you would create a new graph.
% Please modify this file if data storage is needed.
%-----------------------------------------------------------
% --------------Some important variables----------------------
%N_node:Node number 
%matrix_graph: the adjacent matrix of the test graph 
%text_graph: MATLAB graph for the MATLAB shortestpath function 
%-------------------------------------------------------------
warning('off');
%% Generate a Random undirected graph:
N_node=50;%N_node:Node number 
Node_start=2; %starting node 
Node_end=4;% ending node 
value_max=100; %edge's weight ranging from 1 to value_max
matrix_graph_tmp=fix(rand(N_node).*value_max)+1;%create a random matrix
matrix_graph_tmp(Node_end,Node_start)=0;  % cancel the edge connecting the starting node and the ending node directly
matrix_graph_tmp(Node_start,Node_end)=0;  % cancel the edge connecting the starting node and the ending node directly
% If we do not cancel this edge, this edge is highly likely to be the --
% shortest path, which makes the testing less significant.
matrix_graph_tmp2=tril(matrix_graph_tmp)+(tril(matrix_graph_tmp))';%Create a symmetric matrix with diagonal elements all 0
matrix_graph=matrix_graph_tmp2-2.*diag(diag(matrix_graph_tmp));%Create a symmetric matrix with diagonal elements all 0
% so the "matrix_graph" is the adjacent matrix of the test graph 
text_graph=graph(matrix_graph); % MATLAB graph for the MATLAB shortestpath function 

%% Call the MATLAB shortestpath function 
tic; %for executing time calculation 
[P,d] = shortestpath(text_graph,Node_start,Node_end,'Method','positive');
time=toc %%for executing time calculation 
%P is the link set of the shortest path
%D is the length of the shortest path 

%% Call the original PPA
[Qbasic,Dbasic,iteBasic,time_basic,d_PPA,nodelist]=lyameba_basic(matrix_graph,Node_start,Node_end);
iteBasic %output the number of iterations 
time_basic  %output the executing time
%
d
d_PPA



