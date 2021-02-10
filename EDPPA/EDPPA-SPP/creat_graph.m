%% a small program to play with different algorithms 
%% created by Yusheng Huang on 02/02 2021 
warning('off');
%% randomly generate a graph
% randomly generate a fully connected un directed graph 
% with the link connecting the sink and source node disconnected
%N_node number of nodes
%matrix_graph  the adjacent matrix 
N_node=1000;
Node_start=2;
Node_end=4;
value_max=100;
matrix_graph_tmp=fix(rand(N_node).*value_max)+1;
matrix_graph_tmp(Node_end,Node_start)=0;  
matrix_graph_tmp(Node_start,Node_end)=0;  
matrix_graph_tmp2=tril(matrix_graph_tmp)+(tril(matrix_graph_tmp))';
matrix_graph=matrix_graph_tmp2-2.*diag(diag(matrix_graph_tmp));
text_graph=graph(matrix_graph);
%% You could add more algorithms here 
tic;
[P,d] = shortestpath(text_graph,Node_start,Node_end,'Method','positive');
time_DJ=toc

[Q0,itechange2,time_effect,SP]=lyameba_test_Deffective(matrix_graph,Node_start,Node_end);
time_effect
%% output the length of the path
d
SP


