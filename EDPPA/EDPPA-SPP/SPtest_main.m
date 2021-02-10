%% SPtest_main
%% Created by Yusheng Huang on 02/02 2021
% with the following files:
% lyameba_test_Deffective.m the proposed algorithm
% lyameba_IPPA.m the baseline algorithm IPPA
% lyameba_basic2.m the baseline algorithm original PPA
% FastPhysarumSolver_forCompare.m the baseline algorithm APS (created by other authors)
% EHPA_forCompare.m fun_findRoute.m the baseline algorithm EHPA (created by other authors)
% bayPPA.m the baseline algorithm BPPA
% ACO_main.m antrouting.m the baseline algorithm ACO
% findspin_PPA.m used to find the shortest path when the flow matrix
% provided from the PPA-algorithms is given
%% The new test code for PPA_ED-Experiments of shortest path finding 
% We compare the running time and the degree of accuracy between different
% algorithms.
warning('off');
%% settings of experiments
% we randomly generate the tested graph
% (fully connected undirected graph with the link connecting the source node and the sink disconnected )
NNnode=[10 100 500 1000 1500 2000 3000];% the graph sizes of the tested graph
nodenum=length(NNnode);% the number of different graph size
numm=10;% for each graph size, "numm" graph are randomly generated to simulate different topologies
NNN=8; %the number of tested algorithms
times_inner=3;% each algorithm run for "times_inner" times in the same graph
%% the variables for recording
datasets=cell(nodenum,numm);%record the datasets
result_time=zeros(NNN*nodenum,numm+1); %record
result_correct=zeros(NNN*nodenum,numm+1);  % record the degree of accuracy
result_ite=zeros(NNN*nodenum,numm+1);
result_time_innermiddle=zeros(numm*NNN,times_inner+1); % recordthe running time with the same graph size
result_ite_innermiddle=zeros(numm*NNN,times_inner+1); % recordthe iterations with the same graph size
result_correct_innermiddle=zeros(numm*NNN,times_inner+1); % recordthe degree of accuracy  with the same graph size
%% The outer loop--for different graph size
ijk=0; %the counter for the main iteration
for N_node=NNnode(1) % the number of nodes
    ijk=ijk+1 %loop counter
    %% set some parameters of the tested graph
    Node_start=1;% the starting node
    Node_end=N_node;% the ending node
    value_max=100;% the maximum value of the length of the graph
    counttt=0;   %the counter for the loop
    %% The middle loop--for different graphs with the same graph size
    while counttt<numm
        counttt=counttt+1
        %% generate a graph
        matrix_graph_tmp=fix(rand(N_node).*value_max)+1;
        matrix_graph_tmp(Node_end,Node_start)=0;
        matrix_graph_tmp(Node_start,Node_end)=0;
        matrix_graph_tmp2=tril(matrix_graph_tmp)+(tril(matrix_graph_tmp))';
        matrix_graph=matrix_graph_tmp2-2.*diag(diag(matrix_graph_tmp));
        text_graph=graph(matrix_graph);% for the function 'shortestpath'
        inG=matrix_graph;%the adjacent matrix of the graph
        datasets{ijk,counttt}=inG;%record the datasets
        %% The inner loop--each algorithm run for several times in the same graph
        counterer=0; %the counter for the loop
        while counterer<times_inner
            counterer=counterer+1;
            %% execute different algorithms
            tic;
            [P,d] = shortestpath(text_graph,Node_start,Node_end,'Method','positive');
            timeout_DJ=toc;
            [Q_basic,ite_basic,timeout_basic,SP_basic]=lyameba_basic2(inG,Node_start,Node_end);
            [Q_effect,ite_effect,timeout_effect,SP_effect]=lyameba_test_Deffective(inG,Node_start,Node_end);
            [Q_bay,ite_bay,timeout_bay,SP_bay]=bayPPA(inG,Node_start,Node_end);
            [Q_IPPA,ite_IPPA,timeout_IPPA,SP_IPPA]=lyameba_IPPA(inG,Node_start,Node_end,ite_basic);
            [ite_EHPA,time_EHPA, SP_EHPA] = EHPA_forCompare( inG, Node_start, Node_end );
            [ite_APS,time_APS,SP_APS] = FastPhysarumSolver_forCompare(inG,Node_end, Node_start, Node_end);
            [L_ACO,timeout_ACO]=ACO_main(inG);
            %% the degree of accuracy
            if SP_basic==d
                c_basic=1;
            else
                c_basic=0;
            end
            if SP_effect==d
                c_effect=1;
            else
                c_effect=0;
            end
            if SP_bay==d
                c_bay=1;
            else
                c_bay=0;
            end
            if SP_IPPA==d
                c_IPPA=1;
            else
                c_IPPA=0;
            end
            if SP_EHPA==d
                c_EHAP=1;
            else
                c_EHAP=0;
            end
            if SP_APS==d
                c_APS=1;
            else
                c_APS=0;
            end
            if L_ACO==d
                c_ACO=1;
            else
                c_ACO=0;
            end
            result_time_innermiddle((counttt-1)*NNN+1:counttt*NNN,counterer)=[timeout_DJ;timeout_basic;timeout_effect;timeout_bay;timeout_IPPA;time_EHPA;time_APS;timeout_ACO];
            result_correct_innermiddle((counttt-1)*NNN+1:counttt*NNN,counterer)=[1;c_basic;c_effect;c_bay;c_IPPA;c_EHAP;c_APS;c_ACO];
            result_ite_innermiddle((counttt-1)*NNN+1:counttt*NNN,counterer)=[0;ite_basic;ite_effect;ite_bay;ite_IPPA;ite_EHPA;ite_APS;0];
        end
        result_time_innermiddle(counttt*NNN,times_inner+1)=counttt;
        result_correct_innermiddle(counttt*NNN,times_inner+1)=counttt;
        result_ite_innermiddle(counttt*NNN,times_inner+1)=counttt;
        result_time((ijk-1)*NNN+1:ijk*NNN,counttt)=mean(result_time_innermiddle((counttt-1)*NNN+1:counttt*NNN,1:times_inner),2);% the average value of each running time
        result_correct((ijk-1)*NNN+1:ijk*NNN,counttt)=mean(result_correct_innermiddle((counttt-1)*NNN+1:counttt*NNN,1:times_inner),2);% the average value of each running time
        result_ite((ijk-1)*NNN+1:ijk*NNN,counttt)=mean(result_ite_innermiddle((counttt-1)*NNN+1:counttt*NNN,1:times_inner),2);% the average value of each running time
    end
    result_time(ijk*NNN,numm+1)=N_node;
    result_correct(ijk*NNN,numm+1)=N_node;
    result_ite(ijk*NNN,numm+1)=N_node;
    %% recording
    % the detailed results for each running times at each graph with different
    % graph size
    xlswrite('SP_detailed_result.xlsx',result_time_innermiddle,['time',num2str(NNnode(ijk))]);
    xlswrite('SP_detailed_result.xlsx',result_correct_innermiddle,['degree of accuracy',num2str(NNnode(ijk))]);
    xlswrite('SP_detailed_result.xlsx',result_ite_innermiddle,['iterations',num2str(NNnode(ijk))]);
    % the average results for different graphs with different graph size
    xlswrite('SP_average_result.xlsx',result_time,'time');
    xlswrite('SP_average_result.xlsx',result_correct,'degree of accuracy');
    xlswrite('SP_average_result.xlsx',result_ite,'iterations');
end
%% recording
% save the dataset
save datasets
% the detailed results for each running times at each graph with different
% graph size
xlswrite('SP_detailed_result.xlsx',result_time_innermiddle,['time',num2str(NNnode(ijk))]);
xlswrite('SP_detailed_result.xlsx',result_correct_innermiddle,['degree of accuracy',num2str(NNnode(ijk))]);
xlswrite('SP_detailed_result.xlsx',result_ite_innermiddle,['iterations',num2str(NNnode(ijk))]);
% the average results for different graphs with different graph size
xlswrite('SP_average_result.xlsx',result_time,'time');
xlswrite('SP_average_result.xlsx',result_correct,'degree of accuracy');
xlswrite('SP_average_result.xlsx',result_ite,'iterations');



