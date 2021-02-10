%% the main UE experiment testing program in graph SF
% created by Yusheng Huang on 02/04 2021
%reference1:A modified Physarum-inspired model for the user equilibrium traffic assignment problem
%reference2:A Bio-Inspired Approach to Traffic Network Equilibrium Assignment Problem
%% experiment settings
repeat_times=5;%each algorithm runs for repeat_times times
stopping_criteria=1e-7;
% Input the dataset
rawnet_trips = xlsread('SFnetwork.xlsx','demand'); % Input the graph
rawnet_network = xlsread('SFnetwork.xlsx','network');  % Input the graph
nodenum=24; %the node number of the graph
linknum=length(rawnet_network(:,1)); %%the link number of the graph
list_link= rawnet_network; % the graph  format:starting node of link/ending node of link/link's alp/link's cap
alp=zeros(nodenum); %free-flow travel time
cap=zeros(nodenum); %capacity, disconnected link has the value 0
for iijj=1:linknum %construct the alp and cap matrix
    alp(list_link(iijj,1),list_link(iijj,2))=list_link(iijj,3);
    cap(list_link(iijj,1),list_link(iijj,2))=list_link(iijj,4);
end
list_source=1:1:nodenum;
num_source=length(list_source);
list_sinks=cell(num_source,1);
list_demand=cell(num_source,1);
for i_network=1:nodenum
    list_sinks{i_network,1}=[];
    list_demand{i_network,1}=[];
    for i_24=1:24
        if mod(i_24,5)==0
            tmp_count=5;
        else
            tmp_count=mod(i_24,5);
        end
        if rawnet_trips(1+(i_network-1)*7+ceil(i_24/5),2*tmp_count)~=0
            list_sinks{i_network,1}=[list_sinks{i_network,1},rawnet_trips(1+(i_network-1)*7+ceil(i_24/5),(tmp_count-1)*2+1)];
            list_demand{i_network,1}=[list_demand{i_network,1},rawnet_trips(1+(i_network-1)*7+ceil(i_24/5),2*tmp_count)];
        end
    end
end
%% the testing program
for i_times=1:repeat_times
    %% the original UE
    RGAP=1;
    currentresult3=zeros(10000,2);
    ite=0;
    Q_all=0.5.*ones(nodenum);
    L=alp.*(1+0.15.*(Q_all./cap).^4);
    Q_loop=cell(num_source,1);
    D_loop=cell(num_source,1);
    for ii=1:num_source
        D_loop{ii,1}=0.5.*ones(nodenum);
    end
    %% the outter loop
    while RGAP>=stopping_criteria
        ite=ite+1;
        tic;
        Q_all_tmp=zeros(nodenum);
        %% the inner loop
        for ij=1:length(list_source)
            [Q_loop{ij,1},D_loop{ij,1}]= lyameba_basic_trafficUE_mutisink(L,D_loop{ij,1},list_source(ij),list_sinks{ij,1},list_demand{ij,1});
            Q_all_tmp=Q_all_tmp+Q_loop{ij,1};
        end
        Q_all=Q_all_tmp;
        L=(L+alp.*(1+0.15.*(Q_all./cap).^4))./2;
        ite_time=toc;
        %% updateing variables
        L_temper=L;
        L_temper(isnan(L_temper))=0;
        RGAP=abs(cal_RGAP(Q_all,L_temper,num_source,list_source,list_sinks,list_demand,rawnet_network));%calculate the RGAP
        %% store the time and the total cost of each iteration
        currentresult3(ite,:)=[ite_time,RGAP];
        disp(['i_times=',num2str(i_times),' originalUE RGAP=',num2str(RGAP)]);
    end
    % store the time and the total cost of the whole algorithm
    currentresult3(:,1)=cumsum(currentresult3(:,1));
    xlRange=sprintf("%c1",65+(i_times-1)*2);
    xlswrite('result_SF.xlsx',currentresult3,'originalUE_RGAP',xlRange);
    %% the UE-ED(ours)
    currentresult3=zeros(10000,2);
    RGAP=1;
    ite=0;
    Q_all=0.5.*ones(nodenum);
    L=alp.*(1+0.15.*(Q_all./cap).^4);
    Q_loop=cell(num_source,1);
    D_loop=cell(num_source,1);
    for ii=1:num_source
        D_loop{ii,1}=0.5.*ones(nodenum);
    end
    Leffective=L;
    %% the outter loop
    while RGAP>=stopping_criteria
        ite=ite+1;
        tic;
        Q_all_tmp=zeros(nodenum);
        %% the inner loop
        for ij=1:length(list_source)
            [Q_loop{ij,1},D_loop{ij,1}]= lyameba_basic_trafficUE_mutisink(Leffective,D_loop{ij,1},list_source(ij),list_sinks{ij,1},list_demand{ij,1});
            Q_all_tmp=Q_all_tmp+Q_loop{ij,1};
        end
        %% updateing variables
        Q_all=Q_all_tmp;
        Leffective=(1-RGAP).*(L+alp.*(1+0.15.*(Q_all./cap).^4))+RGAP.*((1-log(Q_all./(sum(Q_all,2))))+alp.*(1+0.15.*(Q_all./cap).^4));  %й╫вс2
        Leffective(L==0)=inf;
        L=(L+alp.*(1+0.15.*(Q_all./cap).^4))./2;
        ite_time=toc;
        %% calculate the RGAP
        L_temper=L;
        L_temper(isnan(L_temper))=0;
        RGAP=abs(cal_RGAP(Q_all,L_temper,num_source,list_source,list_sinks,list_demand,rawnet_network));%calculate the RGAP
        %% store the time and the total cost of each iteration
        currentresult3(ite,:)=[ite_time,RGAP];
        disp(['i_times=',num2str(i_times),' UE-ED RGAP=',num2str(RGAP)]);
    end
    % store the time and the total cost of the whole algorithm
    currentresult3(:,1)=cumsum(currentresult3(:,1));
    xlRange=sprintf("%c1",65+(i_times-1)*2);
    xlswrite('result_SF.xlsx',currentresult3,'UE_ED_RGAP',xlRange);
end
