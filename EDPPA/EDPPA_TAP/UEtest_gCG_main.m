%% the main UE experiment testing program in graph CG
% created by Yusheng Huang on 02/04 2021
%reference1:A modified Physarum-inspired model for the user equilibrium traffic assignment problem
%reference2:A Bio-Inspired Approach to Traffic Network Equilibrium Assignment Problem
%% experiment settings
repeat_times=5;%each algorithm runs for repeat_times times
stopping_criteria=1e-5;
% Input the dataset
rawnet_network = xlsread('CGnetwork.xlsx','network');  % Input the graph
nodenum=933;
linknum=length(rawnet_network(:,1)); 
list_link= rawnet_network;
alp=110011.*ones(nodenum); %free-flow travel time %...
cap=zeros(nodenum); %capacity 
explink=zeros(nodenum);
for iijj=1:linknum 
    alp(list_link(iijj,1),list_link(iijj,2))=list_link(iijj,3);
    cap(list_link(iijj,1),list_link(iijj,2))=list_link(iijj,4);
end
explink(find(cap~=0))=4;
explink(find(alp==0))=0;
alp(find(alp==0))=0.034506800000000004/1.15;
alp(find(alp==110011))=0;
num_source=12;
list_source=[1,2,3,12,14,16,18,21,24,30,45,56];
list_sinks=cell(num_source,1); 
list_sinks{1,1}=782;%...
list_sinks{2,1}=915;%...
list_sinks{3,1}=933;%...
list_sinks{4,1}=821;%...
list_sinks{5,1}=912;%...
list_sinks{6,1}=912;%...
list_sinks{7,1}=933;%...
list_sinks{8,1}=908;%...
list_sinks{9,1}=745;%...
list_sinks{10,1}=921;%...
list_sinks{11,1}=904;%...
list_sinks{12,1}=879;%...
list_demand=cell(num_source,1);
list_demand{1,1}=600;%...
list_demand{2,1}=300;%...
list_demand{3,1}=4000;%...
list_demand{4,1}=300;%...
list_demand{5,1}=450;%...
list_demand{6,1}=580;%...
list_demand{7,1}=680;%...
list_demand{8,1}=4000;%...
list_demand{9,1}=400;%...
list_demand{10,1}=500;%...
list_demand{11,1}=2000;%...
list_demand{12,1}=600;%...
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
    xlswrite('result_CG.xlsx',currentresult3,'originalUE_RGAP',xlRange);
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
    xlswrite('result_CG.xlsx',currentresult3,'UE_ED_RGAP',xlRange);
end
