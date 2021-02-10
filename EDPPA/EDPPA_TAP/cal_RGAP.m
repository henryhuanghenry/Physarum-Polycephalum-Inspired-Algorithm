%% calculating the RGAP
% Created by Yusheng Huang in 2021/02
function out=cal_RGAP(Q_all,L_temper,num_source,list_source,list_sinks,list_demand,rawnet_network)
source=rawnet_network(:,1)'; 
sink=rawnet_network(:,2)';  
demand=zeros(1,length(rawnet_network(:,1))); 

for ijk=1:length(rawnet_network(:,1))   
    demand(ijk)=L_temper(rawnet_network(ijk,1),rawnet_network(ijk,2));
end
text_graph=digraph(source,sink,demand); 
%% º∆À„RGAP
sum_cost=0;
for i_source=1:num_source  
    for i_sink=1:length(list_sinks{i_source,1})
        Node_start=list_source(i_source);
        Node_end=list_sinks{i_source,1}(i_sink);
        [P,L_shortestp] = shortestpath(text_graph,Node_start,Node_end,'Method','positive');
        sum_cost=sum_cost+L_shortestp*list_demand{i_source,1}(i_sink);
    end
end
out=1-sum_cost/sum(sum(Q_all.*L_temper)); 
end
