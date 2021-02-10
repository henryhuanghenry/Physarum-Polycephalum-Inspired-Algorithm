%% Created by Yusheng Huang on 02/02 2021 
% used to find the shortest path
function [sum_SP]=findspin_PPA(Node_start,Node_end,Q0,matrix_graph)
graph_flux=Q0;
nodelist=Node_start;
endnode=1;
sum_SP=0;
while nodelist(endnode)~=Node_end
    graph_flux(nodelist(endnode),:)=0;
    endnode=endnode+1;
    node_tmp=find(graph_flux(:,nodelist(endnode-1))==max(graph_flux(:,nodelist(endnode-1))));
    nodelist=[nodelist;node_tmp];
    if size(node_tmp)==0
        sum_SP=NaN;
        break;
    end
    sum_SP=sum_SP+matrix_graph(nodelist(endnode-1),nodelist(endnode));
    if isnan(Q0(nodelist(endnode-1),nodelist(endnode))) || isnan(Q0(nodelist(endnode-1),Node_end))
        sum_SP=NaN;
        break;
    end
end
end