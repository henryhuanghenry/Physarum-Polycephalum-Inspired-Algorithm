%% PPA-ED
%% Created by Yusheng Huang on 02/02 2021 
function [Q,ite,time,SP_effect]= lyameba_test_Deffective(L,sourceNode,sinkNode)
%% the proposed PPA-ED bayPPA
% L:adjacent matrix£»sourceNode£ºstarting node£»sinkNode£ºending node
%% code
% n:number of nodes
tic;
n=length(L(:,1));
L1=L;
L(L==0)=inf;
D=0.5*ones(n);
D=D-diag(diag(D));
A=zeros(n,1);
A(sourceNode)=1;
A(sinkNode)=[];
temp_Q=zeros(n);
Q=ones(n);
ite=0;
Leffective=L;
while sum(sum(abs(temp_Q-Q)))>0.0001
     ite=ite+1;
     temp_Q=Q;
    B=D./Leffective;
    B=diag(sum(B))-B; 
    B(sinkNode,sinkNode)=0;
    B(sourceNode,:)=B(sourceNode,:)-B(sinkNode,:);
    B(sinkNode,:)=[];
    B(:,sinkNode)=[];
    P1=B\A;
    if sinkNode<n&& sinkNode>1
        P=[P1(1:sinkNode-1);0;P1(sinkNode:end)];
    elseif sinkNode==n
        P=[P1;0];
    else
        P=[0;P1];
    end
    tempP=repmat(P,1,n)-repmat(P',n,1);
    Q=(D./Leffective).*tempP;
    Q=abs(Q);
    Leffective=((1-Q).*L./D+Q.*(1-log(4*Q./(sum(Q)+sum(Q,2))))); 
    Leffective(isnan(Leffective))=inf;
    Leffective(L==0)=inf;
     D=(D+Q)/2;
end
[SP_effect]=findspin_PPA(sourceNode,sinkNode,Q,L1);
time=toc;





