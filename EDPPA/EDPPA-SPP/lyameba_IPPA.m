%% IPPA
%% Created by Yusheng Huang on 02/02 2021 according to 
% X. Zhang, Q. Wang, A. Adamatzky, F. T. Chan, S. Mahadevan, and
% Y. Deng, ¡°An improved physarum polycephalum algorithm for the
% shortest path problem,¡± The Scientific World Journal, vol. 2014, 2014.
function [Q,ite,time,SP]= lyameba_IPPA(L,sourceNode,sinkNode,itebasic)
%% The reproduction of IPPA
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
ite=0;
temp_Q=ones(n);
Q=zeros(n);
%% 
while sum(sum(abs(temp_Q-Q)))>0.0001
    temp_Q=Q;
    B=D./L;
    B=diag(sum(B))-B;
    B(sinkNode,sinkNode)=0;
    B(sourceNode,:)=B(sourceNode,:)-B(sinkNode,:);
    B(sinkNode,:)=[];
    B(:,sinkNode)=[];

    P=B\A;
    if sinkNode<n&& sinkNode>1
        P=[P(1:sinkNode-1);0;P(sinkNode:end)];
    elseif sinkNode==n
        P=[P;0];
    else
        P=[0;P];
    end
    tempP=repmat(P,1,n)-repmat(P',n,1);
    Q=(D./L).*tempP;
    Q=abs(Q);
    D=1/2.*((Q.*tempP./(L.*abs(P(sinkNode)-P(sourceNode)))+D));
    D(L==inf)=0;
    D(isnan(D))=0;
    D(D==inf)=0;
    % sometimes, the IPPA could not find the shortest path due to the
    % design of the algorithm, we stop the algorithm when we find it not
    % able to detect the SP or it takes too long
    if (~isempty(find(Q)==Inf) || ~isempty(inan(Q)) || ite>2*itebasic)
        ite=ite+1;
        break;
    end
end
[SP]=findspin_PPA(sourceNode,sinkNode,Q,L1);
time=toc;
