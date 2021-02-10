%% BPPA
%% Created by Yusheng Huang on 02/02 2021 according to 
% Q. Cai and Y. Deng, ¡°A fast bayesian iterative rule in amoeba algorithm.¡±
% International Journal of Unconventional Computing, vol. 14, 2019.
function [Q,ite,time,SP]= bayPPA(f,sourceNode,sinkNode)
%% the reproduction of bayPPA
% L:adjacent matrix£»sourceNode£ºstarting node£»sinkNode£ºending node
%% code
% n:number of nodes
tic;
L=f;
L1=f;
n=length(L(:,1));
L(L==0)=inf;
D=zeros(n,n);
for i=1:n
    for j=1:n
        if L(i,j)~=0
            D(i,j)=1;
        end
    end
end
YU=min(min(abs(L)));
for i=1:n
    for j=1:n
        if L(i,j)<0
            L(i,j)=YU/(1+abs(L(i,j)));
        end
    end
end
D=ones(n);
D=D-diag(diag(D));
A=zeros(n,1);
A(sourceNode)=1;
A(sinkNode)=[];
temp_Q=zeros(n);
Q=ones(n);
ite=0;
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
    M=repmat(sum(Q),n,1).*repmat(sum(Q),n,1)';
    Pryx=2.*Q.*Q.*Q./M;
    D=Pryx;
    D(sourceNode,sinkNode)=0.5.* D(sourceNode,sinkNode);
    D(sinkNode,sourceNode)=0.5.* D(sinkNode,sourceNode);
    ite=ite+1;
end
[SP]=findspin_PPA(sourceNode,sinkNode,Q,L1);
time=toc;
end

