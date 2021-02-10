%% The core of the OPPA-TAP 
%% Created by Yusheng Huang in 2021/02 refering to  
%S. Xu, W. Jiang, X. Deng, Y. Shou, A modified physarum-inspired model
%for the user equilibrium traffic assignment problem, Applied Mathematical
%Modelling 55 (2018) 340¨C353.
function [Q,D]= lyameba_basic_trafficUE_mutisink(L,D,sourceNode,sinkNode,demand)
n=length(L(:,1));
L(L==0)=inf;  
L(isnan(L))=inf;
D(find(L==inf))=0; 
D(find(D==inf))=0;  
D(isnan(D))=0;
D=D-diag(diag(D));
A=zeros(n,1);
A(sourceNode)=sum(demand); 
A(sinkNode)=-demand; 


B=(D./L+(D')./(L')); 
B=diag(sum(B))-B; 
B(:,sourceNode)=[];

P=B\A;

if sourceNode<n&&  sourceNode>1
    P=[P(1: sourceNode-1);0;P( sourceNode:end)];
elseif  sourceNode==n
    P=[P;0];
else
    P=[0;P];
end
tempP=repmat(P,1,n)-repmat(P',n,1);
Q=(D./L).*tempP;
Q(find(Q<0))=0; 
D=(Q+D)./2;
end
