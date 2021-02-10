%% OPPA
%% Created by Yusheng Huang on 02/02 2021 according to 
%  A. Tero, S. Takagi, T. Saigusa, K. Ito, D. P. Bebber, M. D. Fricker,
% K. Yumiki, R. Kobayashi, and T. Nakagaki, ¡°Rules for biologically
%inspired adaptive network design,¡± Science, vol. 327, no. 5964, pp. 439¨C
% 442, 2010.
% T. Nakagaki, M. Iima, T. Ueda, Y. Nishiura, T. Saigusa, A. Tero,
% R. Kobayashi, and K. Showalter, ¡°Minimum-risk path finding by an
% adaptive amoebal network,¡± Physical review letters, vol. 99, no. 6, p.
% 068104, 2007.
function [Q,ite,time,SP]= lyameba_basic2(L,sourceNode,sinkNode)
%% the reproduction of original PPA
% L:adjacent matrix£»sourceNode£ºstarting node£»sinkNode£ºending node
%% code
% n:number of nodes
ite=0;
n=length(L(:,1));
L1=L;
L(L==0)=inf;
D=0.5*ones(n);
D=D-diag(diag(D));
count_sink=length(sinkNode);
A=zeros(n,1);
A(sourceNode)=count_sink;
A(sinkNode)=-1;
temp_Q=zeros(n);
Q=ones(n);
while sum(sum(abs(temp_Q-Q)))>0.0001
     ite=ite+1;
     temp_Q=Q;
    B=D./L;
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
    Q=abs(Q);
    D=(D+Q)/2;
end
[SP]=findspin_PPA(sourceNode,sinkNode,Q,L1);
time=toc;



