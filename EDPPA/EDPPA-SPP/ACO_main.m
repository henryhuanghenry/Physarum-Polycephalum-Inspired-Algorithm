
function [L,timeout]=ACO_main(inG)
tic;

itemax=200;
routemax=10;
if length(inG(1,:))>100
    p_ant=200;
else
    p_ant=30;
end
p_Q=100;
p_beta=1;
p_decay=0.2;
p_gama=1;
p_q0=0.5;
p_mu=2;

best_L=1000.*ones(itemax,1);
best_ant=zeros(itemax,routemax);

size_inG=length(inG(1,:));
insu=ones(size_inG);
insu(find(inG==0))=0;

for i_ite=1:itemax

    ant_matrix=zeros(p_ant,routemax); 
    tmp_bestL=100000;
    tmp_bestant=zeros(1,routemax);
    for i_ant=1:p_ant
        [ant_matrix(i_ant,:),insu,Lant]=antrouting(inG,insu,p_q0,routemax,p_beta,p_decay,p_Q);  %进行每只蚂蚁选路
        if Lant<tmp_bestL 
            tmp_bestL=Lant;
            tmp_bestant=ant_matrix(i_ant,:);
        end
    end
    best_L(i_ite,1)=tmp_bestL;
    best_ant(i_ite,:)=tmp_bestant;
    L_ave=sum(best_L(1:i_ite,1))/i_ite;
    no_L_golbal=find(best_L(1:i_ite,1)==min(best_L(1:i_ite,1)));
    L_golbal=best_L(no_L_golbal(1,1),1);
    if isnan((tmp_bestL-L_golbal)/(L_ave-L_golbal))
        p_xiga=0.5;
    else
        p_xiga=1/pi*atan(p_gama*(tmp_bestL-L_golbal)/(L_ave-L_golbal))+0.5;
    end
    for iii=2:routemax
        if tmp_bestant(iii)==0
            break;
        end
        insu(tmp_bestant(iii-1),tmp_bestant(iii))=insu(tmp_bestant(iii-1),tmp_bestant(iii))+p_mu*p_xiga/tmp_bestL;
        insu(tmp_bestant(iii),tmp_bestant(iii-1))=insu(tmp_bestant(iii),tmp_bestant(iii-1))+p_mu*p_xiga/tmp_bestL;
    end
    if length(no_L_golbal)>=10
        break;
    end
end
nobest=find(best_L==min(best_L));
L=best_L(nobest(1,1));
timeout=toc;
end