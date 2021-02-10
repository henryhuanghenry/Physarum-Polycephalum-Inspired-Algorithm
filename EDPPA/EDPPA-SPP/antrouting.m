function [antroute,insuout,Lant]=antrouting(inG,insu,p_q0,routemax,p_beta,p_decay,p_Q)

endnode=length(inG);
antroute=zeros(1,routemax);
antroute(1)=1; 
tmp_inG=inG;
tmp_inG(find(inG==0))=inf;

tmp_insu=insu;
for ijk=2:routemax
    if ijk==routemax  
        antroute(ijk)=endnode;
        break;
    end

    qq=rand(); 
    if qq>p_q0 
        sum_insu=sum(tmp_insu(antroute(ijk-1),:).*(1./tmp_inG(antroute(ijk-1),:)).^p_beta);
        P_choose=tmp_insu(antroute(ijk-1),:).*(1./tmp_inG(antroute(ijk-1),:)).^p_beta./sum_insu;  
        P_choose=cumsum(P_choose);
        PP=rand();
        no_choose=find(P_choose>=PP,1);
        antroute(ijk)=no_choose(1); 
        tmp_insu(:,no_choose(1))=0; 
    else 
        tmp_insu_matrix=tmp_insu(antroute(ijk-1),:).*(1./tmp_inG(antroute(ijk-1),:)).^p_beta;  
        no_choose=find(tmp_insu_matrix==max(tmp_insu_matrix));  
        antroute(ijk)=no_choose(1,1);  
        tmp_insu(:,no_choose(1,1))=0;
    end
    if antroute(ijk)==endnode  
        break;
    end
end
insuout=insu;

Lant=0;
for iii=2:routemax
    if antroute(iii)==0
        break;
    end
    Lant=Lant+inG(antroute(iii-1),antroute(iii));
end

for iii=2:routemax
    if antroute(iii)==0
        break;
    end
    insuout(antroute(iii-1),antroute(iii))=(1-p_decay)*insu(antroute(iii-1),antroute(iii))+p_decay*p_Q/Lant;
    insuout(antroute(iii),antroute(iii-1))=(1-p_decay)*insu(antroute(iii),antroute(iii-1))+p_decay*p_Q/Lant;
end
end