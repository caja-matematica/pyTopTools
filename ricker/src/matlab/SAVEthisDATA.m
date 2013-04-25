% This function saves data passed to it from the reAP3.m function.

function SAVEthisDATA(n,d,r,T,thr,eth,s,IC,Betti_storage,abundance)
    
%print betti time series
% fname = sprintf('../../BiomathTalk/rand1/Cubical/Bettis/bettiPd_n%g_d%g_r%g_T%g_thr%g_eth%g_s%g_IC%g.txt',n,d*100,r,T,thr,eth,s,IC);
fname = sprintf('../../BiomathTalk/rand1/Cubical/Stoch/Bettis/bettiPd_n%g_d%g_r%g_T%g_thr%g_eth%g_s%g_IC%g.txt',n,d*100,r,T,thr,eth,s,IC);

fid = fopen(fname,'w'); 
fprintf(fid, Betti_storage);
fclose(fid);

% Print total abundance.
% fname = sprintf('../../BiomathTalk/rand1/Cubical/Totals/total_pop_pers_det_n%g_d%g_r%g_T%g_thr%g_eth%g_s%g_IC%g.txt',n,d*100,r,T,thr,eth,s,IC);
fname = sprintf('../../BiomathTalk/rand1/Cubical/Stoch/Totals/total_pop_pers_det_n%g_d%g_r%g_T%g_thr%g_eth%g_s%g_IC%g.txt',n,d*100,r,T,thr,eth,s,IC);
fid = fopen(fname,'w'); 
fprintf(fid,'%g \n', abundance);
fclose(fid);

end


