
function SAVEData(n,dispersal,coef,time,thresh,ethresh,IC,Betti_storage,abundance)

%print betti time series
fname = sprintf('bettiPd_%g_%g_%g_%g_%g_%g_%g.txt',n,dispersal*100,coef,time,thresh,ethresh,IC);
fid = fopen(fname,'w'); 
fprintf(fid, Betti_storage);
fclose(fid);

%print total abundance
fname = sprintf('totalPd_%g_%g_%g_%g_%g_%g_%g.txt',n,dispersal*100,coef,time,thresh,ethresh,IC);
fid = fopen(fname,'w'); 
fprintf(fid,'%g \n', abundance);
fclose(fid);

end