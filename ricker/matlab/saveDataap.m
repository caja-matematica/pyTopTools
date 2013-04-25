% This function saves data passed to it from the reAP3.m function.

function saveDataap(n,d,r,final_time,scale,IC,abundance)
% TODO The paths need to be updated here! 
    
% 'Totals/Coef22/totalPD_%g_%g_%g_%g_%g.txt'

% Print total abundance.
fname = sprintf('total_pop_pers_det_n%g_d%g_r%g_T%g_s%g_IC%g.txt',n,d*100,r,final_time,scale,IC);
fid = fopen(fname,'w'); 
fprintf(fid,'%g \n', abundance);
fclose(fid);

end
