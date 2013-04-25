% Code for the Ricker map code. MATLAB function used to
% Compute Ricker simulations
%
% $$$          n :  The size of the matrix is nxn.
% $$$          d :  Dispersal parameter (percentage of a patch pop to disperse.). 
% $$$          r :  Environemntal fitness parameter.
% $$$          T :  Length of time to run the simulation.
% $$$          s :  For integer conversion in perseus computation. spatial scale.
% $$$         IC :  Number that designates which initial condition to use in simulation.      
% $$$          N :  Population matrix.
% $$$         
% $$$         Example output files:
% $$$         
% $$$         Matrix in Perseus format: 
% $$$             tmpPers_10_10_5_3_1000_1_2_pers.txt
% $$$         
% $$$         Perseus output for 0th homology:
% $$$             tmpPers_10_10_5_3_1000_1_2_pers_0.txt
% $$$         
% $$$         For both of these, the integers, in order, denote the followg, ...
% $$$         with the value above in parentheses:
% $$$        
% $$$         matrix dimension (10)
% $$$         dispersal*100 (10)
% $$$         fitness (r value) (5)
% $$$         final time (T value) (3)
% $$$         scale factor for perseus (1000)
% $$$         initial condition (1)
% $$$         time step (2)
% $$$         pers_* is the homology level (0 above)
%

function reAP3_perseus_det( n, d, r, T, s, IC, N, tmpPerseus, tmpPopMat, tmpPopMovie )

abundance = zeros(1,T); % Initialize abundance array.

D = diag(diag(ones(n),1),1) + diag(diag(ones(n),-1),-1); % Dispersal matrix.

tp = strcat( tmpPerseus, '_0' );
process_perseus_sub_super( N, tp, s ); % Perseus on IC matrix.

tmpPopMat_0 = strcat(tmpPopMat,'_0');  
dlmwrite(tmpPopMat_0, N) % Save IC patch abundances.

% fig = figure;
% 
% name = sprintf('Mov_n%g_d%g_r%g_T%g_s%g_IC%g.avi',n,d*100,r,T,s,IC);
% 
% aviobj = avifile(name,'compression','None');

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for t = 1:T     % Simulation code.
    ti = num2str( t );
    
    abundance(t) = sum(sum(N)); % Sums all elts of N (total pop on grid).
    
    if abundance(t) == 0    % Check for Extinction
        % (If global extinction => stop simulation before reaching T).
        saveDataap(n, d, r, T, s, IC, abundance);
        return;
    end
    
    N = r*N.*exp(-N/s);    % Ricker growth map (each patch, each iteration).
    N = (1 - d)*N + (d/4)*(N * D) + (d/4)*(D * N);      % Dispersal phase.

    
    tp = strcat( tmpPerseus, '_', ti );
    process_perseus_sub_super( N, tp, s );
    
    tmpPopMat_t = strcat(tmpPopMat,'_',ti);
    dlmwrite(tmpPopMat_t,N)
    
    
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Population Matrix Movie Simulation.
    

%     MatName = strcat(tmpPopMovie,'_', ti);
%     PopMat = load(MatName);
% 
%     minPop = min(min(PopMat));
%     maxPop = max(max(PopMat));
% 
%     im = imagesc(PopMat);
% 
%     colormap summer
%     colormap(flipud(colormap))
%     colorbar
% 
%     caxis([0 2])
% 
%     rows = size(PopMat,1);
%     cols = size(PopMat,2);
% 
%     for i = 1:cols
%         matdata = ones(6*rows, 4);
%         k = 1;
% 
%         for j = 0:6:(6*rows - 6)
%             matdata(j + 1:j + 6, :) = PopMat(k,i);
%             k = k + 1;
%         end
% %         set(im(i),'Cdata',matdata)
%     end
% 
%     if(mod(t,2) == 0) % To account for the period 2 flipping (color flips).
% 
%         aviobj = addframe(aviobj, fig);
%     end
    

end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% aviobj = close(aviobj);

saveDataap(n,d,r,T,s,IC,abundance);

end