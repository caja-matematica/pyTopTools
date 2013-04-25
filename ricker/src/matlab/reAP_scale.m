% Code for the Ricker map code. MATLAB function used to
% Compute Ricker simulations
%{
         n :  the size of the matrix is nxn.
         d :  an input dispersal parameter that is percentage of a patch's 
              population that is dispersed. For our paper 0 \leq d \leq 0.40. 
         T :  length of time to run the simulation. 
       thr :  value at-which to threshold entries of the population matrix 
              (N in code and P in the paper). Fixed at 10^(-4) for our paper. 
              Because  ethres>0 this is effectively thresholding at zero.
       eth :  local extinction threshold. 
        IC :  number that signifies which initial condition to use for the 
              simulation.
%}

function reAP_scale(n,d,r,T,thr,eth,s,IC)
dd = num2str(d*100);
ss = num2str(s);
% Create a temporary file that is passed to CHomP:
tmpChompFile = sprintf('../../BiomathTalk/rand1/Cubical/tmpChomp_n%g_d%g_r%g_T%g_thr%g_eth%g_s%g_IC%g.txt',n,d*100,r,T,thr,eth,s,IC);
tmp_string = sprintf(tmpChompFile); % Format data into string.
% tmpPopMat = strcat('../../BiomathTalk/rand1/Cubical/PopMats/s',ss,'/d',dd,'/PopMat'); % Patch population files.
tmpPopMat = strcat('../../BiomathTalk/rand1/Cubical/PopMats/s',ss,'/d',dd,'/PopMat'); % Patch population files.

fig = figure;

Betti_storage = ''; % Empty string for later concatenation of Betti numbers.
abundance = zeros(1,T); % Array of zeros T long. Sum of elements in N into array.

matfile = sprintf('../../BiomathTalk/rand1/rand_%g.txt',IC); % MUST CHANGE when different.
N = importdata(matfile); % Initial Condition matrix. 
N = N * s; % Rescale spatially.

tmpPopMat_0 = strcat(tmpPopMat,'_0');  
dlmwrite(tmpPopMat_0, N)
    
D = diag(diag(ones(n),1),1) + diag(diag(ones(n),-1),-1);

% name = sprintf('../../BiomathTalk/rand1/Cubical/Movies/Mov_n%g_d%g_r%g_T%g_thr%g_eth%g_s%g_IC%g.avi',n,d*100,r,T,thr,eth*100,s,IC);
name = sprintf('../../BiomathTalk/rand1/Cubical/Movies/Mov_n%g_d%g_r%g_T%g_thr%g_eth%g_s%g_IC%g.avi',n,d*100,r,T,thr,eth*100,s,IC);

aviobj = avifile(name,'compression','None');

for t = 1:T % -------------- Simulation code -------------------
    ti = num2str( t );
    
    N(N < eth) = 0; % Extinction.
    
    % ------- Betti Calculations ------- 
    [I,J] = find(N > thr); 
    
    if (isempty(I) == 0)
        tmp_fileID = fopen(tmp_string,'w'); % Open tn file with a temporary ID so that it's writable.
        fprintf(tmp_fileID,'(%g,%g) \n',[I,J]'); % Write data to a text file.
        fclose(tmp_fileID);
        
    % Send temporary file to CHomP, call might need to be changed depending
    % on how CHomP is set up on your computer (i.e. './chomp %s' instead).
        command = sprintf('chomp %s',tmp_string);
        [status, betti] = unix(command);
        Betti_storage = [Betti_storage betti]; % Concatenate Betti number to string.
        
    else
        Betti_storage = [Betti_storage '0 0 0\n']; % Concatenate '0 0 0\n' to Betti str?
    end
    
    clear I J; % Clear these so we can use them for the next t.
    
    % ------- Abundance data ---------------------------------------			
    abundance(t) = sum(sum(N)); 
    
    if abundance(t) == 0 % If global extinction => stop simulation before reaching final T.
        SAVEthisDATA(n,d,r,T,thr,eth,s,IC,Betti_storage,abundance);
        command = sprintf('rm %s',tmpChompFile);
        unix(command);
        return;
    end 

    % ------- Iterate simulation ---------------------------------------
    N = r*N.*exp(-N/s); % Ricker growth function.
    N = (1 - d)*N + (d/4)*(N * D) + (d/4)*(D * N); % Dispersal phase.
%     
%     f = r*exp(-N/s);                         % Ricker growth per patch.
%     N = poissrnd(N.*f);                      % Stoch. Ricker growth all patches.
%     N = (1 - d)*N + (d/4)*(N*D + D*N);       % Det. dispersal.
%     N = floor(N);   
    
        
        
    tmpPopMat_t = strcat(tmpPopMat,'_',ti);
    dlmwrite(tmpPopMat_t,N)
    
    % ------- MOVIE --------------------------------------- 
    if(mod(t,2)==0) % To account for the period 2 flipping (color flips).
        M = zeros(n,n,3);
        M(:,:,1) = 204/255;
        M(:,:,2) = 255/255;
        M(:,:,3) = 204/255;

        [I J] = find(N > thr); % Find treats the matrix as a long vector (row then row).

            for k = 1:length(I)
                M(I(k),J(k),1:3) = [102/255 102/255 51/255];
            end
            
        im = image(M);
        aviobj = addframe(aviobj, fig);

    end
%     
%         M1 = zeros(n,n,3);
%         M1(:,:,1) = 204/255;
%         M1(:,:,2) = 255/255;
%         M1(:,:,3) = 204/255;
% 
%         [K L] = find(N > thr); % Find treats the matrix as a long vector (row then row).
% 
%             for k = 1:length(K)
%                 M1(K(k),L(k),1:3) = [102/255 102/255 51/255];
%             end
%             
%         im1 = image(M1);
%         fn1 = strcat('../../BiomathTalk/rand1/shot00_',ti,'.png');
%         saveas(im1, fn1)
    
end

aviobj = close(aviobj);

SAVEthisDATA(n,d,r,T,thr,eth,s,IC,Betti_storage,abundance);
command = sprintf('rm %s',tmpChompFile);
unix(command);

end