% Code for the Ricker map code. MATLAB function used to
% Compute Ricker simulations
%{
         n :  the size of the matrix is nxn.
 dispersal :  an input dispersal parameter that is percentage of a patch's 
              population that is dispersed. For our paper 0 \leq d \leq 0.40. 
      time :  length of time to run the simulation. 
    thresh :  value at-which to threshold entries of the population matrix 
              (N in code and P in the paper). Fixed at 10^(-4) for our paper. 
              Because  ethres>0 this is effectively thresholding at zero.
   ethresh :  local extinction threshold. 
        IC :  number that signifies which initial condition to use for the 
              simulation.
%}

function reAP3(n,dispersal,coef,time,thresh,ethresh,IC)

% Create a temporary file that is passed to CHomP:
tmpChompFile = sprintf('tmpChomp_%g_%g_%g_%g_%g_%g_%g_.txt',n,dispersal*100,coef,time,thresh,ethresh,IC);
                % Modified. Format data into string.
tmp_string = sprintf(tmpChompFile); % Format data into string.

fig = figure;

% Parameters:
%%coef = 22; % Used in Ricker function. It's the growth rate coefficient.
Betti_storage = ''; % Empty string for later concatenation of Betti numbers.
abundance = zeros(1,time); % Array of zeros T long. Later, for t = 1:T, we'll
                        % input the sum of elements of matrix N into 
                        % position t in array.

% Initial Condition matrix.                        
    matfile = sprintf('../n50/BestParams/rand_%g.txt',IC); % MUST CHANGE when different.
    N = importdata(matfile); 
    
D = diag(diag(ones(n),1),1) + diag(diag(ones(n),-1),-1);
% Dispersal matrix for uniform dispersal. With ones on the super and
% inferior diagonals.

% Simulation code.
name = sprintf('../n50/BestParams/Mov_n%g_d%g_r%g_t%g_thr%g_eth%g_IC%g.avi',n,dispersal*100,coef,time,thresh,ethresh,IC);
% 'ICrands/Mov_percDisp%g_coef%g_time%g'

aviobj = avifile(name,'compression','None');
% -------------- Simulation code -------------------
for t = 1:time
    
    % ---------- Extinction -------------
    % Find all entries less than ethres & set them to zero. We don't do
    % negative populations.
    N(N < ethresh) = 0;
    
    % ------- Betti Calculations ------- 
    % This is not the fastest way to compute Betti numbers from CHomP but  
    % it works and is easy.
    [I,J] = find(N > thresh); % I indicates row position and J indicates 
                              % column position.
    
    if (isempty(I) == 0)
        tmp_fileID = fopen(tmp_string,'w'); % Open tn file with a temporary  
                               % ID so that it's writable.
        fprintf(tmp_fileID,'(%g,%g) \n',[I,J]'); % Write data to a text file.
        fclose(tmp_fileID);
        
    % Send temporary file to CHomP, call might need to be changed depending
    % on how CHomP is set up on your computer (i.e. './chomp %s' instead).
        command = sprintf('chomp %s',tmp_string);
        [status, betti] = unix(command);
        Betti_storage = [Betti_storage betti]; % Concatenate Betti number to string.
        
    else
        Betti_storage = [Betti_storage '0 0 0\n']; % Concatenate '0 0 0\n' to Betti str?
       
    end % End print if.
    
    clear I J; % Clear these so we can use them for the next t.
    
    % ------- Abundance data ---------------------------------------			
    abundance(t) = sum(sum(N)); 
    % sum(N) sums the columns, sum(sum(N)) adds up all the elements of 
    % N. All negative numbers have been set to zero, remember.
    % Gives the total population of the whole model.
    
    % ------- Check for Extinction ---------------------------------
    % If extinction occurred, then stop the simulation before T is reached.
    % No need to go on with simulation is popualtion has died out.
    % Also, recall that all elements of N are >= 0, so the sum is zero if
    % and only if each element is zero.
    if abundance(t) == 0
        SAVEData(n,dispersal,coef,time,thresh,ethresh,IC,Betti_storage,abundance);
        command = sprintf('rm %s',tmpChompFile);
        unix(command);
        return;
    end % End if statement.

    % ------- Iterate simulation ----------------------
    N = coef*N.*exp(-N); % Ricker growth function (new N comes from old N).
    % N = R.*N.*(1-N); % Logistic growth function used in appendix of the 
                       % paper (remember to change R).
                 % N.*exp(-N) is the Haddamard product of N and exp(-N).
                 % (element by element multiplication).
                 % exp(-N) takes the exponent of each element of matrix N.

    % Dispersal phase
    N = (1-dispersal)*N + (dispersal/4)*(N * D) + (dispersal/4)*(D * N);
    
    % MOVIE--------begin
    if(mod(t,2)==0) % To account for the period 2 flipping (color flips).
        %M = N; % want this to be a cell array .. a 3-tuple (RGB value).
        M = zeros(n,n,3);
        M(:,:,1) = 204/255;
        M(:,:,2) = 255/255;
        M(:,:,3) = 204/255;

        [I J] = find(N > thresh); % find treats the matrix as a long vector (row then row).

            for k = 1:length(I)
                M(I(k),J(k),1:3) = [102/255 102/255 51/255];
                %M(I,J,2) = R2;
                %M(I,J,1) = R2;
                %j = find(N <= thresh);
                %M(i) = 255; % Colors... (Could look into RGB scale.).
                %M(j) = 25; % lie 45 and 75 and 25.
            end
            
        im = image(M);
        aviobj = addframe(aviobj, fig);
        %MovieMatrix(count)= getframe();
        %count = count + 1;
    end
    %MOVIE-------end
    
    
end % End the for loop.

aviobj = close(aviobj);

SAVEData(n,dispersal,coef,time,thresh,ethresh,IC,Betti_storage,abundance);
command = sprintf('rm %s',tmpChompFile);
unix(command);


%movie2avi(MovieMatrix,name);


end