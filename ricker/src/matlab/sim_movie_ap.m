
% This program allows the user plug in inputs and then it runs reAP3.m and  
% simulates a movie with the same parameters. 

clear;
% This is a matlab version modified from Ben's use of SpatialChaos.R from 
% Schreiber. 
fig = figure;
% Input your parameters here:
thresh = 0.0001; % All matrix entries above this value will be saved to a txt file.
ethresh = 0.06; % Local time extinction threshold.
dispersal = 0.15; % This is the dispersal rate.
n = 51; % This is the grid size (assume a square grid nXn).
time = 100; % Number of iterations (number time steps).
coef = 22; % Ricker coefficient (sometimes noted r or R).
IC = 1; % Just plugging this in for now. CHANGE LATER WITH REAL IC's. 
        % (Currently a randn matrix).

% Run the program to produce Betti numbers and abundance totals.
%reAP3(dispersal,time,thresh,ethresh,IC)

%skip = 2; % How often to save the data for the movie. Dr. Schreiber often 
% chooses 2 as there often is a period 2 like nature to the dynamics. 
% Skips every other frames (better looking with changing colors).

% Initial Condition
%N = zeros(n); % The matrix of abundnace on the spatial grid.
%N(ceil(n/2), ceil(n/2)) = 50; % Initial abundance in the center of the grid.
matfile = sprintf('ICrands/rand_%g.txt',IC); % MUST CHANGE when different.
N = importdata(matfile); 
% Construct dispersal matrix.
D = diag(diag(ones(n-1)),1) + diag(diag(ones(n-1)),-1); 

%count = 1;
d = dispersal*100;
eth = ethresh*100;
% Simulation code.
name = sprintf('Mov_n%g_d%g_r%g_t%g_thr%g_eth%g_IC%g.avi',n,d,coef,time,thresh,eth,IC);
% 'Movies/Mov_percDisp%g_coef%g_time%g'
aviobj = avifile(name,'compression','None');
for t = 0:time
    N = coef * N.*exp(-N);
    N = (1 - dispersal)*N + (1/4)*dispersal*(D*N + N*D);

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
        %MovieMatrix(count)= getframe(fig);
        aviobj = addframe(aviobj, fig);

        %count = count + 1;
    end

end

aviobj = close(aviobj);

%movie2avi(MovieMatrix,name); 

%close();
