% =======================================================================
%                    Sub-Level Set Movie Simulation.
% =======================================================================



% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                         Movie simulation of sub-level sets progression.
%                         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% 
% Load in the populations matrix at the desired time step.
PopMat = load('Trial10x10/PopMat_1');

n = length(PopMat); % Number of coloumns of square matrix.
T = 0; % Only on one time step. Before any iterations.

% Values that designate where to start and stop your thresholding.
minPop = min(min(PopMat));
maxPop = max(max(PopMat));

fig = figure;

% Title your movie.
name = sprintf('Mov_SubLevel_n%g_T%g.avi',n,T);
aviobj = avifile(name,'compression','None');

for k = minPop - 0.01:0.01:maxPop + 0.1
    
    %M = N; % want this to be a cell array .. a 3-tuple (RGB value).
    M = zeros(n,n,3);
    M(:,:,1) = 250/255;
    M(:,:,2) = 250/255;
    M(:,:,3) = 250/255;

    % < find > treats the matrix as a long vector (row then row).
    [I J] = find(PopMat <= k); 

    for j = 1:length(I)
        M(I(j),J(j),1:3) = [50/255 50/255 50/255];
    end
	
    im = image(M);
    aviobj = addframe(aviobj, fig);
    
end

aviobj = close(aviobj); 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



% =======================================================================
%                               End.
% =======================================================================