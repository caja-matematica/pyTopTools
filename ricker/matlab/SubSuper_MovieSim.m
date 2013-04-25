% =======================================================================
%                  Sub/Super-Level Set Movie Simulation.
% =======================================================================

% INCOMPLETE-- WORK ON LATER.

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                                                               Function.
%                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                   Movie simulation of sub/super-level sets progression.
%                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% 
% 
function SubSuperMovie(N,timestep,ThresStep,s)
            % N: Population abundance matrix. (string format)
                    % ex: N = 'Trial10x10/PopMat_1'
            % timestep: Single time step for thresholding.
            % ThreshStep: Thresholding step values.
                
    % Check sub (= 0) vs. super (= 1) -level set input is correct.
    if (s == 0)     % Sub-level sets.
        level = 'Sub';
    elseif (s == 1) % Super-level sets.
        level = 'Super'; 
    else            % Incorrect input.
        error( 'Input for s must be 0 for SUB-level set or 1 for SUPER level set, NO other options are accepted!' )

    % PopMat = load('Trial10x10/PopMat_1');
    PopMat = load(N); % Load in population matrix at one time step.
    
    minPop = min(min(PopMat)); % Designates where to start threshold.
    maxPop = max(max(PopMat)); % Designates where to stop threshold.
    
    n = length(PopMat); % Number of coloumns of square matrix.
    t = timestep; % Only on one time step. Before any iterations.
    
    dim = num2str(length(N)); % For movie filename.
    time = num2str(timestep); % For movie filename.

    fig = figure; % Need for images as frames of movie.
    
    % name = sprintf('Mov_SubLevel_n%g_T%g.avi',n,T); % Movie filename.
        
    name = strcat('Mov_',level,'_n',dim,'_t',time,'.avi');
    
    
    aviobj = avifile(name,'compression','None');

    for k = minPop - 0.01:ThreshStep:maxPop + 0.1
    
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

end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



% =======================================================================
%                                  End.
% =======================================================================