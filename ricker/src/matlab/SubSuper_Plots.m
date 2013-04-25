% =======================================================================
%                     Plot sub or super-level sets.
% =======================================================================



% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%         Binary plots of sub-level sets - one time step, many threholds.
%         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Load in population matrix.
% PopMat = load('../IC200rands/rand1/OverTime/Pers_n200_d0_r5_T100_s1000_IC1_0_sub.txt');
% PopMat = load('IC200rands/rand1/rand_1.txt');
PopMat = load('../IC200rands/rand1/rand_1.txt');

% % PopMat = [2 2 2 ; 2 3 2 ; 2 2 2 ];
% % PopMat = ( PopMat * 1000 );

PopMat = floor(1000*PopMat); % Make elements integer valued.


% Find the minimum and maximum of the population matrix.
minPop = min(min(PopMat));
maxPop = max(max(PopMat));

% % PopMat = -PopMat + maxPop + 1 ;

% Create 3D figure where color & bar heights indicate patch pop. abundance.
figure()
h = bar3(PopMat);
colormap summer
colorbar;
numBars = size(PopMat,1);
numSets = size(PopMat,2);
for i = 1:numSets
    matdata = ones(6*numBars, 4);
    k = 1;
    for j = 0:6:(6*numBars - 6)
        matdata(j + 1:j + 6, :) = PopMat(k,i);
        k = k + 1;
    end
    set(h(i),'Cdata',matdata)
end
% 
% % Create (binary) plot of sub-level sets for 'k' thresholds.
% for k = minPop:100:maxPop
%     
%     P = PopMat;    % Need to reset population matrix each time.
%     P(P <= k) = 0; % Sub-level sets for threshold 'k'.
%     P(P > k) = 1;  % Patches with abundance too high to be sub-level sets.
%                    % (i.e., super-level sets).
%     figure();    
%     imagesc(P)
%     colormap pink  % Other color scheme options:
%                    % copper, winter, summer, spring, jet, hot, cool, 
%                    % bone, copper, pink, lines, HSV, autumn, gray.
%     % colormap(flipud(colormap)) % To reverse color of numerical values.
%     colorbar; % Show colorbar for associated numerical values.
% end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% k = 500;
% P = PopMat;    % Need to reset population matrix each time.
% P(P <= k) = 0; % Sub-level sets for threshold 'k'.
% P(P > k) = 1;  % Patches with abundance too high to be sub-level sets.
%                % (i.e., super-level sets).
% figure();    
% imagesc(P)
% colormap pink  % Other color scheme options:
%                % copper, winter, summer, spring, jet, hot, cool, 
%                % bone, copper, pink, lines, HSV, autumn, gray.
% % colormap(flipud(colormap)) % To reverse color of numerical values.
% colorbar; % Show colorbar for associated numerical values.


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                               Other plotting options for such matrices.
%                               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Create a simple matrix to play with.
% A = [ 0 1 0; 1 0 1; 0 1 1];
% A = load('Trial10x10/PopMat_1');
% 
% Some plotting options are:
% imagesc(A)
% stem3(A)
% contour(A)
% contour3(A)
% contourf(A)
% mesh(A)
% ngrid(A)
% surf(A)
% pcolor(A)
% meshz(A)
% surfc(A)
% waterfall(A)
% ribbon(A)
% surface(A)
% spy(A)
% 
% colorbar % Always show color bar to see associated numerical values.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



% =======================================================================
%                               End.
% =======================================================================
