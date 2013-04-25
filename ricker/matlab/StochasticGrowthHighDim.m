% =======================================================================
% ~~~~~~~~~~~~~~~~ Converting to stochastic growth. ~~~~~~~~~~~~~~~~~~~~~
% =======================================================================

N = load('ICTinyrand/rand_1.txt'); % Import IC matrix.
n = length(N); % Dimension of square grid is n^2.

r = 13;        % Fitness parameter.
time = 1000;   % Number of iterations.
s = 1000;      % Spatial rescaling. 

trajN = zeros(n,n,time); % Initialize trajectory array.
trajN(:,:,1) = N; % Plug in IC as first entries of trajectories.

D = diag(diag(ones(n),1),1) + diag(diag(ones(n),-1),-1);

f = @(x) r*exp(-x/s); % Ricker growth function. (rescaled: N' = a*N).

for t = 2:time % Third dimension of matrix (trajectory length).
    
    N = poissrnd(N.*f(N)); % Stochatsic growth.
    N = (1 - d)*N + (d/4)*(N*D) + (d/4)*(D*N); % Deterministic dispersal.
    N = floor(N); % Convert all elements in matrix to integers.
    
    trajN(:,:,t) = N; % Record each time step.
    
end

disp(trajN) % Display the matrix that records the trajectiories.
% disp(N)   % Display final  matrix.

tmp1 = zeros(1,time); tmp1(:) = trajN(1,1,:);
tmp2 = zeros(1,time); tmp2(:) = trajN(1,2,:);
tmp3 = zeros(1,time); tmp3(:) = trajN(1,3,:);
tmp4 = zeros(1,time); tmp4(:) = trajN(2,1,:);
tmp5 = zeros(1,time); tmp5(:) = trajN(2,2,:);
tmp6 = zeros(1,time); tmp6(:) = trajN(2,3,:);
tmp7 = zeros(1,time); tmp7(:) = trajN(3,1,:);
tmp8 = zeros(1,time); tmp8(:) = trajN(3,2,:);
tmp9 = zeros(1,time); tmp9(:) = trajN(3,3,:);

figure(1)
p1 = plot(1:time,tmp1,'o'); hold all;
p2 = plot(1:time,tmp2,'o'); hold all;
p3 = plot(1:time,tmp3,'o'); hold all;
p4 = plot(1:time,tmp4,'+'); hold all;
p5 = plot(1:time,tmp5,'+'); hold all;
p6 = plot(1:time,tmp6,'+'); hold all;
p7 = plot(1:time,tmp7,'*'); hold all;
p8 = plot(1:time,tmp8,'*'); hold all;
p9 = plot(1:time,tmp9,'*'); hold all;

% =======================================================================
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% =======================================================================
