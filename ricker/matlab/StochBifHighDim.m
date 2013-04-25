
% ~~~~~~~~~~~~~~~~~~~ Bifurcation diagram.
% ~~~~~~~~~~~~~~~~~~~ Single patch.
% ~~~~~~~~~~~~~~~~~~~ Deterministic and/or Stochastic fitness.

% file = strcat('../../../Desktop/school/Data_2013March/BifDiags10x10/rand_1.txt');
file = strcat('../../../Desktop/school/Data_2013March/11x11/rand_1.txt');
M = load(file);                    % Load in an initial condition matrix.

s = 1000000;       S = num2str(s);     % Spatial rescaling. 
d = 0.30;      D = num2str(d*100); % Dispersal.
n = length(M); N = num2str(n);     % Dim. of grid size.
iter = 1000;                       % Number of iterations.

ic = s*M;                          % Scaled initial condition.
% (This scaling makes the IC an appropriate choice based on the spatial scale factor, s.)

Rmin = 0.5;    Rmax = 70;          % Vary fitness, r, from Rmin to Rmax.
numR = 300;                        % Number of r values to use.
step = (Rmax - Rmin)/numR;         % Distance between r values.

abun = zeros(iter + 1,1);          % Abundance array.
abun(1) = sum(sum(M));             % Initial total population.

centertraj = zeros(numR + 1,250);  % Center patch trajectory array.
center = ceil(n/2);

edgetraj = zeros(numR + 1,250);    % Edge patch trajectory array. (Right middle patch).

x = zeros(n,n,iter + 1);           % Trajectory array for the summ of all elts. in M.
x(:,:,1) = M;                      % Initial individual patch populations.

t = zeros(numR + 1,1);             % Array of r values.
z = zeros(numR + 1,250);           % Long term behavior array.

t2 = zeros(numR + 1,1);            % Array of r values to ...
w = zeros(numR + 1,250);           % ... plot zero values in blue.

D1 = diag(diag(ones(n),1),1);      % Ones on the super-diagonal.
D2 = diag(diag(ones(n),-1),-1);    % Ones on the sub-diagonal.
Di = D1 + D2;                      % Dispersal matrix.

for j = 1:numR + 1
    t(j) = (j - 1)*step + Rmin;    % Input r values into array.
    R = t(j);                      % Current r value.
    
    M = ic*log(R);                 % Scale IC based on the current fitness value.
    
    for i = 1:iter                 % Apply Ricker growth.
        
%         M = R*M.*exp(-M/s);                      % Det. Ricker growth.
%         M = (1 - d)*M + (d/4)*(M*Di + Di*M);     % Det. dispersal phase.

        f = R*exp(-M/s);                         % Ricker growth.
        M = poissrnd(M.*f);                      % Stoch. Ricker growth.
        M = (1 - d)*M + (d/4)*(M*Di + Di*M);     % Det. dispersal.
        M = floor(M);                            % Convert elts to integers.
%         M = ceil(M);                            % Convert elts to integers.

        x(:,:,i + 1) = M;                        % Record trajcetory of each patch.
        abun(i + 1) = sum(sum(M));               % Sum all patch pop.'s (total abundance).
        centerpatch = M(center, center);         % Define center patch value.
        edgepatch = M(center, n);         % Define center patch value.
        
        if (i > 750)                             % Only show last 250 iterations.
            z(j,i - 750) = abun(i + 1);          % Record traj. of total abundance.
            centertraj(j,i - 750) = centerpatch; % Record traj. of center patch.
            edgetraj(j,i - 750) = edgepatch;     % Record traj. of edge patch.

            if (z(j,i - 750) == 0)               % Zeros plotted in blue.
                tmp = R;
                w(j,i - 750) = z(j,i - 750);
            end
            
            t2(j) = tmp;
        end
    end
end

% pathname = strcat('../../../Desktop/school/Data_2013March/11x11/rand1');
% 
% figure()
% plot(t,z,'r.','MarkerSize',4)  % Plot long term stable behavior.
% hold all
% plot(t2,w,'b.','MarkerSize',4) % Plot zeros in blue.
% axis([0 70 0 20*s])
% xlabel(' r ','FontSize',12); 
% ylabel(' Population density ','FontSize',12);
% str = sprintf(' Ricker Bif diagram for ABUNDANCE, s%g   ',s);
% title(str,'FontSize',16,'FontWeight','bold');
% 
% fn = strcat(pathname,'/AbundanceBifs/StochGrowthBif_n',N,'_d',D,'_s',S,'_Abundance.png'); 
% saveas(gcf, fn)                % Save the figure.



% Run this simulation separately.

figure()
plot(t,centertraj,'r.','MarkerSize',4)  % Plot long term stable behavior.
hold all
plot(t2,w,'b.','MarkerSize',4) % Plot zeros in blue.
axis([0 70 0 30*s])
xlabel(' r ','FontSize',12); 
ylabel(' Population density ','FontSize',12);
str = sprintf(' Ricker Bif diagram of CENTER PATCH, s%g ',s);
title(str,'FontSize',16,'FontWeight','bold');

fn1 = strcat(pathname,'/CenterPatchBifs/StochGrowthBif_n',N,'_d',D,'_s',S,'_CenterPatch.png'); 
saveas(gcf, fn1)                % Save the figure.



% Run this simulation separately.

% figure()
% plot(t,edgetraj,'r.','MarkerSize',4)  % Plot long term stable behavior.
% hold all
% plot(t2,w,'b.','MarkerSize',4) % Plot zeros in blue.
% axis([0 70 0 3*10*s])
% xlabel(' r ','FontSize',12); 
% ylabel(' Population density ','FontSize',12);
% str = sprintf(' Ricker Bif diagram of EDGE PATCH, s%g ',s);
% title(str,'FontSize',16,'FontWeight','bold');
% 
% fn2 = strcat(pathname,'/EdgePatchBifs/DetGrowthBif_n',N,'_d',D,'_s',S,'_EdgePatch.png'); 
% saveas(gcf, fn2)                % Save the figure.





