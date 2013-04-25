
% ~~~~~~~~~~~~~~~~~~~ Bifurcation diagram.
% ~~~~~~~~~~~~~~~~~~~ Single patch.
% ~~~~~~~~~~~~~~~~~~~ Deterministic and/or Stochastic fitness.

IC = 0.5;

iter = 1000; % Number of iterations.

s = 100000;    % Spatial rescaling. 

IC = IC * s; % Rescale IC.

x = zeros(iter + 1,1); % Trajectory array.
% x(1) = IC;        % Initial population value of patch.

Rmin = 0.5;  % Vary fitness, r, from Rmin to Rmax.
% Rmax = 30;
Rmax = 70;
% numR = 200;  % Number of r values to use.
numR = 300;
step = (Rmax - Rmin)/numR; % Distance between r values.

t = zeros(numR + 1,1);     % Array of r values.
z = zeros(numR + 1,250);   % Long term behavior array.

t2 = zeros(numR + 1,1);
w = zeros(numR + 1,250);

for j = 1:numR + 1
    t(j) = (j - 1)*step + Rmin; % Input r values into array.
    R = t(j);                   % Current r value.
    
    x(1) = IC*log(R);           % Rescale IC. and reset the next traj.
    
    for i = 1:iter % Apply stochastic Ricker growth.
        
%         x(i + 1) = x(i) * (R*exp(-x(i)/s)); % Spatial Ricker, NO stochasticity.
        
        f = R*exp(-x(i)/s); % Stochstic Ricker growth, step 1.
        x(i + 1) = poissrnd(x(i) * f); % Stochastic Ricker growth, step 2.
        
        if (i > 750)            % Only show last 250 iterations.
            z(j,i - 750) = x(i + 1);
            
%             tmp = (z == 0);
            if (z(j,i - 750) == 0)
                tmp = R;
                w(j,i - 750) = z(j,i - 750);
            end
            
            t2(j) = tmp;
        end
    end
end

% ~~~~~~~~~~~~~~~~~~~ Plot long term stable behavior.
figure()
plot(t,z,'r.','MarkerSize',4)
hold all
plot(t2,w,'b.','MarkerSize',4)
% axis([-1 31 -2 10*s])
axis([-1 71 -2 10*s])
xlabel(' r ','FontSize',12); 
ylabel(' Population density ','FontSize',12);
str = sprintf(' Ricker map Bifurcation diagram, s%g, IC%g ',s,IC);
title(str,'FontSize',16,'FontWeight','bold');

% ~~~~~~~~~~~~~~~~~~~ Save the figure. 
ic = num2str(IC);
S = num2str(s);

fn = strcat('../../../Desktop/school/Data_2013March/BifDiags_1D/StochGrowthBif_s',S,'_IC.5.png'); 
saveas(gcf, fn)


