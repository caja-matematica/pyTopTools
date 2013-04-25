% ~~~~~~~~~~~ Converting to stochastic growth.
TRY = num2str(1);

n = 1; % Grid size -> 1D Ricker map.
r = 16.75; % Fitness parameter.
time = 1000; % Number of iterations.

s = 1; % Spatial rescaling. 

IC = log(r); % Initial condition of population on patch. 
IC = IC*s;
N = zeros(1,time); % Trajectory array.

f = @(x) r*exp(-x/s); % Ricker growth. rescaled.

N(1) = IC;

for t = 1:time - 1
%     N(t + 1) = N(t)*(r*exp(-N(t)/s)); % Deterministic growth.
    N(t + 1) = poissrnd( N(t)*f(N(t)) ); % Stochastic growth.
end

% disp(N)

plot(1:time,N,'mo')

% nn = num2str(n);
% rr = num2str(r);
% tt = num2str(time);
% aa = num2str(s);
% 
% fname = strcat('../ResearchBlogs/StochGrowth_1D/n',nn,'_r',rr,'_t',tt,'_a',aa,'.png');
% fn = sprintf(fname);
% saveas(gcf, fname)


%  Just saving some practice runs, should delete.
R = num2str(r);
ic = num2str(IC);
S = num2str(s);

fn = strcat('../../../Desktop/school/Data_2013March/Trajs/StochGrowth/StochGrowth_r',R,'_s',S,'_IClogr_try',TRY,'.png'); 
saveas(gcf, fn)



