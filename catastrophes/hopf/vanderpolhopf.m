% Define initial conditions.
t0 = 0;
tfinal = 100;
stepsize=0.1;
steps=(tfinal-t0)/stepsize;
x0 = .1; 
y0 = .1;
sigma=0.1;
mu0   =  0;
mu1   =  -1;
mu2   = 1;
% Euler method for stochastic ODE 
x(1) = x0; 
y(1) = y0;
t(1) = t0;
mu(1)= mu1; 
for k = 1:1:steps, 
x(k+1)=x(k)+stepsize*(y(k))+sigma*randn*sqrt(stepsize);
y(k+1)=y(k)+stepsize*(-x(k)+(2*mu(k)-x(k)^2)*y(k))+sigma*randn*sqrt(stepsize);  
t(k+1)=t(k)+stepsize;
mu(k+1)=((t(k+1)-t0)/(tfinal-t0))*(mu2-mu1)+mu1;
end
subplot(1,2,1)
% The plot coordinates are x and y
plot (x(:),y(:),'r')
title ('Phase space')
subplot(1,2,2)
% The plot coordinates are mu (the parameter) and x (the x-variable) 
plot(mu(:),x(:),'r')
title('Bifurcation diagram')
z=[mu(:), x(:)];
% Write data in the file hopf.txt 
dlmwrite('vanderpolhopf.txt',z)