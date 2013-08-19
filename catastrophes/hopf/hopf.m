% mu1,mu2 = the initial and final parameters values to be explored
% x0,y0 = the initial conditions
% sigma = the intensity of the noise
% M. Gidea 2011 
t0 = 0;
tfinal = 100;
stepsize=0.1;
steps=floor((tfinal-t0)/stepsize);
x0 = 0.1; 
y0 = 0.1;
sigma=0.1;
mu1 = -3;
mu2 = 0;
%mu0 = 0;
% Euler method for stochastic ODE 
x=zeros(steps,1);
y=zeros(steps,1);
t=zeros(steps,1);
mu=zeros(steps,1);
x(1) = x0; 
y(1) = y0;
t(1) = t0; 
mu(1)=mu1;
for k = 1:1:steps,
x(k+1)=x(k)+stepsize*(mu(k)*x(k)-y(k)-x(k)*y(k)*y(k))+sigma*randn*sqrt(stepsize);  
y(k+1)=y(k)+stepsize*(x(k)+mu(k)*y(k)-y(k)*y(k)*y(k))+sigma*randn*sqrt(stepsize);
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
dlmwrite('hopf.txt',z)