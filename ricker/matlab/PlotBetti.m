% Converts Betti numbers from .txt into a matrix and plots them. 

%a = load('SmallGrid/bettiPD_10_0.2_500_2_0_1.txt');
%IC = num2str(4);
n = 21;
r = 22;
s = 1;
d = 0.15;
thr = log(r);
eth = 0.00;

dd = num2str(d*100);
%thr = num2str(thresh);
% eth = num2str(ethresh*100);

name = sprintf('../../BiomathTalk/rand1/Cubical/Stoch/Bettis/bettiPd_n21_d%g_r%g_T50_thr%g_eth%g_s%g_IC1.txt',100*d,r,s*thr,s*eth,s);
a = load(name);

% b_range = 1:1:length(a);
b_range = 1:2:length(a);
b_range = b_range.';

b0 = a(:,1);
beta0 = b0(1:2:length(a));
avgb0 = mean(b0);

b1 = a(:,2);
beta1 = b1(1:2:length(a));
avgb1 = mean(b1);

p = figure();
% scatter(b_range,b0,'m')
scatter(b_range,beta0,'m')
hold on
% scatter(b_range,b1,'c')
scatter(b_range,beta1,'c')
str = sprintf('Betti numbers \n r = %g, d = %1.2f, s = %g, thr = %1.4f eth = %1.3f',r,d,s,thr,eth);
title(str);
xlabel('Time (1,...,1000)');
ylabel('Betti number');
loc = 'EastOutside';
legend('\beta_0','\beta_1','location',loc);



% fn = strcat('../Prospectus/Bettis/bettiPD_n51_d',dd,'_r22_t1000_thr',thre,'_eth',eth,'_IC1.png');%%%%%%%
%fn = strcat('../Prospectus/Bettis/bettiPD_n51_d',dd,'_r22_t1000_thr',thre,'_eth5half_IC1.png');%%%%%%%
fn = sprintf('../../BiomathTalk/rand1/Cubical/Stoch/bettiPd_n21_d%g_r%g_T50_thr%g_eth%g_s%g_IC1.png',100*d,r,s*thr,s*eth,s);
saveas(p, fn)

%disp('Average of the zeroth Betti number is')
%disp(avgb0)
%disp('Average of the first Betti number is')
%disp(avgb1)