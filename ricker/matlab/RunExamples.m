% =======================================================================
%                           Trials for codes.
% =======================================================================


% % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % % 
s = 1;
ss = num2str(s);

d = 0.15; 
dd = num2str(100*d);

r = 22; 
rr = num2str(r);

n = 21; 
nn = num2str(n);

T = 50; 
TT = num2str(T);

thr = s * log(r); % OR ZERO.
% thr = 0;

eth = s * 0.00; 

IC = 1;
ic = num2str(IC);

reAP_scale(n,d,r,T,thr,eth,s,IC)


% N = load('../../BiomathTalk/rand1/rand_1.txt'); % Import IC matrix.
% 
% path2 = strcat('../../BiomathTalk/rand1/Persistent');
% 
% tmpPerseus = strcat(path2,'/Pers_n',nn,'_d',dd,'_r',rr,'_T',TT,'_s',ss,'_IC',ic);
% 
% tmpPopMat = strcat(path2,'/PopMat'); 
% 
% tmpPopMovie = strcat(path2);
% 
% 
% % reAP3_perseus_det( n, d, r, T, s, IC, N, tmpPerseus, tmpPopMat, tmpPopMovie )
% 
% 
% for beta = 0:1              % Create full persistent diagrams.
%     bb = num2str(beta);
%     
%     for j = 0:T
%         jj = num2str(j);
% 
%         Pers = strcat(tmpPerseus,'_',jj,'_full_',bb,'.pdia');
% 
%         persdia_sub_super(Pers)
% 
%         fn = strcat(tmpPerseus,'_',jj,'_full_',bb,'.png'); % Save persistence diagrams.
%         saveas(gcf, fn)
%         
%         close all
%     end
% end


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                                          Cubical homology computations. 
%                                          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
%   [[[ Computes Betti #'s through CHoMP and makes movie simulations. ]]]
% 
% reAP3(n, dispersal, coef, time, thresh, ethresh, IC) % Fill in variables.
% reAP3(50,0,5,50,log(5),0.00,1) % CHoMP and movies.
% 
% {{{ [implay] Plays movie simulations from path in MATLAB. }}}
% implay('../../Desktop/School/Data_2012Prospectus/ProspectusMovies/ValueMovies/Mov_n51_d15_r22_t1000_thr0001_eth20_IC1.avi')
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                       Persistent homology computations (ONE time step).
%                       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
%       [[[ Compute persistent homologies on one time step. ]]]
% 
% N = load('Trial10x10/rand_1.txt'); % Import IC matrix.
% process_perseus( N, tmpPerseus, scale, pers_type )
% process_perseus(N,'rand1',1000)
% 
% 
%               [[[ Create persistence diagram. ]]]
% 
% persdia('rand1_scale1000000_pers_1.txt')
% 
% r = 22;
% d = 0.15;
% D = diag(diag(ones(n),1),1) + diag(diag(ones(n),-1),-1);
% N = load('Trial10x10/rand_1.txt');
% 
% % Second iteration of map.
% N = r*N.*exp(-N); 
% N = (1 - d)*N + (d/4)*(N * D) + (d/4)*(D * N); 
% 
% process_perseus(N,'rand1second',100)
% persdia('Trial10x10/rand1second_pers_1.txt')
%
% Create 3D plot showing abundances on patches.
% bar3(N)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                     Persistent homology computations - MANY time steps.
%                     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
%          [[[ Compute persistent homologies over time. ]]]
% 
%   *** MAKE SURE CORRECT N IS INPUTED INTO reAP3_perseus_det.m !!!! ***
% 
% reAP3_perseus_det( dimN, dispersal, fitness, final_time, scale, IC )
% reAP3_perseus_det( n, d, r, final_time, scale, IC )
% reAP3_perseus_det( 10, 0.15, 5, 50, 1000, 1 ) 
% persdia('Trial10x10/5TimeSteps/tmpPers_n10_d15_r5_T5_s1000_IC1_0_pers_0.txt')
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                        Create snapshot image of thresholded simulation.
%                        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
%             [[[ Used prospectus defense presentation. ]]]

% d = 0.15;
% thr = s*log(r);
% eth = 0.75;
% IC = 1;
% 
% N = load('../../BiomathTalk/rand1/rand_1.txt'); % MUST CHANGE when different.
% % N = importdata(matfile); 
% n = length(N);
% D = diag(diag(ones(n),1),1) + diag(diag(ones(n),-1),-1);
% 
% 
% % f = r*exp(-N/s);                         % Ricker growth per patch.
% % N = poissrnd(N.*f);                      % Stoch. Ricker growth all patches.
% % N = (1 - d)*N + (d/4)*(N*D + D*N);       % Det. dispersal.
% % N = floor(N);   
% 
% % N(N < eth) = 0;
% % N = 22*N.*exp(-N/s);
% % N = (1 - d)*N + (d/4)*(N * D) + (d/4)*(D * N);
% % 
% 
% 
% M = zeros(n,n,3);
% M(:,:,1) = 204/255;
% M(:,:,2) = 255/255;
% M(:,:,3) = 204/255;
% 
% N(N < eth) = 0;
% 
% % Find treats the matrix as a long vector (row then row).
% [I J] = find(N > thr); 
% 
% for k = 1:length(I)
%     M(I(k),J(k),1:3) = [102/255 102/255 51/255];
% end
% 
% im = image(M);
% 
% fn = strcat('../../BiomathTalk/rand1/shot70_0.png');
% saveas(im, fn)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                              Create bar3 plot to show matrix abundance.
%                              ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
%        [[[ Used in Prospectus defense (future work slide). ]]]
% 
% dispersal = 0.15;
% ethresh = 0.60;
% IC = 1;
% 
% matfile = sprintf('SmallGrid/rand_%g.txt',IC);
% %matfile = sprintf('ICrands/rand_%g.txt',IC);
% N = importdata(matfile); 
% n = length(N);
% D = diag(diag(ones(n),1),1) + diag(diag(ones(n),-1),-1);
% 
% N(N < ethresh) = 0;
% N = 22*N.*exp(-N);
% N = (1 - dispersal)*N + (dispersal/4)*(N * D) + (dispersal/4)*(D * N);
% 
% bar3(N)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



% =======================================================================
%                                 End. 
% =======================================================================