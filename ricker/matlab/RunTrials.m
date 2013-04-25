% =======================================================================
%                                Trial runs.
% =======================================================================


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                   |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
%                   | Compute persistence and plot persistence diagrams.|
%                   |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

%      ~~~~~~~~ Parameters and initial population matrix. ~~~~~~~~~

% % % IC = 1;   ic = num2str(IC);              % Index of IC matrix.
% % % s = 1;    ss = num2str(s);               % Spatial rescale, to integers.
% % % r = 22;   rr = num2str(r);               % Fitness parameter.
% % % d = 0.00; dd = num2str(d*100);           % Dispersal parameter.
% % % T = 50;   TT = num2str(T);               % Final time step for Perseus. 
% % % 
% % % % path1 = strcat('../IC200rands/rand',ic); % Path to IC matrix file.
% % % % path1 = strcat('../../../Desktop/school/Data_2013March/BestParams'); 
% % % path1 = strcat('../../BiomathTalk/rand1'); 
% % % 
% % % % path2 = strcat(path1,'/r',rr,'_d',dd);       % Path to save files to.
% % % % path2 = strcat(path1,'/d',dd,'_r',rr,'_s',ss); 
% % % path2 = strcat(path1,'/Persistent'); 
% % % 
% % % Mat = strcat(path1,'/rand_',ic,'.txt');  % IC population matrix.
% % % % Mat = strcat(path1,'/IC_dip.txt');  
% % % 
% % % N = load(Mat);                           % Load in population matrix.
% % % n = length(N); nn = num2str(n);          % Dimension of population matrix.
% % % 
% % %                                          % To save persistence data.
% % % tmpPerseus = strcat(path2,'/Pers_n',nn,'_d',dd,'_r',rr,'_T',TT,'_s',ss,'_IC',ic);
% % % 
% % %                                    % To save abundances at each time step.
% % % tmpPopMat = strcat(path2,'/PopMat'); 
% % % 
% % % tmpPopMovie = strcat(path2);
% % % 
% % % 
% % % %           ~~~~~~~~~~~~~ Run functions! ~~~~~~~~~~~~
% % % 
% % % % process_perseus_sub_super(N, tmpPerseus, s) % RUN PERSEUS one time step.
% % % 
% % % reAP3_perseus_det(n,d,r,T,s,IC,N,tmpPerseus,tmpPopMat,tmpPopMat) % RUN PERSEUS many time steps.



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
% 
%  
% for k = 1:T           % Plot population abundances in 3D.
%     kk = num2str(k)
% 
%     tmpMat = strcat(tmpPopMat,'_',kk);
%     bar3(tmpMat)
% end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




% ~~~~~~~~~~~~~~~~ Population matrix movie simulation ~~~~~~~~~~~~~~~~~~~
% (Abundance Movies)

for d = 0.15:0.15:0.30
% for d =  0.50:0.25:0.75
    
%     for eth = 0.00:0.20:0.60
        
%         for s = 1:99:100
    
n = 21;
% d = 0.00;
r = 22;
T = 50;
s = 10000;
IC = 1;
eth = 0.60;


fig = figure;
name = sprintf('../../../Desktop/school/Data_2013_5-April/r%g_eth(s*%g)/MoviesAbundance/Mov_n%g_d%g_r%g_T%g_s%g_IC%g_abundance.avi',r,eth*100,n,d*100,r,T,s,IC);
aviobj = avifile(name,'compression','None');

path2 = sprintf('../../../Desktop/school/Data_2013_5-April/r%g_eth(s*%g)/s%g/d%g',r,eth*100,s,d*100);

for t = 0:T
    
    tt = num2str(t);
    
    MatName = strcat(path2,'/PopMat_', tt);
    PopMat = load(MatName);
    
    minPop = min(min(PopMat));
    maxPop = max(max(PopMat));
    
    im = imagesc(PopMat);

    colormap summer
    colormap(flipud(colormap))
    colorbar
    
    caxis([0 2])

    rows = size(PopMat,1);
    cols = size(PopMat,2);
    
    for i = 1:cols
        matdata = ones(6*rows, 4);
        k = 1;

        for j = 0:6:(6*rows - 6)
            matdata(j + 1:j + 6, :) = PopMat(k,i);
            k = k + 1;
        end
%         set(im(i),'Cdata',matdata)
    end
    
    if(mod(t,2) == 0) % To account for the period 2 flipping (color flips).
        
        aviobj = addframe(aviobj, fig);
    end
end
    
aviobj = close(aviobj);

%         end
%     end
end


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% PlotName = strcat(path2,'/PopMat_50');
% M = load(PlotName);
% figure(1);
% imagesc(M)
% figure(2);
% bar3(M)






% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% for i = 0:2
%     
%     step = num2str(i);
%     
%     dim = num2str(1);
%     
%     fname1 = strcat('../IC200rands/rand1/OverTime/Pers_n',nn,'_d',dd,'_r',rr,'_T',TT,'_s',ss,'_IC',ic,'_',step,'_full_',dim,'.pdia' );
%     persdia_sub_super( fname1 )
%     
%     fn = strcat('../IC200rands/rand1/OverTime/Pers_n',nn,'_d',dd,'_r',rr,'_T',TT,'_s',ss,'_IC',ic,'_',step,'_full_',dim,'.png' );
% 
%     
%     saveas(gcf, fn)
%     
%     close
%     
% end

% 
% fname2 = strcat('../ResearchBlogs/StochGrowth_1D/n',nn,'_r',rr,'_t',tt,'_a',aa,'.png');
% fn = sprintf(fname2);
% saveas(gcf, fname2)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% =======================================================================
%                                   End.
% =======================================================================