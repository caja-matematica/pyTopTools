function n = persdia_sub_super(varargin)
% n = number of persistence intervals
% inputs: mandatory filename = name of file containing 
% birth and death info in column form, and a number plot-type.
% if plot type = 0, we plot ALL intervals, otherwise only the
% ones that actually die.

if nargin==0
    n = -1;
    return
else
    filename = varargin{1};
end    

plot_type = 0;
if nargin == 2 
    plot_type = varargin{2};
end

% set a flag in case the only death is infinite; default is false
maxd_is_inf = 0; 
normal_gens = 0;

% extract birth and death indices
ints = load(filename);
births = ints(:,1)
deaths = ints(:,2)

% extract indices of those intervals which persist throughout
sub_ints = find(deaths == -1);
super_ints = find(deaths == -2);

% extract indices of those intervals which die
try
    normal_ints = find(deaths ~= -1 & deaths ~= -2);
    normal_gens = 1;
catch
    disp( 'Only infinite generators on this homology level.' );
end

% extract maximum death time encountered and min birth times
minb = min( births ) % should always be >= 1
maxbirth = max( births )
maxd = max( deaths )
if (maxd < 0) 
    maxd_is_inf = 1;
    maxd = max( births )
end

% start figure
figure;

% we always plot these
plot(births(normal_ints),deaths(normal_ints),'bo','MarkerFaceColor','b', ...
     'MarkerSize', 2 );
hold on;

max_axis = max( [ maxd, maxbirth ] )
axis([-0.1*maxd,(1.1)*(max_axis+1),-0.1*maxd, (1.1)*(max_axis+1)]);

% plot the diagonal
%diag = -0.1*maxd : (1.1)*maxd;

title(filename);
xlabel('birth');
ylabel('death');

%% set super-level infinite gens to min - 1 (so if min val in
%% matrix is 2, set inf to 1 )
infsub_vec = (max_axis + 1) * ones(size(sub_ints));
%infsuper_vec = (minb -1 ) * ones(size(super_ints));
infsuper_vec = zeros(size(super_ints));
diag = 0:(1.1)*(max_axis+2)/20:(1.1)*(max_axis+2);
%diag = -0.1*maxd : (1.1)*(maxd+1) 
npts = length( diag );
 
% plot diag and horizontal lines
plot(diag,diag,'g-');
plot( zeros( npts ), diag, 'k--' ); % horizontal line
plot( diag, zeros( npts ), 'k--' ); % vertical line

% plot infinite stuff
plot((births(sub_ints)), infsub_vec , 'rd', 'LineWidth', 1 ); %,'MarkerFaceColor','r' );
plot((births(super_ints)), infsuper_vec , 'rd'); %,'MarkerFaceColor','r');    
    
n = size(births);

hold off;
