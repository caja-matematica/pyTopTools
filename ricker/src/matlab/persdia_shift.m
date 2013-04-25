function n = persdia(varargin)
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


% extract birth and death indices
ints = load(filename);
births = ints(:,1);
deaths = ints(:,2);

infdeaths = find(deaths == -1);

minbirth = min(births);

b2 = births - minbirth;
d2 = deaths - minbirth;

d2(infdeaths) = -1;

% extract maximum death time encountered:
maxd = max(deaths);
maxd2 = max(d2);
if (maxd < 0) 
    maxd = max(births) + 1;
end

% extract indices of those intervals which die
normal_ints = find(d2 ~= -1);

figure;

% we always plot these:
plot(b2(normal_ints),d2(normal_ints),'bo','MarkerFaceColor','b');
hold on;
axis([0,maxd2 + 2,0,maxd2 + 2]);

% plot the diagonal
diag = 0:(maxd2 + 2)/20:maxd2 + 2;
plot(diag,diag,'g-');

title(filename);
xlabel('birth shifted');
ylabel('death shifted');

% extract indices of those intervals which persist throughout
if plot_type == 0 
  inf_ints = find(d2 == -1);
  infd_vec = (maxd2)*ones(size(inf_ints));
  plot((b2(inf_ints)), infd_vec , 'rd','MarkerFaceColor','r');
end
    
n = size(births);

hold off;
