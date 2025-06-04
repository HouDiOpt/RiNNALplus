function perf(name,T,logplot)
%PERF    Performace profiles
%
% PERF(T,logplot)-- produces a performace profile as described in
%   Benchmarking optimization software with performance profiles,
%   E.D. Dolan and J.J. More', 
%   Mathematical Programming, 91 (2002), 201--213.
% Each column of the matrix T defines the performance data for a solver.
% Failures on a given problem are represented by a NaN.
% The optional argument logplot is used to produce a 
% log (base 2) performance plot.
%
% This function is based on the perl script of Liz Dolan.
%
% Jorge J. More', June 2004

close all

if (nargin < 2) logplot = 0; end

colors  = ['k' 'b' 'r' 'g' 'c' 'm' 'y'];
lines   = [':' '-' '-.' '--'];
markers = ['v' 's' 'o' '*' 'v' '^' 'o'];
method_names = {'SDP-RLT', 'DNN', 'COMP','SDPNAL+'};
line_width = 2.0;



[np,ns] = size(T);

% Minimal performance per solver

minperf = min(T,[],2);

% Compute ratios and divide by smallest element in each row.

r = zeros(np,ns);
for p = 1: np
  r(p,:) = T(p,:)/minperf(p);
end

if (logplot) r = log2(r); end

max_ratio = max(max(r));

% Replace all NaN's with twice the max_ratio and sort.

r(find(isnan(r))) = 2*max_ratio;
r = sort(r);

% Plot stair graphs with markers.

clf;
for s = 1: ns
 [xs,ys] = stairs(r(:,s),[1:np]/np);
 option = ['-' colors(s) markers(s) ];
 h1=plot(xs,ys,option,'LineWidth', line_width,'DisplayName', method_names{s},'MarkerSize',4);
 set(h1, 'markerfacecolor', get(h1, 'color')); % Use same color to fill in markers
 hold on;
end

% Axis properties are set so that failures are not shown,
% but with the max_ratio data points shown. This highlights
% the "flatline" effect.

axis([ 0 1.1*max_ratio 0 1 ]);

% Legends and title should be added.
xlabel('lg(\tau)')
ylabel('fraction of solver within \tau of the best')
title(['number of tested problems: ',num2str(np)]);

% Create legend and set font size
lgd = legend('Location', 'northeast');  % Start with the northeast location
lgd.Position = lgd.Position + [-0.02, -0.05, 0, 0];  % Move it slightly downward

% lgd = legend('Location', 'southeast');  % Move legend to the bottom right
% lgd = legend('Location', 'northeastoutside');  % Or 'southoutside' for bottom placement
set(lgd, 'FontSize', 14);  % Increase legend font size


% Set figure width and height
% fig = gcf;
% fig.Position = [100, 100, 500, 400];

print(gcf, strcat([name,'.eps']), '-depsc');
