clear all; clc 
%%% Load data
local = [pwd '\..\Datasets\'];    
path = dir([local 'HUTUBS\pp*_HRIRs_simulated.sofa']);
Obj = SOFAload([path(1).folder, '\',path(1).name], 'nochecks');


%% 

% generate spatial sampling grid
g = Obj.SourcePosition;

% generate random data
rng(1)
data = rand(size(g,1), 1);

data = shiftdim(Obj.Data.IR(:,1,:),2).';


% set interpolation parameters
m          = [1 2 3];       % spline order m={1, 2, 3}
lambda     = [0 .01 0.02];  % smoothing factor
do_plot    = 1;             % 0: no plot, 1: planar plot, 2: spherical plot


% generate spatial sampling grid for interpolation
g_ref = AKgreatCircleGrid(-90:5:90, 5, 90);

%% interpolate with different spline orders and smoothing factors
%  (the colored dots in the plots show the interpolated data, the crosses
%  the position of the original data.)
AKf
for mm = 1:numel(m)
    for ll = 1:numel(lambda)
        
        subtightplot(3,3, (m(mm)-1)*3+ll, [0 0])
        
        data_interp = AKsphSplineInterp(g(:,1), g(:,2), data, g_ref(:,1), g_ref(:,2), m(mm), lambda(ll), 'deg', do_plot);
        
        set(gca, 'xTick', 0:90:270, 'yTick', -45:45:45)
        title(''); box on; grid on
        if ll ~= 1; ylabel ''; set(gca, 'yTickLabel', []); end
        if mm ~= 3 ; xlabel ''; set(gca, 'xTickLabel', []); end
        
        text(5, 80, ['m=' num2str(m(mm)) ', lambda=' num2str(lambda(ll))], 'backgroundColor', 'w')
        
    end
end
