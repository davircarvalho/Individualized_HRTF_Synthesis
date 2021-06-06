function Obj = SOFA_barycentric_interp(Obj, out_pos, varargin)
%HRTF_INTERPOLATION_IN_3D Matlab demonstration of an algorithm for 
% interpolating distance-dependent (near-field) head-related transfer 
% functions (HRTFs) in azimuth, elevation, and distance
% 
%       Inputs: 
%               Obj: SOFA SimpleFreeFieldHRIR object
%               out_pos: desired output positions.
%       Outputs: 
%               Obj: SOFA SimpleFreeFieldHRIR object with interpolated
%                    HRIRs for the desired position
% 
% The interpolation is performed using barycentric weights of 4 HRTF 
% measurements forming a tetrahedron that encloses the desired source 
% position.
%
% The main steps of the interpolation algorithm are:
% 1) Initialisation: organise HRTF measurement positions in tetrahedral 
% mesh via Delaunay triangulation and generate an octree;
%
% 2) Source update: find a tetrahedron that encloses the desired source 
% position (via brute-force search, or via adjacency walk with/without 
% octree lookup);
%
% 3) Interpolation: calculate barycentric weights for linear interpolation 
% of HRTF measurements at the vertices of the tetrahedron selected in 2)
%
% Reference:
% Gamper, H. (2013). "Head-related transfer function interpolation in azimuth,
% elevation, and distance", Journal of the Acoustical Society of America
% 134(6), JASA EL547-EL554.
% http://link.aip.org/link/?JAS/134/EL547/1
%
%
% Created by Hannes Gamper.
% 1.0 - Oct 2013 Initial release
% Adapted by Davi R. Carvalho - JUN/2021
% 
% Bug fixes:
% - Jan 2015: removed tilde placeholders in function output calls for 
%   backward compatibility with Matlab versions earlier than 2009b
% - Jan 2015: fixed call to Delaunay to support earlier Matlab versions
%
% Please post comments to the FEX page for this entry if you have any
% bugs or feature requests:
% http://www.mathworks.com/matlabcentral/fileexchange/43809

% Copyright 2015 Hannes Gamper
% Flags
USE_ADJACENCY_WALK = true;  % true:  use adjacency walk
                            %        (quickly find the best tetrahedron
                            %         enclosing the desired position)
                            % false: use brute-force search
                            
USE_OCTREE = false;          % true:  use Octree query to find starting 
                            %        tetrahedron for adjacency walk 
                            %        (speedup initial guess for adjacency walk)
                            % false: use random starting tetrahedron for
                            %        adjacency walk

% constants
binCapacity = 5;            % maximum number of points stored in each leaf 
                            % node of the octree

%% LOAD SOFA
inpt_pos = Obj.SourcePosition;
inpt_pos(:,3) = inpt_pos(:,3)*100; % source position (azimuth [deg], elevation [deg], radius [cm])
out_pos(:,3) = inpt_pos(1,3);% now it's expected you use a single radius HRTF set (you can exclude this, but make sure you have the distance in [cm])

fs = Obj.Data.SamplingRate;    % sampling rate [Hz]

% Generate tetrahedral mesh via Delaunay triangulation
if str2double(datestr(datenum(version('-date')),'YYYY'))>=2013
    dt = delaunayTriangulation(inpt_pos);
    T = dt.ConnectivityList;    % contains the tetrahedral mesh
    X = dt.Points;              % contains the vertices of the tetrahedra
    N = neighbors(dt);          % contains the adjacency list
else
    % for older Matlab releases < 2013
    T = delaunay(inpt_pos(:,1), inpt_pos(:,2), inpt_pos(:,3));
    X = inpt_pos;
    N = dtNeighbours(T);
end


for k_pos = 1:size(out_pos,1)
    posi = out_pos(k_pos,:);
    if USE_OCTREE
        if ~USE_ADJACENCY_WALK
            warning('Octree search can only be used in conjunction with adjacency walk, otherwise the selection might not terminate...');
            USE_ADJACENCY_WALK = true;
        end

        % Generate octree
        OT = getOT(inpt_pos, binCapacity);

        % get centres and one tetrahedron per octree leaf node
        OT = getOTcentres(OT, T);

        % find close tetrahedron
        ti = queryOT(posi, OT);
    else
        ti = 1;
    end

    % variable initialisations
    Niter = 0;
    tetra_indices = zeros(size(T,1),1);
    MIN_GAIN = -0.00001;        % due to numerical issues, very small negative 
                                % gains are considered 0

    %%%% MAIN LOOP %%%%----------------------------------------------------
    % iterate through mesh to find tetrahedron for interpolation
    for t = 1:size(T,1)
        % count iterations
        Niter = Niter+1;

        % vertices of tetrahedron
        HM = X(T(ti,:),:);
        v4 = HM(4,:);
        H = HM(1:3,:) - repmat(v4,3,1);
        tetra_indices(t) = ti;

        % calculate barycentric coordinates
        bary_gains = [(posi-v4)*(pinv(H)), 0];
        bary_gains(4) = 1-sum(bary_gains);

        % check barycentric coordinates
        if all(bary_gains>=MIN_GAIN)
            bary_gains = max(bary_gains, 0);
            tetra_indices = tetra_indices(1:Niter);
            break;
        end

        if USE_ADJACENCY_WALK
            % move to adjacent tetrahedron
            [tmp, bi] = min(bary_gains);
            ti = N(ti,bi);
        else
            % BRUTE-FORCE: move to next tetrahedron in list
            ti = ti+1;
        end
    end

    if tetra_indices(end)==0
        error('No tetrahedron found. Exiting...');
    end


    %%% Interpolation -----------------------------------------------------
    % get hrtfs
    idx_tetra_pos = dsearchn(inpt_pos, HM);
    hrir = Obj.Data.IR(idx_tetra_pos,:,:);
    hrir = shiftdim(hrir,2);
    HRTF = fft(hrir);
    NFFT = size(HRTF,1);
    for k =1:size(hrir,3)
        Hint(k_pos,:,k) = bary_gains*abs(HRTF(1:NFFT/2,:,k)');
    end
end
disp('Barycentric interpolation complete!');

%% Output SOFA Obj
Obj.Data.IR = shiftdim(Hint,1);
Obj.SourcePosition = out_pos;
Obj = SOFAupgradeConventions(Obj);
% Obj = SOFAupdateDimensions(Obj);


%% PLOTS %%%%

if ~isempty(varargin) && any(strcmp(varargin, 'plot'))
    figure('units','normalized','outerposition',[0 0 1 1])
    cols = [1,0,0;0,0.5,0;0,0,1;1,0,1;0,1,1;0,0,0];
    subplot(2,2,[1,3]);
    tetramesh(T(tetra_indices,:), X, 'FaceAlpha', 0.1);
    hold on;
    plot3(posi(1),posi(2),posi(3), 'rx', 'LineWidth', 3, 'MarkerSize', 20);
    grid on;
    xlabel('Azimuth [deg]'); ylabel('Elevation [deg]'); zlabel('Distance [cm]');
    title(['Tetrahedral search: ', num2str(Niter), ' iterations']);

    subplot(2,2,2);
    tetramesh(T(tetra_indices(end),:), X, 'FaceAlpha', 0.1);
    hold on;grid on;
    plot3(posi(1),posi(2),posi(3), 'rx', 'LineWidth', 3, 'MarkerSize', 20);
    ht = zeros(1,4);
    legstr = cell(1,4);
    for vi = 1:4
        ht(vi) = text(HM(vi,1), HM(vi,2), HM(vi,3), ['g', num2str(vi)], 'FontSize', 17, 'Color', cols(vi,:));
        legstr{vi} = ['g', num2str(vi), ' = ', num2str(bary_gains(vi))];
    end
    xlabel('Azimuth [deg]'); ylabel('Elevation [deg]'); zlabel('Distance [cm]');
    title('Barycentric interpolation weights');


    %%% hrtf plot
    subplot(2,2,4);
    cols = [1,0,0;0,0.5,0;0,0,1;1,0,1;0,1,1;0,0,0];
    fft_vect = linspace(0,fs/2,NFFT/2);
    for hi = 1:4
        plot(fft_vect, 20*log10(abs(HRTF(1:NFFT/2,hi,1))), 'Color', cols(hi,:));
        hold on;
    end
    grid on;
    xlim([0,fs/2]);
    hold on;
    plot(fft_vect, 20*log10(squeeze(Hint(end,:,1))), 'k', 'LineWidth', 2);
    legend('h1', 'h2', 'h3', 'h4', 'h interpolated', 'location', 'best');
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [dB]');
    title('Interpolation of (modelled) HRTF data');
end




end





















%% get OCTREE
% Organise 3-D points pts into an octree data structure
% 
% binCapacity: maximum number of points per leaf node
%
% References:
% Samet, H. (1989). "Implementing ray tracing with octrees and neighbor 
% finding", Computers & Graphics 13, pp. 445?460.
function OT = getOT(pts, binCapacity)
% inits
OT = struct;

% start with root
OT.Points = pts;
OT.PointBins = ones(size(pts,1),1);
OT.BinBoundaries = [min(pts), max(pts)];
OT.BinParents = 0;
OT.BinDepths = 0;

% MAIN LOOP
OT = getChildren(OT,1,binCapacity);
OT.BinCount = size(OT.BinParents,1);

end

% Recursive call to split bin into octants
function OT = getChildren(OT, parent, binCapacity)
parentBoundaries = OT.BinBoundaries(parent,:);
for bini = 1:8
    binBoundaries = zeros(1, 6);
    rng = range(reshape(parentBoundaries,3,2),2)';
    binBoundaries(1:3) = bitget(bini-1,1:3).*rng/2 + parentBoundaries(1:3);
    binBoundaries(4:6) = binBoundaries(1:3) + rng/2;
    
    % add bin to list of bins
    OT.BinBoundaries = [OT.BinBoundaries; binBoundaries];
    OT.BinParents = [OT.BinParents; parent];
    OT.BinDepths = [OT.BinDepths; OT.BinDepths(parent)+1];
    currBin = size(OT.BinParents,1);
    
    % check if bin contains any points
    ptinds = (1:size(OT.Points,1))';
    ptinds = ptinds(OT.PointBins==parent);
    
    binPoints = inBin(OT.Points(ptinds,:), binBoundaries);
    binPoints = ptinds(binPoints);
    if ~isempty(binPoints)
        OT.PointBins(binPoints) = currBin;
        if length(binPoints)>binCapacity
            % keep dividing
            OT = getChildren(OT, currBin, binCapacity);
        end
    end
end
end

% Determine which points pts are contained in bin
function inds = inBin(pts, binBoundaries)
N = size(pts,1);
lpts = pts-repmat(binBoundaries(1:3),N,1);
upts = pts-repmat(binBoundaries(4:6),N,1);
inds = find(sum(lpts>=0 & upts<=0,2)==3);
end

% get OT bin centres
function OT = getOTcentres(OT, T)
NPoints = size(OT.Points,1);
NBins = OT.BinCount;
OT.BinCentres = zeros(NBins,3);
OT.BinChildren = zeros(NBins,8);
for bi = 1:NBins
    bin_children = find(OT.BinParents==bi);
    if isempty(bin_children)
        % leaf node: insert a random point as a "child"
        pt = find(OT.PointBins==bi,1);
        if isempty(pt)
            % bin contains no point -> find closest point to bin centre
            binb = OT.BinBoundaries(bi,:);
            binc = (binb(1:3)+binb(4:6))./2;
            binn = sqrt(sum( (OT.Points-repmat(binc, NPoints,1)).^2,2 ));
            [tmp, pt] = min(binn);
        end
        
        % find a tetrahedron containing pt
        Tind = find(T==pt,1);
        Tind = bmod(Tind, size(T,1));
        
        OT.BinChildren(bi,2) = Tind;
        continue;
    end
    boundaries = OT.BinBoundaries(bin_children,:);
    
    % find centre
    x = unique(boundaries(:,[1,4]));
    y = unique(boundaries(:,[2,5]));
    z = unique(boundaries(:,[3,6]));
    OT.BinCentres(bi,:) = [x(2),y(2),z(2)];
    
    % store children in order
    for bch = 1:8
        chb = OT.BinBoundaries(bin_children(bch),:);
        chc = (chb(1:3)+chb(4:6))./2;
        chind = 0;
        for j = 1:3
            if chc(j) >= OT.BinCentres(bi,j)
                chind = chind + 2^(j-1);
            end
        end
        OT.BinChildren(bi,chind+1) = bin_children(bch);
    end
end
end

% query octree
function ti = queryOT(pos, OT)
binind = 1;
maxbiniter = 50;
biniter = 0;
while (1)
    if (biniter > maxbiniter)
        ti = 1;
        fprintf('No bin found...\n');
        break;
    end
    
    biniter = biniter + 1;
    currbin = binind;
    octind = 0;
    binc = OT.BinCentres(binind,:);
    
    if (pos(1)>=binc(1)); octind = bitor(octind,1); end
    if (pos(2)>=binc(2)); octind = bitor(octind,2); end
    if (pos(3)>=binc(3)); octind = bitor(octind,4); end
    
    has_children = OT.BinChildren(binind,1)>0;
    
    if (~has_children)
        ti = OT.BinChildren(currbin,2);
        break;
    end
    
    binind = OT.BinChildren(binind,octind+1);
end
end

%% Get adjacency list for delaunay triangulation
% T: connectivity list (tetrahedral mesh) of delaunay triangulation 
function N = dtNeighbours(T)

% create matrix containing all edges/faces
NT = size(T,1);
Tdim = size(T,2);
face = zeros(Tdim*NT,Tdim-1);
vinds = zeros(Tdim*NT,1);
Tinds = zeros(Tdim*NT,1);
ptr1 = 1;
ptr2 = NT;

for fi = 1:Tdim
    inds = 1:Tdim;
    inds(fi) = [];
    face(ptr1:ptr2,:) = T(:,inds);
    face(ptr1:ptr2,:) = sort(face(ptr1:ptr2,:), 2);
    vinds(ptr1:ptr2) = fi;
    Tinds(ptr1:ptr2) = 1:NT;
    ptr1 = ptr1+NT;
    ptr2 = ptr2+NT;
end
[faceS, indS] = sortrows(face);
vindsS = vinds(indS);
TindsS = Tinds(indS);

N = -ones(size(T));

n = 0;
while n<size(faceS,1)-1
    n = n+1;
    f1 = faceS(n,:);
    f2 = faceS(n+1,:);
    if sum(f1==f2)==(Tdim-1)
        % found matching edge/face!
        t1 = TindsS(n);
        t2 = TindsS(n+1);
        vi1 = vindsS(n);
        vi2 = vindsS(n+1);
        N(t1,vi1) = t2;
        N(t2,vi2) = t1;
        n = n+1;
    end
end

end

%% alternative modulus function
% mod(N,N) returns N (instead of 0)
% useful for indexing in Matlab
function y = bmod(x,a)
y = mod(x-1,a)+1;
end