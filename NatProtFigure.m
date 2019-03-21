% Figure 2 drawing for Nature Protocols paper.
% Want to demonstrate spatial analysis methods on a variety of SMLM
% methods.
% Random, clustered, and maybe ordered (dispersed?) data sets
% - Nearest-neighbor analysis
% - Ripley's K
% - Pair correlation
% - DBSCAN
% - Voronoi/triangulation
%
% Requires ELKI package, available at:
% https://elki-project.github.io/releases/

% Make a dataset of random points + Gaussian clusters
FieldSize = 5000; % in nm
RandFieldPoints = 1e3;
NPeaks = 40;
PeakSigma = [50 80]; % nm
PtsPerPeak = [20 40]; 

% Generate data or load?
dataSource = 'generate'; % 'generate' or 'load'

% Nearest neighbor analysisK = RipleysK(pts,RK.histRange,[minX maxX minY maxY],method)
NN.FieldSamplings = 5; % Times to generate random data
NN.histRange = 0:4:250; % Histogram range, in nm
NN.resampleNum = 100; % Times to resample for random data

% Ripley's K analysis
RK.FieldSizeExpand = 2;
RK.histRange = 0:5:1000;

% Pair correlation
PCF.Resolution_corr = 10; % Resolution of correlated image, in nm.
PCF.rmax = 300; % Maximum distance of correlation to create

% DBSCAN
DB.FieldShrink = 2; % Factor by which to shrink FOV
DB.Epsilon = 100; % Search radius (nm)
DB.minPts = 10; % minimum points/cluster
% DB.BinPixel = 50; % Pixel size of points-binned image
DB.plotPixel = 10;

% Voronoi + Delaunay
VD.FieldShrink = 2; % Factor by which to shrink FOV
VD.AlphaParameter = 20;

% How-they-work images
Demo.FieldShrink = 5;

% Save flag
saveImagesToDisk = false;

%% Generate some initial points

switch dataSource
    
    case 'generate'
    
        Masterpts = GenerateNatProtData(FieldSize, RandFieldPoints, PeakSigma, NPeaks, PtsPerPeak);

        % Points field image
        fig = figure(1);
        clf(fig);
        ax = axes('parent', fig);

        plot(ax, Masterpts(:,1), Masterpts(:,2), '.', 'color', [.5 .5 .5]);
        axis image
        set(ax, 'xlim', [0 FieldSize], 'ylim', [0 FieldSize], 'xtick', [],'ytick', []);
        % xlabel('X Position (nm)'); ylabel('Y Position (nm)');

        patch( [100 600 600 100 100], [100 100 200 200 100], 'k', 'parent', ax);
        set(fig, 'color', [1 1 1]);
        set(ax, 'box', 'on');
        
    case 'load'

        %% -OR- Load dataset from disk

        Masterpts = dlmread('D:\MATLAB\NatProtSMLMAnalysis\MasterPointsNatureProtocol.txt');

        % Points field image
        fig = figure(1);
        clf(fig);
        ax = axes('parent', fig);

        plot(ax, Masterpts(:,1), Masterpts(:,2), '.', 'color', [.5 .5 .5]);
        axis image
        set(ax, 'xlim', [0 FieldSize], 'ylim', [0 FieldSize], 'xtick', [],'ytick', []);
        % xlabel('X Position (nm)'); ylabel('Y Position (nm)');

        patch( [100 600 600 100 100], [100 100 200 200 100], 'k', 'parent', ax);
        set(fig, 'color', [1 1 1]);
        set(ax, 'box', 'on');
        
end


%% Nearest-neighbor analysis

NNhistData = zeros(numel(NN.histRange), NN.FieldSamplings);

for k = 1:NN.FieldSamplings

    pts = GenerateNatProtData(FieldSize, RandFieldPoints, PeakSigma, NPeaks, PtsPerPeak);

    [~, dG] = knnsearch(pts, pts, 'K', 2);
    % For a auto-nearest-neighbor calc, the nearest point is always itself.
    % This finds the two nearest points.  Throw away the first and get the
    % auto-nearest-neighbor without the point itself at 0 distance included.
    dG(:,1) = [];

    % Make a histogram for plotting
    NNhistData(:,k) = histc(dG, NN.histRange);
    
end

NNhist = mean(NNhistData, 2);
NNhist = NNhist/sum(NNhist(:));
NNhistCI = prctile(NNhistData, [2.5 97.5], 2);
NNhistCI(:,1) = NNhistCI(:,1)/sum(NNhistCI(:,1));
NNhistCI(:,2) = NNhistCI(:,2)/sum(NNhistCI(:,2));

% Resample resampleNum number of times to get what 'random' field would
% generate.
% Field here is constrained to rectangular, so generating random spots 
% within this FOV is straightforward.  For polygonal ROIs a triangulation 
% followed by a re-sampling within the triangles (with n points/triangle 
% weighted by triangle area) is a convenient route.
% Also assuming that all data ROIs are identical in size, shape.
% Individual re-sampling can be done and results combined, if needed, 
% but is not necessary here.

NNhistResample = zeros(numel(NN.histRange), NN.resampleNum);

for k = 1:NN.resampleNum
    
    resamplePts = FieldSize*rand(size(pts, 1), 2);
    [~, dGresample] = knnsearch(resamplePts, resamplePts, 'K', 2);
    dGresample(:,1) = [];
    NNhistResample(:,k) = histc(dGresample, NN.histRange);
    
end

NNresampleMean = mean(NNhistResample, 2);
NNresampleMean = NNresampleMean/sum(NNresampleMean(:));
NNresampleCI = prctile(NNhistResample, [2.5 97.5], 2);
NNresampleCI(:,1) = NNresampleCI(:,1)/sum(NNresampleCI(:,1));
NNresampleCI(:,2) = NNresampleCI(:,2)/sum(NNresampleCI(:,2));


% NNdist plot
NNfig = figure(2);
clf(NNfig);
NNax = axes('parent', NNfig);


set(NNax, 'NextPlot', 'add');
% plot(NNax, histRange, NNhistCI(:,1), 'r--');
% plot(NNax, histRange, NNhistCI(:,2), 'r--');
plot(NNax, NN.histRange, NNresampleMean, 'color', rgb(155, 89, 182), 'linewidth', 2);
plot(NNax, NN.histRange, NNresampleCI(:,1), '--', 'color', rgb(155, 89, 182));
plot(NNax, NN.histRange, NNresampleCI(:,2), '--', 'color', rgb(155, 89, 182));

plot(NNax, NN.histRange, NNhist, 'color', rgb(243, 156, 18), 'linewidth', 2);
set(NNax, 'NextPlot', 'replace');
xlabel('Nearest-Neighbor Distance (nm)'); ylabel('PDF');
set(NNax, 'fontsize', 12);
set(NNfig, 'color', [1 1 1]);

%%  Ripley-K

pts = GenerateNatProtData(FieldSize*RK.FieldSizeExpand, RandFieldPoints*RK.FieldSizeExpand, PeakSigma, NPeaks*RK.FieldSizeExpand, PtsPerPeak);

K = RipleysK(pts, RK.histRange, [FieldSize/2 3*FieldSize/2 FieldSize/2 3*FieldSize/2], 1);
Lr_r = 2*sqrt(K/pi) - RK.histRange'; % Factor of 2 is unexpected.

rkFig = figure(3);
clf(rkFig);
rkAx = axes('parent', rkFig);
plot(rkAx, RK.histRange, Lr_r, 'color', rgb(52, 152, 219), 'linewidth', 2);
xlabel(rkAx, 'Distance (nm)'); ylabel(rkAx, 'L(r) - r');
set(rkFig, 'color', [1 1 1]);
set(rkAx, 'Fontsize', 12, 'ytick', 0:50:200);

%% Pair correlation function
% Using the Veatch method, the sparse data set has to be converted to an
% image first. This doesn't *have* to be done, but to keep consistent with
% published methods, this approach is used here.

pts = GenerateNatProtData(FieldSize, RandFieldPoints, PeakSigma, NPeaks, PtsPerPeak);

Edges_x = 0:PCF.Resolution_corr:(FieldSize); 
Edges_y = 0:PCF.Resolution_corr:(FieldSize);

histmat = hist2(pts(:,1), pts(:,2), Edges_x, Edges_y);

Bin_img = (histmat((1:(numel(Edges_x)-1)), (1:(numel(Edges_y)-1))));

[G, r, g, dg, mask] = get_autocorr(Bin_img, ones(size(Bin_img)), PCF.rmax);

[pcFits] = lsqcurvefit(@GaussianPairCorr, [1 1], r(2:end), g(2:end));

detectSize = pcFits(1)*PCF.Resolution_corr;
pcDomain = r(2:end)*PCF.Resolution_corr;
pcRange = g(2:end);

pcFig = figure(4);
clf(pcFig);
pcAx = axes('parent', pcFig);
semilogx(pcAx, pcDomain, pcRange, 'o', 'markersize', 8, 'color', rgb(44, 62, 80));
set(pcAx, 'NextPlot', 'add');
plot(pcAx, linspace(0, PCF.rmax, 500)*PCF.Resolution_corr, GaussianPairCorr(pcFits, linspace(0, PCF.rmax, 500)), 'linewidth', 2, 'color', rgb(231, 76, 60));
set(pcAx, 'NextPlot', 'replace');
xlabel(pcAx, 'Distance (nm)'); ylabel(pcAx, 'PCF');
set(pcAx, 'xlim', [10 2500], 'fontsize', 12, 'ytick', 0.5:1:5.5);
set(pcFig, 'color', [1 1 1]);

%% DBSCAN!
% Going via open-source ELKI package.  Be sure to change path to ELKI
% install in dependent function dbscanViaElki.m so that code can execute.


pts = Masterpts; % take points in range 0 -> round(FieldSize/DB.FieldShrink) in middle of image
ptsRange = FieldSize/2 + [-FieldSize/(2*DB.FieldShrink) FieldSize/(2*DB.FieldShrink)];
pts((pts(:,1) > ptsRange(2)) | (pts(:,1) < ptsRange(1)), :) = [];
pts((pts(:,2) > ptsRange(2)) | (pts(:,2) < ptsRange(1)), :) = [];

clusterIDs = dbscanViaELKI(pts, DB.Epsilon, DB.minPts, 'XY');

dbFig = figure(5);
clf(dbFig);
dbax = axes('parent', dbFig);
set(dbax, 'NextPlot', 'add');

clustList = unique(clusterIDs);
clustColors = [rgb(26, 188, 156); rgb(52, 152, 219); rgb(230, 126, 34); ...
    rgb(243, 156, 18); rgb(192, 57, 43); rgb(41, 128, 185); rgb(22, 160, 133); ...
    rgb(231, 76, 60); rgb(46, 204, 113); rgb(241, 196, 15); rgb(142, 68, 173); rgb(155, 89, 182)];

% Find 'noise' cluster and set that to gray
% Picking this one out as the one that has the greatest variance in
% position along the x axis.  

clustXVar = zeros(numel(clustList), 1);
for k = 1:numel(clustList)
    clustXVar(k) = var(pts(clusterIDs == clustList(k), 1));
end
clustColors(clustXVar == max(clustXVar(:)), :) = [0.4 0.4 0.4];

convHullDia = zeros(numel(clustColors), 1);
gaussImage = zeros(numel(0:DB.plotPixel:round(FieldSize/DB.FieldShrink)), numel(0:DB.plotPixel:round(FieldSize/DB.FieldShrink)), 3);
for k = 1:numel(clustList)
   
    pointsList = pts(clusterIDs == clustList(k),:);
    plot(pointsList(:,1), pointsList(:,2), '.', 'color', clustColors(k,:));
    
    % Segment clusters with convex hull, then determine equivalent
    % diameter'
    if clustXVar(k) ~= max(clustXVar(:))
        
        [hull, vol] = convhull(pointsList(:,1), pointsList(:,2));
        
        plot(pointsList(hull, 1), pointsList(hull, 2), 'color', clustColors(k,:), 'linewidth', 2);
        convHullDia(k) = 2*sqrt(vol/pi);
   
    

%     
%     fitDomX = min(pointsList(:,1)):DB.BinPixel:max(pointsList(:,1));
%     fitDomY = min(pointsList(:,2)):DB.BinPixel:max(pointsList(:,2));
%     if numel(fitDomX) ~= numel(fitDomY)
%         while numel(fitDomX) < numel(fitDomY)
%             fitDomX = [fitDomX, max(fitDomX)+DB.BinPixel];
%         end
%         while numel(fitDomX) > numel(fitDomY)
%             fitDomY = [fitDomY, max(fitDomY)+DB.BinPixel];
%         end
%     end
%     
%         % Fit each set of points to a 2D gaussian?
%     subGauss = hist2(pointsList(:,1), pointsList(:,2), ...
%         fitDomX, fitDomY);
%     
%     [gfit] = lsqcurvefit(@gauss_2D, [max(subGauss(:)) mean(pointsList(:,1)) mean(pointsList(:,2)) 50 50 min(subGauss(:))], ...
%         [fitDomX; fitDomY]', ...
%         subGauss);
%     
%     % Add gaussian to the plot
%     
%     invColor = rgb2hsv(clustColors(k,:));
%     gImg = gauss_2D(gfit, [0:10:round(FieldSize/DB.FieldShrink); 0:10:round(FieldSize/DB.FieldShrink)]');
%     gI = cat(3, ones(size(gImg, 1), size(gImg, 2))*(1-invColor(:,1)), ones(size(gImg, 1), size(gImg, 2)), ...
%         gImg);
%     gaussImage = gaussImage + hsv2rgb(gI./max(gI(:)));
%     
%     gaussImage = 1-gaussImage;
    
    
    
    end
    
end

set(dbax, 'xlim', ptsRange, 'ylim', ptsRange);
set(dbax, 'xtick', [], 'ytick', [], 'box', 'on');
set(dbFig, 'color', [1 1 1]);

axis square

% Scale bar
patch( [100 600 600 100 100] + ptsRange(1), [100 100 150 150 100] + ptsRange(1), 'k', 'parent', dbax); 


%% Voronoi + Delaunay

pts = Masterpts; % take points in range 0 -> round(FieldSize/DB.FieldShrink) in middle of image
ptsRange = FieldSize/2 + [-FieldSize/(2*VD.FieldShrink) FieldSize/(2*VD.FieldShrink)];
pts((pts(:,1) > ptsRange(2)) | (pts(:,1) < ptsRange(1)), :) = [];
pts((pts(:,2) > ptsRange(2)) | (pts(:,2) < ptsRange(1)), :) = [];

% Triangulate and Voronoi
dT = delaunayTriangulation(pts(:,1), pts(:,2));

shp = alphaShape(pts(:,1), pts(:,2), 50);
shp.HoleThreshold = (max(pts(:,1)) - min(pts(:,1)))*(max(pts(:,2)) - min(pts(:,2)));
shp.RegionThreshold = ((max(pts(:,1)) - min(pts(:,1)))*(max(pts(:,2)) - min(pts(:,2))))/1000;


% edgelist = shp.boundaryFacets();
% borderPolygon = zeros(size(edgelist, 1)+1, 2);
% 
% edgeNext = edgelist(1,:);
% borderPolygon(1, :) = [shp.Points(edgeNext(1), 1), shp.Points(edgeNext(1), 2)];
% 
% for k = 2:size(edgelist, 1)
%     borderPolygon(k, :) = [shp.Points(edgeNext(1), 1), shp.Points(edgeNext(1), 2)];
%     edgeNext = edgelist((edgelist(:,1) == edgeNext(2)), :);
% %     disp(edgeNext)
%     
% end
% 
% borderPolygon(end,:) = borderPolygon(1,:);


vFig = figure(7);
clf(vFig);
vAx = axes('parent', vFig);
v = voronoi(vAx, dT);
set(v(2), 'color', [52, 152, 219]/255);
set(v(1), 'color', 'none');
hold on
% triplot(dT, 'color', [230, 126, 34]/255);
set(findobj('marker', '.', 'parent', gca), 'color', 'none');


plot(shp, 'facecolor', [46, 204, 113]/255, 'facealpha', 0.6, 'edgecolor', 'none');
plot(vAx, pts(:,1), pts(:,2), '.', 'color', [.5 .5 .5]);

hold off
set(vAx, 'xlim', ptsRange, 'ylim', ptsRange);
set(vAx, 'xtick', [], 'ytick', []);
set(vAx, 'color', [1 1 1]);

patch( [100 600 600 100 100] + ptsRange(1), [100 100 150 150 100] + ptsRange(1), 'k', 'parent', vAx); 


%%
% Demo plots
% Small plots zooming into center, with additional annotation

pts = Masterpts; % take points in range 0 -> round(FieldSize/DB.FieldShrink) in middle of image
ptsRange = FieldSize/2 + [-FieldSize/(2*Demo.FieldShrink) FieldSize/(2*Demo.FieldShrink)];
pts((pts(:,1) > ptsRange(2)) | (pts(:,1) < ptsRange(1)), :) = [];
pts((pts(:,2) > ptsRange(2)) | (pts(:,2) < ptsRange(1)), :) = [];

%%
% Nearest Neighbor pairs
[IDX, D] = knnsearch(pts, pts, 'K', 2);
IDX(:,1) = [];
D(:,1) = [];

nnDemoFig = figure(8);
clf(nnDemoFig);
nnDemoAx = axes('parent', nnDemoFig);
hold on
for k = 1:numel(IDX)
    plot(nnDemoAx, [pts(k,1) pts(IDX(k), 1)], [pts(k,2) pts(IDX(k), 2)], 'color', rgb(231, 76, 60), 'linewidth', 2);
end

plot(nnDemoAx, pts(:,1), pts(:,2), '.', 'markersize', 12, 'color', rgb(127, 140, 141));

hold off
patch( [50 550 550 50 50] + ptsRange(1), [30 30 50 50 30] + ptsRange(1), 'k', 'parent', nnDemoAx); 
set(nnDemoFig, 'position', [100 100 250 250], 'color', [1 1 1]);
set(nnDemoAx, 'xtick', [], 'ytick', [], 'box', 'on');

%%
% RipleyK
rkDemoFig = figure(9);
clf(rkDemoFig);
rkDemoAx = axes('parent', rkDemoFig);

% Choose peak closest to middle of image
midPkDist = sum(abs((pts - mean(ptsRange))), 2);
midPkID = find(midPkDist == min(midPkDist));

[cx, cy] = pol2cart((0:.01:2*pi)', repmat(100, numel(0:0.01:2*pi), 1));
cx = cx + pts(midPkID, 1);
cy = cy + pts(midPkID, 2);
set(rkDemoAx, 'nextplot', 'add');
patch(cx, cy, rgb(46, 204, 113), 'facealpha', 0.6, 'linewidth', 2, 'edgecolor', rgb(46, 204, 113));
plot(rkDemoAx, [pts(midPkID, 1), cx(200)], [pts(midPkID, 2), cy(200)], 'k', 'linewidth', 2);
text(2.3792e3,    2.6728e3, 'r = 100 nm');

[cx, cy] = pol2cart((0:.01:2*pi)', repmat(250, numel(0:0.01:2*pi), 1));
cx = cx + pts(midPkID, 1);
cy = cy + pts(midPkID, 2);
set(rkDemoAx, 'nextplot', 'add');
patch(cx, cy, rgb(46, 204, 113), 'facealpha', 0.4, 'linewidth', 2, 'edgecolor', rgb(46, 204, 113));
plot(rkDemoAx, [pts(midPkID, 1), cx(300)], [pts(midPkID, 2), cy(300)], 'k', 'linewidth', 2);
text(2.0169e3,    2.6066e3, 'r = 250 nm');

[cx, cy] = pol2cart((0:.01:2*pi)', repmat(400, numel(0:0.01:2*pi), 1));
cx = cx + pts(midPkID, 1);
cy = cy + pts(midPkID, 2);
set(rkDemoAx, 'nextplot', 'add');
patch(cx, cy, rgb(46, 204, 113), 'facealpha', 0.2, 'linewidth', 2, 'edgecolor', rgb(46, 204, 113));
plot(rkDemoAx, [pts(midPkID, 1), cx(400)], [pts(midPkID, 2), cy(400)], 'k', 'linewidth', 2);
text(2.0653e3,    2.1323e3, 'r = 400 nm');

uistack(findobj('color', 'k', 'parent', rkDemoAx), 'top');
plot(rkDemoAx, pts(:,1), pts(:,2), '.', 'markersize', 12, 'color', rgb(127, 140, 141));

set(rkDemoFig, 'position', [100 100 250 250], 'color', [1 1 1]);
set(rkDemoAx, 'xtick', [], 'ytick', [], 'box', 'on');
patch( [50 550 550 50 50] + ptsRange(1), [30 30 50 50 30] + ptsRange(1), 'k', 'parent', rkDemoAx); 

%% 
% Pair correlation
pcDemoFig = figure(10);
clf(pcDemoFig);
pcDemoAx = axes('parent', pcDemoFig);

% Choose peak closest to middle of image
midPkDist = sum(abs((pts - mean(ptsRange))), 2);
midPkID = find(midPkDist == min(midPkDist));

[cx, cy] = pol2cart((0:.01:2*pi)', repmat(410, numel(0:0.01:2*pi), 1));
cx = cx + pts(midPkID, 1);
cy = cy + pts(midPkID, 2);
set(pcDemoAx, 'nextplot', 'add');
patch(cx, cy, rgb(41, 128, 185), 'facealpha', 0.6, 'linewidth', 1, 'edgecolor', rgb(41, 128, 185));
plot(pcDemoAx, [pts(midPkID, 1), cx(400)], [pts(midPkID, 2), cy(400)], 'k', 'linewidth', 2);
[cx, cy] = pol2cart((0:.01:2*pi)', repmat(390, numel(0:0.01:2*pi), 1));
cx = cx + pts(midPkID, 1);
cy = cy + pts(midPkID, 2);
set(pcDemoAx, 'nextplot', 'add');
patch(cx, cy, 'w', 'linewidth', 1, 'edgecolor', rgb(41, 128, 185));
text(2.0653e3,    2.1323e3, 'r = 400 nm');
text(2.8996e3,    2.5756e3 , '\delta r');
plot(pcDemoAx, [cx(1)-5, cx(1)+25], [cy(1), cy(1)], 'k', 'linewidth', 2);
plot(pcDemoAx, [cx(1)-12, cx(1)-12], [cy(1)-22, cy(1)+22], 'k', 'linewidth', 2);
plot(pcDemoAx, [cx(1)+32, cx(1)+32], [cy(1)-22, cy(1)+22], 'k', 'linewidth', 2);

[cx, cy] = pol2cart((0:.01:2*pi)', repmat(260, numel(0:0.01:2*pi), 1));
cx = cx + pts(midPkID, 1);
cy = cy + pts(midPkID, 2);
set(pcDemoAx, 'nextplot', 'add');
patch(cx, cy, rgb(41, 128, 185), 'facealpha', 0.6, 'linewidth', 1, 'edgecolor', rgb(41, 128, 185));
plot(pcDemoAx, [pts(midPkID, 1), cx(300)], [pts(midPkID, 2), cy(300)], 'k', 'linewidth', 2);
[cx, cy] = pol2cart((0:.01:2*pi)', repmat(240, numel(0:0.01:2*pi), 1));
cx = cx + pts(midPkID, 1);
cy = cy + pts(midPkID, 2);
set(pcDemoAx, 'nextplot', 'add');
patch(cx, cy, 'w', 'linewidth', 1, 'edgecolor', rgb(41, 128, 185));
text(2.0169e3,    2.6066e3, 'r = 250 nm');

[cx, cy] = pol2cart((0:.01:2*pi)', repmat(110, numel(0:0.01:2*pi), 1));
cx = cx + pts(midPkID, 1);
cy = cy + pts(midPkID, 2);
set(pcDemoAx, 'nextplot', 'add');
patch(cx, cy, rgb(41, 128, 185), 'facealpha', 0.6, 'linewidth', 1, 'edgecolor', rgb(41, 128, 185));
plot(pcDemoAx, [pts(midPkID, 1), cx(200)], [pts(midPkID, 2), cy(200)], 'k', 'linewidth', 2);
[cx, cy] = pol2cart((0:.01:2*pi)', repmat(90, numel(0:0.01:2*pi), 1));
cx = cx + pts(midPkID, 1);
cy = cy + pts(midPkID, 2);
set(pcDemoAx, 'nextplot', 'add');
patch(cx, cy, 'w', 'linewidth', 1, 'edgecolor', rgb(41, 128, 185));
text(2.3792e3,    2.6728e3, 'r = 100 nm');
% 

% 
uistack(findobj('color', 'k', 'parent', pcDemoAx), 'top');
plot(pcDemoAx, pts(:,1), pts(:,2), '.', 'markersize', 12, 'color', rgb(127, 140, 141));

set(pcDemoFig, 'position', [100 100 250 250], 'color', [1 1 1]);
set(pcDemoAx, 'xtick', [], 'ytick', [], 'box', 'on');
patch( [50 550 550 50 50] + ptsRange(1), [30 30 50 50 30] + ptsRange(1), 'k', 'parent', pcDemoAx); 

%%
% DBSCAN!



pts = Masterpts; % take points in range 0 -> round(FieldSize/DB.FieldShrink) in middle of image
ptsRange = FieldSize/2 + [-FieldSize/(1.3*Demo.FieldShrink) FieldSize/(1.3*Demo.FieldShrink)];
pts((pts(:,1) > ptsRange(2)) | (pts(:,1) < ptsRange(1)), :) = [];
pts((pts(:,2) > ptsRange(2)) | (pts(:,2) < ptsRange(1)), :) = [];

dbDemoFig = figure(11);
clf(dbDemoFig);
dbDemoAx = axes('parent', dbDemoFig, 'nextplot', 'add');

% Choose peak closest to middle of image
midPkDist = sum(abs((pts - 2800)), 2);
midPkID = find(midPkDist == min(midPkDist));

clusterIDs = dbscanViaELKI(pts, DB.Epsilon, DB.minPts, 'XY');
clustList = unique(clusterIDs);
clustXVar = zeros(numel(clustList), 1);
for k = 1:numel(clustList)
    clustXVar(k) = var(pts(clusterIDs == clustList(k), 1)) + var(pts(clusterIDs == clustList(k), 2));
end
noiseID = clustList(clustXVar == max(clustXVar(:)));

[cx, cy] = pol2cart((0:.01:2*pi)', repmat(DB.Epsilon, numel(0:0.01:2*pi), 1));

clustColors = [rgb(41, 128, 185); rgb(46, 204, 113); rgb(241, 196, 15);      rgb(26, 188, 156);  0.6 0.6 0.6; 0 0 1];

for k = 1:numel(clustList)
    
    if clustList(k) ~= noiseID
        
        ptsHere = pts(clusterIDs == clustList(k), :);
        
        for m = 1:size(ptsHere, 1)
            
            plot(dbDemoAx, cx + ptsHere(m, 1), cy + ptsHere(m, 2), 'color', clustColors(k, :));
            
        end
    end
end
        
        
for k = 1:numel(clustList)
    plot(dbDemoAx, pts(clusterIDs == clustList(k), 1), pts(clusterIDs == clustList(k),2), '.', 'markersize', 12, 'color', clustColors(k, :)*0.8);
end

ptsRange = FieldSize/2 + [-FieldSize/(2*Demo.FieldShrink) FieldSize/(2*Demo.FieldShrink)];
set(dbDemoFig, 'position', [100 100 250 250], 'color', [1 1 1]);
set(dbDemoAx, 'xtick', [], 'ytick', [], 'box', 'on', 'xlim', ptsRange, 'ylim', ptsRange);
patch( [50 550 550 50 50] + ptsRange(1), [30 30 50 50 30] + ptsRange(1), 'k', 'parent', dbDemoAx); 

%%
% Voronoi + alpha shape

pts = Masterpts; % take points in range 0 -> round(FieldSize/DB.FieldShrink) in middle of image
ptsRange = FieldSize/2 + [-FieldSize/(2*VD.FieldShrink) FieldSize/(2*VD.FieldShrink)];
pts((pts(:,1) > ptsRange(2)) | (pts(:,1) < ptsRange(1)), :) = [];
pts((pts(:,2) > ptsRange(2)) | (pts(:,2) < ptsRange(1)), :) = [];

% Triangulate and Voronoi
dT = delaunayTriangulation(pts(:,1), pts(:,2));

shp = alphaShape(pts(:,1), pts(:,2), 50);
shp.HoleThreshold = (max(pts(:,1)) - min(pts(:,1)))*(max(pts(:,2)) - min(pts(:,2)));
shp.RegionThreshold = ((max(pts(:,1)) - min(pts(:,1)))*(max(pts(:,2)) - min(pts(:,2))))/1000;

vDemoFig = figure(12);
clf(vDemoFig);
vDemoAx = axes('parent', vDemoFig);
v = voronoi(vDemoAx, dT);
set(v(2), 'color', [52, 152, 219]/255);
set(v(1), 'color', 'none');
hold on
triplot(dT, 'color', ([250, 180, 150]/255));
set(findobj('marker', '.', 'parent', gca), 'color', 'none');


plot(shp, 'facecolor', [46, 204, 113]/255, 'facealpha', 0.7, 'edgecolor', 'none');
plot(vDemoAx, pts(:,1), pts(:,2), '.', 'color', [.6 .6 .6], 'markersize', 12);

hold off
ptsRange = FieldSize/2 + [-FieldSize/(2*Demo.FieldShrink) FieldSize/(2*Demo.FieldShrink)];
set(vDemoFig, 'position', [100 100 250 250], 'color', [1 1 1]);
set(vDemoAx, 'xtick', [], 'ytick', [], 'box', 'on', 'xlim', ptsRange, 'ylim', ptsRange);

patch( [50 550 550 50 50] + ptsRange(1), [30 30 50 50 30] + ptsRange(1), 'k', 'parent', vDemoAx); 


%% Output all figures as SVG and PNGs

if saveImagesToDisk
    figList = findobj('type', 'figure');
    for k = 1:numel(figList)

        set(figList(k), 'PaperPositionMode', 'auto');
    %     print(figList(k), sprintf('%s\\FigFile_%.0d.svg', pwd, k), '-dsvg');
        print(figList(k), sprintf('%s\\KG01Figures\\FigFile_%.0d.eps', pwd, k), '-depsc');
        print(figList(k), sprintf('%s\\KG01Figures\\FigFile_%.0d.png', pwd, k), '-dpng', '-r900');

    end  
end






