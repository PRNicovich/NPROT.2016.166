% Draw timeline figure for Nature Protocols paper
% Want some filamentous stuff with some clusters around/in it (at
% junctions?) that can then be sampled and measured a few different ways

FieldSize = 1000; % in nm, for simulations
RandFieldPoints = 5e2;
PeakSigma = [40 60]; % nm
PtsPerPeak = [60 80]; 
NFilaments = 40;
FilamentCoverage = 0.7;
OversampleImage = 4;
PtsPerFilSpot = [1 5];
FilSpotSigma = [10 20]; % nm
SMLMImageSize = 5000; % in nm
PSFWidth = 100; % Widefield PSF width
widefieldImgSize = 75; % Pixels

PCF.Resolution_corr = 5;
PCF.rmax = 500;

DB.Epsilon = 50;
DB.minPts = 20;

% Make filament image
% Draw a bunch of polynomial curves, then rotate a few of them slightly
% Overlay should look something like microtubules

baseFig = false(FieldSize*OversampleImage, FieldSize*OversampleImage);

figure(1)
clf(1)
hold on

polyDom = linspace(0, pi, FieldSize*OversampleImage) - pi/2;

for k = 1:NFilaments
   
    poly = rand(1,1)*sin((pi + rand(1,1))+5*polyDom*(rand(1,1))) + polyDom + pi*rand(1,1);
    
    rotTheta = pi*rand(1,1);
    rotMat = [cos(rotTheta) -sin(rotTheta); sin(rotTheta), cos(rotTheta)];
    
    poly = rotMat*[poly' polyDom']';
    poly = poly';
    
    domPix = OversampleImage*FieldSize*(poly(:,1) - min(poly(:,1)))/(max(poly(:,1)) - min(poly(:,1)));
    domPoly = OversampleImage*FieldSize*(poly(:,2) - min(poly(:,2)))/(max(poly(:,2)) - min(poly(:,2)));
    domPix(1) = [];
    domPoly(1) = [];
    
    domPoly(domPix < 1) = [];
    domPix(domPix < 1) = [];

    domPix(domPoly < 1) = [];
    domPoly(domPoly < 1) = [];
    
    for m = 1:numel(domPix)
        baseFig(round(domPix(m)), round(domPoly(m))) = baseFig(round(domPix(m)), round(domPoly(m))) + 1;
    end
    
end

baseFig = baseFig(0.5*size(baseFig,1):0.6*size(baseFig,1), 0.5*size(baseFig,2):0.6*size(baseFig,2));


% % Bin down to OversampleImage size
% % Downsample baseImg to dsImg
% dsImg = zeros(FieldSize, FieldSize);
% for p = 2:FieldSize
%     for m = 2:FieldSize
%         
%         kernHere = baseFig(((p-1)*OversampleImage+1):((p*OversampleImage)), ...
%             ((m-1)*OversampleImage+1):((m*OversampleImage)));
%         dsImg(p, m) = sum(kernHere(:));
%         
%     end
% end
snapnow;





%%
skelImg = bwmorph(bwmorph(baseFig, 'dilate', 2), 'erode', 1);
skelImg = bwmorph(skelImg, 'thin', Inf);

crossPoints = bwmorph(skelImg, 'branchpoints');
[cpY, cpX] = find(crossPoints);
% 
imagesc(baseFig);
hold on
for k = 1:numel(cpX)
    plot(cpX, cpY, 'rx');
end
hold off

% From idealized image above, make a SMLM image + widefield image

[c1x, c1y] = find(baseFig);
% Stretch to SMLMImageSize 
c1x = c1x*SMLMImageSize/size(baseFig,1) + rand(1,1);
c1y = c1y*SMLMImageSize/size(baseFig,2) + rand(1,1);

cwhich = randsample(1:numel(c1x), round(numel(c1x)*FilamentCoverage));

chan1 = [c1x(cwhich), c1y(cwhich)];

cpY = cpY*SMLMImageSize/size(baseFig,1) + rand(1,1);
cpX = cpX*SMLMImageSize/size(baseFig,2) + rand(1,1);
chan2 = [cpY, cpX];

% Widefield image from chan1 and chan2 data points
wfImgDom = linspace(0, SMLMImageSize, widefieldImgSize);
wfImgChan1 = hist2(chan1(:,1), chan1(:,2), wfImgDom, wfImgDom);
wfImgChan2 = hist2(chan2(:,1), chan2(:,2), wfImgDom, wfImgDom);

PSFInPixels = PSFWidth/(SMLMImageSize/widefieldImgSize);
PSFMatrix = fspecial('gaussian', 71, PSFInPixels); % Ensure second value here is odd

wfImgChan1 = imfilter(wfImgChan1, PSFMatrix,  'symmetric', 'same');
wfImgChan2 = imfilter(wfImgChan2, PSFMatrix,  'symmetric', 'same');

wfImgChan1 = 4095*(wfImgChan1 - min(wfImgChan1(:)))/(max(wfImgChan1(:)) - min(wfImgChan1(:)));
wfImgChan1 = wfImgChan1 + 200;
wfImgChan1 = poissrnd(wfImgChan1);

wfImgChan2 = 4095*(wfImgChan2 - min(wfImgChan2(:)))/(max(wfImgChan2(:)) - min(wfImgChan2(:)));
wfImgChan2 = wfImgChan2 + 200;
wfImgChan2 = poissrnd(wfImgChan2);


figure(2)
wfImg = imfuseSpecColors((wfImgChan1/max(wfImgChan1(:))), (wfImgChan2/max(wfImgChan2(:))), rgb(39, 174, 96), rgb(142, 68, 173), 'white');
image(wfImg);
axis image
set(gca, 'xtick', [], 'ytick', []);
set(gcf, 'color', [1 1 1]);

% SMLM image from same data set

chan1Pts = [];
for k = 1:size(chan1, 1)
   
    pS = diff(FilSpotSigma)*rand(1)+min(FilSpotSigma);
    pN = round(diff(PtsPerFilSpot)*rand(1)+min(PtsPerFilSpot));
    chan1Pts = [chan1Pts; pS*randn(pN, 2)+repmat(chan1(k,:), pN, 1)];
    
    
end

randPts = SMLMImageSize*rand(RandFieldPoints, 2);
chan1Pts = [chan1Pts; randPts];

chan2Pts = [];
for k = 1:size(chan2, 1)
   
    pS = diff(PeakSigma)*rand(1)+min(PeakSigma);
    pN = round(diff(PtsPerPeak)*rand(1)+min(PtsPerPeak));
    chan2Pts = [chan2Pts; pS*randn(pN, 2)+repmat(chan2(k,:), pN, 1)];
    
    
end

randPts = SMLMImageSize*rand(RandFieldPoints/2, 2);
chan2Pts = [chan2Pts; randPts];

SRdom = linspace(1, SMLMImageSize, SMLMImageSize);
chan1SRhist = hist2(chan1Pts(:,1), chan1Pts(:,2), SRdom, SRdom);
chan2SRhist = hist2(chan2Pts(:,1), chan2Pts(:,2), SRdom, SRdom);

PSFInPixels = mean(FilSpotSigma);
PSFMatrix = fspecial('gaussian', 71, PSFInPixels); % Ensure second value here is odd
chan1SRImg = imfilter(chan1SRhist, PSFMatrix,  'symmetric', 'same');
chan2SRImg = imfilter(chan2SRhist, PSFMatrix,  'symmetric', 'same');

figure(3)
srImg = imfuseSpecColors((chan1SRImg/max(chan1SRImg(:))), (chan2SRImg/max(chan2SRImg(:))), rgb(39, 174, 96), rgb(142, 68, 173), 'white');
image(srImg);
axis image
set(gca, 'xtick', [], 'ytick', []);
set(gcf, 'color', [1 1 1]);

%%
% %% Analysis images
% clusterIDs = dbscanViaELKI(chan2Pts, 60, DB.minPts, 'XY');
% 
% dbFig = figure(5);
% clf(dbFig);
% dbax = axes('parent', dbFig);
% set(dbax, 'NextPlot', 'add');
% 
% clustList = unique(clusterIDs);
% 
% clustColors = jet(numel(clustList));
% 
% % clustXVar = zeros(numel(clustList), 1);
% % for k = 1:numel(clustList)
% %     clustXVar(k) = var(chan2Pts(clusterIDs == clustList(k), 1));
% % end
% % clustColors(clustXVar == max(clustXVar(:)), :) = [0.4 0.4 0.4];
% 
% convHullDia = zeros(numel(clustColors), 1);
% for k = 1:numel(clustList)
%    
%     pointsList = chan2Pts(clusterIDs == clustList(k),:);
%     plot(pointsList(:,1), pointsList(:,2), '.', 'color', clustColors(k,:));
%     
%     % Segment clusters with convex hull, then determine equivalent
%     % diameter'
%     if clustXVar(k) ~= max(clustXVar(:))
% %         
%         [hull, vol] = convhull(pointsList(:,1), pointsList(:,2));
%         
% %         plot(pointsList(hull, 1), pointsList(hull, 2), 'color', clustColors(k,:), 'linewidth', 2);
%         convHullDia(k) = 2*sqrt(vol/pi);
% %         
%     end
% end


%%

pts = GenerateNatProtData(SMLMImageSize, RandFieldPoints, PeakSigma, 100, PtsPerPeak);

Edges_x = 0:PCF.Resolution_corr:(SMLMImageSize); 
Edges_y = 0:PCF.Resolution_corr:(SMLMImageSize);

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
set(pcAx, 'xlim', [10 2500], 'fontsize', 12, 'ytick', 0:4:16);
set(pcFig, 'color', [1 1 1]);
%% SMLM stack images

FracPerImg = 0.0008;

for m = 1:5;
singleImg = zeros(widefieldImgSize, widefieldImgSize);
whichPts = randsample(1:size(chan1Pts, 1), round(size(chan1Pts, 1)*FracPerImg));
for k = 1:numel(whichPts)
   
    intensity = rand(1)*1000 + 500;
    Parameters = [intensity, chan1Pts(whichPts(k), 1)/(SMLMImageSize/widefieldImgSize), ...
        chan1Pts(whichPts(k), 2)/(SMLMImageSize/widefieldImgSize),...
        sqrt(PSFInPixels), sqrt(PSFInPixels), 0];
    singleImg = singleImg + gauss_2D(Parameters, (1:widefieldImgSize)');
   

    
end

singleImg = poissrnd(singleImg + 200);

    figure(11 + m)
    imagesc(singleImg)
    axis image
    set(gca, 'xtick', [], 'ytick', []);
    colormap('gray');
    set(gcf, 'color', [1 1 1]);
    drawnow
    
    pause(1)
    
end

%% Localization figures

locImg = zeros(12, 12);

xCent = 0.9*rand(1)+6;
yCent = 0.9*rand(1)+6;

intensity = rand(1)*1000 + 500;
    Parameters = [intensity, xCent, ...
        yCent,...
        sqrt(PSFInPixels)/2, sqrt(PSFInPixels)/2, 0];
    locImg = locImg + gauss_2D(Parameters, (1:12)');
    
    locImg = poissrnd(locImg + 200);
    
locFig = figure(20);
clf(locFig);
locAx = axes('parent', locFig);
imagesc(locImg)
colormap('gray')
set(locAx, 'xtick', [], 'ytick', [], 'nextplot', 'add');
set(locFig, 'color', [1 1 1]);
axis image
xlabel('X Position', 'FontSize', 12); ylabel('Y Position', 'FontSize', 12);

[a, b] = find(locImg == max(locImg(:)));

plot(locAx, 1:12, a*ones(1, 12), '--', 'color', rgb(52, 152, 219), 'linewidth', 2);
plot(locAx, b*ones(1, 12), 1:12, '--', 'color', rgb(230, 126, 34), 'linewidth', 2);

set(locAx, 'xtick', [], 'ytick', [], 'nextplot', 'replace');

% cross-section plots
xCrossFig = figure(21);
clf(xCrossFig);
xCrossAx = axes('parent', xCrossFig);
bar(xCrossAx, 1:12, locImg(a, :)'-min(locImg(:)), 1, 'facecolor', rgb(52, 152, 219), 'edgecolor', rgb(41, 128, 185), 'facealpha', 0.2);
set(xCrossAx, 'xtick', [], 'ytick', [], 'nextplot', 'add');
set(xCrossFig, 'color', [1 1 1]);
plot(xCrossAx, 1:.1:12, gauss_1D([intensity + 200 - min(locImg(:)), xCent, sqrt(PSFInPixels)/2, 0], (1:.1:12)'), ...
    'color', rgb(52, 152, 219), 'linewidth', 2);
set(xCrossAx, 'nextplot', 'replace', 'xlim', [0.5 12.5], 'box', 'off');
xlabel('X Position', 'FontSize', 12); ylabel('Intensity', 'FontSize', 12);


yCrossFig = figure(22);
clf(yCrossFig);
yCrossAx = axes('parent', yCrossFig);
bar(yCrossAx, 1:12, locImg(:,b)'-min(locImg(:)), 1, 'facecolor', rgb(230, 126, 34), 'edgecolor', rgb(211, 84, 0), 'facealpha', 0.2);
set(yCrossAx, 'xtick', [], 'ytick', [], 'nextplot', 'add');
set(yCrossFig, 'color', [1 1 1]);
plot(yCrossAx, 1:.1:12, gauss_1D([intensity + 200 - min(locImg(:)), yCent, sqrt(PSFInPixels)/2, 0], (1:.1:12)'), ...
    'color', rgb(230, 126, 34), 'linewidth', 2);
set(yCrossAx, 'nextplot', 'replace', 'xlim', [0.5 12.5], 'box', 'off');
xlabel('X Position', 'FontSize', 12); ylabel('Intensity', 'FontSize', 12);

%% Draw voronoi for illustration purposes

figure(31)
v = voronoi(chan1Pts(:,2), chan1Pts(:,1));
set(v(2), 'color', rgb(52, 152, 219));
set(v(1), 'color', 'none');

hold on
b = voronoi(chan2Pts(:,2), chan2Pts(:,1));
set(b(2), 'color', rgb(243, 156, 18));
set(b(1), 'color', 'none');
hold off

set(gca, 'ydir', 'reverse', 'xlim', [0 5000], 'ylim', [0 5000], 'xtick', [], 'ytick', []);
set(gcf, 'color', [1 1 1]);


%% Print images to file

figList = findobj('type', 'figure');
for k = 1:numel(figList)
    
    set(figList(k), 'PaperPositionMode', 'auto');
%     print(figList(k), sprintf('%s\\FigFile_%.0d.svg', pwd, k), '-dsvg');
    print(figList(k), sprintf('%s\\KG01Figures\\Fig1_%.0d.eps', pwd, k), '-depsc');
    print(figList(k), sprintf('%s\\KG01Figures\\Fig1_%.0d.png', pwd, k), '-dpng', '-r900');
    
end

