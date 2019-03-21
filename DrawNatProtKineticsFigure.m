% Nature protocols figure
% Clustering w/ kinetics to get absolute number of molecules

% Make a dataset of random points + Gaussian clusters
FieldSize = 200; % in nm
RandFieldPoints = 100;
NPeaks = 5;
PeakSigma = [20 30]; % nm
PtsPerPeak = [20 60];
molDiameter = 30;

DB.Epsilon = 25;
DB.minPts = 20;
DB.plotPixel = 10;
DB.FieldShrink = 2;

%%%%%%%%%
randPts = FieldSize*rand(RandFieldPoints, 2);

molPoints = linspace(0, 2*pi, NPeaks + 1);
molPoints(end) = [];
molPoints = molPoints + 0.1;
[px, py] = pol2cart(molPoints, repmat(molDiameter/2, numel(molPoints), 1)');

peakPost = [px+FieldSize/2; py+FieldSize/2]'; % add some rotation so it's not straight up
peakPts = [];
ptsCell = cell(NPeaks+1, 1);
ptsCell{end} = randPts;
for k = 1:NPeaks;
    
    pS = diff(PeakSigma)*rand(1)+min(PeakSigma);
    pN = round(diff(PtsPerPeak)*rand(1)+min(PtsPerPeak));
    
    ptsCell{k} = pS*randn(pN, 2)+repmat(peakPost(k,:), pN, 1);
    
    peakPts = [peakPts; ptsCell{k}];
    
end
pts = [randPts; peakPts];

%%
% Cluster peaks as a first operation
clusterIDs = dbscanViaELKI(pts, DB.Epsilon, DB.minPts, 'XY');

dbFig = figure(5);
clf(dbFig);
dbax = axes('parent', dbFig);
set(dbax, 'NextPlot', 'add');

clustList = unique(clusterIDs);
clustColors = [rgb(52, 152, 219); rgb(52, 152, 219)];

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
    plot(pointsList(:,1), pointsList(:,2), '.', 'color', hsv2rgb(rgb2hsv(clustColors(k,:)).*[1 0.5 1]), 'markersize', 12);
    
        if clustXVar(k) ~= max(clustXVar(:))
        
        [hull, vol] = convhull(pointsList(:,1), pointsList(:,2));
        
        plot(pointsList(hull, 1), pointsList(hull, 2), 'color', clustColors(k,:), 'linewidth', 2);
        
        end
    
end

% Clean up figure
set(dbax, 'xtick', [], 'ytick', [], 'box', 'on');
axis image
set(dbFig, 'color', [1 1 1]);
patch( [10 60 60 10 10], [10 10 13 13 10], 'k', 'parent', dbax); 


    
    
%% Carry forward those points identified within a cluster

noiseID = find(clustXVar == max(clustXVar(:)));
inCluster = pts(clusterIDs ~= clustList(noiseID), :);

%% Those peaks incorrectly segmented into the cluster need to get stuck into one of the underlying clusters
% or points will be lost in next step
clustAccount = false(numel(inCluster), 1);
for k = 1:(length(ptsCell)-1)
    clustAccount(ismember(inCluster, ptsCell{k})) = true;
end

needAssignment = inCluster(~clustAccount);
assignedTo = randsample(NPeaks, numel(needAssignment), true);

for k = 1:NPeaks

    ptsCell{k} = [ptsCell{k}; inCluster(assignedTo == k, :)];
    
end


%%

kinFig = figure(6);
clf(kinFig);
kax = axes('parent', kinFig);
set(kax, 'NextPlot', 'add');

colorList = [rgb(155, 89, 182); rgb(26, 188, 156); rgb(241, 196, 15); rgb(231, 76, 60); rgb(52, 152, 219)];

plot(kax, pts(clusterIDs == clustList(noiseID), 1), pts(clusterIDs == clustList(noiseID), 2), '.', 'color', [0.4, 0.4, 0.4], 'markersize', 14);

for k = 1:(length(ptsCell)-1)
    
    goodPoints = ptsCell{k}(ismember(ptsCell{k}(:,1), inCluster(:,1)), :);
    
    plot(kax, goodPoints(:,1), goodPoints(:,2), '.', 'color', hsv2rgb(rgb2hsv(colorList(k,:)).*[1 0.5 1]), 'markersize', 14);
    [hull, vol] = convhull(goodPoints(:,1), goodPoints(:,2));
    plot(kax, goodPoints(hull, 1), goodPoints(hull, 2), 'color', colorList(k,:), 'linewidth', 2);
    
    plot(kax, mean(goodPoints(:,1)), mean(goodPoints(:,2)), 'x', 'markersize', 20, 'color', colorList(k,:), 'linewidth', 3);
    
end

uistack(findobj('marker', 'x', 'parent', kax), 'top');

patch( [10 60 60 10 10], [10 10 13 13 10], 'k', 'parent', kax);

% Clean up figure
set(kax, 'xtick', [], 'ytick', [], 'box', 'on');
axis image
set(kinFig, 'color', [1 1 1]);

%% Kinetics figure
% Illustrative figure that's basically a kymograph. 
% Inspired by Annebale + Radenovic PLoS One 2012 Fig 1d

% mEos2 parameters from PNAS, 10/23/12, 17437
kill_rate = .1; % kb, 1/seconds
blink_rate = 20; % kd, 1/seconds
return_rate_1 = .1; % kr1, 1/seconds
return_rate_2 = 10; % kr2, 1/seconds
return_alpha = 0.2; % alpha, unitless

Ex_time = 30;
Steps = 10;

N_half_lives = 20;
N_frames = 10000;
N_dyes = NPeaks*10;
Intensity_fraction_threshold = 0.1;
Activation_half_life = (N_frames*Ex_time)/N_half_lives; % seconds



Step_time = Ex_time/Steps;

% Each transition probability stated as I2J, the probability in the next
% time step for a molecule in state I to transition to state J.

Activation_rate = (log(2)/Activation_half_life);
Pre2Pre = exp(-Activation_rate*Step_time);
Pre2On = 1 - exp(-Activation_rate*Step_time);
% Pre2Pre = (1 - 1/(N_frames*Steps))/Step_time;
% Pre2On = (1/(N_frames*Steps))/Step_time;
Pre2Blink1 = 0;
Pre2Blink2 = 0;
Pre2Dead = 0;

On2Pre = 0;
On2On = exp(-(blink_rate + kill_rate)*Step_time);
% On2Blink1 = return_alpha*(1 - exp(-blink_rate*Step_time));
% On2Blink2 = (1-return_alpha)*(1 - exp(-blink_rate*Step_time));
% On2Blink1 = (1 - exp(-((blink_rate/(1+return_alpha))*Step_time)));
% On2Blink2 = (1 - exp(-(blink_rate*return_alpha/(1+return_alpha))*Step_time));
% % On2Blink1 = (1 - exp(-blink_rate*return_alpha*Step_time));
% % On2Blink2 = (1 - exp(-blink_rate*(1 - return_alpha)*Step_time));
% On2Dead = (blink_rate/(blink_rate + kill_rate))*(1 - exp(-kill_rate*Step_time));
On2Blink1 = (blink_rate/((1 + return_alpha)*(blink_rate + kill_rate)))*(1 - exp(-(blink_rate + kill_rate)*Step_time));
On2Blink2 = (return_alpha*blink_rate/((1 + return_alpha)*(blink_rate + kill_rate)))*(1 - exp(-(blink_rate + kill_rate)*Step_time));
On2Dead = (kill_rate/(blink_rate + kill_rate))*(1 - exp(-(kill_rate + blink_rate)*Step_time));


Blink12Pre = 0;
Blink12On = 1 - exp(-return_rate_1*Step_time);
Blink12Blink1 = exp(-return_rate_1*Step_time);
Blink12Blink2 = 0;
Blink12Dead = 0;

Blink22Pre = 0;
Blink22On = 1 - exp(-return_rate_2*Step_time);
Blink22Blink1 = 0;
Blink22Blink2 = exp(-return_rate_2*Step_time);
Blink22Dead = 0;

Dead2Pre = 0;
Dead2On = 0;
Dead2Blink1 = 0;
Dead2Blink2 = 0;
Dead2Dead = 1;

Rate_constants = [Pre2Pre    Pre2On    Pre2Blink1    Pre2Blink2    Pre2Dead;
    On2Pre     On2On     On2Blink1     On2Blink2     On2Dead;
    Blink12Pre Blink12On Blink12Blink1 Blink12Blink2 Blink12Dead;
    Blink22Pre Blink22On Blink22Blink1 Blink22Blink2 Blink22Dead;
    Dead2Pre   Dead2On   Dead2Blink1   Dead2Blink2   Dead2Dead];


Transitions = Rate_constants;
%
% Sum of probabilities in a single line should equal 1.
% If not you've done something wrong above.
% for k = 1:length(Transitions)
%
%     sum_here = 1:length(Transitions);
%     Transitions(k,:) = Transitions(k,:)./(sum(Transitions(k,:)));
%     %Transitions(k, k) = 1 - sum(Transitions(k, sum_here(~ismember(sum_here, k))));
%
% end

Transitions(Rate_constants == 0) = 0;

Emissions = [1 0;
    0 1;
    1 0;
    1 0;
    1 0];

Symbols = [0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Hidden Markov Model traces for each dye

States_matrix = zeros(N_frames*Steps, N_dyes);
Em_matrix = zeros(N_frames*Steps, N_dyes);
disp('Modeling time traces')

for m = 1:N_dyes
    
    
    
    [seq, states] = hmmgenerate(N_frames*Steps, Transitions, Emissions, 'symbols', Symbols);
    
    States_matrix(:,m) = states;
    Em_matrix(:,m) = seq;
    
end


% Collapse Em_matrix by number of Steps to Int_matrix

bin_vect = 1:Steps:(N_frames*Steps);
Int_matrix = zeros(N_frames, size(Em_matrix,2));
States_binned = zeros(N_frames, size(Em_matrix,2));

for k = 1:(numel(bin_vect))
    
    Int_matrix(k,:) = sum(Em_matrix((bin_vect(k)):(bin_vect(k) + Steps - 1), :), 1);
    
    States_hold = States_matrix((bin_vect(k)):(bin_vect(k) + Steps - 1), :);
    
    States_binned(k, ge(sum(States_hold == 1, 1), 1)) = 1; % Pre
    States_binned(k, ge(sum(States_hold == 5, 1), 1)) = 5; % Dead
    States_binned(k, ge(sum(States_hold == 3, 1), 1)) = 3; % Blink 1
    States_binned(k, ge(sum(States_hold == 4, 1), 1)) = 4; % Blink 2
    States_binned(k, ge(sum(States_hold == 2, 1), 1)) = 2; % Active
    
    
end

% Make Em_matrix simple 1 or 0 for on/off during that frame, binned over steps.

Em_matrix = double(Int_matrix > Steps*Intensity_fraction_threshold);

keepEm = Em_matrix(:, find(any(Em_matrix), 5, 'first'));


%% Once a decent trace is made, plot it

traceImg = zeros(21, size(keepEm, 1));
traceImg(3:18, any(keepEm, 2)') = 1;
traceImg(:, (find(any(traceImg), 1, 'last')+10):end) = [];

colorList = [rgb(155, 89, 182); rgb(26, 188, 156); rgb(241, 196, 15); rgb(231, 76, 60); rgb(52, 152, 219)];

figure(15)
clf;
imagesc(traceImg)
set(gca, 'xtick', [], 'ytick', [])
colormap('gray')
xlabel('Time');
hold on
for k = 1:size(keepEm, 2)
    
    stfrm = find(keepEm(:,k), 1, 'first');
    enfrm = find(keepEm(:,k), 1, 'last');
    
    plot([stfrm-5 stfrm-5 enfrm+5 enfrm+5 stfrm-5], [2 19 19 2 2], 'color', colorList(k,:), 'linewidth', 2);
end

set(gcf, 'color', [1 1 1]);

%% Print 'em!

% figList = findobj('type', 'figure');
% for k = 1:numel(figList)
%     
%     set(figList(k), 'PaperPositionMode', 'auto');
%     print(figList(k), sprintf('%s\\KG01Figures\\KineticsFigsFile_%.0d.eps', pwd, k), '-depsc');
%     print(figList(k), sprintf('%s\\KG01Figures\\KineticsFigsFile_%.0d.png', pwd, k), '-dpng', '-r900');
%     
% end

    



