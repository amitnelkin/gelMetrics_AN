close all;   % Close all open figure windows
clearvars;       % Clear workspace
clc;

% Algorithm Parameters
movingAverageWindowSpan = 13;             % Moving average window span for signal denoising
stringencyParameter = 14;                 % Parameter for burst events detection
frameInterval_min = 10/60;                % frame interval in minutes - set to your acquisition interval
timeVecTotal = (0:359) * frameInterval_min; % Time steps vector (was linspace(0,60,360): dt = 60/359 min)

% Load signal and background data
% StringencyParameter (SP) and truncation percentile (GOF) for replicating the exact plots: 
% [InVitro 100:1/10:1 - SP=12,GOF=10] ; [InVitro 1:1 - SP=12,GOF=15] ; [PP7/PP7-FUS 4x - SP=10,GOF=15]; [PP7/PP7-FUS 30x - SP=12,GOF=10] ; [Qb - SP=14,GOF=15]
% To recreate the figures in the paper please uncomment any of the following signals-backgrounds pairs:

% RNA-protein granules formed with 1:1 protein to RNA ratio
% signals = struct2cell(load('./in_vitro_experiments/protein_RNA_1_1_data.mat'));
% backgrounds = struct2cell(load('./in_vitro_experiments/protein_RNA_1_1_background.mat'));

% RNA-protein granules formed with 10:1 protein to RNA ratio
% signals = struct2cell(load('./in_vitro_experiments/protein_RNA_10_1_data.mat'));
% backgrounds = struct2cell(load('./in_vitro_experiments/protein_RNA_10_1_background.mat'));

% RNA-protein granules formed with 100:1 protein to RNA ratio
% signals = struct2cell(load('./in_vitro_experiments/protein_RNA_100_1_data.mat'));
% backgrounds = struct2cell(load('./in_vitro_experiments/protein_RNA_100_1_background.mat'));

% In vivo granules formed with Qb-5x
% signals = struct2cell(load('./in_vivo_experiments/5Qb_data.mat'));
% backgrounds = struct2cell(load('./in_vivo_experiments/5Qb_background.mat'));

% In vivo granules formed with Qb-10x
% signals = struct2cell(load('./in_vivo_experiments/10Qb_data.mat'));
% backgrounds = struct2cell(load('./in_vivo_experiments/10Qb_background.mat'));

% In vivo granules formed with PP7 & QCP-5x/PCP-4x
% signals = struct2cell(load('./in_vivo_experiments/4PP7_data.mat'));
% backgrounds = struct2cell(load('./in_vivo_experiments/4PP7_background.mat'));

% In vivo granules formed with PP7-FUS & QCP-5x/PCP-4x
% signals = struct2cell(load('./in_vivo_experiments/4PP7-FUS_data.mat'));
% backgrounds = struct2cell(load('./in_vivo_experiments/4PP7-FUS_background.mat'));

% In vivo granules formed with PP7 & PCP-30x
% signals = struct2cell(load('./in_vivo_experiments/30PP7_data.mat'));
% backgrounds = struct2cell(load('./in_vivo_experiments/30PP7_background.mat'));

% In vivo granules formed with PP7-FUS & PCP-30x
% signals = struct2cell(load('./in_vivo_experiments/30PP7-FUS_data.mat'));
% backgrounds = struct2cell(load('./in_vivo_experiments/30PP7-FUS_background.mat'));


% Signal plotting controls
% To plot a few sample signals you can change these variables
% plotMin is the first signal plotted
% plotMax is the last signal plotted
% Example: to plot the first 5 signals change plotMax to 5
plotMin = 0;
plotMax = 0;

signals = signals{1};
backgrounds = backgrounds{1};

if ~iscolumn(signals)
    signals = signals(2:end)';
else
    signals = signals(2:end);
end

if ~iscolumn(backgrounds)
    backgrounds = backgrounds(2:end)';
else
    backgrounds = backgrounds(2:end);
end

% Perform analysis
events_d = dictionary;
k_d = dictionary;
lam_d = dictionary;
MSE_d = dictionary;
[timeBetweenPos, timeBetweenNeg, ampPos, ampNeg, ampVar, eventsPerSignal] = ...
        analysis(signals, backgrounds, timeVecTotal, ...
        movingAverageWindowSpan, stringencyParameter, plotMin, plotMax);

% Remove the -Inf placeholder used to initialize the variable
ampPos = ampPos(2:end);
ampNeg = ampNeg(2:end);
ampVar = ampVar(2:end);

fprintf('  • Found %d positive burst events: \n', length(ampPos));
fprintf('  • Found %d unclassified events: \n', length(ampVar));
fprintf('  • Found %d negative burst events: \n', length(ampNeg));

timeBetweenPos = timeBetweenPos(2:end);
timeBetweenNeg = timeBetweenNeg(2:end);    

figure;
histogram(ampVar(abs(ampVar)<1500),50);
ylabel('Counts');
xlabel('Intesity [A.U]');
title('Quiescent Amplitudes');


% Plot histogram of all amplitudes 
figure; 
histogram(ampVar(abs(ampVar)<3000),'FaceColor','b','BinWidth',100);
hold on;
histogram(ampPos(ampPos<3000),'FaceColor','g','BinWidth',100); 
histogram(ampNeg(abs(ampNeg)<3000),'FaceColor','r','BinWidth',100);
ylabel('Counts'); xlim([-3000,3000])
title('Segment amplitudes','FontSize',12);
xlabel('Amplitudes [A.U]','FontSize',12);

% Plot boxplots of durations between burst events 
figure;
plot_var = [timeBetweenPos;timeBetweenNeg];
group_var = zeros(length(plot_var),1);
group_var(1:length(timeBetweenPos))=2;
group_var(length(timeBetweenPos)+1:end)=3;
boxplot(plot_var,group_var);
xticklabels({'Positive','Negative'});
ylabel('Time [min]');
title('Time between burst events');


lambda = 1:5;
K0Vector = 1:500; % Experimental data is noisy by nature, using standard errors will lead to unreasonably large grids
truncPercent = 15; % left-truncation percentile, passed explicitly to modifiedPoissonFit

% Poisson analysis on positive amplitudes
[~, ~, ~, ~, ~, w_score] = modifiedPoissonFit(rmoutliers(ampPos,"percentiles",[0,90]),lambda,K0Vector,10,'Positive amplitudes',truncPercent);

% Final ranking of lambdas
for ll = 1:length(lambda)
    [min_w_vec{ll}, min_w_k0{ll}] = min(w_score(ll,:));
    min_w_k0{ll} = K0Vector(min_w_k0{ll}); % convert index to k0 value (only equal when K0Vector = 1:N)
end
if any(cell2mat(min_w_k0) == K0Vector(1) | cell2mat(min_w_k0) == K0Vector(end))
    warning('Best k0 is at the edge of K0Vector - extend the grid');
end

[min_w_scores, min_w_idx] = sort(cell2mat(min_w_vec));
min_w_lambda = lambda(min_w_idx);
min_w_k0 = min_w_k0(min_w_idx);

fprintf('Positive amplitudes results:\n')
for ww = 1:length(lambda)
    fprintf('  #%d: lambda=%d, k0=%d, WS=%g\n', [ww,min_w_lambda(ww),min_w_k0{ww},round(min_w_scores(ww),2,"significant")]);
end

% Poisson analysis on negative amplitudes
[~, ~, ~, ~, ~, w_score] = modifiedPoissonFit(rmoutliers(abs(ampNeg),"percentiles",[0,90]),lambda,K0Vector,10,'Negative amplitudes',truncPercent);

% Final ranking of lambdas
for ll = 1:length(lambda)
    [min_w_vec{ll}, min_w_k0{ll}] = min(w_score(ll,:));
    min_w_k0{ll} = K0Vector(min_w_k0{ll}); % convert index to k0 value (only equal when K0Vector = 1:N)
end
if any(cell2mat(min_w_k0) == K0Vector(1) | cell2mat(min_w_k0) == K0Vector(end))
    warning('Best k0 is at the edge of K0Vector - extend the grid');
end

[min_w_scores, min_w_idx] = sort(cell2mat(min_w_vec));
min_w_lambda = lambda(min_w_idx);
min_w_k0 = min_w_k0(min_w_idx);

fprintf('Negative amplitudes results:\n')
for ww = 1:length(lambda)
    fprintf('  #%d: lambda=%d, k0=%d, WS=%g\n', [ww,min_w_lambda(ww),min_w_k0{ww},round(min_w_scores(ww),2,"significant")]);
end

