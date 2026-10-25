% VALENCY_REAL  Valency-linked joint fits (analysis/gelmetricsJoint.m) of the
% real data, per protein. Reviewer points R2 4b, R1 2; Qbeta vs PP7 claim.
% Groups (datasets ordered by valency / protein amount; ratioMax = maximum
% brightness relative to the first dataset = valency ratio):
%   Qb       5Qb, 10Qb                     [1 2]       R2 re-extraction, SP 14
%   PP7      4PP7, 30PP7                   [1 7.5]     R2 re-extraction, SP 12
%   PP7-FUS  4PP7-FUS, 30PP7-FUS           [1 7.5]     R2 re-extraction, SP 12
%   invitro  protein:RNA 1:1, 10:1, 100:1  [1 10 14]   stored data (PCP-14x RNA;
%            at most 1, 10, 14 proteins per RNA), SP 12
% One SP per group so brightness is compared under one detection setting;
% truncation 15%, top trim 10%, window 13; negative bursts.
% Fits: separate MLE per dataset; joint 'fixed' (full occupancy); joint
% 'bounded' (brightness non-decreasing, at most full occupancy).
% CIs: 50 granule-bootstrap resamples per joint fit (traces resampled within
% each dataset; seeded). Writes analysis_outputs/valency/real_fits.csv, real_boot.mat.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'valency');
if ~exist(outDir, 'dir'), mkdir(outDir); end
R2 = 'in_vivo_reextracted/after_R2_offset'; IV = 'in_vitro_experiments';
G = {'Qb',      {'5Qb','10Qb'},                  R2, 14, [1 2];
     'PP7',     {'4PP7','30PP7'},                R2, 12, [1 7.5];
     'PP7-FUS', {'4PP7-FUS','30PP7-FUS'},        R2, 12, [1 7.5];
     'invitro', {'protein_RNA_1_1','protein_RNA_10_1','protein_RNA_100_1'}, IV, 12, [1 10 14]};
tVec = (0:360) * 10/60; nBoot = 50; rows = {}; boot = struct();
for g = 1:size(G,1)
    names = G{g,2}; m = numel(names); T = cell(1,m); A = cell(1,m); sg = zeros(1,m);
    for i = 1:m
        [S, B] = gm_loadDataset(repoDir, G{g,3}, names{i});
        Ti = gm_amplitudes_by_trace(S, B, tVec, 13, G{g,4});
        T{i} = Ti([Ti.len] >= 180);
        A{i} = abs(vertcat(T{i}.ampNeg)); sg(i) = gm_sigma_from_ampvar(vertcat(T{i}.ampVar));
        Rs = gelmetricsMLE(A{i}, sg(i), struct('profile', false));
        rows(end+1,:) = {G{g,1}, 'separate', names{i}, numel(T{i}), Rs.n, sg(i), Rs.lambda, NaN, NaN, Rs.k0, NaN, NaN, NaN, Rs.nll}; %#ok<SAGROW>
    end
    for md = {'fixed', 'bounded'}
        Rj = gelmetricsJoint(A, sg, G{g,5}, md{1});
        rng(100 + g);
        bl = nan(nBoot, m); bk = nan(nBoot, m);
        for b = 1:nBoot
            Ab = cell(1,m); sb = zeros(1,m);
            for i = 1:m
                Tb = T{i}(randi(numel(T{i}), numel(T{i}), 1));
                Ab{i} = abs(vertcat(Tb.ampNeg)); sb(i) = gm_sigma_from_ampvar(vertcat(Tb.ampVar));
            end
            Rb = gelmetricsJoint(Ab, sb, G{g,5}, md{1}, struct('start', Rj));
            bl(b,:) = Rb.lambda; bk(b,:) = Rb.k0;
        end
        boot.(matlab.lang.makeValidName([G{g,1} '_' md{1}])) = struct('lambda', bl, 'k0', bk);
        for i = 1:m
            rows(end+1,:) = {G{g,1}, md{1}, names{i}, numel(T{i}), Rj.n(i), sg(i), Rj.lambda(i), ...
                prctile(bl(:,i), 2.5), prctile(bl(:,i), 97.5), Rj.k0(i), prctile(bk(:,i), 2.5), prctile(bk(:,i), 97.5), ...
                Rj.ratio(i), Rj.nll}; %#ok<SAGROW>
        end
        fprintf('%s %s: lambda %s k0 %s ratio %s nll %.1f\n', G{g,1}, md{1}, mat2str(Rj.lambda, 3), ...
            mat2str(Rj.k0, 3), mat2str(Rj.ratio, 3), Rj.nll);
    end
end
writetable(cell2table(rows, 'VariableNames', {'group','fit','dataset','nTraces180','nFit','sigma', ...
    'lambda','lamLo','lamHi','k0','k0Lo','k0Hi','ratio','nll'}), fullfile(outDir, 'real_fits.csv'));
save(fullfile(outDir, 'real_boot.mat'), '-struct', 'boot');
disp('VALENCY_REAL_DONE');
