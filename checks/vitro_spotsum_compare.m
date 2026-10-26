function vitro_spotsum_compare()
% VITRO_SPOTSUM_COMPARE  In vitro results with the spot mean (the paper) vs the spot
% sum (as in vivo), RESULTS.md section 53. Same granules (same_granules_R2), same
% settings: GelMetricsReport with the paper seeds and nCopies 1, and with nCopies 20;
% GelMetricsCompare (10000 movie shuffles) for 100:1 vs 1:1 and 10:1 vs 1:1.
% With the spot sum, burst fractions, turnover, frame noise and net change use the
% integrated granule light above the dilute phase; brightness and burst amplitudes
% in a.u. still use the spot mean.
% Writes analysis_outputs/vitro_spotsum/{reports.csv, tests.csv}.
repo = fileparts(fileparts(mfilename('fullpath'))); addpath(repo, fullfile(repo, 'analysis'), fullfile(repo, 'legacy'));
out = fullfile(repo, 'analysis_outputs', 'vitro_spotsum'); if ~exist(out, 'dir'), mkdir(out); end
V = {'vitro_1_1', 'vitro_10_1', 'vitro_100_1'};
F = {'frameNoise', 'sizeGranules', 'sizeCopies', 'sizeBeyond', 'turnover', 'noiseShare', 'net', 'turnover1', 'turnover2', 'net1', 'net2', 'medianBurst', 'meanBurst', 'brightness'};
rows = {}; trow = {};
for meas = {'spotmean', 'spotsum'}
    for K = [1 20]
        R = cell(1, 3);
        for k = 1:3
            D = gm_paper_dataset(V{k}, meas{1});
            R{k} = GelMetricsReport(D, struct('nCopies', K, 'seed', D.paperSeed, 'verbose', false));
            r = {meas{1}, K, V{k}, R{k}.granules, R{k}.movies};
            for q = 1:numel(F), r = [r, {R{k}.(F{q}), R{k}.([F{q} '_lo']), R{k}.([F{q} '_hi'])}]; end %#ok<AGROW>
            r = [r, {R{k}.pHalvesTurnover, R{k}.pHalvesNet}]; %#ok<AGROW>
            rows(end+1, :) = r; %#ok<AGROW>
            fprintf('%s nCopies %d %s: beyond %.1f [%.1f, %.1f], turnover %.1f, net %.1f [%.1f, %.1f], noise %.1f%%\n', meas{1}, K, V{k}, ...
                R{k}.sizeBeyond, R{k}.sizeBeyond_lo, R{k}.sizeBeyond_hi, R{k}.turnover, R{k}.net, R{k}.net_lo, R{k}.net_hi, R{k}.frameNoise);
        end
        for pr = [3 1; 2 1]'
            T = GelMetricsCompare(R{pr(1)}, R{pr(2)});
            for m = {'sizeBeyond', 'turnover', 'net', 'sizeGranules', 'medianBurst'}
                i = strcmp(T.measure, m{1});
                trow(end+1, :) = {meas{1}, K, [V{pr(1)} ' vs ' V{pr(2)}], m{1}, T{i, 2}, T{i, 3}, T.difference(i), T.p(i), T.p_AgtB(i)}; %#ok<AGROW>
            end
        end
    end
end
names = {'measure', 'nCopies', 'dataset', 'granules', 'movies'};
for q = 1:numel(F), names = [names, {F{q}, [F{q} '_lo'], [F{q} '_hi']}]; end %#ok<AGROW>
names = [names, {'pHalvesTurnover', 'pHalvesNet'}];
writetable(cell2table(rows, 'VariableNames', names), fullfile(out, 'reports.csv'));
writetable(cell2table(trow, 'VariableNames', {'measure', 'nCopies', 'comparison', 'quantity', 'A', 'B', 'difference', 'p', 'p_AgtB'}), fullfile(out, 'tests.csv'));
disp('SPOTSUM_DONE');
end
