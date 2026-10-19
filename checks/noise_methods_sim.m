function noise_methods_sim()
% NOISE_METHODS_SIM  Test bench for noise-robust scores (RESULTS.md section 19).
% Reviewer points R2 6 (noise floor), R2 4a/4b. Seeded simulations with KNOWN
% granule dynamics and noise varied independently:
%   relative burst size s in {0.03, 0.06, 0.12} of the current granule level
%   (+-20% jitter), burst rate r in {0, 3, 10} per hour (Poisson; half losses,
%   half gains, multiplicative steps), noise SD n in {0.005, 0.01, 0.02} x the
%   initial granule level; granule level 3000-8000, bleaching 0.0026/min on
%   granule and background (1200 + SD 10); 100 traces x 361 frames, seeds 1:2.
% True negative turnover = 100 * (r/2) * s  (% of the granule shed per hour).
% For every simulated set, exactly the real-data code path:
%   - burst detection (analysis.m, window 13, SP 10/12/14) on the signal
%     (granule + background, as for the spot sum) and on its burst-free
%     surrogate (cubic trend + white noise at the trace's own level, as
%     score_export2 'poly');
%   - structure function of the relative granule intensity
%     y = g / cubicTrend(g) - 1:  D(tau) = mean (y(t+tau) - y(t))^2 for
%     tau = 1 2 3 6 12 18 30 45 60 frames, real and surrogate.
% Writes analysis_outputs/noise_methods/sim_bursts.csv, sim_traces.csv, sim_sf.csv.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'noise_methods');
if ~exist(outDir, 'dir'), mkdir(outDir); end
taus = [1 2 3 6 12 18 30 45 60]; Lf = 361; nT = 100; tVec = (0:360) * 10/60; tMin = (0:Lf-1) * 10/60;
bu = cell(0, 11); trs = cell(0, 6); sf = cell(0, 8);
for s = [0.03 0.06 0.12]
    for r = [0 3 10]
        for n = [0.005 0.01 0.02]
            for seed = 1:2
                rng(round(1e4*s) + 100*r + round(1e4*n) + 7*seed);
                tag = sprintf('s%.2f_r%d_n%.3f_%d', s, r, n, seed);
                bl = exp(-0.0026 * tMin);
                G = zeros(nT, Lf); Bk = zeros(nT, Lf);
                for i = 1:nT
                    L0 = 3000 + 5000 * rand; lvl = L0 * ones(1, Lf);
                    ne = poissrnd(r * Lf * 10 / 3600);
                    ev = sort(randi([2 Lf], 1, ne)); sg = sign(rand(1, ne) - 0.5);
                    for e = 1:ne
                        f = 1 + sg(e) * s * (1 + 0.2 * (2 * rand - 1));
                        lvl(ev(e):end) = lvl(ev(e):end) * f;
                    end
                    G(i,:) = lvl .* bl + n * L0 * randn(1, Lf);
                    Bk(i,:) = 1200 * bl + 10 * randn(1, Lf);
                end
                S = G + Bk;
                Ssur = S; Gsur = G;
                for i = 1:nT
                    Ssur(i,:) = surrogate(S(i,:)); Gsur(i,:) = surrogate(G(i,:));
                end
                for i = 1:nT, trs(end+1,:) = {tag, s, r, n, i, Lf}; end %#ok<AGROW>
                for kind = {'real', 'surr'}
                    if strcmp(kind{1}, 'real'), X = S; Y = G; else, X = Ssur; Y = Gsur; end
                    for SP = [10 12 14]
                        ws = warning('off', 'all');
                        [~, ~, ~, ~, ~, ~, seg] = analysis(X, Bk, tVec, 13, SP, 0, 0);
                        warning(ws);
                        for q = 1:size(seg,1)
                            bu(end+1,:) = {tag, s, r, n, kind{1}, SP, seg(q,1), seg(q,2), seg(q,5), seg(q,6), seg(q,3)}; %#ok<AGROW>
                        end
                    end
                    for i = 1:nT
                        D = structFun(Y(i,:), taus);
                        sf(end+1,:) = {tag, s, r, n, kind{1}, i, mat2str(taus), mat2str(D, 6)}; %#ok<AGROW>
                    end
                end
                fprintf('%s done\n', tag);
            end
        end
    end
end
writetable(cell2table(bu, 'VariableNames', {'set','s','r','n','kind','SP','trace','sign','amp','levelBefore','startFrame'}), ...
    fullfile(outDir, 'sim_bursts.csv'));
writetable(cell2table(trs, 'VariableNames', {'set','s','r','n','trace','lengthFrames'}), fullfile(outDir, 'sim_traces.csv'));
writetable(cell2table(sf, 'VariableNames', {'set','s','r','n','kind','trace','taus','D'}), fullfile(outDir, 'sim_sf.csv'));
disp('NOISE_SIM_DONE');
end

function y = surrogate(x)
x = x(:)'; t = 1:numel(x); trd = polyval(polyfit(t, x, 3), t);
sg = 1.4826 * median(abs(diff(x) - median(diff(x)))) / sqrt(2);
y = trd + sg * randn(size(x));
end

function D = structFun(g, taus)
g = g(:)'; t = 1:numel(g); trd = polyval(polyfit(t, g, 3), t);
y = g ./ trd - 1;
D = arrayfun(@(k) mean((y(1+k:end) - y(1:end-k)).^2), taus);
end
