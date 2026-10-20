function ms_sim_test()
% MS_SIM_TEST  The manuscript's own approach on simulated data (RESULTS.md
% section 30): simulated burst signals as MS section III (generateSim type 3,
% 1000 signals x 361 frames, k0 = 60, white noise 40 a.u., bleaching
% 0.0026 /min), the manuscript's burst detection (analysis.m, window 13,
% SP 15 as for the MS simulations), negative bursts, no data removed.
% lambda in {0.5 1 2 3 4 5}; noise 40 (MS) and 4 (nearly noise-free, best case);
% 5 independent datasets per setting (seeded).
% Also saved per dataset: the true negative jump sizes of the noise-free twin
% (same seed, noise 0: generateSim draws the same random numbers), i.e. what
% Eq. (1.9) is meant to describe; and for one dataset (lambda 2, noise 40) three
% example traces with their processed signal and detected bursts.
% Writes analysis_outputs/ms_sim_test/{amplitudes.csv, true_jumps.csv,
% example_traces.csv, example_bursts.csv}.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'ms_sim_test'); if ~exist(outDir, 'dir'), mkdir(outDir); end
lams = [0.5 1 2 3 4 5]; noises = [40 4]; nRep = 5; nT = 1000; nF = 361; k = 0.0026; SP = 15; w = 13;
tVec = (0:360) * 10/60; A = cell(0, 6); J = cell(0, 5); ex = cell(0, 5); eb = cell(0, 6); id = 0;
for nz = noises
    for lam = lams
        for rep = 1:nRep
            id = id + 1; seed = 900 + id;
            rng(seed); [S, B] = generateSim(nT, nF, 3, k, nz, lam, 60);
            rng(seed); S0 = generateSim(nT, nF, 3, k, 0, lam, 60);
            ws = warning('off', 'all');
            [~, ~, ~, ~, ~, ~, seg] = analysis(num2cell(S, 2), num2cell(B, 2), tVec, w, SP, 0, 0);
            warning(ws);
            neg = seg(seg(:,2) < 0, :);
            for q = 1:size(neg, 1), A(end+1,:) = {id, lam, nz, rep, neg(q,1), abs(neg(q,5))}; end %#ok<AGROW>
            u = S0 ./ exp(-k * linspace(0, 60, nF)); du = diff(u, 1, 2);
            [ti, ~] = find(du < -1e-9); v = -du(du < -1e-9);
            for q = 1:numel(v), J(end+1,:) = {id, lam, nz, rep, v(q)}; end %#ok<AGROW>
            if lam == 2 && nz == 40 && rep == 1
                for i = 1:3
                    [proc, sg] = detectOne(S(i,:)', B(i,:)', tVec, w, SP);
                    for t = 1:nF, ex(end+1,:) = {i, t, S(i,t), proc(t), S0(i,t)}; end %#ok<AGROW>
                    for q = 1:size(sg, 1), eb(end+1,:) = {i, sg(q,2), sg(q,3), sg(q,4), sg(q,6), sg(q,7)}; end %#ok<AGROW>
                end
            end
            fprintf('lambda %g noise %d rep %d: %d bursts\n', lam, nz, rep, size(neg, 1));
        end
    end
end
writetable(cell2table(A, 'VariableNames', {'id','lambda','noise','rep','trace','amp'}), fullfile(outDir, 'amplitudes.csv'));
writetable(cell2table(J, 'VariableNames', {'id','lambda','noise','rep','jump'}), fullfile(outDir, 'true_jumps.csv'));
writetable(cell2table(ex, 'VariableNames', {'trace','frame','raw','processed','noiseFree'}), fullfile(outDir, 'example_traces.csv'));
writetable(cell2table(eb, 'VariableNames', {'trace','sign','firstFrame','lastFrame','levelStart','levelEnd'}), fullfile(outDir, 'example_bursts.csv'));
disp('MS_SIM_DONE');
end

function [data, seg] = detectOne(s, b, tVec, w, SP)
% the processed signal exactly as analysis.m builds it, plus its bursts
ws = warning('off', 'all');
tv = tVec(1:numel(s));
f = fit(tv', b, 'poly3'); bf = feval(f, tv);
x = filloutliers(s', 'nearest', 'movmedian', 3 * w);
x = smoothdata(x, 'movmean', w);
data = ((x ./ bf' - 1)) * bf(1);
[~, ~, ~, ~, ~, ~, seg] = analysis({s'}, {b'}, tVec, w, SP, 0, 0);
warning(ws);
end
