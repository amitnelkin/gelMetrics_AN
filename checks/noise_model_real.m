function noise_model_real(part)
% NOISE_MODEL_REAL  Real data with the noise model chosen on simulations
% (checks/noise_model_validation.m): 'mixH+free' = burst scatter fitted from the
% bursts + a false-burst component whose width is fixed at h0K*sigma
% (calibrated on burst-free simulations; a free width swallowed real bursts in vitro). R2 4b (lambda), R2 6 (noise floor), Qbeta
% vs PP7, valency.
% part 'A': each dataset separately, paper SP and truncation; mix+free with
%   100-resample granule bootstrap (negative bursts); point estimates for
%   floor+free (noise floor), mix+free (free width) and the old fixed-sigma fit.
% part 'B': valency-linked joint fits (gelmetricsJoint, opts.mixFree) per
%   protein group at one SP per group (as checks/valency_real.m), mixH+free, 'fixed'
%   (full occupancy) and 'bounded' (at most full); 30 bootstrap resamples.
% Data: in vitro stored; PP7, PP7-FUS, Qbeta from the R2 re-extraction.
% Writes analysis_outputs/noise_models/real_separate.csv / real_joint.csv.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'noise_models');
C = readtable(fullfile(outDir, 'calibration.csv'));
R2 = 'in_vivo_reextracted/after_R2_offset'; IV = 'in_vitro_experiments';
tVec = (0:360) * 10/60; rows = {};
if strcmp(part, 'A')
    D = {'protein_RNA_1_1',IV,12,15; 'protein_RNA_10_1',IV,12,10; 'protein_RNA_100_1',IV,12,10;
         '4PP7',R2,10,15; '4PP7-FUS',R2,10,15; '30PP7',R2,12,10; '30PP7-FUS',R2,12,10; '5Qb',R2,14,15; '10Qb',R2,14,15};
    for d = 1:size(D,1)
        [S, B] = gm_loadDataset(repoDir, D{d,2}, D{d,1});
        T = gm_amplitudes_by_trace(S, B, tVec, 13, D{d,3}); T = T([T.len] >= 180);
        sg = gm_sigma_from_ampvar(vertcat(T.ampVar)); A = abs(vertcat(T.ampNeg));
        fk = C.floorK(C.SP == D{d,3}); hk = C.h0K(C.SP == D{d,3});
        o = gm_fit_dataset(T, struct('sign','Neg','truncPercent',D{d,4},'nBoot',100,'seed',d,'grid',false, ...
            'noise','free','falseBursts','mix','h0K',hk));
        m = o.mle;
        rows(end+1,:) = {D{d,1}, 'mixH+free', m.n, sg, m.lambda, o.boot.lamCI(1), o.boot.lamCI(2), m.k0, ...
            o.boot.k0CI(1), o.boot.k0CI(2), m.noiseParam, m.w0, m.h0}; %#ok<AGROW>
        others = {'floor+free', struct('noise','free','floorK',fk);
                  'mix+free', struct('noise','free','falseBursts','mix');
                  'fixed (old)', struct('noise','fixed')};
        for k = 1:size(others,1)
            oo = others{k,2}; oo.truncPercent = D{d,4}; oo.profile = false;
            Rk = gelmetricsMLE(A, sg, oo);
            w = NaN; h = NaN; if isfield(Rk, 'w0'), w = Rk.w0; h = Rk.h0; end
            rows(end+1,:) = {D{d,1}, others{k,1}, Rk.n, sg, Rk.lambda, NaN, NaN, Rk.k0, NaN, NaN, Rk.noiseParam, w, h}; %#ok<AGROW>
        end
        fprintf('%s: mixH+free lam %.2f [%.2f-%.2f] k0 %.0f [%.0f-%.0f] sigmaFit %.0f (ampVar %.0f) w0 %.2f\n', D{d,1}, ...
            m.lambda, o.boot.lamCI, m.k0, o.boot.k0CI, m.noiseParam, sg, m.w0);
    end
    writetable(cell2table(rows, 'VariableNames', {'dataset','model','nFit','sigmaAmpVar','lambda','lamLo','lamHi', ...
        'k0','k0Lo','k0Hi','noiseParam','w0','h0'}), fullfile(outDir, 'real_separate.csv'));
else
    G = {'Qb', {'5Qb','10Qb'}, R2, 14, [1 2];
         'PP7', {'4PP7','30PP7'}, R2, 12, [1 7.5];
         'PP7-FUS', {'4PP7-FUS','30PP7-FUS'}, R2, 12, [1 7.5];
         'invitro', {'protein_RNA_1_1','protein_RNA_10_1','protein_RNA_100_1'}, IV, 12, [1 10 14]};
    nBoot = 30;
    for g = 1:size(G,1)
        names = G{g,2}; m = numel(names); T = cell(1,m); A = cell(1,m); sg = zeros(1,m);
        sep = struct('lambda', {}, 'k0', {}, 'sigmaFit', {}, 'w0', {}, 'h0', {});
        for i = 1:m
            [S, B] = gm_loadDataset(repoDir, G{g,3}, names{i});
            Ti = gm_amplitudes_by_trace(S, B, tVec, 13, G{g,4}); T{i} = Ti([Ti.len] >= 180);
            A{i} = abs(vertcat(T{i}.ampNeg)); sg(i) = gm_sigma_from_ampvar(vertcat(T{i}.ampVar));
            hk = C.h0K(C.SP == G{g,4});
            Rs = gelmetricsMLE(A{i}, sg(i), struct('profile', false, 'noise', 'free', 'falseBursts', 'mix', 'h0K', hk));
            sep(i) = struct('lambda', Rs.lambda, 'k0', Rs.k0, 'sigmaFit', Rs.noiseParam, 'w0', Rs.w0, 'h0', Rs.h0);
            rows(end+1,:) = {G{g,1}, 'separate', names{i}, Rs.n, Rs.lambda, NaN, NaN, Rs.k0, NaN, NaN, NaN, Rs.nll}; %#ok<AGROW>
        end
        for md = {'fixed', 'bounded'}
            st = struct('lambda', [sep.lambda], 'k0', [sep.k0], ...
                'ratio', min(max([sep.k0] / sep(1).k0, 1), G{g,5}), ...
                'extra', reshape([log([sep.sigmaFit]); log(min([sep.w0], 0.49) ./ (1 - min([sep.w0], 0.49)))], 1, []));
            Rj = gelmetricsJoint(A, sg, G{g,5}, md{1}, struct('mixFree', true, 'h0K', hk, 'sep', sep, 'start', st));
            rng(200 + g); bl = nan(nBoot, m); bk = nan(nBoot, m);
            for b = 1:nBoot
                Ab = cell(1,m); sb = zeros(1,m);
                for i = 1:m
                    Tb = T{i}(randi(numel(T{i}), numel(T{i}), 1));
                    Ab{i} = abs(vertcat(Tb.ampNeg)); sb(i) = gm_sigma_from_ampvar(vertcat(Tb.ampVar));
                end
                try
                    Rb = gelmetricsJoint(Ab, sb, G{g,5}, md{1}, struct('mixFree', true, 'h0K', hk, 'sep', sep, 'start', Rj));
                    bl(b,:) = Rb.lambda; bk(b,:) = Rb.k0;
                catch err
                    fprintf('bootstrap %d failed: %s [%s line %d]\n', b, err.message, err.stack(1).name, err.stack(1).line);
                end
            end
            for i = 1:m
                rows(end+1,:) = {G{g,1}, md{1}, names{i}, Rj.n(i), Rj.lambda(i), prctile(bl(:,i),2.5), prctile(bl(:,i),97.5), ...
                    Rj.k0(i), prctile(bk(:,i),2.5), prctile(bk(:,i),97.5), Rj.ratio(i), Rj.nll}; %#ok<AGROW>
            end
            fprintf('%s %s: lambda %s k0 %s ratio %s nll %.1f\n', G{g,1}, md{1}, mat2str(Rj.lambda,3), mat2str(Rj.k0,3), mat2str(Rj.ratio,3), Rj.nll);
        end
    end
    writetable(cell2table(rows, 'VariableNames', {'group','fit','dataset','nFit','lambda','lamLo','lamHi', ...
        'k0','k0Lo','k0Hi','ratio','nll'}), fullfile(outDir, 'real_joint.csv'));
end
disp(['NOISE_REAL_' part '_DONE']);
end
