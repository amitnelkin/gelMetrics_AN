% LAMBDA_K0_PRODUCT_RECOVERY  How often are lambda, and lambda*k0, recovered? (RESULTS.md section 89)
% Same simulation as the Fig. 16 / Sec. II.D numbers (AppendixFigures.m Appendix G settings, same seeds):
% 200 tracks, k0 = 60, SP 15, window 13, true lambda 1-5, 20 sets each, WS fit (PoissonFit2WS, lambda 1-5).
% Noise 40 a.u. (one unit / noise = 1.5, the setting of the MS sentence), plus 12 (SNR 5) and 120 (SNR 0.5).
% "Recovered" = within +-20% of the true value (for lambda on the 1-5 grid this is an exact hit).
% Estimates of lambda*k0: (1) the WS fit's lambda times its k0; (2) the mean detected burst size.
% Output: analysis_outputs/lambda_k0_product/recovery.csv
timeVecTotal = linspace(0, 60, 361);
lambdaSim = [1 2 3 4 5]; noiseList = [12 40 120]; k0 = 60;
nSets = 20; lambdaGrid = 1:5; K0Vector = 1:500; expConst = 0.0026;
rows = {};
for ni = 1:numel(noiseList)
    for l = 1:numel(lambdaSim)
        for r = 1:nSets
            rng(5000 + 100 * (l + 1) + r);      % the same seeds as the 31% run (its lambda list started at 0.5)
            [signals, backgrounds] = generateSim(200, 361, 3, expConst, noiseList(ni), lambdaSim(l), k0);
            [~, ~, ~, ampNeg] = analysis(signals, backgrounds, timeVecTotal, 13, 15, 0, 0);
            amps = abs(ampNeg(2:end));
            [lamWS, k0WS] = PoissonFit2WS(amps, lambdaGrid, K0Vector, 10);
            rows(end+1, :) = {noiseList(ni), k0 / noiseList(ni), lambdaSim(l), r, lamWS, k0WS, lamWS * k0WS, mean(amps), numel(amps)}; %#ok<SAGROW>
        end
        fprintf('noise %g lambda %g done\n', noiseList(ni), lambdaSim(l));
    end
end
T = cell2table(rows, 'VariableNames', {'noise', 'snr', 'lambdaTrue', 'set', 'lambdaWS', 'k0WS', 'productWS', 'meanBurst', 'nBursts'});
if ~exist('analysis_outputs/lambda_k0_product', 'dir'), mkdir('analysis_outputs/lambda_k0_product'); end
writetable(T, 'analysis_outputs/lambda_k0_product/recovery.csv');
