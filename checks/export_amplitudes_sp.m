% EXPORT_AMPLITUDES_SP  Negative and positive burst amplitudes of the PP7
% datasets at stringency parameters SP = 10..15, no fitting. For recreating
% MS Fig. 16 (PP7 vs PP7-FUS histograms vs SP) with the re-extracted data.
% Reviewer points R2 4a (threshold dependence) and the FUS claim.
% Data: R2 re-extraction (in_vivo_reextracted/after_R2_offset: red channel,
% multithresh spot, offset 200) and, for reference, the stored MS data.
% Same code path as the runners: analysis.m per trace (gm_amplitudes_by_trace),
% window 13, traces >= 180 frames.
% Writes analysis_outputs/fig16_sp/amplitudes_sp.csv (data, dataset, SP, sign, trace, amp).
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'fig16_sp');
if ~exist(outDir, 'dir'), mkdir(outDir); end
src = {'R2_reextracted', 'in_vivo_reextracted/after_R2_offset'; 'stored_MS', 'in_vivo_experiments-old'};
names = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS'};
tVec = (0:360) * 10/60; win = 13;
C = cell(0, 6);
for v = 1:size(src,1)
    for d = 1:numel(names)
        [S, B] = gm_loadDataset(repoDir, src{v,2}, names{d});
        for SP = 10:15
            T = gm_amplitudes_by_trace(S, B, tVec, win, SP);
            T = T([T.len] >= 180);
            for i = 1:numel(T)
                for sg = {'Neg', 'Pos'}
                    a = abs(T(i).(['amp' sg{1}]));
                    for k = 1:numel(a)
                        C(end+1,:) = {src{v,1}, names{d}, SP, sg{1}, i, a(k)}; %#ok<SAGROW>
                    end
                end
            end
            fprintf('%s %s SP=%d done\n', src{v,1}, names{d}, SP);
        end
    end
end
writetable(cell2table(C, 'VariableNames', {'data','dataset','SP','sign','trace','amp'}), ...
    fullfile(outDir, 'amplitudes_sp.csv'));
disp('EXPORT_SP_DONE');
