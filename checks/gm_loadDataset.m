function [signals, backgrounds] = gm_loadDataset(repoDir, folder, name, minLen)
% GM_LOADDATASET  Loads a <name>_data.mat / <name>_background.mat pair the
% same way ExperimentalAnalysis.m does (drops the empty first cell).
% If minLen > 0, keeps only traces with at least minLen frames.

if nargin < 4, minLen = 0; end
S = struct2cell(load(fullfile(repoDir, folder, [name '_data.mat'])));
B = struct2cell(load(fullfile(repoDir, folder, [name '_background.mat'])));
S = S{1}; B = B{1};
S = S(2:end); B = B(2:end);
S = S(:); B = B(:);
if minLen > 0
    keep = cellfun(@numel, S) >= minLen;
    S = S(keep); B = B(keep);
end
signals = S; backgrounds = B;
end
