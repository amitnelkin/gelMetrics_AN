function T = compareConditions_ref(MA, MB, nameA, nameB, nPerm, seed)
% COMPARECONDITIONS_REF  (the function version of compareConditions.m, kept for checks/replicate_test.m and PaperFigures_revised.m)  Is condition A different from condition B?
%   T = compareConditions(MA, MB, nameA, nameB)              MA, MB from burstMeasures.m
%   T = compareConditions(MA, MB, nameA, nameB, nPerm, seed) (defaults 10000 shuffles, seed 1)
% For each measure: A, B, the difference A - B, and p values from shuffling whole
% movies between the two conditions (movies are the independent replicates):
%   p        two-sided
%   p_AgtB   one-sided, for A > B (use it only when that direction was stated in advance)
% Prints the table and returns it.
if nargin < 5, nPerm = 10000; end
if nargin < 6, seed = 1; end
rng(seed);
measures = {'sizeGranules', 'sizeCopies', 'sizeBeyond', 'turnover', 'net', 'medianBurst', 'meanBurst', 'brightness', 'brightnessPooled'};
nM = numel(measures);

% put both conditions into one list of granules, keeping each movie's granules together
M = joinConditions(MA, MB);
moviesA = unique(MA.movie); moviesB = unique(MB.movie);
nA = numel(moviesA); nTotal = nA + numel(moviesB);

% observed values
VA = summarizeGranules(MA, 1:MA.nGranules);
VB = summarizeGranules(MB, 1:MB.nGranules);
valueA = zeros(nM, 1); valueB = zeros(nM, 1);
for j = 1:nM
    valueA(j) = VA.(measures{j}); valueB(j) = VB.(measures{j});
end
observed = valueA - valueB;

% shuffle whole movies between the two conditions
countTwo = zeros(nM, 1); countOne = zeros(nM, 1);
for p = 1:nPerm
    order = randperm(nTotal);
    picksA = find(ismember(M.movie, order(1:nA)));          % granules of the movies now called "A"
    picksB = find(ismember(M.movie, order(nA+1:end)));
    V1 = summarizeGranules(M, picksA); V2 = summarizeGranules(M, picksB);
    for j = 1:nM
        d = V1.(measures{j}) - V2.(measures{j});
        if abs(d) >= abs(observed(j)) - 1e-12, countTwo(j) = countTwo(j) + 1; end
        if d >= observed(j) - 1e-12, countOne(j) = countOne(j) + 1; end
    end
end
p = (countTwo + 1) / (nPerm + 1);
p_AgtB = (countOne + 1) / (nPerm + 1);
T = table(measures', valueA, valueB, observed, p, p_AgtB, ...
    'VariableNames', {'measure', 'A', 'B', 'difference', 'p', 'p_AgtB'});
fprintf('\n%s (A) vs %s (B): %d vs %d movies, p from %d movie shuffles\n', nameA, nameB, nA, nTotal - nA, nPerm);
disp(T);
end


function M = joinConditions(MA, MB)
% one combined list of granules; movies renumbered 1..nA (condition A) and nA+1.. (condition B)
[~, ~, movieA] = unique(MA.movie); [~, ~, movieB] = unique(MB.movie);
nA = max(movieA);
M.movie = [movieA(:); movieB(:) + nA];
M.nGranules = MA.nGranules + MB.nGranules;
M.nCopies = MA.nCopies;
if MA.nCopies ~= MB.nCopies, error('Use the same number of copies for both conditions.'); end
perGranule = {'level', 'noise', 'hours', 'hours1', 'hours2', 'granuleNet', 'granuleNet1', 'granuleNet2'};
for j = 1:numel(perGranule)
    M.(perGranule{j}) = [MA.(perGranule{j}); MB.(perGranule{j})];
end
lists = {'ampAU', 'fracReal', 'startReal', 'fracCopy', 'startCopy'};    % one list of bursts per granule
for j = 1:numel(lists)
    M.(lists{j}) = [MA.(lists{j}); MB.(lists{j})];
end
end
