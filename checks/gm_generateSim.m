function [signals, backgroundSignals, t] = gm_generateSim(numTracks, type, expConst, noiseStd, lambda, fixBleach, negFrac, k0, refractory, pBurst)
% GM_GENERATESIM  generateSim.m with switches for the issues being tested.
%
%   fixBleach  - false (default): baseline NOT bleached, as in generateSim.m
%                true: whole signal bleached, (base + bursts).*exp(-k t),
%                which matches the model assumed in analysis.m
%   negFrac    - fraction of negative bursts (generateSim.m uses 0.6)
%   k0         - intensity per unit (generateSim.m uses 60)
%   refractory - minimum frames between bursts (generateSim.m uses 18)
%   pBurst     - burst probability per frame (generateSim.m uses 0.025)
%
% type: 1 = constant, 3 = bursts (type 2 is in the original generateSim.m)

if nargin < 6  || isempty(fixBleach),  fixBleach  = false; end
if nargin < 7  || isempty(negFrac),    negFrac    = 0.6;   end
if nargin < 8  || isempty(k0),         k0         = 60;    end
if nargin < 9  || isempty(refractory), refractory = 18;    end
if nargin < 10 || isempty(pBurst),     pBurst     = 0.025; end

T = 360;
t = linspace(0, 60, T);

backgroundSignals = randi([1000,1500], numTracks, 1) .* ones(numTracks, T);
backgroundSignals = backgroundSignals .* exp(-expConst*t);
backgroundSignals = backgroundSignals + noiseStd*randn(numTracks, T);

signals = zeros(numTracks, T);
switch type
    case 1
        for i = 1:numTracks
            signals(i,:) = randi([2000,10000]) .* exp(-expConst*t) + noiseStd*randn(1, T);
        end
    case 3
        for i = 1:numTracks
            pos = -10; level = 0; trace = zeros(1, T);
            for j = 1:T
                if rand < pBurst && abs(j - pos) > refractory
                    pos = j;
                    change = k0 * poissrnd(lambda);
                    if rand < negFrac
                        change = -change;
                    end
                    level = level + change;
                end
                trace(j) = level;
            end
            base = randi([2000,10000]);
            if fixBleach
                signals(i,:) = (base + trace) .* exp(-expConst*t) + noiseStd*randn(1, T);
            else
                signals(i,:) = base + trace .* exp(-expConst*t) + noiseStd*randn(1, T);
            end
        end
    otherwise
        error('type must be 1 or 3');
end
end
