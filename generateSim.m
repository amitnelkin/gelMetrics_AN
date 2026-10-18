% generateSim function generates simulated trajectories
% Function recieces as input: 
% numTracks - number of tracks to generate
% trackLength - the number of points in a track
% type - 1 for constant, 2 for slope, 3 for bursts
% expConst - constant for photobleaching exponential component
% noiseStd - standard deviation of the gaussian noise component
% lambda - Poisson parameter for the bursty signals

function [signals, backgroundSignals] = generateSim(numTracks, trackLength, type, expConst, noiseStd, lambda)
    % Initialize the time vector
    t = linspace(0, 60, trackLength)';
    
    backgroundSignals = randi([1000,1500],numTracks,1) .* ones(numTracks,trackLength); % Base background signal
    backgroundSignals = backgroundSignals .* exp(-expConst*t'); % Add exponential component
    backgroundSignals = backgroundSignals + noiseStd*randn(numTracks,trackLength); % Add Gaussian noise component

    % Initialize the signals matrix
    signals = zeros(numTracks, trackLength);
    
    % Generate signals based on the specified type
    switch type
        case 1
            % Constant signals with exponential component and gaussian noise
            for i = 1:numTracks
                signals(i,:) = randi([2000,10000]) .* exp(-expConst * t') + noiseStd * randn(1, trackLength);
            end
        case 2
            % Sloped signals with exponential component and gaussian noise
            for i = 1:numTracks
                slope = randi(100);
                intercept = randi([2000,10000]);
                signals(i,:) = slope * t + intercept;
                signals(i,:) = signals(i,:) .* exp(-expConst * t') + noiseStd * randn(1, trackLength);
            end 
        case 3
            % Intermittent bursts signals with exponential component and gaussian noise
            k = 60; % Assuming a constant k
            for i = 1:numTracks
                currentBurstPos = -10;
                currentLevel = 0;
                for j = 1:trackLength
                    distFromBurst = abs(j-currentBurstPos);
                    if rand < 0.025 && distFromBurst > 18 % 2.5% chance of a burst at each point
                        currentBurstPos = j; 
                        change = k * poissrnd(lambda);
                        % We noticed a positive bias in our simulated signals
                        % we counter this by giving a higher chance for a
                        % negative burst which reduces signal intensity
                        if rand < 0.5
                            change = -change; % symmetric: 50% chance of being negative
                        end
                        currentLevel = currentLevel + change;
                    end
                    signals(i, j) = currentLevel;
                end
                signals(i, :) = (randi([2000,10000]) + signals(i, :)) .* exp(-expConst * t') + noiseStd * randn(1, trackLength); % whole signal bleaches, matching analysis.m
                
            end
            
        otherwise
            error('Invalid type. Type must be 1, 2, or 3.');
    end
end
