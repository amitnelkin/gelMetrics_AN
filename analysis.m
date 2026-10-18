% Analysis function constucts the signal from the measured spot signal and
% measured background signal. It then locates activity bursts in the signal
% and classifies them to positive, negative or unclassified.

function [timeBetweenPos,timeBetweenNeg,ampPos,ampNeg,ampVar,eventsPerSignal] = ...
    analysis(signals,background,timeVecTotal,movingAverageWindowSpan,...
    stringencyParameter,plotMin,plotMax)


% Initialize variables for data aggregation
timeBetweenPos = -Inf;
timeBetweenNeg = -Inf;

ampPos = -Inf;
ampNeg = -Inf;
ampVar = -Inf;

eventsPerSignal = [];

for i = 1:size(signals,1)
    % Experimental data is in a cell array due to signals might not having
    % the same length (e.g., a spot could be tracked for less than 60 min)
    if isa(signals,'cell')
        currentSignal = signals{i,:};
        currentBackground = background{i,:};
    else
        currentSignal = signals(i,:);
        currentBackground = background(i,:);
    end

    if iscolumn(currentSignal)
        currentSignal = currentSignal';
    end
    if iscolumn(currentBackground)
        currentBackground = currentBackground';
    end
    
    % Discard experimental signals with less than 180 frames (30 min at 10 s/frame)
    if length(currentSignal) < 180
        continue;
    end
    
    
    timeVec = timeVecTotal(1:length(currentSignal));
    
    posSegmentStart = [];
    negSegmentStart = [];
    
    % Fit background signal to a 3rd degree polynomial to capture behavior
    [f,gof,output] = fit(timeVec',currentBackground','poly3');
    currentBackground = feval(f,timeVec);
    % Deal with outliers (mostly results of experimental noise)
    currentSignal = filloutliers(currentSignal,'nearest','movmedian',3*movingAverageWindowSpan); % local outliers only (global mean could clip real plateaus)
    % Apply moving average filter on the signal
    currentSignal = smoothdata(currentSignal, 'movmean', movingAverageWindowSpan);
    % Remove background and denoise from photobleaching component
    data = ((currentSignal./currentBackground'-ones(size(currentSignal)))) * currentBackground(1);

    [startPosValid,lengthsValid,thresholdPos,thresholdNeg] = getPosition(data',stringencyParameter);
    
    if ~isempty(startPosValid)
        eventsPerSignal = [eventsPerSignal,length(startPosValid)];
    end
    
    % Locate bursts in the signal
    % The concept is to fit the signal segments to linear lines and extract
    % the relevant observables
    for tt = 1:length(startPosValid)
        x = linspace(startPosValid(tt),startPosValid(tt)+lengthsValid(tt),...
            lengthsValid(tt)+1);
        if x(end)>length(data)
            x = x(x<length(data));
            p{tt} = polyfit(x',data(startPosValid(tt):length(data)-1),1);
        else
            p{tt} = polyfit(x,data(startPosValid(tt):startPosValid(tt)+lengthsValid(tt)),1); 
        end
        
        y = polyval(p{tt},x);
        if y(end)-y(1) >= 0 && length(x) >= thresholdPos
            ampPos = [ampPos; y(end)-y(1)];
            posSegmentStart = [posSegmentStart; x(1)];
        end
        if y(end)-y(1) < 0 && length(x) >= thresholdNeg
            ampNeg = [ampNeg; y(end)-y(1)];
            negSegmentStart = [negSegmentStart; x(1)];
        end
    end
        
    % find segments that cannot be classified as bursts
    for tt = 1:length(startPosValid)-1
        xVar = linspace(startPosValid(tt)+lengthsValid(tt)+1,startPosValid(tt+1)-1,...
            startPosValid(tt+1) - startPosValid(tt) - lengthsValid(tt)-1);
        if isempty(xVar)
            continue;
        end
        if xVar(end)>length(data)
            xVar = xVar(1:end-1);
            pVar{tt} = polyfit(xVar,data(startPosValid(tt)+lengthsValid(tt)+1:startPosValid(tt+1)-2),1);
        else
            pVar{tt} = polyfit(xVar,data(startPosValid(tt)+lengthsValid(tt)+1:startPosValid(tt+1)-1),1);
        end
            yVar = polyval(pVar{tt},xVar);
            ampVar = [ampVar; (yVar(end)-yVar(1))];
    end
        
        timeBetweenPos = [timeBetweenPos; (timeVec(2)-timeVec(1))*diff(posSegmentStart)];
        timeBetweenNeg = [timeBetweenNeg; (timeVec(2)-timeVec(1))*diff(negSegmentStart)];
        
        % segment classification over 
        if ~exist('p','var')
            p = [];
        end
        if ~exist('pVar','var')
            pVar = [];
        end
        
        % Plot signals according to input variables
        if plotMin < plotMax
            plotTracksCustom(timeVec, data, plotMin, startPosValid,...
                lengthsValid, p, pVar, [])
        end
        plotMin = plotMin + 1;
end