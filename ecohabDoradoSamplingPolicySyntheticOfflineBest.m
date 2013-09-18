function gulpIndices = ecohabDoradoSamplingPolicySyntheticOfflineBest(dummyDivisorUnused, auvdataRaw,score)
%gulps = [1051 3842 6592 10086 12181 14938 17690 20468 23257 26153]
DEPTH_THRESH = 2;
%DIVISOR = 2;
TOTALTICKS = length(auvdataRaw);
NO_WINDOWS = 9;
GULP_WINDOW_SIZE = floor(TOTALTICKS/NO_WINDOWS);

indices = 1:GULP_WINDOW_SIZE:GULP_WINDOW_SIZE*(NO_WINDOWS+1);
indices(end) = min(TOTALTICKS, indices(end));

gulpIndices = [];
validInd = [];
for i=1:length(indices)-1
    startInd = indices(i);
    endInd = indices(i+1);
    indInWindow = startInd:endInd;
    inDepthInd = indInWindow(find(auvdataRaw(indInWindow,8) > DEPTH_THRESH));
    [val,indMax] = max(score(inDepthInd));
    gulpIndices = [gulpIndices ;inDepthInd(indMax) 1];
end
% if(length(gulpIndices) <9)
%     keyboard
% end
end
