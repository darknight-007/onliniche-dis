function gulpIndices = ecohabDoradoSamplingPolicySyntheticRandom(dummyDivisorUnused,auvdataRaw,score)
%gulps = [1051 3842 6592 10086 12181 14938 17690 20468 23257 26153]

%auvdataRaw = analyzeAUVctdDataForEnvFeatures('~/Downloads/Dorado389_2013_074_02_074_02.mat');
%auvdataRaw = analyzeAUVctdDataForEnvFeatures('~/Downloads/Dorado389_2013_075_05_075_06.mat');
%auvdataRaw = analyzeAUVctdDataForEnvFeatures('~/Downloads/Dorado389_2013_076_01_076_02.mat');

NO_GULPERS = 9;
TOTALTICKS = length(auvdataRaw);
WINDOWSIZE = floor(TOTALTICKS/9);
gulpIndices = [];
for i = 1:NO_GULPERS
    ind = randi(WINDOWSIZE,1);
    ind = ind + (i-1)*WINDOWSIZE;
    gulpIndices = [gulpIndices; ind 1];
end



end
