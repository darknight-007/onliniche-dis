function gulpIndices = ecohabDoradoSamplingPolicySyntheticNbyCustom(DIVISOR, auvdataRaw,LgridVec)
%gulps = [1051 3842 6592 10086 12181 14938 17690 20468 23257 26153]

%auvdataRaw = analyzeAUVctdDataForEnvFeatures('~/Downloads/Dorado389_2013_074_02_074_02.mat');
%auvdataRaw = analyzeAUVctdDataForEnvFeatures('~/Downloads/Dorado389_2013_075_05_075_06.mat');
%auvdataRaw = analyzeAUVctdDataForEnvFeatures('~/Downloads/Dorado389_2013_076_01_076_02.mat');

NO_GULPERS = 9;
shouldPlot  = 0;
DEPTH_THRESH = 2;
highThresh = 1;
TOTALTICKS = length(auvdataRaw);
NO_WINDOWS = 9;
GULPS_PER_WINDOW =  floor(NO_GULPERS/NO_WINDOWS)
GULP_WINDOW_SIZE = floor(TOTALTICKS/NO_WINDOWS)
eWINDOW = floor(GULP_WINDOW_SIZE/(DIVISOR))



fl = auvdataRaw(:,6);
tempVec = auvdataRaw(:,2);
salinityVec  = auvdataRaw(:,3);
nitrateVec = auvdataRaw(:,4);
xVec  = auvdataRaw(:,10);
yVec  = auvdataRaw(:,9);
zVec  = -auvdataRaw(:,8);

chlLim = [0 30]
load hs2fit;
load jetplus;



mSP = LgridVec;

if(shouldPlot)
    subplot(122)
    hist(mSP)
    % a = vline(0.4,'r')
    % set(a,'LineWidth',2)
    xlabel('score','FontSize',13);
    ylabel('Num. samples','FontSize',13);
    % get these numbers from ROC analysis.
    % One analysis for bloom (high-thresh), another for non-bloom (low thresh).
    % With threshold gap being the margin.
    
    subplot(121)
    hist(chlVec)
    % a = vline(0.4,'r')
    % set(a,'LineWidth',2)
    xlabel('Chl. conc','FontSize',13);
    ylabel('Num. samples','FontSize',13);
end

% figure;
% subplot(131)
% histfit(salinityVec(LgridVec>0.95),10)
%
% subplot(132)
% histfit(tempVec(LgridVec>0.95),10)
%
% subplot(133)
% histfit(chlVec(LgridVec>0.95),10)

gulpDetails = [];
gulpIndices = [];
gulpDetailsGlobal = [];
candidates = [];
wCtr = 1;
gulpCtr = 0;
maxScore = 0;
okToGulp = 1;
gulpInWindowCtr = 0;
corr = 0;
haveGulpsForWindow = 1;

for i=1:length(mSP)
    if(gulpCtr < NO_GULPERS)
        testVector = [auvdataRaw(i,[2,3]) chlVec(i)];
        
        
        % 0.0105
        if(wCtr < eWINDOW)
            if(mSP(i) > maxScore)
                maxScore = min(mSP(i),highThresh);
            end
            okToGulp = 0;
        else
            okToGulp =1;
        end
        
        if (gulpInWindowCtr < GULPS_PER_WINDOW)
            haveGulpsForWindow = 1;
        else
            haveGulpsForWindow = 0;
        end
        
        if(wCtr < GULP_WINDOW_SIZE) % within gulpable window
            wCtr = wCtr+1;
            if(gulpCtr == 9 && wCtr == 5552)
                ['break here']
            end
        else % end of gulpable window
            wCtr = 1;
            maxScore = 0;
            if(haveGulpsForWindow == 1 )
                gulpDetailsGlobal = [gulpDetailsGlobal; testVector  mSP(i) corr 0];
                gulpIndices = [gulpIndices ; i 0];
                [i 0]
                
                gulpCtr = gulpCtr +1
                gulpDetails = [];
                [0]
            end
            gulpInWindowCtr = 0;
            haveGulpsForWindow = 0;
        end
        
        
  %[isCorrOk, corr] = isDecorrelated(testVector,mSP(i),gulpDetails,CORR_CUTOFF);
        isCorrOk = 1;
        corr = 99;
        isDepthOK = auvdataRaw(i,8) >=DEPTH_THRESH;
        isThresholdOK = (mSP(i) > maxScore);
        %  [i wCtr gulpCtr gulpInWindowCtr  mSP(i) maxScore okToGulp isThresholdOK isCorrOk  isDepthOK]
        candidates = [candidates; isThresholdOK isCorrOk isDepthOK mSP(i) corr];
        
        if(haveGulpsForWindow==1 && okToGulp==1 &&  (isThresholdOK ==1) && (isCorrOk == 1) && (isDepthOK == 1)) % OK to gulp
            gulpDetails = [gulpDetails; testVector  mSP(i) corr];
            gulpDetailsGlobal = [gulpDetailsGlobal; testVector  mSP(i) corr 1];
            gulpIndices = [gulpIndices ; i 1];
            gulpCtr = gulpCtr +1;
            gulpInWindowCtr = gulpInWindowCtr + 1;
            [i 1 gulpInWindowCtr]
        end
    end
end
%
% histogramNormalizedJD(5,NO_GULPERS,gulpDetailsGlobal(:,3),'Phytoplankton abundance (OD)','Gulp count');
% %xlim([0.06     0.15])
% ylim([0 NO_GULPERS])

if(shouldPlot==1)
    
    figure('Position',[167         130        1380         804])
    subplot(611)
    scatter(1:length(xVec),zVec,40,tempVec,'.')
    xlabel('Experiment time (hour)');
    ylim([ -70 1]);
    ylabel('Depth (m)');
    hold on;plot(gulpIndices(:,1),zVec(gulpIndices(:,1)),'k+','markerSize',20)
    plot(gulpIndices(:,1),zVec(gulpIndices(:,1)),'ko','markerSize',20)
    vline(eWINDOW:GULP_WINDOW_SIZE:TOTALTICKS)
    hold on;vline(1:GULP_WINDOW_SIZE:TOTALTICKS,'r');
    c = colorbar;
    xlabel(c,'Temperature')
    
    
    
    subplot(612)
    scatter(1:length(xVec),zVec,40,salinityVec,'.')
    ylim([ -70 1]);
    xlabel('Experiment time (hour)');
    ylabel('Depth (m)');
    hold on;plot(gulpIndices(:,1),zVec(gulpIndices(:,1)),'k+','markerSize',20)
    plot(gulpIndices(:,1),zVec(gulpIndices(:,1)),'ko','markerSize',20)
    vline(eWINDOW:GULP_WINDOW_SIZE:TOTALTICKS)
    hold on;vline(1:GULP_WINDOW_SIZE:TOTALTICKS,'r');
    c = colorbar;
    xlabel(c,'Salinity')
    
    
    
    
    
    
    subplot(613)
    scatter(1:length(xVec),zVec,40,chlVec,'.')
    ylim([ -70 1]);
    xlabel('Experiment time (hour)');
    ylabel('Depth (m)');
    hold on;plot(gulpIndices(:,1),zVec(gulpIndices(:,1)),'k+','markerSize',20)
    plot(gulpIndices(:,1),zVec(gulpIndices(:,1)),'ko','markerSize',20)
    caxis(chlLim)
    vline(eWINDOW:GULP_WINDOW_SIZE:TOTALTICKS)
    hold on;vline(1:GULP_WINDOW_SIZE:TOTALTICKS,'r');
    c = colorbar;
    xlabel(c,'chl')
    
    
    subplot(614)
    scatter(1:length(xVec),zVec,40,LgridVec,'.')
    ylim([ -70 1]);
    xlabel('Experiment time (hour)');
    ylabel('Depth (m)');
    hold on;plot(gulpIndices(:,1),zVec(gulpIndices(:,1)),'k+','markerSize',20)
    plot(gulpIndices(:,1),zVec(gulpIndices(:,1)),'ko','markerSize',20)
    caxis([0.5 1])
    vline(eWINDOW:GULP_WINDOW_SIZE:TOTALTICKS)
    hold on;vline(1:GULP_WINDOW_SIZE:TOTALTICKS,'r');
    c = colorbar;
    xlabel(c,'bloom score')
    
    
    subplot(615)
    plot(1:length(xVec),mSP)
    %zlim([ -70 1]);
    %caxis([0 15])
    xlabel('Experiment time (hour)');
    ylabel('Bloom score');
    hold on;plot(gulpIndices(:,1),mSP(gulpIndices(:,1)),'k+','markerSize',20,'displayName','niche based')
    plot(gulpIndices(:,1),mSP(gulpIndices(:,1)),'ko','markerSize',20)
    ylim([0 0.6]);
    vline(eWINDOW:GULP_WINDOW_SIZE:TOTALTICKS)
    hold on;vline(1:GULP_WINDOW_SIZE:TOTALTICKS,'r');
    colorbar
    
    subplot(616)
    plot(1:length(xVec),chlVec)
    %zlim([ -70 1]);
    %caxis([0 15])
    xlabel('Experiment time (hour)');
    ylabel('Chl ');
    hold on;plot(gulpIndices(:,1),chlVec(gulpIndices(:,1)),'k+','markerSize',20,'displayName','niche based')
    plot(gulpIndices(:,1),chlVec(gulpIndices(:,1)),'ko','markerSize',20)
    ylim(chlLim);
    vline(eWINDOW:GULP_WINDOW_SIZE:TOTALTICKS)
    hold on;vline(1:GULP_WINDOW_SIZE:TOTALTICKS,'r');
    colorbar
    
    figure;
    scatter3(salinityVec,tempVec, chlVec,20,LgridVec,'.')
    hold on;plot3(salinityVec(gulpIndices(:,1)),tempVec(gulpIndices(:,1)),chlVec(gulpIndices(:,1)),'kp-','LineWidth',2)
    ind = find(gulpIndices(:,2)<1);
    %hold on;plot3(salinityVec(gulpIndices(ind,1)),tempVec(gulpIndices(ind,1)),chlVec(gulpIndices(ind,1)),'rp','MarkerSize',20)
    
    hold on;plot3(salinityVec(gulpIndices(ind,1)),tempVec(gulpIndices(ind,1)),chlVec(gulpIndices(ind,1)),'k*','MarkerSize',20)
    
    xlim([min(salinityVec) max(salinityVec)])
    ylim([min(tempVec) max(tempVec)])
    zlim([min(chlVec) max(chlVec)])
    zlim(chlLim)
    text(salinityVec(gulpIndices(1,1)),tempVec(gulpIndices(1,1)),chlVec(gulpIndices(1,1)),'Gulp 1','FontSize',15,'FontWeight','bold');
    %drawSphere(figGulps, CORR_CUTOFF, salinityVec(gulpIndices),tempVec(gulpIndices),chlVec(gulpIndices))
    
end
end

function [isDecorrelated_,minDistance] = isDecorrelated(auvdataRaw,mSP,gulpDetails,distThresh)
isDecorrelated_ = 0;
minDistance = minDistanceFromExistingDataPoints([auvdataRaw],gulpDetails);
if(isempty(gulpDetails) || minDistance > distThresh)
    isDecorrelated_ = 1;
end
end


function [minDistance] = minDistanceFromExistingDataPoints(testVector,gulpDetails)
minDistance = 9999;
for i=1:size(gulpDetails,1)
    dist = envDistance(testVector,gulpDetails(i,1:3));
    if(dist < minDistance)
        minDistance = dist;
    end
end
end

function dist_ = envDistance(x,z)
%sq_dist(normalizeEnvDataDoradoTSC([10 33.4 0]'), normalizeEnvDataDoradoTSC([14 33.9 15]'))
dist_ = sq_dist(normalizeEnvDataDoradoTSC(x'),normalizeEnvDataDoradoTSC(z'));
end

function drawSphere(fig, radius, xCent,yCent,zCent)
figure(fig);
% function [result] = normalizeEnvDataDoradoTSC(x)
% normVec(1) = (x(1)-10)/8;
% normVec(2) = (x(2)-32)/2;
% normVec(3) = (x(1))/20;
% result = normVec';
% end

hold on;
[x,y,z] = sphere;
for i=1:length(xCent)
    surf((x*radius*2)+xCent(i),y*radius*8+yCent(i),(z*radius*20)+zCent(i));
    shading flat;
    alpha(0.2);
end
end
