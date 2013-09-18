function nicheModelRealFluoroSamplingComparisonsCustomMultiCore(howMany, rootFolderName, dataFolderName)
if(isempty(rootFolderName))
    rootFolderName =  '/Users/jdas/Dropbox/Research/MBARI/Documents/dissertalkPitch/';
end
if(isempty(dataFolderName))
    dataFolderName = '/Users/jdas/Research/doradoMatFiles/';
end

% x1 = 32.6:.02:34; x2 = 12:.02:20;
% [X1_,X2_] = meshgrid(x1,x2);
% X1 = reshape(X1_,length(x1)*length(x2),1);
% X2 = reshape(X2_,length(x1)*length(x2),1);
%tTSONFBDLonLat

rankVec = [];
shouldPlot = 0;
NUM_HOTSPOTS = 1;
NO_GULPERS = 9;
CLIPSIZE = 0;
maxVec = [];
StatsAll = [];

if(shouldPlot)
    statFig = figure;
end

for howManyCtr = 1:howMany
    
    try
        Yvecs = [];
        MAX_PROBABILITY = 0.3;
        THRESH_FLUORO = 0.001;
        DEPTH_MAX = 200;
        DEPTH_MIN = 2;
        INITIAL_N = NO_GULPERS*1;
        samplingMethods = {@ecohabDoradoSamplingPolicySyntheticRandom, @ecohabDoradoSamplingPolicySyntheticOfflineBest, @ecohabDoradoSamplingPolicySyntheticNbyCustom, @ecohabDoradoSamplingPolicySyntheticNbyCustom,@ecohabDoradoSamplingPolicySyntheticNbyCustom,@ecohabDoradoSamplingPolicySyntheticNbyCustom}
        samplingMethodsParam = [0, 1, exp(1), 4,9,12]
        methodName = {'Random','Offline','N/e ',  'N/4 ', 'N/9','N/12'};
        patterns = {'b', 'c--','r', 'k:', 'm-', 'g-','k--'};
        methodLineWidth = [1,2,2,2,2,2,3]
        
        timenow = datestr(now,'YYYY-mm-DD_HHMMSS');
        newFolderName = ['logs/' timenow  '_' num2str(NUM_HOTSPOTS) '_' num2str(INITIAL_N) '_divisor_' sprintf('%.1f_%d_%d_%d',samplingMethodsParam(3:6)) '_' computer '_' version('-release') '/'];
        mkdir(newFolderName);
        DEFAULT_NICHE_MODEL_FILENAME = [newFolderName 'nicheModelInitial.mat'];
        axis_ = [ 33                     34                   10                     18];
        
        S = [];
        Y = [];
        Stats = [];
        isHandpicked = [];
        
        if(shouldPlot)
            perfFig =  figure;
            transectFig = figure;
            initialFig =  figure;
            h = figure;
        end
        
        X1 = [];
        X2 = [];
        pHeadGrid = [];
        pHeadAll = [];
        
        
        filename = [];
        
        
        fid = fopen([dataFolderName    'doradoLocoListTrimmed.txt']);
        ctr = 1;
        while ~feof(fid)
            fname = fgetl(fid);
            filename{ctr} =  [dataFolderName fname];
            ctr = ctr + 1;
        end
        Y0 = [];
        
        while(isempty(find(Y0>0)))
            ['Drawing random samples with successes']
            auvdataRaw = analyzeAUVctdDataForEnvFeatures(filename{1});
            [x1, x2,X1, X2, pHeadGrid,pHeadAll] =  mixtureofGaussianNiche(MAX_PROBABILITY, [],NUM_HOTSPOTS);
            resultAll = ((auvdataRaw(:,6)>THRESH_FLUORO))*2-1;
            candidateInd = find(auvdataRaw(:,8)<DEPTH_MAX & auvdataRaw(:,8)>DEPTH_MIN);
            candidateIndTrain =  randi(length(candidateInd),1,INITIAL_N);
            S0 = [auvdataRaw(candidateIndTrain,[3,2])];
            
            result = resultAll(candidateIndTrain);
            Y0 = result;
        end
        
        
        NUM_METHODS = length(samplingMethods);
        newportPierModelGPClassificationForSyntheticData(S0,Y0,[],2,1,DEFAULT_NICHE_MODEL_FILENAME)
        
        for methodCtr = 1:NUM_METHODS
            S{methodCtr} = S0;
            Y{methodCtr} = Y0;
            Stats{methodCtr}.proportionTargetSuccess = [];
            isHandpicked{methodCtr} = [];
            copyfile(DEFAULT_NICHE_MODEL_FILENAME,[newFolderName 'nicheModel_' num2str(methodCtr) '_1.mat'])
        end
        
        
        if(shouldPlot)
            figure(initialFig)
            
            X1vec = reshape(X1,numel(X1),1);
            X2vec = reshape(X2,numel(X1),1);
            [a b c d predictedYmesh] = newportPierModelGPClassificationForSyntheticData(S{1},Y{1},[X1vec X2vec],2,0,DEFAULT_NICHE_MODEL_FILENAME);
            
            surf(X1,X2,reshape(predictedYmesh,size(X1,1),size(X1,2)));
            xlabel('x1'); ylabel('x2'); zlabel('Probability Density');
            view(2)
            colorbar
            caxis auto
            shading interp
            
            
            alpha(0.5)
            hold on;plot(S0(result==1,1),S0(result==1,2),'r+','MarkerSize',10)
            hold on;plot(S0(result==-1,1),S0(result==-1,2),'b.')
            
            pause(2)
            
            figure(h)
            set(h,'position',[100 100 1200*1.2 600*1.2]);
            
        end
        
        CAMPAIGN_LENGTH = length(filename)-1-CLIPSIZE;
        for ctr = 2:CAMPAIGN_LENGTH
            
            auvdataRaw = analyzeAUVctdDataForEnvFeatures(filename{ctr});
            
            S_all = [auvdataRaw(:,3) auvdataRaw(:,2)];
            resultAll = ((auvdataRaw(:,6)>THRESH_FLUORO))*2-1;
            
            
            for methodCtr = 1:NUM_METHODS
                
                if(shouldPlot)
                    titleStr = [num2str(ctr) ': FILENAME=' filename{ctr}]
                    figure(h)
                    sp1 = subplot(2,NUM_METHODS,methodCtr)
                    cla
                end
                
                OLD_MODEL_FILE = [newFolderName 'nicheModel_' num2str(methodCtr) '_' num2str(ctr-1) '.mat'];
                NEW_MODEL_FILE = [newFolderName 'nicheModel_' num2str(methodCtr) '_' num2str(ctr) '.mat'];
                
                [a b c d predictedY] = newportPierModelGPClassificationForSyntheticData(S{methodCtr},Y{methodCtr},[auvdataRaw(:,3) auvdataRaw(:,2)],2,0,OLD_MODEL_FILE);
                
                
                score = predictedY;
                
                
                gulpIndices = samplingMethods{methodCtr}(samplingMethodsParam(methodCtr),auvdataRaw, score)
                
                newS = [auvdataRaw(gulpIndices(:,1),3) auvdataRaw(gulpIndices(:,1),2)];
                result = resultAll(gulpIndices(:,1));
                isHandpicked{methodCtr} = [isHandpicked{methodCtr} ; gulpIndices(:,2)];
                
                
                Stats{methodCtr}.proportionTargetSuccess = [Stats{methodCtr}.proportionTargetSuccess ; sum(gulpIndices(:,2) > 0)/NO_GULPERS];
                Y{methodCtr} = [Y{methodCtr} ; result];
                S{methodCtr} = [S{methodCtr} ; newS];
                thisS = S{methodCtr};
                thisY = Y{methodCtr};
                thisS = thisS(INITIAL_N+1:end,:);
                thisY = thisY(INITIAL_N+1:end,:);
                newportPierModelGPClassificationForSyntheticData(S{methodCtr},Y{methodCtr},[],2,1,NEW_MODEL_FILE)
                
                X1vec = reshape(X1,numel(X1),1);
                X2vec = reshape(X2,numel(X1),1);
                [a b c d predictedYmesh] = newportPierModelGPClassificationForSyntheticData(S{methodCtr},Y{methodCtr},[X1vec X2vec],2,0,NEW_MODEL_FILE);
                
                if(shouldPlot)
                    pcolor(X1,X2,reshape(predictedYmesh,size(X1,1),size(X1,2)));
                    xlabel('x1'); ylabel('x2'); zlabel('Probability Density');
                    view(2)
                    colorbar
                    caxis auto
                    shading interp
                    alpha(0.5)
                    hold on
                    caxis auto
                    xlabel('x1'); ylabel('x2'); zlabel('Probability Density');
                    view(2)
                    shading interp
                    colorbar
                    axis(axis_)
                    hold on;plot(thisS(thisY==1,1),thisS(thisY==1,2),'r+','MarkerSize',15,'LineWidth',3)
                    hold on;scatter(thisS(:,1),thisS(:,2),(isHandpicked{methodCtr}.*50)+20,isHandpicked{methodCtr},'ko')
                    
                    
                    subplot(2,NUM_METHODS,methodCtr+NUM_METHODS)
                    cla
                    pcolor(X1,X2,reshape(d,size(X1,1),size(X1,2)));
                    xlabel('x1'); ylabel('x2'); zlabel('Probability Density');
                    view(2)
                    colorbar
                    caxis auto
                    shading interp
                    alpha(0.5)
                    axis(axis_)
                    hold on;plot(thisS(thisY==1,1),thisS(thisY==1,2),'r+','MarkerSize',15)
                    hold on;plot(thisS(thisY==1,1),thisS(thisY==1,2),'ro','MarkerSize',10)
                    hold on;plot(thisS(thisY==-1,1),thisS(thisY==-1,2),'k.','MarkerSize',10)
                    title([num2str(ctr) '/' num2str(CAMPAIGN_LENGTH)], 'FontSize',15,'FontWeight','bold');
                end
                pause(0.1)
            end
            
            
            if(shouldPlot)
                figure(perfFig)
                subplot(121)
                for ll = 1:NUM_METHODS
                    plot(cumsum((Y{ll}+1)/2),patterns{ll},'displayName',methodName{ll},'LineWidth',methodLineWidth(ll))
                    hold on
                end
                xlabel('Total sample count ','FontSize',15','FontWeight','bold')
                ylabel('Positive sample count ','FontSize',15','FontWeight','bold')
                ll = legend('show')
                set(ll,'FontSize',14);
                set(ll, 'Location','northwest')
                set(ll,'FontWeight','bold');
                xlhand = get(gca,'xlabel')
                set(xlhand,'fontsize',15)
                ylhand = get(gca,'ylabel')
                set(ylhand,'fontsize',15)
                yvec = zeros(1,NUM_METHODS);
                optRatio =  zeros(1,NUM_METHODS);
                ylim([0 length(Y{1})]);
                
                subplot(122)
                hold on
                for ll = 1:NUM_METHODS
                    plot(Stats{ll}.proportionTargetSuccess  ,patterns{ll},'displayName',methodName{ll},'LineWidth',methodLineWidth(ll))
                    hold on
                end
                hold on
                
                xlabel('Total sample count ','FontSize',15','FontWeight','bold')
                ylabel('Targeting success rate ','FontSize',15','FontWeight','bold')
                ll = legend('show')
                set(ll,'FontSize',14);
                set(ll, 'Location','northwest')
                set(ll,'FontWeight','bold');
                xlhand = get(gca,'xlabel')
                set(xlhand,'fontsize',15)
                ylhand = get(gca,'ylabel')
                set(ylhand,'fontsize',15)
                
            end
            
            for i = 1:NUM_METHODS
                a = cumsum((Y{i}+1)/2)
                [sortA, ind] = sort(1./a)
                yvec(i) = a(end);
                optRatio(i) = sum(Y{i}(21:end)>0)/sum(Y{2}(21:end)>0);
                hold on;
            end
            Yvecs = [Yvecs ; yvec]
            
           
            
            if(shouldPlot)
            figure(h)
            end
            
        end
        
        maxVec = [maxVec ;Yvecs(end,:)];
        if(shouldPlot)
            if(size(maxVec,1)>1)
                figure(statFig);clf;
                plotBoxPlotMethodComparison(methodName,maxVec);
            end
            hgsave(h,[newFolderName 'nicheSamplesDistributionAndUncertinty.fig'])
            saveas(h, [newFolderName 'nicheSamplesDistributionAndUncertinty.png'],'png');
            hgsave(perfFig,[newFolderName 'nicheSamplesPerformanceComparison.fig'])
            saveas(perfFig, [newFolderName 'nicheSamplesPerformanceComparison.png'],'png');
            save([newFolderName 'workspaceSnapshot.mat'])
            close(perfFig)
            close(initialFig)
            close(h)
        end
        StatsAll{ctr} = Stats;
        
    catch err
        [err.message]
        rethrow(err)
    end
end
keyboard

end


