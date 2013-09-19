function [d predictedY] = nicheModelClassificationLogReg(trainDataX_,trainDataY_,trainDataL_,testData,nDim,shouldTrainOnly,filename)

trainDataX = trainDataX_;
trainDataY = trainDataL_;

trainDataX_ = trainDataX;
testData_ = testData;


a = []; 
b = []
c = []
d = []
predictedY = [];


X = trainDataX_;
Xtest = testData_;
Y = trainDataY;



if(shouldTrainOnly ==1)
    bHat = glmfit(X,(Y+1)/2,'binomial')
    save(filename,'bHat','X','Y');
else
    load(filename)
    yHat = 1./(1+exp( -[ones( size(testData_,1),1 ), testData_] *bHat));
    predictedY = yHat;
    d = zeros(length(yHat),1)+0.5;
end
end