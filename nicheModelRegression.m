function [d predictedY] = nicheModelRegression(trainDataX_,trainDataY_,trainDataL_,testData,nDim,shouldTrainOnly,filename)

trainDataX = trainDataX_;
trainDataY = trainDataY_;

trainDataX_ = trainDataX;
testData_ = testData;


scales = [0.1 0.001];



d = []
predictedY = [];


X = trainDataX_;
Xtest = testData_;
Y = trainDataY;

par = [scales 1];

covfunc = @covSEard;   hyp.cov = log(par); hyp.lik = log(0.1);
likfunc = @likGauss;


if(shouldTrainOnly ==1)
    hyp = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, X, Y);
    save(filename,'hyp','covfunc','likfunc','X','Y');
else
    load(filename)
    [mSP, s2SP] = gp(hyp, @infExact, [], covfunc, likfunc, X, Y, Xtest);
    predictedY = mSP;
    d = s2SP;
    
end
end



