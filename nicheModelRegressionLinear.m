function [d predictedY] = nicheModelRegressionLinear(trainDataX_,trainDataY_,trainDataL_,testData,nDim,shouldTrainOnly,filename)

trainDataX = trainDataX_;
trainDataY = trainDataY_;

trainDataX_ = trainDataX;
testData_ = testData;


scales = [0.1 1];



d = []
predictedY = [];


X = trainDataX_;
Xtest = testData_;
Y = trainDataY;

par = [scales 1];
% 
% covfunc = @covSEiso; likfunc = @likGauss; sn = 0.01; hyp.lik = log(sn);
% hyp2.cov = [1 ; 1];
% hyp2.lik = log(0.1);
% hyp2 = minimize(hyp2, @gp, -100, @infExact, [], covfunc, likfunc, Xtrain, Ytrain);
% exp(hyp2.lik)
% [mSP s2SP] = gp(hyp2, @infExact, [], covfunc, likfunc, Xtrain, Ytrain,Xtest_sigma_points);


covfunc = @covSEard;   hyp.cov = log(par); hyp.lik = log(0.1);
likfunc = @likGauss;


if(shouldTrainOnly ==1)
    w = ((X'*X)*X')*Y;
    save(filename,'w','X','Y');
else
    load(filename)
    predictedY = Xtest*w;
    d = zeros(length(predictedY),1)+0.5;
    
end
end



