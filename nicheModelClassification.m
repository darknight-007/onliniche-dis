function [d predictedY] = nicheModelClassification(trainDataX_,trainDataY_,trainDataL_,testData,nDim,shouldTrainOnly,filename)

trainDataX = trainDataX_;
trainDataY = trainDataL_;

trainDataX_ = trainDataX;
testData_ = testData;


scales = [0.1 1];


a = []; 
b = []
c = []
d = []
predictedY = [];


X = trainDataX_;
Xtest = testData_;
Y = trainDataY;

par = [scales 1];

meanfunc = @meanConst; hyp.mean = 0;
covfunc = @covSEard;   hyp.cov = log(par);
likfunc = @likErf;


if(shouldTrainOnly ==1)
    hyp = minimize(hyp, @gp, -100, @infEP, meanfunc, covfunc, likfunc, X, Y);
    save(filename,'hyp','meanfunc','covfunc','likfunc','X','Y');
else
    load(filename)
    [a b c d lp] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, X, Y, Xtest,ones(length(Xtest),1));
    predictedY = exp(lp);
    
end
end