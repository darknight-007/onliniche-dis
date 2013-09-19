function score =  predictionTimesUncert(pred, uncert)
score = zeros(length(pred),1);

parfor i = 1:length(pred)
    score(i) = pred(i)*uncert(i);
end
end