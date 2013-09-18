function [x1, x2,X1, X2, pHeadGrid,pHeadAll] =  mixtureofGaussianNiche(MAX_PROBABILITY, S_all,numHotspots)

pHeadAll = [];
x1 = 33:0.01:34; x2 = 10:0.1:18;

[X1,X2] = meshgrid(x1,x2);
F = zeros(length(X1(:)),1);
X1vec = X1(:);
X2vec = X2(:);

S_grid = [X1vec, X2vec];
mu1 = [33.59   13.76];

%mu1 = [33.5929 15.7071]; % geo localized to north nonterey bay
Sigma1 =[ 0.00035	0; 0	0.1];

mu2 = [33.83   10.5];
Sigma2 =[ 0.00035*10	0; 0	0.1];



pHeadAll = unnormProbFunction(S_grid,mu1,Sigma1,mu2,Sigma2,numHotspots);
max_pHead = max(pHeadAll);
pHeadGrid = reshape((pHeadAll/max_pHead)*MAX_PROBABILITY, size(X1));
%pcolor(X1,X2,pHeadGrid);


if(~isempty(S_all))
    pHeadUnnorm = unnormProbFunction(S_all,mu1,Sigma1,mu2,Sigma2,numHotspots);
    pHeadAll = (pHeadUnnorm/max_pHead)*MAX_PROBABILITY;
end

end


function pHead = unnormProbFunction(S_all,mu1,Sigma1,mu2,Sigma2,numHotspots)
pHead = [];
if(numHotspots == 2)
    pHead = mvnpdf(S_all,mu1,Sigma1) + mvnpdf(S_all,mu2,Sigma2);
else
    pHead = mvnpdf(S_all,mu1,Sigma1);
end
end