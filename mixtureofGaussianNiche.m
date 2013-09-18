function [x1, x2,X1, X2, pHeadGrid,pHeadAll] =  mixtureofGaussianNiche(MAX_PROBABILITY, S_all,numHotspots)

pHeadAll = [];
pHeadGrid = [];
x1 = 33:0.01:34; x2 = 10:0.1:18;
[X1,X2] = meshgrid(x1,x2);
end