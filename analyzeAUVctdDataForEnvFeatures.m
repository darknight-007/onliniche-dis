function auvdata =  analyzeAUVctdDataForEnvFeatures(filename)
%returns tTSONFBDLonLat, 10 dim

load(filename,'time');
mvec = structfun(@(st)min(st),time);
t1 = min(mvec(~isnan(mvec)));
mvec = structfun(@(st)max(st),time);
tend = max(mvec(~isnan(mvec)));

% CTD
load(filename,'ctd1');
filt = ~(isnan(ctd1.T) | isnan(ctd1.S) | isnan(ctd1.O2) | isnan(time.ctd1));
Tf = ctd1.T(filt);
Sf = ctd1.S(filt);
O2f = ctd1.O2(filt);
ctdtf = time.ctd1(filt);
t1 = max([t1,ctdtf(1)]);
tend = min([tend,ctdtf(end)]);

% depth and profile count
load(filename,'depth');
filt = ~(isnan(time.nav) | isnan(depth.nav));
navdf = depth.nav(filt);
navtf = time.nav(filt);
t1 = max([t1,navtf(1)]);
tend = min([tend,navtf(end)]);

% HS2
load(filename,'hs2');
% hs2 is a little more complicated because the wavelengths are
% not fixed.
filt = ~isnan(time.hs2);
hs2vecs = zeros(length(time.hs2),3);
fn = sort(fieldnames(hs2));
for ii=1:3
    hs2vecs(:,ii) = hs2.(fn{ii});
    filt = filt & ~isnan(hs2vecs(:,ii));
end
BBSf = hs2vecs(filt,1); % short wavelength
BBLf = hs2vecs(filt,2); % long wavelength
FLLf = hs2vecs(filt,3); % data product
hs2tf = time.hs2(filt);
t1 = max([t1,hs2tf(1)]);
tend = min([tend,hs2tf(end)]);

% Isus
load(filename,'isus');
filt = ~(isnan(isus.nitrate) | isnan(time.isus));
Nf = isus.nitrate(filt);
isustf = time.isus(filt);
t1 = max([t1,isustf(1)]);
tend = min([tend,isustf(end)]);

tvec = (t1:.5:tend)';

load(filename,'auvnav');
%TSONFB
% tT
[~,uniqueHS2F,~]=unique(hs2tf);
hs2tf =  hs2tf(uniqueHS2F);
FLLf = FLLf(uniqueHS2F);
BBSf = BBSf(uniqueHS2F);

[~,uniqueHS2F,~]=unique(hs2tf);
hs2tf =  hs2tf(uniqueHS2F);
FLLf = FLLf(uniqueHS2F);
BBSf = BBSf(uniqueHS2F);
[~,uniqueHS2F,~]=unique(hs2tf);
hs2tf =  hs2tf(uniqueHS2F);
FLLf = FLLf(uniqueHS2F);
BBSf = BBSf(uniqueHS2F);
[~,uniqueHS2F,~]=unique(hs2tf);
hs2tf =  hs2tf(uniqueHS2F);
FLLf = FLLf(uniqueHS2F);
BBSf = BBSf(uniqueHS2F);
[~,uniqueHS2F,~]=unique(hs2tf);
hs2tf =  hs2tf(uniqueHS2F);
FLLf = FLLf(uniqueHS2F);
BBSf = BBSf(uniqueHS2F);

[~,uniqueIsustf,~] = unique(isustf);
isustf = isustf(uniqueIsustf);
Nf = Nf(uniqueIsustf);

auvdata = [tvec interp1(ctdtf,Tf,tvec,'nearest'),interp1(ctdtf,Sf,tvec,'nearest'), interp1(ctdtf,O2f,tvec,'nearest'), interp1(isustf,Nf,tvec,'nearest') , interp1(hs2tf,FLLf,tvec,'nearest'),  interp1(hs2tf,BBSf,tvec,'nearest') interp1(navtf,navdf,tvec) interp1(time.nav,auvnav.dvlLon,tvec) interp1(time.nav,auvnav.dvlLat,tvec) ];
end