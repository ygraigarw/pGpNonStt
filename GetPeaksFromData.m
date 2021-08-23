%% Clear everything
clear all;

%% Load the data file you want to process
load historical.mat

%% Find peaks over threshold
XTim=data(:,1); %this variable must be the time of the event
XVal=data(:,3); %this variable must be the size of the event
Thr=5;          %the threshold to use
[PkTim, PkVal]=pPot(XTim, XVal, Thr);

%% Make a plot
clf; hold on;
plot(XTim,XVal,'k.');
plot(PkTim,PkVal,'ro');

%% Save data read for pGpNonStt
%X.nYr ; %  1   x 1 number of years
%X.nT  ; %  1   x 1 number of occurrences
%X.Tim ; % nT   x 1 years on [0,1] (so that floor((X.Tim*X.nYr)+1) gives the year number
%X.Dat ; % nT   x 1 data

X.nYr=range(PkTim);
X.nT=size(PkVal,1);
X.Tim=(PkTim-PkTim(1))./range(PkTim);
X.Dat=PkVal;

save historical_peak.mat X;

