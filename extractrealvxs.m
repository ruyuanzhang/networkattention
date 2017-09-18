%% extract 50 neurons 

close all;clc;clear all;
datadir = '/Users/ruyuan/Dropbox/stonesync/attentionprf';

%% ============ read in the data ====================
% load data;
data = matfile(sprintf('%s/fMRIdata/dataset01.mat',datadir));
tmp1 = data.betas;tmp1 = tmp1(:,2); % we only analyze 2nd task
data = matfile(sprintf('%s/fMRIdata/dataset02.mat',datadir));
tmp2 = data.betas;tmp2 = tmp2(:,2);
data = matfile(sprintf('%s/fMRIdata/dataset03.mat',datadir));
tmp3 = data.betas;tmp3 = tmp3(:,2);
%concatenate across three subject. we pool voxels from thee subjects into one.
betas = cellfun(@vertcat,tmp1,tmp2,tmp3,'UniformOutput',0);
fprintf('finish download data ....\n'); clear tmp1 tmp2 tmp3;

% Pool voxels across hemisphere
for iRoi = 1:7
    betas{iRoi} = vertcat(betas{iRoi},betas{iRoi+7});
end
betas = betas(1:7);

%% read in pRF fit results
digitfit = load(sprintf('%s/digittaskCSSfit.mat',datadir));
facefit = load(sprintf('%s/facetaskCSSfit.mat',datadir));
for iRoi =1:14 %loop roi
    tmp = digitfit.paramsacrossroi_digit{iRoi}.params;
    tmp(:,end+1,:)=sqrt((squeeze(tmp(:,1,:))-100).^2+(squeeze(tmp(:,2,:))-100).^2); % calculate eccentricity in pixels
    digitfit.paramsacrossroi_digit{iRoi}.params = tmp;
    clear tmp;
    
    tmp = facefit.paramsacrossroi_face{iRoi}.params;
    tmp(:,end+1,:)=sqrt((squeeze(tmp(:,1,:))-100).^2+(squeeze(tmp(:,2,:))-100).^2); % calculate eccentricity in pixels
    facefit.paramsacrossroi_face{iRoi}.params = tmp;
    clear tmp;
end
%%
for iRoi = 1:7 
    digitfit.paramsacrossroi_digit{iRoi}.params = cat(3, digitfit.paramsacrossroi_digit{iRoi}.params, digitfit.paramsacrossroi_digit{iRoi+7}.params);
    facefit.paramsacrossroi_face{iRoi}.params = cat(3, facefit.paramsacrossroi_face{iRoi}.params, facefit.paramsacrossroi_face{iRoi+7}.params);
    digitfit.paramsacrossroi_digit{iRoi}.trainperformance = cat(2, digitfit.paramsacrossroi_digit{iRoi}.trainperformance, digitfit.paramsacrossroi_digit{iRoi+7}.trainperformance);
    facefit.paramsacrossroi_face{iRoi}.trainperformance = cat(2, facefit.paramsacrossroi_face{iRoi}.trainperformance, facefit.paramsacrossroi_face{iRoi+7}.trainperformance);
end
% calculate the x,y as the eccentricity

%% let's figure out the distribution of all 
tRoi = 7; % we focus on mFUS
cpsfigure(1,3);
for i=1:3
    ax(i)=subplot(1,3,i);
    histogram(digitfit.paramsacrossroi_digit{tRoi}.params(:,2+i,:));hold on;
    histogram(facefit.paramsacrossroi_face{tRoi}.params(:,2+i,:));
end

% plot relationship between eccentricity and pRF size
ecc = squeeze(facefit.paramsacrossroi_face{tRoi}.params(:,5,:))';
width = squeeze(facefit.paramsacrossroi_face{tRoi}.params(:,3,:))';

% choose good vxs
ind = find(digitfit.paramsacrossroi_digit{tRoi}.trainperformance>75&facefit.paramsacrossroi_face{tRoi}.trainperformance>75);
ind = ind(1:50);

% save the results
beta = betas{tRoi}(ind,[1:25 51:75],:);
betamn = median(beta,3);
pRFparams=cell(1,2);
pRFparams{1} = squeeze(digitfit.paramsacrossroi_digit{tRoi}.params(:,:,ind));
pRFparams{2} = squeeze(facefit.paramsacrossroi_face{tRoi}.params(:,:,ind));
save('nnattenddata.mat','beta','betamn','pRFparams','image');



