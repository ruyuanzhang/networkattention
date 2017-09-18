% This is a simulation program to attempt to develop a dynamic network model for
% attention-induced receptive-field shift
% 20170915 by RZ

clear all;close all;clc;

nNeuron = 50;  % # of numbers to simulate
sigma_r = 40; % sigma for the hyper gaussian prior for radius of the RF positions
sigma_s = 10; % sigma of the hyper gaussian prior for sizes of the RFs
sigma_g = 1.3; % % sigma of the hyper gaussian prior for sizes of the RFs
nDim = 200;  % pixel of spatial scale of visual space.
lim = [-nDim/2+1, nDim/2];
% We assumes screen size is [nDim, nDim]. So we consider the visual space as a
% coordinate with [-nDim/2+1, nDim/2] range

% Let's first simulate neurons. To describe a receptive field, we need
% five parameters. The center (x,y), size(s), gain(g) and nonlinearity(n). This
% is according the population receptive field model of a fMRI voxel. To ease the
% simulation, we fix nonlinearity to 0.2. Then other parameters of all neurons
% are drawn from a hyper distribution. These settings for the hyper
% distributions are consistent with the fMRI data. To drive the RF center(x,y),
% we first sample it in polar coordinate and then convert them to cartesian
% coordiante

n = 0.2; % we fix nonlinearity to 0.2
g = abs(sigma_g * randn(1, nNeuron)); % gain of the RFs,use positive values
r = abs(sigma_r * randn(1, nNeuron)); % radius of the center of RFs, most neurons' RF center around fixation
s = 0.3*r + 14; % we assume size of RF linearly increases as eccentricity increases.
angd = ceil(360*rand(1, nNeuron)); % angles of the center of RFs
% convert polar to cartesian coordinate 
x = r.*cosd(angd); y = r.*sind(angd);

% let's visualize 5 neurons 
iNeuron = [3, 10, 22, 30, 40, 45];  % index of the neurons to visualize.Randomly pick
close all;
figure;
for i=1:length(iNeuron)
    plotcircle(x(iNeuron(i)), y(iNeuron(i)), 2*s(iNeuron(i))); hold on;
end
ylim(lim);xlim(lim);


%% now we create 25 stimuli.
constimages = load('nnattenddata.mat','image');% load the stimuli
constimages = constimages.image; % 800x800x25. The third dimension indicate stimuli number
constimages = processmulti(@imresize,constimages,[nDim nDim]); % resize all images to 200x200
% visualize the 25 stimuli
close all; figure;
imagesc(makeimagestack(constimages));colormap(gray); 
% 25 stimuli follow the column order-the first 5 are the 1st column.

% convert it to 2D
constimages = reshape(constimages,size(constimages,1)*size(constimages,2),size(constimages,3))'; % 25x200*200;
% all stimuli have been converted to binary masks
%derive the center of each stimulus;
stimloc = [3.2500    4.7500    6.2500    7.7500    9.2500]/12.5*nDim;% center location of each image
[centerx,centery]=meshgrid(1:5,1:5);
stimx=flatten(stimloc(centerx));
stimy=flatten(stimloc(centery));
stimecc = sqrt((stimx-100).^2+(stimy-100).^2)/nDim*12.5; %eccentricties(deg) of 25 positions
clear centerloc centerx centery;
[~,xx,yy]= makegaussian2d(nDim,nDim/2,nDim/2,10,10); % obtain xx,yy to speed up;


%% ======== compute neuronal responses of a single neuron ============
% Given the 25 stimuli and a neuron's RF, we can predict 25 responses of this
% neuron. The stimulus-response model first takes dot product between the RF and
% the stimulus, and then passes result through nonlinearity and adds gain onto
% it.

% let's first use one neuron as an example
iNeuron = 30; % let use this neuron as an example
% create the RF of this neuron
RF = makegaussian2d(nDim, x(iNeuron)+100, y(iNeuron)+100, s(iNeuron), s(iNeuron), xx, yy)/(2*pi*s(iNeuron)^2); % nDimxnDim
RF = vflatten(RF); % convert it to a column vector. nDim*nDim x 1
resp = g(iNeuron)*(constimages*RF).^n; % resp is a 25x1 vector
% visualize the RF of this neuron
close all;figure;
ax(1)=subplot(1,2,1);
plotcircle(x(iNeuron), y(iNeuron), 2*s(iNeuron)); hold on;
xlim(lim);ylim(lim);
ax(2)=subplot(1,2,2); % plot the 25 response
bar(1:25, resp);
ylabel('Response');xlabel('Stimuli');
xlim([0 26]);
%% ======== compute neuronal responses of all neurons ============
% now compute all neurons responses. We have nNeuron, and each neuron has 25
% responses towards 25 stimuli.
RF_all = zeros(nDim*nDim, nNeuron);
for i=1:nNeuron % loop neurons
    RF_all(:, i) = vflatten(makegaussian2d(nDim, x(i)+99, y(i)+99, s(i), s(i), xx, yy)/(2*pi*s(i)^2));
end
g_all = repmat(g, [25, 1]); % we convert g to a 25xnNeuron matrix
resp_all = g_all.*(constimages*RF_all).^n; %resp_all is a 25xnNeuron matrix
% resp_all is a 25xnNeuron matrix, each column is 25 responses of a neuron.
% we can plot 5 neurons
iNeuron = [3, 10, 22, 30, 40, 45];  % index of the neurons to visualize,randomly pick
close all;
cpsfigure(2,length(iNeuron));
for i=1:length(iNeuron)
    ax(i) = subplot(2,length(iNeuron), i);
    plotcircle(x(iNeuron(i)), y(iNeuron(i)), 2*s(iNeuron(i))); hold on;
    xlim(lim);ylim(lim);
    ax(i+length(iNeuron)) = subplot(2,length(iNeuron),i+length(iNeuron));
    bar(1:25, resp_all(:,iNeuron(i)));
    xlim([0 26]);
end
% we can also derive the tuning similarity matrix by computing the correlation
% between RFs
figure;
imagesc(corr(RF_all));
title('RF correlation matrix')

%% =============== deal with real voxel data ======================
% all analyses above are based on simulation data. I tried to make the
% simulation very similar to real dataset. Here, I provide 50 real neurons(voxels)
% load the data
data = load('nnattenddata.mat');
% data has four fields:
%   beta: 50(nNenuron) x 50(25 stim x 2 tasks) x 100 bootstrapsamples
%   betamn: take the median(beta,3), to obtain the median of responses
%   image:800x800x25 stimuli images, used as before
%   pRFparams: a cell array with RF parameters for the fixation task in element{1},
%       and attention task in element{2}. Each task has a 5(params)x50(nNeurons)
%       matrix. 5 params are (x,y),size,gain,eccentriciy. All params are in
%       pixel unit, except gain.

% we can inspect the attention effect
fixparams = median(data.pRFparams{1}, 2);
attendparams = median(data.pRFparams{2}, 2);
figure;
ax(1)=subplot(1,3,1); % plot RF size
bar([fixparams(3), attendparams(3)]);
ylabel('RF size');
ax(2)=subplot(1,3,2); % plot gain
bar([fixparams(4), attendparams(4)]);
ylabel('Gain');
ax(3)=subplot(1,3,3); % plot gain
bar([fixparams(5), attendparams(5)]);
ylabel('Eccentricity');
set(ax,'XTickLabels',{'fix','att'});
% As you might see, in the attention task, Neuron's RF size, gain and eccentricity
% increases, compared to the fixation task.

      


