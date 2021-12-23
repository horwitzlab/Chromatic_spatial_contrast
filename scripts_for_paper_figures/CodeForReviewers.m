% For the Reviewers of "Coding of chromatic spatial contrast by macaque V1
% neurons". 

% Simulation 1
% A neuron whose receptive field consists of two linear subunits that are 
% individually half-wave rectified and then multiplied.
%
% For this hypothetical neuron, the PC1 is not larger than expected by
% chance (because the nonlinearity manifests as a stimulus dimension of low
% variance), but white noise NLI is large.

n = 1000000; % number of stimuli
nbins = 20; % For plotting 2-D firing rate map
stim = normrnd(0,1,n,6);
kernel1 = [.5 -.5 .5];
kernel2 = -kernel1;

% Simulation 1
lingen = stim*[[kernel1';0;0;0], [0;0;0;kernel2']];
resp = prod(max(lingen+1,0),2);

%resp = sum(lingen,2);
resp = (resp-min(resp))./(max(resp)-min(resp));
resp = resp>unifrnd(zeros(n,1),ones(n,1));

STA = resp'*stim;
% Projecting stim orthgonal to STA
P = eye(length(STA))-STA'*inv(STA*STA')*STA;
[v,d] = eig(P*cov(stim(resp,:))*P');
maxeigval = max(diag(d));
% Done projecting stimuli orthogonal to STA

% Permutation test on PC1
niter = 2000;
data = nan*ones(niter,1);
for i = 1:niter
    i
    tmpresp = resp(randperm(length(resp)));
    tmpSTA = tmpresp'*stim;
    P = eye(length(tmpSTA))-tmpSTA'*inv(tmpSTA*tmpSTA')*tmpSTA;
    [~,d] = eig(P*cov(stim(tmpresp,:))*P');
    data(i) = max(diag(d));
end

figure; subplot(2,1,1); hold on;
hist(data,linspace(min(data),max([data; maxeigval]),20))
plot(maxeigval,250,'rv');
axis square;

clear spikestim allstim
allstim(:,1)=stim*[STA(1:3),0 0 0]';
allstim(:,2)=stim*[0 0 0 STA(4:6)]';
spikestim(:,1)=stim(resp,:)*[STA(1:3),0 0 0]';
spikestim(:,2)=stim(resp,:)*[0 0 0 STA(4:6)]';

bins = [nbins prctile(allstim(:),5) prctile(allstim(:),95)]';

[allstim_hist,~]=hist2(allstim,[bins bins]);
[spikestim_hist,~]=hist2(spikestim,[bins bins]);
subplot(2,1,2);
imagesc(spikestim_hist./allstim_hist*255); colormap(gray(255)); hold on;
contour(spikestim_hist./allstim_hist*255,'LineColor','w'); axis square;
axis square;
set(gca,'Xtick',[],'Ytick',[]);

% Fit a GLM & GQM after the frames have been projected onto the subunits
mdllin1 =  fitglm(lingen,resp,'linear','Distribution','binomial','Link','logit');
mdlquad1 =  fitglm(lingen,resp,'quadratic','Distribution','binomial','Link','logit');
predlin1 = predict(mdllin1,lingen); % perdiction from GLM
predquad1 = predict(mdlquad1,lingen); % perdiction from GQM
Error_lin1 = 1-rocN(predlin1(resp),predlin1(~resp));
Error_quad1 = 1- rocN(predquad1(resp),predquad1(~resp));
WhiteNoise_NLI_1 = log10(Error_lin1/Error_quad1)
%%
% Simulation 2
% Two subunits that have a half-wave rectified response to one color
% channel and a full-wave rectified response to another, which are then
% added linearly. This cell has a significant PC1 because the full-wave 
% rectified responses  manifest as a stimulus direction of high variance.
% The NLI is low, however, because the GLM an GQM fitting is in the 2D
% space of projections onto the STA, and the output of the two (nonlinear)
% subfields combine linearly.

lingen = [];
lingen(:,1) = stim(:,1).*kernel1(1)+(stim(:,2).*kernel1(2)).^2+stim(:,3).*kernel1(3); % nonlinear subfield
lingen(:,2) = stim(:,[4:6])*kernel2'; % linear subfield
resp = sum(lingen,2);
resp = (resp-min(resp))./(max(resp)-min(resp));
resp = resp>unifrnd(zeros(n,1),ones(n,1));

STA = resp'*stim;

% Projecting stim orthgonal to STA
P = eye(length(STA))-STA'*inv(STA*STA')*STA;
[v,d] = eig(P*cov(stim(resp,:))*P');
% Done projecting stimuli orthogonal to STA

% Permutation test on PC1
niter = 2000;
data = nan*ones(niter,1);
for i = 1:niter
    i
    tmpresp = resp(randperm(length(resp)));
    tmpSTA = tmpresp'*stim;
    P = eye(length(tmpSTA))-tmpSTA'*inv(tmpSTA*tmpSTA')*tmpSTA;
    [~,d] = eig(P*cov(stim(tmpresp,:))*P');
    data(i) = max(diag(d));
end

figure; subplot(2,1,1); hold on;
hist(data,linspace(min(data),max([data; maxeigval]),20))
plot(maxeigval,250,'rv');
axis square;

clear allstim spikestim
allstim(:,1)=stim*[STA(1:3),0 0 0]';
allstim(:,2)=stim*[0 0 0 STA(4:6)]';
spikestim(:,1)=stim(resp,:)*[STA(1:3),0 0 0]';
spikestim(:,2)=stim(resp,:)*[0 0 0 STA(4:6)]';

bins = [nbins prctile(allstim(:),5) prctile(allstim(:),95)]';

[allstim_hist,~]=hist2(allstim,[bins bins]);
[spikestim_hist,~]=hist2(spikestim,[bins bins]);
subplot(2,1,2);
imagesc(spikestim_hist./allstim_hist*255); colormap(gray(255)); hold on;
contour(spikestim_hist./allstim_hist*255,'LineColor','w'); axis square;
axis square;
set(gca,'Xtick',[],'Ytick',[]);

% Fit a GLM & GQM after the frames have been projected onto the subunits
mdllin2 =  fitglm(lingen,resp,'linear','Distribution','binomial','Link','logit');
mdlquad2 =  fitglm(lingen,resp,'quadratic','Distribution','binomial','Link','logit');
predlin2 = predict(mdllin2,lingen); % perdiction from GLM
predquad2 = predict(mdlquad2,lingen); % perdiction from GQM
Error_lin2 = 1-rocN(predlin2(resp),predlin2(~resp));
Error_quad2 = 1- rocN(predquad2(resp),predquad2(~resp)); 
WhiteNoise_NLI_2 = log10(Error_lin2/Error_quad2)