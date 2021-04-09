% spectral super resolution, sai ravela (C) Copyright Sai Ravela, ESSG, MIT
% 2020.
% Produces synthetic fields with a convergent spectra (exponential rate)
% first run randmaps_datageneration.m
N=256; 
[x0,y0]=meshgrid(1:N,1:N);
filename='sresanim.gif';

%% Load data
load zrec;

%% Select one example and move square

% Just select one
sel_exp=3;
zhi = zrechi(:,:,sel_exp);
zlo = zreclo(:,:,sel_exp);

zout = zeros(size(zhi));
maz = max(zhi(:)); miz = min(zhi(:));
set(gcf,'color','w')

%Move square across the grid and create input and output pairs in zout
for i = 1:32: 256-32
  for j = 1:32:256-32
      
      zout(j:j+64-1,i:i+64-1)=zhi(j:j+64-1,i:i+64-1);
      subplot(121); imagesc(zlo);
      title('Low Resolution Input');
      axis('image'); axis('off');
      rectangle('Position',[i j 64 64],'LineWidth',2)
      subplot(122);
      imagesc(zout,[miz,maz]);
      title('High Resolution Output');
      axis('image');axis('off');
      drawnow;
          frame = getframe(gcf);
         im = frame2im(frame);
     [A,map] = rgb2ind(im,256);
    if i == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.2);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.2);
    end

  end
end

%%
ntrain = 512;
nsz = size(zreclo,3); % p
ztrainhi = zeros(64*64,ntrain);ztrainlo =ztrainhi;

% get ntrain randomly picked samples of 64x64 from random experiments at
% random locations (64*65=4096,512)
for i = 1:ntrain
    idx = randi(nsz); % pick random experiment
    loc = randi([33 256-32],2); % pick random location on the grid
    tmp =  squeeze(zreclo(loc(1)-31:loc(1)+32,loc(2)-31:loc(2)+32,idx));
    ztrainlo(:,i) = tmp(:); % get the low res secctio
    tmp = squeeze(zrechi(loc(1)-31:loc(1)+32,loc(2)-31:loc(2)+32,idx));
    ztrainhi(:,i)=tmp(:);
    disp(i);
end

%% Define models
%get mean and offset of the mean of trainset
mlo = mean(ztrainlo,2); dlo = ztrainlo - mlo;
mhi = mean(ztrainhi,2); dhi = ztrainhi - mhi;

% define models:
modl = (dhi*dlo')*pinv(dlo*dlo'); %Cov est
modl2 = (ztrainhi*ztrainlo')*pinv(ztrainlo*ztrainlo'); % Corr. est
modl3 = ztrainhi*pinv(ztrainlo); % Cheap

%% Test models

% pick a random experiment and random location
idx = randi(nsz)
loc = randi([33 256-32],2);

%get the lo res and hi res test square
ztestlo =  squeeze(zreclo(loc(1)-31:loc(1)+32,loc(2)-31:loc(2)+32,idx));
ztesthi = squeeze(zrechi(loc(1)-31:loc(1)+32,loc(2)-31:loc(2)+32,idx));

% reconstruct high resolution with models
yy  = modl*ztestlo(:);
yy2 = modl2*ztestlo(:);
yy3 = modl3*ztestlo(:);

subplot(325); imagesc(reshape(yy3,[64 64]));colorbar;axis('image');title('est cheap');
subplot(324); imagesc(reshape(yy2,[64 64]));colorbar;axis('image');title('Est Corr');
subplot(323); imagesc(reshape(yy,[64 64]));colorbar;axis('image');title('Est Cov');
subplot(321); imagesc(ztestlo);colorbar;axis('image');title('Input');
subplot(322); imagesc(ztesthi);colorbar;title('Truth');
axis('image');

%% Apply models to many test examples to get enough data for svd
for i = 1:ntrain
    idx = randi(nsz);
    loc = randi([33 256-32],2);
    
    ztestlo =  squeeze(zreclo(loc(1)-31:loc(1)+32,loc(2)-31:loc(2)+32,idx));
    ztesthi = squeeze(zrechi(loc(1)-31:loc(1)+32,loc(2)-31:loc(2)+32,idx));
    
    yy=modl*ztestlo(:);
    yy3=modl3*ztestlo(:);
    ztestenshi(:,i) = yy;
    ztestenshi3(:,i)= yy3;
end
dtest = ztestenshi-mean(ztestenshi,2);

subplot(326); 
%svd gives just the singluar values (diagonal of S matrix) of the svd
plot(svd(dlo)); hold on; 
plot(svd(dhi));hold off;
plot(svd(dtest));
plot(svd(ztestenshi3));hold off;
axis([0 600 0 50]);
legend('lo','true','est cov','cheap');