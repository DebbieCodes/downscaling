% spectral super resolution, sai ravela (C) Copyright Sai Ravela, ESSG, MIT
% 2020.
% Produces synthetic fields with a convergent spectra (exponential rate)
filename='sresanim.gif';
N=256; 
[x0,y0]=meshgrid(1:N,1:N);

wd = 0;
for exp = 1:1000
ph = rand(N+2*wd);

th = rand*2*pi;
R = [cos(th) -sin(th); sin(th) cos(th)];
S = [1 0;...
     0 1+(rand-1/2)]; 
p1 = R*S*[x0(:)' ; y0(:)'];

x = reshape(p1(1,:)',[N N]); mX = max(abs(x(:)));
y = reshape(p1(2,:)',[N N]); mY = max(abs(y(:)));
zlo =zeros(size(x));
for wx = 2:N/32
    for wy=2:N/32
        m=sqrt(wx.*wx+wy.*wy);
zlo=zlo+m^(-2)*sin(2*pi*(wy*y/mY+...
                                 ph(wy,wx)+ wx*x/mX));
    end
end
zhi=zlo;
for wx = N/32:N/4
    for wy=N/32:N/4
        m=sqrt(wx.*wx+wy.*wy);
zhi=zhi+m^(-5/3)*sin(2*pi*(wx*x/mX+...
                         ph(wy,wx) + wy*y/mY));
    end
end
%zlo=zlo./var(zlo(:));
%zhi = zhi./var(zhi(:));

% colormap('jet');
% subplot(121);imagesc(zlo);
% axis('image');title('LoRes');
% subplot(122);imagesc(zhi);
% axis('image');title('HiRes');drawnow;
zreclo(:,:,exp) = zlo;
zrechi(:,:,exp) = zhi;
disp(exp)
end

save zrec zreclo zrechi;
%%
load zrec;
zout = zeros(size(zhi));
maz = max(zhi(:)); miz = min(zhi(:));
set(gcf,'color','w')
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
nsz = size(zreclo,3);
ztrainhi = zeros(64*64,ntrain);ztrainlo =ztrainhi;
for i = 1:ntrain
    idx = randi(nsz);
    loc = randi([33 256-32],2);
    tmp =  squeeze(zreclo(loc(1)-31:loc(1)+32,loc(2)-31:loc(2)+32,idx));
    ztrainlo(:,i) = tmp(:);
    tmp = squeeze(zrechi(loc(1)-31:loc(1)+32,loc(2)-31:loc(2)+32,idx));
    ztrainhi(:,i)=tmp(:);
    disp(i);
end
mlo = mean(ztrainlo,2); dlo = ztrainlo - mlo;
mhi = mean(ztrainhi,2); dhi = ztrainhi - mhi;
modl = (dhi*dlo')*pinv(dlo*dlo');

modl2 = (ztrainhi*ztrainlo')*pinv(ztrainlo*ztrainlo');
modl3 = ztrainhi*pinv(ztrainlo);
%%
idx = randi(nsz)
    loc = randi([33 256-32],2);
    ztestlo =  squeeze(zreclo(loc(1)-31:loc(1)+32,loc(2)-31:loc(2)+32,idx));
    ztesthi = squeeze(zrechi(loc(1)-31:loc(1)+32,loc(2)-31:loc(2)+32,idx));
    
yy=mhi+modl*ztestlo(:)-mlo;
yy2 = modl2*ztestlo(:);
yy3 = modl3*ztestlo(:);
subplot(326); plot(svd(dlo)); hold on; plot(svd(dhi));hold off;
subplot(325); imagesc(reshape(yy3,[64 64]));colorbar;axis('image');title('est cheap');
subplot(324); imagesc(reshape(yy2,[64 64]));colorbar;axis('image');title('Est Corr');
subplot(323); imagesc(reshape(yy,[64 64]));colorbar;axis('image');title('Est Cov');
subplot(321); imagesc(ztestlo);colorbar;axis('image');title('Input');
subplot(322); imagesc(ztesthi);colorbar;title('Truth');
axis('image');

%%
for i = 1:ntrain
idx = randi(nsz);
    loc = randi([33 256-32],2);
    ztestlo =  squeeze(zreclo(loc(1)-31:loc(1)+32,loc(2)-31:loc(2)+32,idx));
    ztesthi = squeeze(zrechi(loc(1)-31:loc(1)+32,loc(2)-31:loc(2)+32,idx));
yy=mhi+modl*ztestlo(:)-mlo;
yy3=modl3*ztestlo(:);
ztestenshi(:,i) = yy;
ztestenshi3(:,i)=yy3;
end
dtest = ztestenshi-mean(ztestenshi,2);
subplot(326); hold on; plot(svd(dtest));plot(svd(ztestenshi3));hold off;axis([0 600 0 50]);
legend('lo','true','est cov','cheap');