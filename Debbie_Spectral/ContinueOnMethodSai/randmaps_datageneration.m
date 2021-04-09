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