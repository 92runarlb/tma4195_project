clc; clear; close all;

g = 9.81;


xstart=-1000;
xend = 1000;

T = 1;
dt = 0.00001;
dx = 0.1;

epsilon = 1e2;
seaLevel = 0;
u_0 = @(x) 0*ones(size(x,1),1);
%z_0 = @(x) min([(x-2),ones(size(x,1),1)],[],2);
z_0 = @(x) 500/xend*(x-xend);
h_0 = @(x) seaLevel*ones(size(x,1),1)+ 1*exp(-x.^2/epsilon)-min([z_0(x),seaLevel*ones(size(x))],[],2)+ 0.001;

fluxW = @(h1,u1,h2,u2,lambda) 0.5*(h1.*u1+h2.*u2) %-1./(2*lambda).*(h1-h2);
fluxZ = @(h1,u1,h2,u2,lambda) 0.5*(h1.*u1.^2 + g*h1.^2/2 + h2.*u2.^2 + g*h2.^2/2) ...
                            %-1./(20*lambda).*(h1.*u1-h2.*u2);

x = (xstart:dx:xend)';
z = z_0(x);
u = u_0(x);
h = h_0(x);
lambda = dt./dx;
W = h;
Z = h.*u;

F = [];
f = figure('units','normalized','outerposition',[0 0 1 1]);
takePic = 100;
i = 0;
for t= 0:dt:T
    zExt = [z(1); z; z(end)];
    hExt = [h(1); h; h(end)];
    uExt = [u(1); u; u(end)];
    zFace = max([zExt(1:end-1),zExt(2:end)],[],2);
    hFace_l = max([zeros(size(hExt,1)-1,1), hExt(1:end-1) + zExt(1:end-1) - zFace], [], 2);
    hFace_r = max([zeros(size(hExt,1)-1,1), hExt(2:end) + zExt(2:end) - zFace], [], 2); 
    Wsource = zeros(size(x,1),1);
    Zsource = g/2*(hFace_l(2:end).^2 - hFace_r(1:end-1).^2);
    Wface_l = hFace_l(2:end); % do not think I need this
    Wface_r = hFace_r(2:end);
    Zface_l = hFace_l(2:end).*u;
    Zface_r = hFace_r(2:end).*u;

    
    FWface = fluxW(hFace_l, uExt(1:end-1), hFace_r, uExt(2:end),lambda);
    FZface = fluxZ(hFace_l, uExt(1:end-1), hFace_r, uExt(2:end),lambda);
    
    Wrhs = Wsource - (FWface(2:end) - FWface(1:end-1));
    Zrhs = Zsource - (FZface(2:end) - FZface(1:end-1));
    
    W = W + lambda*Wrhs;
    Z = Z + lambda*Zrhs;
    
    h = W;
    u = Z./h;
    
    clf
    plot(x, h+z)
    hold on
    plot(x,z);
    axis([0-50,xend,seaLevel-1,seaLevel+1])
    title(t)
    if mod(i,takePic)==0
        F = [F, getframe(f)];    
    end
end

close all

name = strcat('shallowWaterChangingBottom',num2str(xstart),'_',num2str(xend),'_',num2str(T),'_','.avi');
myVid = VideoWriter(name);
myVid.FrameRate = round(24/dt);
open(myVid);
writeVideo(myVid, F);
close(myVid)
