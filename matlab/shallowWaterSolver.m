clc; clear; close all

nx = 2000;
xstart=-10000;
xend = 20000;
plotMin = -750;
plotMax = 750;
T = 10*60;
x_0 = linspace(xstart, xend, nx)';
a = 0.1;
epsilon = 12000e1;
seaLevel = 300;
u_0 = @(x) 0.0*ones(size(x,1),1);
eta_0 = @(x) seaLevel*ones(size(x,1),1) + 50*exp(-x.^2/epsilon);

x = x_0;
u = u_0(x);
sqrtEta = sqrt(eta_0(x));
t = zeros(size(x,1),1);
V = u + 2*sqrtEta;
W = u - 2*sqrtEta;
[~, plotFollow] = min(abs(x));

F = [];
f = figure('units','normalized','outerposition',[0 0 1 1]);
while max(t)<T
    slopeV = u(1:end-1)  + sqrtEta(1:end-1);
    slopeW = u(2:end) - sqrtEta(2:end);
    assert(all(slopeW <= slopeV));
    dx = diff(x);
    dudx = -diff(u);
    etaPlus = sqrtEta(1:end-1) + sqrtEta(2:end);
    dt1 = (dx + (u(2:end)-sqrtEta(2:end)).*(t(1:end-1)-t(2:end)))./(dudx + etaPlus);
    dt2 = (dx + (u(1:end-1)+sqrtEta(1:end-1)).*(t(1:end-1)-t(2:end)))./(dudx + etaPlus);
    t1 = t(1:end-1) + dt1;
    t2 = t(2:end) + dt2;
    t = t1;
    assert(all(t1-t2<1e-10));
    xV = x(1:end-1) + (u(1:end-1) + sqrtEta(1:end-1)).*dt1;
    xW = x(2:end) + (u(2:end) - sqrtEta(2:end)).*dt2;
    assert(all(xV-xW<1e-10));
    x = xV;
    V = V(1:end-1);
    W = W(2:end);
    sqrtEta = 0.25*(V - W);
    uV = V - 2*sqrtEta;
    uW = W + 2*sqrtEta;
    assert(all(uV-uW<1e-10));
    u = uV;
    plot(x,sqrtEta.^2,'linewidth', 4);
    plotCenter = x(plotFollow);
    axis([plotMin+plotCenter,plotMax+plotCenter,seaLevel-1,seaLevel+50])
    xlabel('Distance from Rock Fall [-]')
    ylabel('Distance from Bottom [-]')
    title(strcat('Time: ', num2str(max(t))))
    F = [F, getframe(f)];
    clf
    
    if size(x,1)<2
        break
    end
    if plotCenter>10000
        break
    end
end

close all

name = strcat('ShallowWater_char_',num2str(nx),'_',num2str(epsilon),'.avi');
myVid = VideoWriter(name);
myVid.FrameRate = 80;
open(myVid);
writeVideo(myVid, F);
close(myVid)

