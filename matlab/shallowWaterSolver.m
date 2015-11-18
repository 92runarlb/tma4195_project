clc; clear; close all

nx = 800;
xstart=-100;
xend = 100;
plotFollow = ceil(nx/2);
plotMin = -5;
plotMax = 5;
T = 20;
x_0 = linspace(xstart, xend, nx)';
a = 0.1;
epsilon = 1e0;
seaLevel = 10;
u_0 = @(x) 0.0*ones(size(x,1),1);
eta_0 = @(x) seaLevel*ones(size(x,1),1) + 1*exp(-x.^2/epsilon);

x = x_0;
u = u_0(x);
sqrtEta = sqrt(eta_0(x));
t = zeros(size(x,1),1);
V = u + 2*sqrtEta;
W = u - 2*sqrtEta;


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
    assert(all(xV-xW<1e-12));
    x = xV;
    V = V(1:end-1);
    W = W(2:end);
    sqrtEta = 0.25*(V - W);
    uV = V - 2*sqrtEta;
    uW = W + 2*sqrtEta;
    assert(all(uV-uW<1e-10));
    u = uV;
    plot(x,sqrtEta.^2);
    plotCenter = x(plotFollow);
    axis([plotMin+plotCenter,plotMax+plotCenter,seaLevel-1,seaLevel+1])
    pause(0.1)
    if size(x,1)<2
        break
    end
end


