clear; close all

mrstModule add mimetic
g = 9.81;

% Set up a Cartesian grid.
nx = 400;
ny = 10;
xMax = 46000;
yMax = 300;
gridLimits = [xMax, yMax];
dx = 1/nx;
G = cartGrid([nx, ny], [xMax, yMax]);
G = computeGeometry(G);

T = 700;
k = 2;

% Define initial values

add = 22000;
%hStart = 1; hEnd = 0.2; x1 = 0.75; x2 = 1.3;
%% HELLESYLT
%hStart = yMax; hEnd = 0.1*yMax; x1 = 7000 + xMax/2; x2 = 12000 + xMax/2; %0.8*xMax; x2 = 0.9*xMax;
%% GEIRANGER
hStart = yMax; hEnd = 0.1*yMax; x1 = 7000 + add; x2 = 23000 + add; %0.8*xMax; x2 = 0.9*xMax;


xCenter = 1.3;

bottomType = 'Geir';

%h = @(x) -0.25*atan(60*(x-xCenter)/xMax)./pi +0.145;

%h = @(x) ones(size(x,1),1);
%h = @(x) max(hStart,min((hEnd-hStart)/(x2-x1)*(x-x1) + hStart,hEnd));%ones(size(x,1),1);
h = @(x) min(hStart,max((hEnd-hStart)/(x2-x1)*(x-x1) + hStart,hEnd));%ones(size(x,1),1);
hPlot = @(x) h(x); % 0.1*min(hStart,max((hEnd-hStart)/(x2-x1)*(x-x1) + hStart,hEnd));
epsilon = 120000;
eta_0 = @(x) 50*exp(-(x-add).^2/epsilon);%zeros(size(x,1),1);
phi_top_0 = @(x) zeros(size(x,1),1);%3700*exp(-(x-0.2*xMax).^2/epsilon); %zeros(size(x,1),1);

top_faces = (1 : G.faces.num)';
top_centroids = G.faces.centroids(G.faces.centroids(:, 2) == gridLimits(2));

eta_old = eta_0(G.nodes.coords(:,1));
phi_top_old = phi_top_0(top_centroids);
phi_top = phi_top_old;
eta = eta_old;

F = [];
f = figure('units','normalized','outerposition',[0 0 1 1]);

for t=0:k:T
    % Calculate phi^n+1
    [phi, gradPhi_top, Gnew] = poissonMimetic2D(G,h,phi_top, eta, gridLimits, 0);
    % Remove extra etas
    
    eta = eta(1:nx+1);  
    
    % Calculate the derivative
    etax = diff(eta)./dx;
    % Center the etas at the face centers
        
    node_coords = Gnew.nodes.coords(G.nodes.coords(:,2)==yMax, 1);
    centroids = Gnew.faces.centroids(G.faces.centroids(:,2)==yMax, 1);
    
    eta = 0.5*(eta(2:end)+eta(1:end-1));
    
    
    % Find eta^n+1
    %eta_centroid = eta_centroids - k*sum(gradPhi_top.*[etax,-ones(size(etax,1),1)],2);
    %eta_centroids = eta_centroids - k*sum(gradPhi_top.*[gradEta_centroids,-ones(size(gradEta_centroids,1),1)],2);
    %centroids = G.faces.centroids(top_faces);
    spaceStep = diff(centroids(:,1));
    eta = findEta(eta, gradPhi_top, spaceStep, k);
    %eta = (eta - k*sum(gradPhi_top.*[etax,-ones(size(etax,1),1)],2));%plot(eta)
    
    % Calculate phi^n+2
    phi_top = phi_top-k*(0.5*sum(gradPhi_top.^2,2) - g*eta);
    
    
    % Change eta back to the nodes.
    eta = 0.5*(eta(2:end)+eta(1:end-1));
    deta1 = eta(2)-eta(1);         %Should divide on dx, but cancels:
    deta2 = eta(end)-eta(end-1);
    eta1 = eta(1)- deta1;     % HERE
    eta2 = eta(end) + deta2;
    eta = [eta1; eta; eta2];

    x = 0:xMax/nx:xMax;
    %plot(x,eta,x,-hPlot(x));
    %axis([0.4, xMax, -0.3, 0.030]);
    %axis off;
    bl = [0 0.4470 0.7410];
    plotGrid(Gnew,'facecolor', bl, 'edgecolor', 'blue');
    etaMax = max(eta(centroids > add));
    etaMax = round(etaMax(1));
    %axis([xMax/2, xMax, -0.2*yMax, 0.2*yMax]);
    axis([add, xMax, -yMax, 0.2*yMax]);
    curtick = get(gca, 'XTick');
    set(gca, 'XTickLabel', cellstr(num2str(curtick(:)-add)));
    
    titleString = {strcat('Time: ', num2str(t), ' s') ...
                   strcat('\eta: ', num2str(etaMax), ' m')};
    dim = [.5 .6 .3 .3];
    annotation('textbox', dim, 'String', titleString, 'FitBoxToText', 'on');
    title('Tsunami approaching Geiranger')
%     tbx = uicontrol('style','text')
%     set(tbx,'String',titleString)
    xlabel('x [m]');
    ylabel('z [m]');
    F = [F, getframe(f)];
%     x = 0:xMax/(nx):xMax;
%     x = x(1:end-1);
%     figure(2);
%     v = sqrt(sum(gradPhi_top.^2,2));
%     plot(x,v,x,-h(x));
%     
    %plotCellData(Gnew,phi)
    %caxis([-0.01, 0.01])
    %prep eta for next iteration
    eta = repmat(eta,ny+1,1);
    pause(0.0001)
    %pause()
    clf;
end

close all

name = strcat('wave_',num2str(nx),'_',num2str(ny),'_',bottomType,'.avi');
myVid = VideoWriter(name);
myVid.FrameRate = 8;
open(myVid);
writeVideo(myVid, F);
close(myVid)

    %figure
    %plotCellData(Gnew,phi)
    %colorbar