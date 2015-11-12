clear; close all

mrstModule add mimetic
g = 9.81;

% Set up a Cartesian grid.
xmax = 1;
ymax = 3;
zmax = 1;
nx = 40;
ny = 120;
nz = 20;
dx = 1/nx;
G = cartGrid([nx, ny, nz], [xmax, ymax, zmax]);
G = computeGeometry(G);

T = 2;
k = 0.05;

% Define initial values
h = @(x,y) -ones(size(x,1),1);
epsilon = 9e-3;
eta_0 = @(x,y) -0.2*zeros(size(x,1),1);
phi_top_0 = @(x,y) 0.05*exp(-((x-0.5).^2 +(y-0.5).^2)/epsilon); %zeros(size(x,1),1);

top_faces = (1 : G.faces.num)';
top_centroids = G.faces.centroids(G.faces.centroids(:, 3) == 1,:);


eta_old = eta_0(G.nodes.coords(:,1),G.nodes.coords(:,2));

phi_top_old = phi_top_0(top_centroids(:,1),top_centroids(:,2));
phi_top = phi_top_old;
eta = eta_old;

%ylim([-0.1,0.15])
%plot(eta(1:nx),'r')
F = [];
i = 1
f = figure('units','normalized','outerposition',[0 0 1 1])

for t=0:k:T
    % Calculate phi^n+1
    [phi, gradPhi_top, top_cells, top_faces, Gnew] = poissonMimeticDirich3d(G,h,phi_top, eta,[xmax,ymax,zmax], 0);
    
    node_coords = Gnew.nodes.coords(G.nodes.coords(:,3)==1, 1:2);
    centroids = Gnew.faces.centroids(G.faces.centroids(:,3)==1, 1:2);

        %Could use,  eta = node_coords(:,3);
    eta = eta(1:(nx+1)*(ny+1));

    %interpolote eta to face entroids 
    gradEta = centralDiffCartesian(eta,nx+1,ny+1, node_coords);

    etaFunc = scatteredInterpolant(node_coords, eta);
    gradEtaxFunc = scatteredInterpolant(node_coords,gradEta(:,1));
    gradEtayFunc = scatteredInterpolant(node_coords,gradEta(:,2));
    
    gradEta_centroids = [gradEtaxFunc(centroids),gradEtayFunc(centroids)];
    eta_centroids = etaFunc(centroids);
%    gradEta_centroids = zeros(nx*ny,2);

    
%    eta_centroids = zeros(nx*ny,1);
%     for i = 0:ny-1
%         gradEta_centroids(i*nx+1:(i+1)*nx,:) = 0.5*(gradEta(i*(nx+1)+2:(i+1)*(nx+1),:) ... 
%                                                    +gradEta(i*(nx+1)+1:(i+1)*(nx+1)-1,:));
%                                                
%         eta_centroids(i*nx+1:(i+1)*nx) = 0.5*(eta(i*(nx+1)+2:(i+1)*(nx+1)) ... 
%                                              +eta(i*(nx+1)+1:(i+1)*(nx+1)-1));
%     end
    
    
%     x = linspace(0,xmax,nx);
%     y = linspace(0,ymax,ny);
%     %X = top_centroids(:,1:2);
%     X = node_coords(:,1:2); 
%     [XM,YM] = meshgrid(x,x);
%     Z = griddata(X(:,1),X(:,2),sum((gradEta).^2,2),XM,YM);
%  
   %surf(XM, YM, Z)
    % Find eta^n+1
    
     eta_centroids = (eta_centroids - k*sum(gradPhi_top.*[gradEta_centroids,-ones(size(gradEta_centroids,1),1)],2));
     eta_centroid_func = scatteredInterpolant(centroids, eta_centroids);
%     For implicit eta
%     centroids = G.faces.centroids(top_faces);
%     spaceStep = diff(centroids(:,1));
%     
    %eta = findEta3d(eta, gradPhi_top, spaceStep, k, nx, ny);
    

    
    % Calculate phi^n+2
    phi_top = phi_top-k*(0.5*sum(gradPhi_top.^2,2) - g*eta_centroids);

    
%     mainDiag = 1/2*ones(nx,1);
%     mainDiag(1) = 1;
%     endMat = spdiags([mainDiag(end:-1:1),mainDiag],-1:0,nx+1,nx);
%     
    
    % Change eta back to the nodes.
%     etaMat_cent = vec2mat(eta_centroids,nx)';
%     etaMat = zeros(nx+1,ny+1);
%     etaMat(1,:) = [etaMat_cent(1,1), 0.5*(etaMat_cent(1,1:end-1)+etaMat_cent(1,2:end)),etaMat_cent(1, end)];
%     etaMat(end,:) = [etaMat_cent(end,1), 0.5*(etaMat_cent(end,1:end-1)+etaMat_cent(end,2:end)),etaMat_cent(end, end)];
%     etaMat(:,1) = [etaMat_cent(1,1); 0.5*(etaMat_cent(1:end-1,1)+etaMat_cent(2:end,1));etaMat_cent(end,1)];
%     etaMat(:,end) = [etaMat_cent(1,end); 0.5*(etaMat_cent(1:end-1,end)+etaMat_cent(2:end,end));etaMat_cent(end,end)];
%     etaMat(2:end-1,2:end-1) = 1/4*(etaMat_cent(1:end-1,1:end-1) + etaMat_cent(2:end,1:end-1) ...
%                                   +etaMat_cent(1:end-1,2:end) + etaMat_cent(2:end,2:end));
%     etaMat = etaMat;
%     eta = etaMat(:);
% %
    eta = eta_centroid_func(node_coords);
    

    x = linspace(0,xmax,nx+1);
    y = linspace(0,ymax,ny+1);
    %X = top_centroids(:,1:2);
    %X = node_coords(:,1:2); 
    [XM,YM] = meshgrid(x,y);
    etaMat = eta_centroid_func(XM,YM); 
    
    
    surf(XM,YM, etaMat)
    zlim([-0.04, 0.04])
    caxis([-0.01, 0.01])
      view(145,15)
    axis equal off
  
    
    F = [F, getframe(f)];

    
    %figure()
    %Z = griddata(centroids(:,1),centroids(:,2), phi_top,XM, YM);
    %surf(XM,YM,Z);
    %zlim([-0.2, 0.5])
    %eta = etaMat(:);


    %prep eta for next iteration
    eta = repmat(eta,nz+1,1);
    
    pause(0.0001)


end

close all

myVid = VideoWriter('WAWY.avi')
myVid.FrameRate = 8;
open(myVid);
writeVideo(myVid, F);
close(myVid)

%      xlabel('x-axis')
%      ylabel('y-axis')
%      zlabel('z-axis')