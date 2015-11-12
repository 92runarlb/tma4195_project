clear; close all

mrstModule add mimetic
g = 9.81;

% Set up a Cartesian grid.
xmax = 1;
ymax = 1;
zmax = 1;
nx = 40;
ny = 40;
nz = 20;
dx = 1/nx;
G = cartGrid([nx, ny, nz], [xmax, ymax, zmax]);
G = computeGeometry(G);

T = 5;
k = 0.1;

% Define initial values
h = @(x,y) -ones(size(x,1),1);
epsilon = 1e-2;
eta_0 = @(x,y) 0*0.2*sin(pi*x) ;
phi_top_0 = @(x,y) 0.005*exp(-((x-0.5).^2 +(y-0.5).^2)/epsilon); %zeros(size(x,1),1);

top_faces = (1 : G.faces.num)';
top_centroids = G.faces.centroids(G.faces.centroids(:, 3) == 1,:);


eta_old = eta_0(G.nodes.coords(:,1),G.nodes.coords(:,2));

phi_top_old = phi_top_0(top_centroids(:,1),top_centroids(:,2));
phi_top = phi_top_old;
eta = eta_old;
figure()
%ylim([-0.1,0.15])
%plot(eta(1:nx),'r')

for t=0:k:T
    % Calculate phi^n+1
    [phi, gradPhi_top, top_cells, top_faces, Gnew] = poissonMimeticDirich3d(G,h,phi_top, eta, 0);
    
    node_coords = Gnew.nodes.coords(G.nodes.coords(:,3)==1,:);
    eta = node_coords(:,3);
    gradEta = centralDiffCartesian(eta,nx+1,ny+1, node_coords);
    eta_centroids = zeros(nx*ny,1);
    for i = 0:ny-1
        gradEta_centroids(i*nx+1:(i+1)*nx,:) = 0.5*(gradEta(i*(nx+1)+2:(i+1)*(nx+1),:) ... 
                                                   +gradEta(i*(nx+1)+1:(i+1)*(nx+1)-1,:));
                                               
        eta_centroids(i*nx+1:(i+1)*nx) = 0.5*(eta(i*(nx+1)+2:(i+1)*(nx+1)) ... 
                                             +eta(i*(nx+1)+1:(i+1)*(nx+1)-1));
    end
    %interpolote eta to face entroids
    % Could use, eta = eta(1:(nx+1)*(ny+1));
       
    gradPhiMatx = zeros((nx+1),(ny+1));
    gradPhiMat_cent= vec2mat(gradPhi_top(:,1),nx)';
    
    gradPhiMatx(1,:) = [gradPhiMat_cent(1,1), 0.5*(gradPhiMat_cent(1,1:end-1)+gradPhiMat_cent(1,2:end)),gradPhiMat_cent(1, end)];
    gradPhiMatx(end,:) = [gradPhiMat_cent(end,1), 0.5*(gradPhiMat_cent(end,1:end-1)+gradPhiMat_cent(end,2:end)),gradPhiMat_cent(end, end)];
    gradPhiMatx(:,1) = [gradPhiMat_cent(1,1); 0.5*(gradPhiMat_cent(1:end-1,1)+gradPhiMat_cent(2:end,1));gradPhiMat_cent(end,1)];
    gradPhiMatx(:,end) = [gradPhiMat_cent(1,end); 0.5*(gradPhiMat_cent(1:end-1,end)+gradPhiMat_cent(2:end,end));gradPhiMat_cent(end,end)];
    gradPhiMatx(2:end-1,2:end-1) = 1/4*(gradPhiMat_cent(1:end-1,1:end-1) + gradPhiMat_cent(2:end,1:end-1) ...
                                  +gradPhiMat_cent(1:end-1,2:end) + gradPhiMat_cent(2:end,2:end));
                              
    gradPhiMaty = zeros((nx+1),(ny+1));
    gradPhiMat_cent= vec2mat(gradPhi_top(:,2),nx)';
    
    gradPhiMaty(1,:) = [gradPhiMat_cent(1,1), 0.5*(gradPhiMat_cent(1,1:end-1)+gradPhiMat_cent(1,2:end)),gradPhiMat_cent(1, end)];
    gradPhiMaty(end,:) = [gradPhiMat_cent(end,1), 0.5*(gradPhiMat_cent(end,1:end-1)+gradPhiMat_cent(end,2:end)),gradPhiMat_cent(end, end)];
    gradPhiMaty(:,1) = [gradPhiMat_cent(1,1); 0.5*(gradPhiMat_cent(1:end-1,1)+gradPhiMat_cent(2:end,1));gradPhiMat_cent(end,1)];
    gradPhiMaty(:,end) = [gradPhiMat_cent(1,end); 0.5*(gradPhiMat_cent(1:end-1,end)+gradPhiMat_cent(2:end,end));gradPhiMat_cent(end,end)];
    gradPhiMaty(2:end-1,2:end-1) = 1/4*(gradPhiMat_cent(1:end-1,1:end-1) + gradPhiMat_cent(2:end,1:end-1) ...
                                  +gradPhiMat_cent(1:end-1,2:end) + gradPhiMat_cent(2:end,2:end));
                              
    gradPhiMatz = zeros((nx+1),(ny+1));
    gradPhiMat_cent= vec2mat(gradPhi_top(:,3),nx)';
    
    gradPhiMatz(1,:) = [gradPhiMat_cent(1,1), 0.5*(gradPhiMat_cent(1,1:end-1)+gradPhiMat_cent(1,2:end)),gradPhiMat_cent(1, end)];
    gradPhiMatz(end,:) = [gradPhiMat_cent(end,1), 0.5*(gradPhiMat_cent(end,1:end-1)+gradPhiMat_cent(end,2:end)),gradPhiMat_cent(end, end)];
    gradPhiMatz(:,1) = [gradPhiMat_cent(1,1); 0.5*(gradPhiMat_cent(1:end-1,1)+gradPhiMat_cent(2:end,1));gradPhiMat_cent(end,1)];
    gradPhiMatz(:,end) = [gradPhiMat_cent(1,end); 0.5*(gradPhiMat_cent(1:end-1,end)+gradPhiMat_cent(2:end,end));gradPhiMat_cent(end,end)];
    gradPhiMatz(2:end-1,2:end-1) = 1/4*(gradPhiMat_cent(1:end-1,1:end-1) + gradPhiMat_cent(2:end,1:end-1) ...
                                  +gradPhiMat_cent(1:end-1,2:end) + gradPhiMat_cent(2:end,2:end));
    
    gradPhiMatx = gradPhiMatx(:);
    gradPhiMaty = gradPhiMaty(:);
    gradPhiMatz = gradPhiMatz(:);
    gradPhi_nodes = [gradPhiMatx, gradPhiMaty, gradPhiMatz]; 
%     x = linspace(0,xmax,nx);
%     y = linspace(0,ymax,ny);
%     %X = top_centroids(:,1:2);
%     X = node_coords(:,1:2); 
%     [XM,YM] = meshgrid(x,x);
%     Z = griddata(X(:,1),X(:,2),sum((gradEta).^2,2),XM,YM);
%  
   %surf(XM, YM, Z)
    % Find eta^n+1
    
     eta = eta - k*sum(gradPhi_nodes.*[gradEta,-ones(size(gradEta,1),1)],2);

%     For implicit eta
%     centroids = G.faces.centroids(top_faces);
%     spaceStep = diff(centroids(:,1));
%     
    %eta = findEta3d(eta, gradPhi_top, spaceStep, k, nx, ny);
    

    
    % Calculate phi^n+2
    phi_top = phi_top-k*0.5*(sum(gradPhi_top.^2,2) - g*eta_centroids);
    

    x = linspace(0,xmax,nx+1);
    y = linspace(0,ymax,ny+1);
    %X = top_centroids(:,1:2);
    X = node_coords(:,1:2); 
    [XM,YM] = meshgrid(x,y);
    etaMat = griddata(X(:,1),X(:,2),eta,XM,YM)'; %OBS Transpose
    surf(XM,YM,etaMat)
    zlim([-0.2, 0.2])

    %prep eta for next iteration
    eta = repmat(eta,nz+1,1);
    pause(0.0001)

end
figure
plotCellData(Gnew, phi);
xlabel('x')
ylabel('y')
zlabel('z')



%      xlabel('x-axis')
%      ylabel('y-axis')
%      zlabel('z-axis')