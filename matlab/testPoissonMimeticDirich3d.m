clear; close all

mrstModule add mimetic
g = 9.81;

% Set up a Cartesian grid.
nx = 20;
ny = 16;
nz = 20;
dx = 1/nx;
G = cartGrid([nx, ny, nz], [1,1,1]);
G = computeGeometry(G);

T = 5;
k = 0.01;

% Define initial values
h = @(x,y) 0.2*cos(6*pi*x)+ones(size(x,1),1);
epsilon = 5e-3;
eta_0 = @(x,y) 0.2*sin(pi*x) + 0.2*cos(pi*y);
phi_top_0 = @(x,y) 0.05*exp(-((x-0.5).^2 +(y-0.5).^2)/epsilon); %zeros(size(x,1),1);

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
    [phi, gradPhi_top, top_cells, top_faces, Gnew] = poissonMimeticDirich3d(G,h,phi_top, eta,1);
    
    node_coords = Gnew.nodes.coords(G.nodes.coords(:,3)==1,:);
    eta = node_coords(:,3);
    %interpolote eta to face entroids
    % Could use, eta = eta(1:(nx+1)*(ny+1));
    gradEta = centralDiffCartesian(eta,nx+1,ny+1, node_coords);
    gradEta_centroids = zeros(nx*ny,2);
    for i = 0:ny-1
        gradEta_centroids(i*nx+1:(i+1)*nx,:) = 0.5*(gradEta(i*(nx+1)+2:(i+1)*(nx+1),:) ... 
                                                   +gradEta(i*(nx+1)+1:(i+1)*(nx+1)-1,:));
    end
    
    Ngrid = 40;
    x = linspace(0,1,Ngrid);
    %X = top_centroids(:,1:2);
    X = node_coords(:,1:2); 
    [XM,YM] = meshgrid(x,x);
     Z = griddata(X(:,1),X(:,2),sum((gradEta).^2,2),XM,YM);
 
      surf(XM, YM, Z)
    % Find eta^n+1
    eta = eta - k*sum(gradPhi_top.*[gradEta_centroids',-ones(size(etax,1),1)],2);
    centroids = G.faces.centroids(top_faces);
    spaceStep = diff(centroids(:,1));
    eta = findEta3d(eta, gradPhi_top, spaceStep, k, nx, ny);
    
    %plot(eta)
    
    % Calculate phi^n+2
    phi_top = phi_top-k*0.5*(sum(gradPhi_top.^2,2) - g*eta);
    
    % Change eta back to the nodes.
    eta = 0.5*(eta(2:end)+eta(1:end-1));
    deta1 = eta(2)-eta(1);         %Should divide on dx, but cancels:
    deta2 = eta(end)-eta(end-1);
    eta1 = eta(1)-0.5*deta1;     % HERE
    eta2 = eta(end) + 0.5*deta2;
    eta = [eta1; eta; eta2];
    plot(eta)
    ylim([-0.1,0.15]);
    %prep eta for next iteration
    eta = repmat(eta,ny+1,1);
    % Calculate new G
    G = cartGrid([nx, ny], [1, 1]);
    G = computeGeometry(G);
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