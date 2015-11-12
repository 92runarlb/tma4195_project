clear; close all

mrstModule add mimetic
g = 9.81;

% Set up a Cartesian grid.
nx = 60;
ny = 20;
dx = 1/nx;
G = cartGrid([nx, ny], [1, 1]);
G = computeGeometry(G);

T = 5;
k = 0.01;

% Define initial values
h = @(x) ones(size(x,1),1);
epsilon = 5e-3;
eta_0 = @(x) zeros(size(x,1),1);
phi_top_0 = @(x) 0.05*exp(-(x-0.25).^2/epsilon); %zeros(size(x,1),1);

top_faces = (1 : G.faces.num)';
top_centroids = G.faces.centroids(G.faces.centroids(:, 2) == 1);

eta_old = eta_0(G.nodes.coords(:,1));
phi_top_old = phi_top_0(top_centroids);
phi_top = phi_top_old;
eta = eta_old;
figure()
ylim([-0.1,0.15])
plot(eta(1:nx),'r')

for t=0:k:T
    % Calculate phi^n+1
    [phi, gradPhi_top, top_cells, top_faces, G] = poissonMimeticDirich(G,h,phi_top, eta,0);
    % Remove extra etas
    eta = eta(1:nx+1);
    % Calculate the derivative
    etax = diff(eta)./dx;
    % Center the etas at the face centers
    eta = 0.5*(eta(2:end)+eta(1:end-1));
    % Find eta^n+1
    %eta = eta - k*sum(gradPhi_top.*[etax,-ones(size(etax,1),1)],2);
    centroids = G.faces.centroids(top_faces);
    spaceStep = diff(centroids(:,1));
    eta = findEta(eta, gradPhi_top, spaceStep, k);
    
    plot(eta)
    
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
    plotCellData(G,phi)
    colorbar