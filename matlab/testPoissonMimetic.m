clear; close all

mrstModule add mimetic
g = 9.81;

% Set up a Cartesian grid.
nx = 120;
ny = 20;
G = cartGrid([nx, ny], [1, 1]);
G = computeGeometry(G);

T = 0.001;
k = 0.001;

% Define initial values
h = @(x) ones(size(x,1),1);
epsilon = 5e-3;
eta_0 = @(x) k*exp(-(x-0.5).^2/epsilon) %zeros(size(x,1),1);
etat_0 = @(x) 0.01*k*cos(pi*x);

phi_top_old = zeros(nx,1);

top_faces = (1 : G.faces.num)';
top_centroids = G.faces.centroids(G.faces.centroids(:, 2) == 1);
etat = etat_0(top_centroids(:,1));
eta_old = eta_0(G.nodes.coords(:,1));
eta = eta_old;
figure()
ylim([-1,0.5])
hold on
plot(eta(1:nx),'r')

for t=0:k:T
    [phi, gradPhiSqr_top, top_cells, G] = poissonMimetic(G,h,eta,etat);
    phit_top = (phi(top_cells)-phi_top_old)/k;
    
    eta = 1/g*(phit_top+0.5*gradPhiSqr_top);
    eta = [eta(1); eta];
    eta = repmat(eta,ny+1,1);
    etat = (eta(1:nx)-eta_old(1:nx))./k;
    phi_top_old = phi(top_cells);
    plot(eta(1:nx))

    G = cartGrid([nx, ny], [1, 1]);
    G = computeGeometry(G);
    pause()
end

    figure
    plotCellData(G,phi)
    colorbar