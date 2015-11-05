clear all; close all
% mrstModule add mimetic;
run '/global/etc/matlab/mrst/startup.m'
require mimetic;

g = 9.81;

% Set up a Cartesian grid.
nx = 50;
ny = 50;
G = cartGrid([nx, ny], [1, 1]);
G = computeGeometry(G);

T = 1;
k = 0.1;

% Define initial values
h = @(x) 1;
epsilon = 5e-2;
eta_0 = @(x) (0*exp(-(x-0.2).^2/(epsilon^2)));
etat_0 = @(x) (0.01*exp(-(x-0.2).^2/(epsilon^2)));

phi_old = zeros(nx*ny,1);

top_faces = (1 : G.faces.num)';
top_centroids = G.faces.centroids(G.faces.centroids(:, 2) == 1);
etat = etat_0(top_centroids(:,1));
eta = eta_0(G.nodes.coords(:,1));


for t = 0:k:T
    phi = solverPoissonMimetic(G, h, eta, etat);
    gradPhi = centralDiffCartesian(G,phi);
    phit = backwardDiff(G,phi_old,phi,k);
    eta = -1/g*(0.5*sum(gradPhi.^2,2) + phit);
    plot(eta(end-nx+1:end));
    axis([0,1,-1,1]);
    pause()
end
    
