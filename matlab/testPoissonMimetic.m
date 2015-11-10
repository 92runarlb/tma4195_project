clear; close all
g = 9.81;

% Set up a Cartesian grid.
nx = 40;
ny = 40;
G = cartGrid([nx, ny], [1, 1]);
G = computeGeometry(G);

T = 0;
k = 0.1;

% Define initial values
h = @(x) 0.0*cos(6*pi*x);
epsilon = 5e-2;
eta_0 = @(x) ones(size(x,1),1); %+0.3*sin(pi*x);
etat_0 = @(x) ones(size(x,1),1);

phi_old = zeros(nx*ny,1);

top_faces = (1 : G.faces.num)';
top_centroids = G.faces.centroids(G.faces.centroids(:, 2) == 1);
etat = etat_0(top_centroids(:,1));
eta = eta_0(G.nodes.coords(:,1));


[phi, gradPhiSqr] = poissonMimetic(G,h,eta,etat);