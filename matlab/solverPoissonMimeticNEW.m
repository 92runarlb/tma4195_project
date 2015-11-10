function [phi, gradPhiSqr] = poissonMimetic(G, h, eta, etat);
    DIM = G.cartDims;
    nx = DIM(1);
    ny = DIM(2);
    nc = G.cells.num;
    nf = G.faces.num;
    half_faces = G.cells.faces(:, 1);
    nhf = numel(half_faces);