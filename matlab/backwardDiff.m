function phit = backwardDiff(G,phi_old,phi_new, k)
    phit = 1/k*(phi_new - phi_old);    

    
    %         DIM = G.cartDims;
%    Nx = DIM(1);
%    Ny = DIM(2);
%
%     X = G.cells.centroids;
% 
%     x = linspace(0,1,Ngrid);
%     [XM,YM] = meshgrid(x,x);
%     Z = griddata(X(:,1),X(:,2),phit,XM,YM);
% 
%     surf(XM, YM, Z)
%     xlabel('x-axis')
%     ylabel('y-axis')
%     zlabel('z-axis')
end