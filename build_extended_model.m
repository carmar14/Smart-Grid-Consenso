function [Az,Bz,Cz] = build_extended_model(A,B,C)

nx = size(A,1);
ny = size(C,1);

Az = [A zeros(nx,ny);
     -C eye(ny)];

Bz = [B; zeros(ny,size(B,2))];

Cz = [C zeros(ny,ny)];

end