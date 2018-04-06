%% Descrition
% Gives the projector onto the positive subspace of a Hermitian matrix
%   Input: 
%          A:       Hermitian matrix
%   Output
%          posproj: projector onto the positive subspace of A
%   authors: Marco Piani, Tom Bintener
function [posproj] = PosProj(A)
if(issparse(A))
    [V,D] = eigs(A);
else
    [V,D] = eig(A);
end
    posproj = V*((D+abs(D))>0)*V';
end

