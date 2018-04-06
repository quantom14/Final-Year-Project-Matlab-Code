%% Description
%   Computes the tensor product of two states |n> and |m> of dimension d. 
%   Inputs:
%           n and m :   such that n,m is in [0,d-1]
%           d:          dimension of state
%   Output: 
%           tensor product of two states |n> and |m>
% authors: Tom Bintener
%%
function[prod] = nn(n,m,d)
    if(n < 0 || m < 0 || n >= d|| m >= d)
        error('n and m have to be in the interval [0,d-1]');
    else
    %Identity matrix 
    I = eye(d);
    %tensor product
    prod = kron(I(:,n+1),I(:,m+1));
    end
end