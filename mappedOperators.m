%% Description
%   This function computes the action of the inverse of the Choi operator W
%   on a set of states.
%
%   OUTPUTS:
%
%   SIGMA:      A cell containing the density matrices of the mapped states
%               of the corresponding input states.     
%   INPUTS:
%
%   W:          A Choi operator
%   
%   RHO:        A cell containing the density matrices of input states.
%               All input states should be of the same dimension.
%       
%   REQUIRES:   PartialTrace.m 
%
%   AUTHOR:     Tom Bintener
%%
function[sigma] = mappedOperators(W,rho)
    
    n = length(rho);        %Determine the number of states to be mapped.
    sigma = cell(1,n);      %Initialise the output cell.
    
    dIn = length(rho{1});   %Total input dimension of the input states.
    dOut = length(W)/dIn;   %Total output dimension of the output states.
    
    %Compute each mapped state using the inverse of the Choi.
    for i = 1:n
        temp = W*Tensor(eye(dOut), transpose(rho{i}));
        sigma{i} = dIn*PartialTrace(temp,2,[dOut,dIn]);
    end  
end
