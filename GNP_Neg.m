%% Description
%   This function computes the minimal Genuine Multipartite Negativity 
%   (GMN) of the normalised Choi operator of a channel, over all channels  
%   that map each state rho{i} into sigma{i}.
%   
%   OUTPUTS:
%
%   NEG:        The minimum negativity required to map rho{i} into sigma{i}
%   JOUT:       The optimal Choi operator which maps the states with the
%               least negativity.
%   WOUT:       The Choi containing all the different slices that sum up to
%               J.
%   INPUTS:
%   
%   RHO:        A cell containing the density matrices of the corresponding 
%               input states 
%   DIN:        A vector containing all the dimensions of the subsystems of
%               rho{i}. (We assume all states have the same subsystem
%               dimensions.
%   SIMGA:      A cell containing the density matrices of the corresponding
%               output states.
%   DOUT:       A vector containing all the dimensions of the subsystems of
%               sigma{i}. (We assume all states have the same subsystem
%               dimensions.)
%       
%   REQUIRES:   SetPartition.m, Genuine_Negativity.m, PartialTrace.m, 
%               mappedOperators.m, CVX
%   AUTHOR:     Tom Bintener
%%

function [neg,Jout,Wout] = GNP_Neg(rho,dIn,sigma,dOut)

%% Initialising required variables.
    n = length(dIn);                %number of parties
    m = length(rho);                %number of states to be mapped
    l = length(SetPartition(n,2));  %number of NT bipartitions
    dInTot = length(rho{1});    %total input dimension
    dOutTot = length(sigma{1}); %total output dimension
    dim = [dOut,dIn];           %vector containing all sub-dimensions
    
    perm(1:2:2*n) = 1:n;        %permutation vector that contains the 
    perm(2:2:2*n) = n+1:2*n;    %indices, [1,n+1,2,n+2,....n,2*n] 
    
    DIM = zeros(1,n);                      %vector containing products of  
    for i = 1:n; DIM(i)=dIn(i)*dOut(i);end %pairs of inputs and outputs.
    
    %% Error checking 
        if(n ~= length(dOut))
            error('Number of input systems must match number of output'+... 
            'systems dIn and dOut cannot differ in length.');
        end
        if(m ~= length(sigma))
            error('Number of states in rho must match the number of states in sigma.');
        end
        
        for i = 1:length(rho)
            if(length(rho{i}) ~= prod(dIn))
                error('The product of the dimensions in dIn must match the length of rho{i}');
            end
            if(length(sigma{i}) ~= prod(dOut))
                error('The product of the dimensions in dOut must match the length of sigma{i}');
            end
        end
    %% --------------------------------------------------------------------   
  
%% CVX    
cvx_begin sdp quiet

    %We need a SDP variable for the total Choi and one that contains l
    %slices of which we will compute the Negativity across a different
    %bipartition.
    variable J(dInTot*dOutTot,dInTot*dOutTot) sparse hermitian semidefinite;
    variable W(dInTot*dOutTot,dInTot*dOutTot,l) sparse hermitian semidefinite;
    %We want to minimise the GMN of W across all NT bipartitions.
    
    minimise Genuine_Negativity(W,dim);
    
    subject to
        
        %The sum of the slices should equal the total Choi.
        tot = 0;
        for i = 1:l
            tot = tot + W(:,:,i);
        end
            tot == J;
        
        %We want J to be trace-preserving, we trace over the output
        %subsystems.
        PartialTrace(J,1,[dOutTot,dInTot]) == eye(dInTot)/dInTot; 
    
        %Using the inverse of the Choi we impose the mapping of rho{i} into
        %sigma{i}
        map = mappedOperators(J,rho);
        for i = 1:m
            map{i} == sigma{i};  
        end
cvx_end

%% Outputs
    neg = cvx_optval;
    Jout = J;
    Wout = W;
end