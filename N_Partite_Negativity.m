%% Description
%   This function computes the minimum N-Partite Negativity of the 
%   normalised Choi operator of a channel, over all channels that map each 
%   state rho{i} into sigma{i}.
%
%   OUTPUTS:
%
%   NEGATIVITY: The minimum negativity required to map rho{i} into sigma{i}
%   CHOI:       The optimal Choi operator which maps the states with the
%               least negativity.
%   CHOI_PERM:  The Choi operator where the subsystems have been permuted. 
%               (see CVX part of code)
%       
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
%   REQUIRES:   SetPartition.m, Multi_Negativity.m, PartialTrace.m, 
%               mappedOperators.m, CVX
%   AUTHOR:     Tom Bintener
%%
function [negativity,Choi,Choi_perm] = N_Partite_Negativity(rho,dIn,sigma,dOut)

%% Initialising required variables
    n = length(dIn);            %number of parties
    m = length(rho);            %number of states to be mapped
    dInTot = length(rho{1});    %total input dimension
    dOutTot = length(sigma{1}); %total output dimension
    dim = [dOut,dIn];           %vector containing all sub-dimensions
    
    perm(1:2:2*n) = 1:n;        %permutation vector that contains the 
    perm(2:2:2*n) = n+1:2*n;    %indices, [1,n+1,2,n+2,....,n,2*n] 
    
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

    %We need a SDP variable for the Choi
    variable W(dInTot*dOutTot,dInTot*dOutTot) sparse hermitian semidefinite;
 
    %We want to minimise the N-Partite Negativity of W across all NT
    %bipartitions.
  
    minimise Multi_Negativity(PermuteSystems(W,perm,dim),DIM);
    
    subject to
        
        %We want W to be trace-preserving, we trace over the output
        %subsystems.
        PartialTrace(W,1,[dOutTot,dInTot]) == eye(dInTot)/dInTot; 
    
        %Using the inverse of the Choi we impose the mapping of rho{i} into
        %sigma{i}
        map = mappedOperators(W,rho(i));
        for i = 1:m
             map{i} == sigma{i};  
        end
cvx_end

%% Outputs
    negativity = cvx_optval;
    Choi = W;
    
    %Optional representation of the Choi (permuted version)
    Choi_perm = PermuteSystems(W,perm,dim);

end