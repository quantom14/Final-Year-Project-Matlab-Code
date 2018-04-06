%%  Descritption  Computes the genuine negativity of a density matrix
%   
%   This function computes the Genuine Multipartite Negativity of rho. 
%   
%   OUTPUTS: 
%
%   NEGATIVITY:     The genuine multipartite Negativity (GMN) contained in 
%                   a density matrix rho.
%
%   INPUTS:
%
%   RHO:            Either an input state RHO or a Choi operator (to be
%                   used in an SDP programme). 
%                   If rho is a 2D-matrix, an SDP programme computes the
%                   GMN of rho.
%                   If rho is a 3D-matrix, we sum the Negativites of the
%                   slices rho(:,:,i) across all non-trivial bipartitions
%                   of rho. 
%
%   DIM:            A vector containing the dimensions of the subsystems of
%                   rho.
%
%   REQUIRES:       SetPartition.m, PartialTranspose, TraceNorm.m, 
%                   CVX (only for 2D-matrices)
%   AUTHOR:         Tom Bintener
%%
function [negativity,rho_out] = Genuine_Negativity(rho,dim)

%Dimension checking   
    dimTot = length(rho);
    if(dimTot ~= prod(dim))
        error('The dimension of rho should equal the product of the'+... 
              'dimensions of the subsystems');
    end
    
    switch ndims(rho)
        case 2  %determine genuine multipartite negativity of rho
        
        %Determine the number of parties
        n = length(dim);
        %First we store all non-trivial bipartitions "A:B" in a cell 
        partition = SetPartition(n,2);
        %Number of non-trivial bipartitions
        l = length(partition);    
            
        cvx_begin sdp quiet
        
        %initialise 3D SDP variable 
        variable rho_l(dimTot,dimTot,l) hermitian semidefinite;
        
        %we sum the Negativities across all bipartitions in the different
        %slices of rho_l
        neg = 0;
        for i = 1:l
            pt = PartialTranspose(rho_l(:,:,i),partition{i}{1},dim);
            neg = neg +0.5* (TraceNorm(pt)-trace(rho_l(:,:,i) ) );
        end
        
        minimise neg
        
        subject to
            %The slices of the SDP variable have to equal the input state
            tot = 0;
            for i = 1:l
                tot = tot + rho_l(:,:,i);
            end
                tot == rho;
                
        cvx_end
        
        rho_out = rho_l;
        negativity = cvx_optval;  
        
        case 3  %determine genuine multipartite negativity of rho 
                %to be used inside SDP programme. 
        %Determine the number of parties
        n = length(dim);
        %First we store all non-trivial bipartitions "A:B" in a cell 
        partition = SetPartition(n/2,2);
        %Number of non-trivial bipartitions
        l = length(partition);
        
        %In this case we assume dim contains the dimensions of all in- and
        %output systems: dim = [dOut,dIn] and we need to permute the slices.
        perm(1:2:n) = 1:n/2;        %permutation vector that contains the 
        perm(2:2:n) = (n/2)+1:n;    %indices, [1,n+1,2,n+2,....n,2*n] 
        
        dOut = dim(1:n/2);      %vector containing output dimensions.
        dIn= dim(n/2+1:n);      %vector containing input dimensions.
        DIM = zeros(1,n/2);                      %vector containing products of  
        for i = 1:n/2; DIM(i)=dIn(i)*dOut(i);end %pairs of inputs and outputs.        
                
        %It is assumed that rho is a Choi operator, it is a 3D matrix
        %and in the W(:,:,i) slices we compute the negativity across
        %the different bipartitions.
        neg = 0;
        for i = 1:l
            permuted = PermuteSystems(rho(:,:,i),perm,dim);
            pt = PartialTranspose(permuted,partition{i}{1},dim);
            neg = neg +0.5* (TraceNorm(pt)-trace(permuted));
        end
        
        negativity = neg;
    end           
end

