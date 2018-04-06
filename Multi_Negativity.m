%%  Descritption  Computes the negativity of a N-partite density matrix
%   
%   To compute the N-Partite Negativity of rho, this function requires a
%   2D matrix RHO and a vector containing the dimensions of the N
%   subsystems.
%   
%   OUTPUTS:
%
%   NEGATIVITY:     The N-Partite Negativity contained in rho, the 
%                   Negativity is computed across all non trivial 
%                   bipartitions.
%   INPUTS:
%
%   RHO:            A density matrix 
%
%   DIM:            A vector containing the dimensions of the subsystems of
%                   rho.
%
%   REQUIRES:       PartialTranspose.m, TraceNorm.m, SetPartition.m
%   AUTHOR:         Tom Bintener
%%
function negativity = Multi_Negativity(rho,dim)
%Determine the number of parties
    n = length(dim);
%Dimension checking   
    dimTot = length(rho);
    if(dimTot ~= prod(dim))
        error('The dimension of rho should equal the product of the'+... 
              'dimensions of the subsystems');
    end

    %First we store all non-trivial bipartitions "A:B" in a cell 
    partition = SetPartition(n,2);
    %Number of non-trivial bipartitions
    l = length(partition);

    %We go through the permutations to add up the Negativities across
    %all bipartitions
    tr_rho = 0;
    for i = 1:l                      
        %The permutations of the labels are stored in cells in the
        %variable partition. We want to calculate the partial tranpose
        %of rho over the subsystem A( or B) where A and B are the
        %non-trivial bipartitions of rho.

        %Each iteration computes the Negativity across a new 
        %bipartition and adds it to the previous ones.         
        tr_rho=tr_rho+ TraceNorm(PartialTranspose(rho,partition{i}{1},dim));
    end

    %We have only added the tracenorms of the partial transposes for
    %each permutation, now we subtract l times the trace of rho and
    %divide by 2. 
    negativity = 1/l * 0.5 * (tr_rho - l* trace(rho));        
    
end

