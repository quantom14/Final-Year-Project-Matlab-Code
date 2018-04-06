%% Descritpion
%   This function computes an improved lower bound to the Negativity
%   required for physical transformations, by optimizing over new input
%   states on which a map is partially applied via its Choi-operator.
%   Inputs:
%   Choi:       A Choi operator which encodes some physical transformation
%   dim:        A vector of the dimensions [in A, in B, out C, out D]         
%   varargin:   
%       max_iterations: maximum number of iterations before the
%                       algorithm stops. If set to 0, it is assigned
%                       the value 10.
%       error:          when the improvement, i.e. the absolute value
%                       of the difference between two steps is less
%                       than or equal to 'error' the algorithm stops. When
%                       set to 0, the default error is taken to be 10^(-4)
%       vec:            either true or false, default is false.
%                       When true the function returns the negativities
%                       at each step, when false just the final value.
%   Outputs:
%       neg_out:        Improved Negativity, that is obtained by optimizing
%                       over new inputs.
%       rho_improved:   Optimal state which achieves the best Negativity.
%   Additional functions:   PosProj
%   Additional package requirements:
%     + CVX     http://cvxr.com/cvx/
%     + Qetlab  http://www.qetlab.com
%               - PartialTranspose
%               - Negativity
%               - PartialTrace 
%               - Tensor
%               - PermuteSystems
% authors: Marco Piani, Tom Bintener
%%
function [neg_out, rho_improved,r] = improved_bound1(Choi,dim,var_argin)
    %varargin: max_iterations, error, output final or vector
    dim_cell = num2cell(dim); 
    [dA,dB,dC,dD] = dim_cell{:};
    
    dInTot = dA*dB;
    dOutTot = dC*dD;    
    
    argin_cell = num2cell(var_argin);
    [max_iterations, error, vec] = argin_cell{:};
    
    if max_iterations == 0, max_iterations = 10; end
    if error == 0, error = 10^(-4); end 
    
    rhoOut = Choi;
    
    iteration = 1;
    neg(1) = real(Negativity(PermuteSystems(rhoOut,[1,3,2,4],[dC,dD,dA,dB]),dA*dC));

    not_converged = true;
while not_converged    
    %determine optimal M
    M = PosProj(-PartialTranspose(rhoOut,[1,3],[dC dD dA dB]));

    cvx_begin sdp quiet
    cvx_precision best
        variable rho(dInTot*dInTot,dInTot*dInTot) hermitian semidefinite
        %partial transpose of rho and 
        pt_rho = PartialTranspose(rho,[1,2],[dA dB dA dB]);
        
        %partial action of map on rho and its partial transpose
        rhoOut = dInTot * PartialTrace(Tensor(Choi,eye(dInTot))*...
            Tensor(eye(dOutTot),pt_rho),[3,4],[dC,dD,dA,dB,dA,dB]);            
        pt_rhoOut = PartialTranspose(rhoOut,[1,3],[dC,dD,dA,dB]);           
        maximize -real(trace(pt_rhoOut * M))
        subject to           
            PartialTranspose(rho,[1,3],[dA,dB,dA,dB]) >= 0;                 
            trace(rho) == 1;                                                
    cvx_end
    rho_improved = rho;
    r = rhoOut;
    iteration = iteration + 1;
    
    neg(iteration) = real(Negativity(...
                PermuteSystems(rhoOut,[1,3,2,4],[dC,dD,dA,dB]),dA*dC));   
        %Compare current value with the previous one
        if abs(neg(iteration)- neg(iteration-1)) <= error 
            not_converged = false;
            fprintf('Converged at iteration %d \n' ,iteration)
        end 
        if iteration > max_iterations
            not_converged = false;
            fprintf('Stopped at maximal iteration allowed.\n')
        end
end
    if  vec == 0
        neg_out = neg(length(neg));
    else 
        neg_out = neg;
    end
end