%% Descrition
% Computes the minimum negativity N of the normalized Choi operator of a
% channel, over all channels AB -> CD that map each state setrhoAB{i}
% into setsigmaCD{i}; if a value for epsilon is provided an approximate 
% mapping is considered. 
% required arguments: 
% setrhoAB:     Cell array containing the set of input states; all such
%               states should have the same input dimension dInTot
% dA:           Dimension of input system A
% setsigmaCD:   Cell array containing the set of output states; all such 
%               states should have the same output dimension dOutTot 
% dC:           Dimension of output system C
% optional arguments:
% epsilon:      When a value for epsilon is provided an approximate mapping
%               is considered 
% Output:
% neg:          Minimum negativity of a normalised Choi state isomorphic to
%               a bipartite channel AB->CD that maps each input state
%               setrho{j} into the output state setsigma{j}
% W:            an optimal normalised Choi isomorphic state that achieves
%               the minimum negativity
% Additional package requirements:
%     + CVX     http://cvxr.com/cvx/
%     + Qetlab  http://www.qetlab.com
%               - Negativity
%               - Swap
%               - PartialTrace
%               - PartialTranspose 
%               - Tensor
% authors: Marco Piani, Tom Bintener
%% 
function [neg,W_out] = NegMapEns2(setrhoAB,dA,setsigmaCD,dC,varargin)
%% Initialization of values and error checking
% Calculate dimensions that are not given already
%number of input states
n = length(setrhoAB);
%total input dimension (the dimension of each state setrhoAB{j}
dInTot = length(setrhoAB{1});
%total output dimension (the dimension of each state setsigmaAB{j}
dOutTot = length(setsigmaCD{1});
 if n == length(setsigmaCD)
     for i = 2:n
         if length(setrhoAB{i}) ~= dInTot 
             error('All input states should have the same dimension dInTot');
         elseif length(setsigmaCD{i}) ~= dOutTot
             error('All output states should have the same dimension dOutTot');
         end 
     end
 else
     error('Number of input and output systems must be equal.')
 end
%dimension of input system B
dB = dInTot/dA ;
%dimension of output system D
dD = dOutTot/dC ;
% Set optional argument 
    if nargin == 4    
    elseif nargin == 5 
        if length(varargin{1}) == 1
            epsilon = varargin{1};               
        else
            error('epsilon should be a number.')
        end
    else 
        error('Unexpected number of arguments.') 
    end    
%% CVX part    
cvx_begin sdp quiet
    
%W is the Choi SDP variable
    variable W(dInTot*dOutTot,dInTot*dOutTot) sparse hermitian semidefinite;
%X is the negativity SDP variable
    variable X(dInTot*dOutTot,dInTot*dOutTot) sparse hermitian semidefinite;
%E is the SDP variable 
    variable E(dOutTot,dOutTot) hermitian sparse semidefinite;
%target funtion (Negativity) to be minimized
    minimise 0.5 * (trace(X)-1);
%conditions
    subject to
%conditions for the trace norm
        -X <= PartialTranspose(W,[1,3],[dC,dD,dA,dB]) <= X;
%Condition on E
        if nargin == 5
            trace(E) <= epsilon;
        end
%trace-preservation of the channel imposed via condition on partial trace
%of W; trace is done on output systems C and D
        PartialTrace(W,[1,2],[dC,dD,dA,dB]) == eye(dInTot)/dInTot;
%Imposing the mapping using W^(-1) of setrhoAB{i} into setsigmaCD{i}
    switch nargin 
        case 4
             for j=1:n                 
                 dInTot*PartialTrace(W*Tensor(eye(dOutTot),transpose(setrhoAB{j})),...
                 [3,4],[dC,dD,dA,dB]) == setsigmaCD{j};                 
              end
        case 5
            for j=1:n                   
                cond = (dInTot*PartialTrace(W*Tensor(eye(dOutTot),transpose(setrhoAB{j})),...
                   [3,4],[dC,dD,dA,dB]) - setsigmaCD{j});                   
                  -E <= cond+cond' <= E ;
                  -E <= -1i*(cond-cond') <= E ;                  
            end
    end       
%end of CVX
cvx_end
%minimum negativity required
neg = cvx_optval;
%Choi of the (an) optimal map
W_out = W;
end

