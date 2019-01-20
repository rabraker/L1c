% function [eta, LBRes] = l1qc(x0, b, pix_idx, opts)
%
% Solves the Denoising Basis Pursuit problem
%
% min  ||eta||_1   s.t. ||E*M*eta - b||_2 < epsilon                  (1)
%  x
%
% Where: M = inverse DCT
%        b \in R^m is a vector of length m of subsampled data.
%        E \in R^{m x n} is the sampling matrix
%        eta \in R^n is the recovered signal in the DCT domain.
%
%       Typically, m << n.
%
%  Arguments
%  ---------
%    x0:  length n vector, is the initial guess. Aften, x0 = E^T*dct(b) works
%         well.
% 
%     b:  length m vector.
%         pix_idx: index set corresponding to the samples if b. In general, if 
%         I is an n by n identity, the E = I(pix_idx, :). pix_idx must use
%         0-based indexing.
% 
%   opts: a struct (create with CsTools.l1qc_opts) with the following fields:
%     .epsilon : Denoising parameter. l1qc will enforce ||A*eta - b|| < epsilon
% 
%     .mu :   Log-barrier weight is multiplied by mu each LB iteration
%             Default is  10.
%  
%     .lbtol : Tolerance for the log-barrier iterations and by
%              default, determines the total number of log-barrier 
%              iterations. Default is 1e-3.
% 
%     .l1_tol:  The newton iterations will stop if 
%                 ||f(x_k)|_1 - |f(x_k-1)|_1 |/|f(x_k)|_1 <l1_tol
%               Default is 0, which will match the Matlab code.
% 
%     .newton_tol : The Newton iterations terminate when the 
%                   newton decrement < newton_tol (try newton_tol = lbtol).
%
%     .newton_max_iter : Maximum number of newton iterations.
%                        Default is 50.
% 
%     .cgtol:   tolerance for the conjugate gradient (CG) solver.
%           Default is 1e-8.
% 
%     .cgmaxiter: maximum iterations for the CG solver. Default is 200.
% 
%     .warm_start_cg : (0|1|2) Warm start options for CG solver. This
%                      seems only marginally succesfull and possibly detrimental.
%                  -If 0, use dx*0.
%                  -If 1, use use previously computed value of dx as warmstart.
%                  -If 2, use dx*step from previous (newton)  iteration.
% 
%     .verbose : How much to print. 
%                - If 0, nothing besides warnings are printed.
%                - If 1, prints info each log-barrier iteration.
%                - If 2, prints same info as 1, but also info each newton 
%                  iteration.
%
%
% Returns
% --------
%   eta: vector of length n, which is the recovered signal, in the DCT domain.
%   To recover the signal in the time/spatial domain, x_recovered = idct(eta).
%   LBRes : a struct with stats from the optimization with the following
%           fields:
%         .l1 : (scalar) is the achieved l1 norm of eta.
%         .total_newton_iter : is the total number of newton iterations (across
%               all log-barrier iterations).
%         .status: 0 if l1qc finished without errors, 1 otherwise.
% 
