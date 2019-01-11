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
%     b:  length m vector.
%         pix_idx: index set corresponding to the samples if b. In general, if 
%         I is an n by n identity, the E = I(pix_idx, :). pix_idx must use
%         0-based indexing.
%   opts: a struct (create with CsTools.l1qc_opts) with the following fields:
%         .epsilon : denoising tolerance (as in (1) above). 
%         .mu :      log-barrier weight is multiplied by mu each LB iteration
%                    (try 10).
%         .cgtol:    tolerance for the conjugate gradient (CG) solver (try 1e-8).
%         .cgmaxiter: maximum iterations for the CG solver (try 200).
%         .lbtol :   Determines the total number of log-barrier iterations
%                    (try  1e-3).
%         .newton_tol : newton iterations terminate when the 
%                       newton decrement < newton_tol (try newton_tol = lbtol).
%         .newton_max_iter : Maximum number of newton iterations (try 50).
%         .verbose : How much to print. Silent if verbose =0, prints info each
%         log-barrier iteration if verbose =1, and also each newton iteration
%         if verbose = 2.
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
