function [x, LBRes] = l1qc_dct(N, M, b, pix_idx, opts)
% [x, LBRes] = l1qc_dct(N, M, b, pix_idx, opts)
%
% Solves the Denoising Basis Pursuit problem with a DCT sparsifying
% basis. Consider given a true signal $z$, and a subsampled vector
% $b$ such that $b=E*z$, then let $x=M^T*z$. $M$ is the IDCT
% transformation matrix and $E$ is an identity with rows removed.
% 
% We then solve the optimization
% 
%
% min  ||x||_1   s.t. ||E*M*x - b||_2 < epsilon                  (1)
%  x
%
% Where: M = inverse DCT
%        $E \in\mathbb{R}^{m x n} is the sampling matrix
%        $b \in\mathbb{R}^m$ is a vector of length m of subsampled
%        data. Given the true signal $z$ in standard coordinates, b=E*z.
%        ${\eta \in \mathbb{R^n}$ is the recovered signal in the DCT domain.
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
%    pix_idx: A vector of indeces of where the measurements in b
%             were taken. Given the true signal, x, then it should
%             hold that b = x(pix_idx) and
%             length(b)=length(pix_idx). That is, pix_idx should
%             NOT be boolean mask. Finally, pix_idx should have
%             1-based indexing: it will be converted to c's 0-based
%             indexing inside the mex-function
% 
%   opts: a struct (l1qc_dct_opts) with the following fields:
%     .epsilon : Denoising parameter. l1qc will enforce ||A*eta - b|| < epsilon
% 
%     .mu :   Log-barrier weight is multiplied by mu each LB iteration
%             Default is  10.
%  
%     .lbtol : Tolerance for the log-barrier iterations and by
%              default, determines the total number of log-barrier 
%              iterations. Default is 1e-3.
% 
%     .l1_tol:  The newton iterations will stop if the relative
%               difference in the l1-norms of $x$ between newton
%               iterations falls below l1_tol. That is, if
%                | ||x_k||_1 - ||x_{k-1}||_1 |
%                ------------------------------ < l1_tol
%                            ||x_k||_1 
%               The default is 0, which means this criteria will
%               never be exercised.
% 
%     .newton_tol : The Newton iterations terminate when the 
%                   newton decrement < newton_tol (try newton_tol = lbtol).
%
%     .newton_max_iter : Maximum number of newton iterations.
%                        Default is 50.
% 
%     .cgtol:   tolerance for the conjugate gradient (CG) solver.
%               Default is 1e-8.
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
%   x: vector of length n, which is the recovered signal, in the DCT domain.
%      To recover the signal in the time/spatial domain,
%      z_recovered = idct(x).
% 
%   LBRes : a struct with stats from the optimization with the following
%           fields:
% 
%         .l1 : (scalar) is the achieved l1 norm of x.
% 
%         .total_newton_iter : is the total number of newton iterations (across
%               all log-barrier iterations).
% 
%         .total_cg_iter : is the total number of conjugate gradiant
%                          iterations, across all log-barrier and all Newton 
%                          iterations.
% 
%         .status: 0 if l1qc finished without errors, 1 otherwise.
% 
