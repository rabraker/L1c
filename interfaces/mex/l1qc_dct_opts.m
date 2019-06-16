
function opts = l1qc_dct_opts(varargin)
% opts = l1qc_opts(varargin)
% 
% Will build an options struct for l1qc. Vargin are name-value
% pairs with the following options.
% 
% Arguments
%   epsilon : Denoising parameter. l1qc will enforce ||A*eta - b|| < epsilon
% 
%   mu :      log-barrier weight is multiplied by mu each LB iteration
%             Default is  10.
%  
%  lbtol :   Tolerance for the log-barrier iterations and by
%            default, determines the total number of log-barrier 
%            iterations. Default is 1e-3.
% 
%  l1_tol:  The newton iterations will stop if 
%                 ||f(x_k)|_1 - |f(x_k-1)|_1 |/|f(x_k)|_1 <l1_tol
%          Default is 0, which will match the Matlab code.
% 
%  newton_tol : The Newton iterations terminate when the 
%                       newton decrement < newton_tol (try newton_tol = lbtol).
%
%  newton_max_iter : Maximum number of newton iterations.
%                   Default is 50.
% 
%  cgtol:   tolerance for the conjugate gradient (CG) solver.
%           Default is 1e-8.
% 
%  cgmaxiter: maximum iterations for the CG solver. Default is 200.
% 
%  warm_start_cg : (0|1|2) Warm start options for CG solver. This
%                  seems only marginally succesfull and possibly detrimental.
%               -If 0, use dx*0.
%               -If 1, use use previously computed value of dx as warmstart.
%               -If 2, use dx*step from previous (newton)
%               iteration.
% 
%  verbose : How much to print. 
%          - If 0, nothing besides warnings are printed.
%          - If 1, prints info each log-barrier iteration.
%          - If 2, prints same info as 1, but also info each newton iteration.
%

  p = inputParser();
  p.addParameter('epsilon', 0.1);
  p.addParameter('mu', 10);
  p.addParameter('lbtol', 1e-3);
  p.addParameter('tau', 0);
  p.addParameter('lbiter', 0);
  p.addParameter('newton_tol', 1e-3);
  p.addParameter('newton_max_iter', 50);
  p.addParameter('verbose', 2);
  p.addParameter('l1_tol', 1e-5);
  p.addParameter('cgtol', 1e-8);
  p.addParameter('cgmaxiter', 200);
  p.addParameter('warm_start_cg', 0);
  p.addParameter('dct_mode', 1);
  p.parse(varargin{:});
  opts = p.Results();
  
end
