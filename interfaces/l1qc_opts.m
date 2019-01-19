
function opts = l1qc_opts(varargin)
% Will build an options struct for l1qc.
  p = inputParser();
  p.addParameter('epsilon', 0.1);
  p.addParameter('mu', 10);
  p.addParameter('cgtol', 1e-8);
  p.addParameter('cgmaxiter', 200);
  p.addParameter('warm_start_cg', 0);
  p.addParameter('lbtol', 1e-3);
  p.addParameter('newton_tol', 1e-3);
  p.addParameter('newton_max_iter', 50);
  p.addParameter('verbose', 2);
  p.addParameter('l1_tol', 1e-5);
  
  p.parse(varargin{:});
  opts = p.Results();
  
end
