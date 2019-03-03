%
%
function [x_out, lb_res, stat] = l1qc_so_wrap(N, b, pix_idx, opts)

%   opts2.epsilon = opts.epsilon;
%   opts2.mu = opts.mu;
%   opts2.lbtol = opts.lbtol;
%   opts2.tau = opts.tau;
%   opts2.lbiter = opts.lbiter;
%   opts2.newton_tol = opts.newton_tol;
%   opts2.newton_max_iter = opts.newton_max_iter;
%   opts2.verbose = opts.verbose;
%   opts2.l1_tol = opts.l1_tol;
%   opts2.cgtol = opts.cgtol;
%   opts2.cgmaxiter = opts.cgmaxiter;
%   opts2.warm_start_cg = opts.warm_start_cg;

  libname = 'libl1qc_dct';
  hfile = 'libl1qc_dct.h';

  if not(libisloaded(libname))
    loadlibrary(libname, hfile)
  end
  %libfunctionsview(libname)

  lb_res = struct('l1', 0, 'total_newton_iter', 0,...
    'total_cg_iter', 0, 'status', 0);
  x_out = zeros(1, N);
  M = length(b);
  [stat, x_out,~, ~, lb_res] = calllib(libname, 'l1qc_dct',...
                                       N, 1, x_out, M, b, int32(pix_idx)-1, opts, lb_res);


end

% if ~libisloaded('shrlibsample')
%    addpath(fullfile(matlabroot,'extern','examples','shrlib'))
%    loadlibrary('shrlibsample')
% end


%
