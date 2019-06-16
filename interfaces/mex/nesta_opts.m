% opts = nesta_opts(varargin)
% Builds an options structure for nesta_dctTV().
% 
% Supply name-value pairs. The possible names are:
% 
%  alpha_v - weight for vertical TV (default 0)
%  alpha_h - weight for horizontal TV (default 0)
% 
%  n_continue - number of continuation steps. The default is 5
% 
%  tol - tolerance for the stopping criteria.
% 
%  verbose - if 0 no output is displayn.
%
%  mode - (string) ['analysis'|'synthesys']
% 
function opts = nesta_opts(varargin)
    
   p = inputParser();
   p.addParameter('sigma', 0.01);
   p.addParameter('mu',1e-5);
   p.addParameter('tol',1e-3);
   p.addParameter('alpha_v', 0);
   p.addParameter('alpha_h', 0);

   p.addParameter('verbose', 5);
   p.addParameter('n_continue',5);
   p.addParameter('bp_mode', 'analysis');
   p.addParameter('dct_mode', 1);
   
   
   p.parse(varargin{:});
   
   opts = p.Results;
   
end
