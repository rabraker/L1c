
fname = 'test_data/f_eval_data.json'
rng(1);

N = 50;

x = randn(N,1);
% u = randn(N,1);
u = (0.95)*abs(x) + (0.10)*max(abs(x));
r = randn(N,1)*0.001;

epsilon = 0.01;
tau = 0.1;


fu1 = x - u;
fu2 = -x - u;
fe = 1/2*(r'*r - epsilon^2)
f = sum(u) - (1/tau)*(sum(log(-fu1)) + sum(log(-fu2)) + log(-fe));

data = struct('x', x(:)', 'u', u(:)', 'r', r(:)', 'tau', tau, 'epsilon', epsilon,...
  'fu1_exp', fu1(:)', 'fu2_exp', fu2(:)', 'fe_exp', fe, 'f_exp', f);

opts.FloatFormat = '%.20f';
json_str = savejson('', data, opts)

fid = fopen(fname, 'w');
fprintf(fid, '%s', json_str);
fclose(fid)