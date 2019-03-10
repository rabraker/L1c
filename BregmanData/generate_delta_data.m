
n = 4;
m = 4;

A = round(rand(n,m)*256, 0);


% [diff(A); zeros(1,4)]
Dx_exp = [diff(A, 1,2), zeros(n,1)];
Dy_exp = [diff(A);, zeros(1,m)]

dx_exp = Dx_exp';
dx_exp = dx_exp(:);

dy_exp = Dy_exp';
dy_exp = dy_exp(:);


a = A';
a = a(:);




dx = Dx(a, n, m);
reshape(dx', m, n)' - Dx_exp

dy = Dy(a, n, m);
reshape(dy', m, n)' - Dy_exp

x = round(randn(1, 32), 2);
gam = 0.25;
x_sh = shrink1(x, gam)

y = round(randn(1, 32), 2);
z = -x + y;

nrm_err = norm(x-y,2);

jopts.FileName = 'test_data/bregman_diff.json';
savejson('', struct('A', a(:)', 'dy', dy_exp(:)', 'dx', dx_exp(:)', 'n', n, 'm', m), jopts);

jopts.FileName = 'test_data/bregman.json';
savejson('', struct('x', x(:)', 'x_shrunk', x_sh(:)', 'gamma', gam,...
  'y', y(:)', 'z', z(:)', 'nrm_err', nrm_err), jopts);

function y = shrink1(x, gamma)
  
  y = sign(x).* max(abs(x) - gamma, 0);
  
end
function dy = Dy(X, n, m)
  
  dy = zeros(1,n*m);
  
  for row=0:n-2
    for col= row*m:(row+1)*m - 1
      i = col +1; %row*m +col + 1
      
      dy(i) = X(i+m) - X( i);
    end
  end
  
end


function dx = Dx(X, n, m)
  
  dx = zeros(1,n*m);
  
  for row=0:n-1
    for col= row*m:(row+1)*m - 2
      [row+1, col+1];
      i = col +1; %row*m +col + 1
      
      dx(i) = X(i+1) - X( i);
    end
  end
  
end
