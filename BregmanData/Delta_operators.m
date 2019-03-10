% Assume the image matrix is concatenated in row major order.

rng(1);
n = 4;
m = 3;

A = round(rand(n,m)*256, 0);
a = A';
a = a(:);

% Thus, the Dx kernel looks like (kernel for signal of length m)
% Dm_kernel = [-1, 1, 0, 0;
%              0,  -1, 1, 0;
%              0,   0, -1 1;
%              0,   0  0  0];


Dm_kernel = -eye(m) + [zeros(m, 1), eye(m, m-1)];
Dm_kernel(end,end) = 0;
Dx_kernel = Dm_kernel;

Dn_kernel = -eye(n) + [zeros(n, 1), eye(n, n-1)];
Dn_kernel(end,end) = 0;
Dy_kernel = Dn_kernel;

% And the full Dx matrix is
In = eye(n);
Im = eye(m);
Dx = kron(In, Dx_kernel);

%And the full Dy matrix is
Dy = kron(Dy_kernel, Im);

% Check:
Dx_a = reshape(Dx*a, m, [])';
Dy_a = reshape(Dy*a, m, [])'

Dx_exp = [diff(A, 1,2), zeros(n,1)];
Dy_exp = [diff(A); zeros(1,m)]

any(any(Dx_exp ~= Dx_a))
any(any(Dy_exp ~= Dy_a))

%%
lambda = 2.0;

lap_y_exp = lambda*(Dy'*Dy)*a
lap_y = lambda*DyTDy(a, n, m)


% Now, we figure out DyT
dyt_exp = lambda*Dy'*a
dyt = lambda*DyT(a, n,m)

function dyt = DyT(a, n, m)
  dyt = a*0;
  for i=1:n*m
    if i<=(m-1)*n
      Ai = a(i);
    else
      Ai = 0;
    end
    if i<= m
      Ai_min_m = 0;
    else
      Ai_min_m = a(i-m);
    end
    dyt(i) = Ai_min_m - Ai;
  end
  
end


%%
function lap_y = DyTDy(a, n, m)
  lap_y = a*0;
  % First, do the diagonal and upper diagonal
  for i = 1:(m-1)*n
    if i <= m
      D_ii = 1;
    else
      D_ii = 2;
    end
    lap_y(i) = a(i)*D_ii -a(i+m);
  end
  % Now, subtract the lower diagonal:
  for i=m+1:m*n
    if i <= (m-1)*n
      D_ii = 0; % Already accounted for.
    else
      D_ii = 1;
    end
    lap_y(i) = lap_y(i) - a(i-m) + D_ii*a(i);
  end
end