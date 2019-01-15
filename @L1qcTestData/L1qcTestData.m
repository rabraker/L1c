classdef L1qcTestData
% The static class provides an interface to build the test data needed by the
% c-code unit tests. To build all the data, call, e.g.,
%
% L1qcTestData.build_all_test_data(test_data_root);
  properties
  end
  
  methods (Static)
    function build_all_test_data(test_data_root)
    % Main method to build all of the test data.
      addpath ~/matlab/afm-cs/matlab-code/functions
      
      L1qcTestData.build_newton_init_data(test_data_root);
      L1qcTestData.build_feval_test_data(test_data_root);
      L1qcTestData.build_cgsolve_test_data(test_data_root);
      L1qcTestData.build_dct_test_data(test_data_root);
      L1qcTestData.build_l1qc_newton_test_data(test_data_root);
      L1qcTestData.build_logbarrier_test_data(test_data_root);
    end
    
    xp = build_logbarrier_test_data(data_root, lbiter);
    build_newton_init_data(test_data_root);
     
    build_l1qc_newton_test_data(test_data_root);
    
    build_cgsolve_test_data(test_data_root);
    
    build_dct_test_data(test_data_root);
    
    build_feval_test_data(data_root) 
    
    [x, res, iter] = cgsolve(A, b, tol, maxiter, verbose, x0);
    [xp, up, niter, cgit_tot] = l1qc_newton(x0, u0, A, At, b, epsilon, tau, newtontol,...
      newtonmaxiter, cgtol, cgmaxiter, Tii, verbose, warm_start_cg);
    xp = l1qc_logbarrier(x0, A, At, b, epsilon, lbtol, mu, ...
                         cgtol, cgmaxiter, lbiter, verbose)  

                       
    function [ b ] = Afun_dct(eta, pix_idx)
    % Given the CS equation
    % b = E * M * eta
    % where E is the subsampling matrix and M is the idct.
      b = idct(eta);
      b = b(pix_idx);
    end

    function [ eta ] = Atfun_dct(b, pix_idx, len_eta)
    % Computes the adjoint of the CS equation
    % b = E * M * eta
    % where E is the subsampling matrix and M is the idct. That is
    %
    % eta = M'*E'*b
    %
    % N2 is the length eta.
      
      eta = zeros(len_eta, 1);
      eta(pix_idx) = b;
      eta = dct(eta);
      
    end
                       
    function vpix = pixmat2vec(Mat)
      % Convert the matrix X to a single vector of pixels, going along the rows of X.
      %
      % In other words,
      % v = [X(1,:), X(2,:), ...]'
      
      [n,m] = size(Mat);
      
      vpix = reshape(Mat', n*m,1);
      
    end
    
    function Mpix = pixvec2mat(vpix, nrows)
      % Converts the (forced to be) column vector vpix into a matrix such that
      % The rows of Mpix are taken as contiguous chunks of vpix. I.e.,
      %
      % Mpix = [ vpix(1:ncol)';
      %          vpix(ncol+1:2*ncol)'
      %              :              ]
      % where ncol = length(vpix)/nrows;
      %
      % See Also: PixelVectorToMatrix, which is a less efficient implementation of
      % the same operation.
      
      % assure vpix is a column vector to start
      vpix = vpix(:);
      
      Mpix = reshape(vpix, [], nrows)';
      
    end
  end
  
end