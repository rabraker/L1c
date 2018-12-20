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
    
    [x, res, iter] = cgsolve(A, b, tol, maxiter, verbose);
    [xp, up, niter] = l1qc_newton(x0, u0, A, At, b, epsilon, tau, newtontol,...
      newtonmaxiter, cgtol, cgmaxiter, Tii);

    function [ output ] = DCTfun( z,E )
      t = zeros(length(E),1);
      t(E>0.5)=z;
      output = dct(t);
    end
    
    function [ output ] = IDCTfun( z,E )
      output = idct(z);
      output = output(E>0.5);
    end
    
  end

 
  
end