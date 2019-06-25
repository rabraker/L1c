classdef TestL1cMex < matlab.unittest.TestCase
  
  properties
    b_small;
    pix_idx_small;
    x_orig_small;
    x0_small;
    
    b_large;
    pix_idx_large;
    x_orig_large;
    x0_large;
  
    opts;
  end
  
  methods
    function self = TestL1cMex()
      rng(1);
      Nsmall = 1024;
      alpha = 0.1;
      t = [0:Nsmall-1]';
      self.x_orig_small = cos(2*pi*t/Nsmall) + 0.00*randn(Nsmall,1) + sin(2*pi*t*5/Nsmall);
      
      self.x_orig_small(300:400) = 1;
      
      mask = rand(Nsmall, 1);
      
      idx_keep = find(mask < alpha);
      idx_drop = find(mask >= alpha);
      mask(idx_keep) = 1;
      mask(idx_drop) = 0;
      
      self.pix_idx_small = find(mask == 1);
      self.b_small = self.x_orig_small(self.pix_idx_small);
      self.x0_small = self.DCTfun_IDX(self.b_small, self.pix_idx_small, Nsmall);
      self.opts = struct('epsilon', 0.01, 'mu', 10, 'cgtol', 1e-8, 'cgmaxiter', 200,...
        'lbtol', 1e-3, 'newton_tol', 1e-4, 'newton_max_iter', 50, 'verbose', 0,...
        'warm_start_cg', 0);

    end
    function [ output ] = IDCTfun_IDX(self, z, pix_mask_IDX)
      % A(x) = E*M*x
      output = idct(z);
      output = output(pix_mask_IDX);
      
    end
    
    function [ output ] = DCTfun_IDX(self, v, pix_mask_idx, N)
      % At(x) = M^T*E^T * x
      t = zeros(N,1);
      t(pix_mask_idx)=v;
      output = dct(t);
    end     

end
  
  
  methods (Test)
    
    function test_matrix_x_invalid(self)
      x = [1,2,3,4;
        1,2,3,4];
      b = [1,2,3,4];
      pix_idx = [1,2,3,4];
      exp_error = 'l1qc:l1qc_log_barrier:notVector';
      self.verifyError(@()l1qc_helper(x,b, pix_idx, self.opts), exp_error)
      
      x = [1,2,3,4, 5,6];
      b = [1,2;3,4];
      self.verifyError(@()l1qc_helper(x,b, pix_idx, self.opts), exp_error)
      
      x = [1,2,3,4, 5,6];
      b = [1,2,3,4];
      pix_idx = [1,2;3,4];
      self.verifyError(@()l1qc_helper(x,b, pix_idx, self.opts), exp_error)
    end
    
    function test_struct_required(self)
      exp_error = 'l1qc:l1qc_log_barrier:notStruct';
      
      x = [1,2,3,4, 5,6];
      b = [1,2,3,4];
      pix_idx = [1,2,3,4];
      self.verifyError(@()l1qc_helper(x,b, pix_idx, 5), exp_error)
      
    end
    
    
    function test_sizes(self)
      
      x = [1,2,3,4];
      b = [1,2,3,4];
      pix_idx = [1,2,3];
      
      exp_error = 'l1qc:l1qc_log_barrier:incompatible_dimensions';
      self.verifyError(@()l1qc_helper(x,b,pix_idx, self.opts), exp_error)
      
    end
    function test_num_in_out(self)
      x = [1,2,3,4];
      b = [1,2,3,4];
      pix_idx = [1,2,3];
      exp_error = 'l1qc:l1qc_log_barrier:nrhs';
      self.verifyError(@()l1qc(x,b, self.opts), exp_error)
      
      exp_error = 'l1qc:l1qc_log_barrier:nlhs';
      self.verifyError(@()l1qc(x,b, pix_idx, self.opts), exp_error)

    end
    function test_inputs_copied(self)
    % Check that the mex function copies the inputs so that our workspace data
    % is not perturbed.

      b_tmp = self.b_small;
      pix_idx_tmp = self.pix_idx_small;
      
      % It seems we have to force matlab to make a copy, ie, it seems that just
      % declaring x0_tmp = self.x0_small leaves the two linked, until we change
      % one here, so our test passes when it shouldn't
      x0_tmp = self.x0_small * 2;
      x0_tmp = x0_tmp/2;
      self.verifyEqual(self.x0_small, x0_tmp, 'AbsTol', 1e-15);
      % -------------------------------------------------
      opts = self.opts;
%       A = @(x)self.IDCTfun_IDX(x, self.pix_idx_small);
%       At = @(x)self.DCTfun_IDX(x, self.pix_idx_small, length(self.x_orig_small));
%       xopt_ml = L1qcTestData.l1qc_logbarrier(x0_tmp, A, At, b_tmp, opts.epsilon, opts.lbtol,...
%         opts.mu, opts.cgtol, opts.cgmaxiter, opts.verbose);
%     
%       figure, plot([xopt_ml, x0_tmp, dct(self.x_orig_small)])
%       figure, plot([idct(x0_tmp), idct(xopt_ml), self.x_orig_small])
% 
      
      xopt = l1qc(x0_tmp, b_tmp, pix_idx_tmp-1, self.opts);
      
      self.verifyEqual(self.x0_small, x0_tmp, 'AbsTol', 1e-15);
      self.verifyEqual(self.x0_small, x0_tmp);
      self.verifyEqual(self.b_small, b_tmp);
      self.verifyEqual(self.pix_idx_small, pix_idx_tmp);
      
    
    end
    
  end
  
  
end

function l1qc_helper(x,b, pix_idx, opts)
  % We need an output argument.
  x = l1qc(x,b,pix_idx, opts)
  
end