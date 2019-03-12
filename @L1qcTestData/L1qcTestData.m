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
      
      % Make us a test image to work with.
      L1qcTestData.make_test_image(test_data_root);
      
      L1qcTestData.build_l1qc_newton_test_data(test_data_root);

      L1qcTestData.build_logbarrier_test_data(test_data_root);

    end
    function make_test_image(test_data_root)
    % A build a 256 x 256 test image of the CS20ng grating.
      x_start = 13;
      y_start = 13;
      npix = 256;
      nholes = 10;
      img_mat = L1qcTestData.make_CS20NG(x_start, y_start, npix, nholes);
      pix_idx = L1qcTestData.mu_path_mask(15, npix, npix, 0.15, false);
      
      xorig = L1qcTestData.pixmat2vec(img_mat);
      save(fullfile(test_data_root, 'test_image_data.mat'), 'xorig', 'pix_idx');
    end
    
    xp = build_logbarrier_test_data(data_root, lbiter);
         
    build_l1qc_newton_test_data(test_data_root);
    
            
    [x, res, iter] = cgsolve(A, b, tol, maxiter, verbose, x0);
    [xp, up, niter, cgit_tot] = l1qc_newton(x0, u0, A, At, b, epsilon, tau, newtontol,...
      newtonmaxiter, cgtol, cgmaxiter, Tii, verbose, warm_start_cg);
    xp = l1qc_logbarrier(x0, A, At, b, epsilon, lbtol, mu, ...
                         cgtol, cgmaxiter, lbiter, verbose);
    %
    [pix_idx, pix_mask_mat] = mu_path_mask(mupath_len,n,m,samplingRatio,repeat_sampling);
    

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
    
    function img_mat = make_CS20NG(x_start, y_start, npix, nholes)
  
      % Make these options later
      if nargin <3
        npix = 512;
      end
      
      if nargin <4
        hole_width = npix/20;
        pitch = npix/10;
      else
        pitch = npix/nholes;
        hole_width = npix/2/nholes;
      end
      
      img_mat = ones(npix, npix);
      m_pitch = ceil(pitch);
      n_hole = ceil(hole_width);
      
      % First lets, create a prototype hole. Then we will insert it into the master
      % image periodically.
      sq = ones(n_hole,n_hole);
      % circle radius
      r = n_hole/2;
      % center
      xc = n_hole/2;
      yc = n_hole/2;
      
      
      for x_pt=1:n_hole
        for y_pt=1:n_hole
          
          rad_pt = sqrt( (x_pt - xc)^2 + (y_pt - yc)^2 );
          
          if rad_pt <= r
            sq(y_pt, x_pt) = 0; % make it black
          end
        end
      end
      
      x_lft = x_start;
      while x_lft+n_hole < npix +n_hole
        y_tp = y_start;
        while y_tp+n_hole <= npix + n_hole
          img_mat(y_tp:y_tp+n_hole-1, x_lft:x_lft+n_hole-1) = sq;
          y_tp = y_tp + m_pitch;
        end
        x_lft = x_lft + m_pitch;
        
      end
      img_mat = img_mat(1:npix, 1:npix);
      
    end
    function print_bt(ME)
      fprintf('%s\n', ME.message);
      for i=1:length(ME.stack)
        fprintf('Error in %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
      end
    end
  end
  
end