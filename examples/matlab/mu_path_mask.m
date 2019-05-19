% [pix_idx, pix_mask_mat] = muPathMaskGen(mupathLength, n, m, sampling_frac)
% Return Mu path sampling mask
%
% Arguments
% ----------
%   mupath_len : size (in pixels) of mu path pattern
%   n,m : image size
%   sampling_frac : total pixels to sample = samplingRatio*n*m
% Returns
% --------
%   pix_idx : a vector of sampled pixel indeces.
%   pix_mask_mat : an n by m mu path pattern mask, which contains 1's in the
%                  mu-path areas, and zeros elsewhere.


function [pix_idx, pix_mask_mat] = mu_path_mask(mupath_len, n, m, sampling_frac)
   
    pix_mask_mat = zeros(n, m);

    while (sum(sum(pix_mask_mat))<sampling_frac*n*m)
        
        rand_i = randi(n);
        rand_j = randi(m - mupath_len + 1);
        if sum(pix_mask_mat(rand_i, rand_j:(rand_j+mupath_len - 1))) < 0.5
            pix_mask_mat(rand_i, rand_j:(rand_j + mupath_len - 1)) = 1;
        end
        
    end
    
    pix_vec = reshape(pix_mask_mat', n*m,1);
    pix_idx = find(pix_vec > 0.5);
   
end



