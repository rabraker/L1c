% Return Mu path sampling mask

% ------------- INPUT -------------------
% mupathLength - size of mu path pattern
% n,m - image size
% samplingRatio(Percentage) - total pixels to sample = samplingRatio*n*m
% RepeatSamplingFlag - true: Repeat Sampling ALLOWED; false: Repeat Sampling NOT ALLOWED
% ------------- OUTPUT -------------------
% pixelifsampled - mu path pattern mask 


function pixelifsampled = muPathMaskGen(mupathLength,n,m,samplingRatio,RepeatSamplingFlag)
   
    if nargin < 5
        RepeatSamplingFlag = false;
    end

    pixelifsampled = zeros(n,m);

    while (sum(sum(pixelifsampled))<samplingRatio*n*m)
        
        if RepeatSamplingFlag
            
            rand_i = randi(n);
            rand_j = randi([2-mupathLength m]);            
            pixelifsampled(rand_i,max(rand_j,1):min(rand_j+mupathLength-1,n)) = 1;            
            
        else            
            
            rand_i = randi(n);
            rand_j = randi(m-mupathLength+1);
            if sum(pixelifsampled(rand_i,rand_j:rand_j+mupathLength-1)) < 0.5
                pixelifsampled(rand_i,rand_j:rand_j+mupathLength-1) = 1;
            end            
            
        end

    end
    
    
    %imshow(pixelifsampled, [0 1]);
   
   
end



