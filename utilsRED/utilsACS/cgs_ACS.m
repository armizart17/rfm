function [ACS, n_ori] = cgs_ACS(A, SLD_term)
% function [ACS, n_ori] = cgs4acs(A, SLD_term)

    u = cgs(A'*A,A'*SLD_term(:)); % IS A VECTOR
    
    % GET ACS and N components
    ACS = u(1:end/2); 
    n_ori = u(end/2+1:end);
    ACS = 8.686*ACS;   % [dB/cm/Hz]
    
    % RESHAPE AS IMAGES
    [H,W,~] =size(SLD_term);
    ACS = reshape(ACS,H,W);
    n_ori = reshape(n_ori,H,W);

end