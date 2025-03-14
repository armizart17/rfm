function iBSC = compute_iBSC(b,n,band)
    band = band(:); % assure col vector @
    band = band'; % make it row vetor @
    F = repmat(band,[length(n),1]);
    N = repmat(n,[1,length(band)]);
    B = repmat(b,[1,length(band)]);
    iBSC = sum(B.*(F.^N),2)/(length(band)*(band(end) - band(1)));
end