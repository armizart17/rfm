function [fLeft,fRight] = findFreqBand(f, y, ratio)
% function [fLeft,fRight] = findFreqBand(f, y, ratio)
% ====================================================================== %
% Function that finds the left and right points in the f axis such that 
% they cut the y axis at ratio*max(y). Used to find the frequency band.
% Inputs:
%       - f         Frequencies vector 
%       - y         Magnitude spectrum vector, must be the same size of f
%       - ratio     Scalar (before do db2pow())
% Outputs:
%       - fLeft     Low Frequency
%       - fRight    High Frequency
% Example
% %% BW from spectrogram
% [pxx,fpxx] = pwelch(sam1-mean(sam1),nz,nz-wz,nz,fs);
% meanSpectrum = mean(pxx,2);
% meanSpectrum(1) = 0;
% % figure, plot(fpxx/1e6,db(meanSpectrum/max(meanSpectrum))),grid on
% if ~fixedBW
%     [freq_L,freq_H] = findFreqBand(fpxx, meanSpectrum, db2mag(-30));
% end
% ====================================================================== %
    y = y / max(y);
    df = f(2)-f(1);
    N = length(y);
    
    [~,imax]=max(y);
    if imax == 1 || imax == N
        disp('Error');
        ix0 = 0; ixf = 0;
        return;
    end
    
    for iLeft = imax:-1:1
        if y(iLeft) < ratio
            fLeft = f(iLeft) + df*(ratio-y(iLeft))/(y(iLeft+1) - y(iLeft));
            break;
        end
        if (iLeft == 1)
            disp('Error');
        end
    end
    
    for iRight = imax:N
        if y(iRight) < ratio
            fRight = f(iRight) - df*(y(iRight) - ratio)/(y(iRight) - y(iRight-1));
            break;
        end
        if (iRight == N)
            disp('Error');
        end
    end
end