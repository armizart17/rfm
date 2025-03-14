
function [out]=my_RF2Bmode(RF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: generateBmode of RF data
% INPUTS:
%       RF: RF data
%       window = [m n] size of window
%       W: weight to be applied on center pixel of each window
% OUTPUT:
%       out: Bmode image
% AUTHOR: EDMUNDO AROM MIRANDA ZARATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%     h = waitbar(0, 'Generando las imagenes de brillo');
    iUS = zeros(size(RF));
    for index=1:size(RF,3)
        US = 20*log10(abs(hilbert(RF(:,:,index)))); 
        iUS(:,:,index) = US - max(US(:));
%         waitbar(index/size(RF,3),h)
    end
%     close(h);

    out = iUS;

end
