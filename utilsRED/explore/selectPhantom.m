function [phantomInfo] = selectPhantom(operator, numPhantom, freqValue)

%%%% MATCH PHANTOM INFORMATION ACCORDING TO SPECS %%% 
% function [phantomInfo] = selectPhantom(operator, numPhantom, freqValue)
% Example: operator = 'AH' | 'ZT', numPhantom = 20 | 22, freqValue = 2.5 | 4.0 | 5.5

%%% PHANTOM 20 %%%
% FREQ = [2.5 4.0 5.5]
% CASE = [130640 130724 130825]
%%% PHANTOM 22 %%%
% FREQ = [2.5 4.0 5.5]
% CASE = [134425 134447 134610]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% OPERATOR ZT %%%

%%% PHANTOM 20 %%%
% FREQ = [2.5 4.0 5.5]
% CASE = [133512 133529 133609]

%%% PHANTOM 22 %%%
% FREQ = [2.5 4.0 5.5]
% CASE = [140408 140427 140503]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FREQ = ["2.5", "4.0", "5.5"]; % [MHz]

lchar_freq = cell2mat(convertStringsToChars(FREQ(:)));
lnum_freq = str2num(lchar_freq);

idx = find(lnum_freq == freqValue);

if strcmp(operator, 'AH')
    if (numPhantom == 20)
        list_case = [130640 130724 130825];
        freqChar = lchar_freq(idx,:);
        numCase = list_case(idx);
    end 
    if (numPhantom == 22)
        list_case = [134425 134447 134610];
        freqChar = lchar_freq(idx,:);
        numCase = list_case(idx);
    end
end

if strcmp(operator, 'ZT')
    if (numPhantom == 20)
        list_case = [133512 133529 133609];
        freqChar = lchar_freq(idx,:);
        numCase = list_case(idx);
    end
    if (numPhantom == 22)
        list_case = [140408 140427 140503];
        freqChar = lchar_freq(idx,:);
        numCase = list_case(idx);
    end 
end

phantomInfo.operator = operator;
phantomInfo.numPhantom = numPhantom;
phantomInfo.freqChar = freqChar;
phantomInfo.numCase = numCase;

end