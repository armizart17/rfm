function window = window_choice(len_window,choice)
% -------------------------------------------------------------------------
% DESCRIPTION 
% Selects a window
% -------------------------------------------------------------------------
% INPUTS
%       len_window      Length of the window
%       choice          Number which defines how the window is generated:
%                       1: Hanning
%                       2: Tuckey with cosine fraction 0.25
%                       3: Hamming
%                       4: Tchebychev
% -------------------------------------------------------------------------
% OUTPUTS
%       window          Generated window
% -------------------------------------------------------------------------
% AUTHOR: Unknown (adapted by Hector Chahuara)
% CONTACT: hector.chahuara@pucp.edu.pe
% -------------------------------------------------------------------------

    switch choice
        case 1
            window = hann(len_window);
        case 2
            window = tukeywin(len_window,0.25);
        case 3
            window = hamming(len_window);
        case 4
            window = chebwin(len_window);
    end

end
