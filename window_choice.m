function [WINDOW] = window_choice(NB_SAMPLE, CHOICE)
% input:    - lenght of the window 
%           - choice, a number which defines the nature of the window
%           function.
%           1. Hanning
%           2. Tuckey 0.25
%           3. Hamming
%           4. Tchebychev
% output:   - the window

switch CHOICE
    case 1
        WINDOW = hann(NB_SAMPLE);
    case 2
        WINDOW = tukeywin(NB_SAMPLE,0.25);
    case 3
        WINDOW = hamming(NB_SAMPLE);
    case 4
        WINDOW = chebwin(NB_SAMPLE);
end
