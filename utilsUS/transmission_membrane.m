function t = transmission_membrane(band)
a = -0.087146134643529;
b = 0.642790369805658;
c = 1.052562238394482;
t = (a*band.^(b)+c).^4;
end