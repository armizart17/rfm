

%% CHANGE OF PATH

cd('C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\Data UIUC\');
%% UICI REPOSITORY

%%%%%% CODE TO READ DATA .RFD %%%%%%
addpath('.\Siemens code\')
addpath('.\Siemens code\URI_code\')

%%%%%% DATA .RFD %%%%%%
addpath('.\Siemens_S3000_6C1HD_THI_OFF\RF_data_acquired_by_operator1_AH\P20_THI_Off_F2.5MHz_130640\')
addpath('.\Siemens_S3000_6C1HD_THI_OFF\RF_data_acquired_by_operator1_AH\P20_THI_Off_F2.5MHz_130640\')
    


[RF_data axial lateral fs parameters] = read_siemens(filename, c, whichSlice)    