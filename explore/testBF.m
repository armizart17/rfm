
clear; clc; close all;
fileName = 'Q:\CAvendanoData\42460445\42460445_IHR_F\42460445_IHR_F.mat';
preSetName = 'Q:\CAvendanoData\42460445\42460445_IHR_F\42460445_IHR_F_preSet.mat';


ax = axes(figure);
generateBMode(fileName, preSetName, ax)