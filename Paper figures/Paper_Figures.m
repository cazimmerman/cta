clearvars; close all; clc
addpath(genpath('figures'))
addpath(genpath('utils'))
warning('off','all')
data_path = uigetdir([],'Select folder with Figshare data');

fig_1(data_path)
fig_2(data_path)
fig_3(data_path)
fig_4(data_path)
fig_5(data_path)

fig_ED1(data_path)
fig_ED2(data_path)
fig_ED3(data_path)
fig_ED4(data_path)
fig_ED5(data_path)
fig_ED6(data_path)
fig_ED8(data_path)
fig_ED9(data_path)
fig_ED10(data_path)