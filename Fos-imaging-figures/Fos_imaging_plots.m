%% Setup
addpath(genpath('Z:\Chris\matlab\violin-plot\'));
addpath(genpath('Z:\Chris\matlab\heatmap\'));
clear all; close all; clc
path.main = 'Z:\Chris\data\clearmap2\';
path.file = '\rawdata\resolution_3.6x\region_summary_statistics_classified_350size_0intensity.csv';
library = 'no_cerebellum_OLD';
GLMMdir = ['GLMM_',path.file(63:end-4),'/'];
load(['regions_allen_',library,'.mat'],'RegionLibrary')
%% Load data
fname  = cell(0,0);
flavor = cell(0,0);
phase  = cell(0,0);
batch  = cell(0,0);
sex    = cell(0,0);
code   = cell(0,0);

%%%%% BATCH 0 %%%%%
% Familiar/Consumption/Male
fname{end+1} = 'zimmerman_05\zimmerman_05-216\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch00'; sex{end+1} = 'Male'; code{end+1} = '00-000';
fname{end+1} = 'zimmerman_05\zimmerman_05-217\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch00'; sex{end+1} = 'Male'; code{end+1} = '00-001';
fname{end+1} = 'zimmerman_05\zimmerman_05-218\imaging_request_2'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch00'; sex{end+1} = 'Male'; code{end+1} = '00-002';
fname{end+1} = 'zimmerman_05\zimmerman_05-220\imaging_request_2'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch00'; sex{end+1} = 'Male'; code{end+1} = '00-003';
fname{end+1} = 'zimmerman_05\zimmerman_05-221\imaging_request_2'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch00'; sex{end+1} = 'Male'; code{end+1} = '00-004';
fname{end+1} = 'zimmerman_05\zimmerman_05-222\imaging_request_2'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch00'; sex{end+1} = 'Male'; code{end+1} = '00-005';
% Novel/Consumption/Male
fname{end+1} = 'zimmerman_05\zimmerman_05-211\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch00'; sex{end+1} = 'Male'; code{end+1} = '00-006';
fname{end+1} = 'zimmerman_05\zimmerman_05-212\imaging_request_2'; flavor{end+1} = 'Novel'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch00'; sex{end+1} = 'Male'; code{end+1} = '00-007';
fname{end+1} = 'zimmerman_05\zimmerman_05-213\imaging_request_2'; flavor{end+1} = 'Novel'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch00'; sex{end+1} = 'Male'; code{end+1} = '00-008';
fname{end+1} = 'zimmerman_05\zimmerman_05-223\imaging_request_2'; flavor{end+1} = 'Novel'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch00'; sex{end+1} = 'Male'; code{end+1} = '00-009';
fname{end+1} = 'zimmerman_05\zimmerman_05-224\imaging_request_3'; flavor{end+1} = 'Novel'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch00'; sex{end+1} = 'Male'; code{end+1} = '00-010';
fname{end+1} = 'zimmerman_05\zimmerman_05-225\imaging_request_2'; flavor{end+1} = 'Novel'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch00'; sex{end+1} = 'Male'; code{end+1} = '00-011';

%%%%% BATCH 1 %%%%%
% Familiar/Consumption/Female
fname{end+1} = 'zimmerman_08\zimmerman_08-976\imaging_request_2'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch01'; sex{end+1} = 'Female'; code{end+1} = '01-012';
fname{end+1} = 'zimmerman_08\zimmerman_08-977\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch01'; sex{end+1} = 'Female'; code{end+1} = '01-013';
fname{end+1} = 'zimmerman_08\zimmerman_08-978\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch01'; sex{end+1} = 'Female'; code{end+1} = '01-014';
fname{end+1} = 'zimmerman_08\zimmerman_08-979\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch01'; sex{end+1} = 'Female'; code{end+1} = '01-015';
fname{end+1} = 'zimmerman_08\zimmerman_08-980\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch01'; sex{end+1} = 'Female'; code{end+1} = '01-016';
fname{end+1} = 'zimmerman_08\zimmerman_08-984\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch01'; sex{end+1} = 'Female'; code{end+1} = '01-017';
% Novel/Consumption/Female
fname{end+1} = 'zimmerman_08\zimmerman_08-981\imaging_request_3'; flavor{end+1} = 'Novel'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch01'; sex{end+1} = 'Female'; code{end+1} = '01-018';
fname{end+1} = 'zimmerman_08\zimmerman_08-982\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch01'; sex{end+1} = 'Female'; code{end+1} = '01-019';
fname{end+1} = 'zimmerman_08\zimmerman_08-983\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch01'; sex{end+1} = 'Female'; code{end+1} = '01-020';
fname{end+1} = 'zimmerman_08\zimmerman_08-985\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch01'; sex{end+1} = 'Female'; code{end+1} = '01-021';
fname{end+1} = 'zimmerman_08\zimmerman_08-986\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch01'; sex{end+1} = 'Female'; code{end+1} = '01-022';
fname{end+1} = 'zimmerman_08\zimmerman_08-987\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Consumption'; batch{end+1} = 'Batch01'; sex{end+1} = 'Female'; code{end+1} = '01-023';

%%%%% BATCH 2 %%%%%
% Familiar/Malaise/Male
fname{end+1} = 'zimmerman_09\zimmerman_09-930\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch02'; sex{end+1} = 'Male'; code{end+1} = '02-024';
fname{end+1} = 'zimmerman_09\zimmerman_09-931\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch02'; sex{end+1} = 'Male'; code{end+1} = '02-025';
fname{end+1} = 'zimmerman_09\zimmerman_09-932\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch02'; sex{end+1} = 'Male'; code{end+1} = '02-026';
fname{end+1} = 'zimmerman_09\zimmerman_09-933\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch02'; sex{end+1} = 'Male'; code{end+1} = '02-027';
fname{end+1} = 'zimmerman_09\zimmerman_09-934\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch02'; sex{end+1} = 'Male'; code{end+1} = '02-028';
fname{end+1} = 'zimmerman_09\zimmerman_09-935\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch02'; sex{end+1} = 'Male'; code{end+1} = '02-029';
% Novel/Malaise/Male
fname{end+1} = 'zimmerman_09\zimmerman_09-926\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch02'; sex{end+1} = 'Male'; code{end+1} = '02-030';
fname{end+1} = 'zimmerman_09\zimmerman_09-927\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch02'; sex{end+1} = 'Male'; code{end+1} = '02-031';
fname{end+1} = 'zimmerman_09\zimmerman_09-928\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch02'; sex{end+1} = 'Male'; code{end+1} = '02-032';
fname{end+1} = 'zimmerman_09\zimmerman_09-929\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch02'; sex{end+1} = 'Male'; code{end+1} = '02-033';
fname{end+1} = 'zimmerman_09\zimmerman_09-936\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch02'; sex{end+1} = 'Male'; code{end+1} = '02-034';
fname{end+1} = 'zimmerman_09\zimmerman_09-937\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch02'; sex{end+1} = 'Male'; code{end+1} = '02-035';

%%%%% BATCH 3 %%%%%
% Familiar/Malaise/Female
fname{end+1} = 'zimmerman_11\zimmerman_11-991\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch03'; sex{end+1} = 'Female'; code{end+1} = '03-036';
fname{end+1} = 'zimmerman_11\zimmerman_11-992\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch03'; sex{end+1} = 'Female'; code{end+1} = '03-037';
fname{end+1} = 'zimmerman_11\zimmerman_11-993\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch03'; sex{end+1} = 'Female'; code{end+1} = '03-038';
fname{end+1} = 'zimmerman_11\zimmerman_11-994\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch03'; sex{end+1} = 'Female'; code{end+1} = '03-039';
fname{end+1} = 'zimmerman_11\zimmerman_11-995\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch03'; sex{end+1} = 'Female'; code{end+1} = '03-040';
fname{end+1} = 'zimmerman_11\zimmerman_11-997\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch03'; sex{end+1} = 'Female'; code{end+1} = '03-041';
% Novel/Malaise/Female
fname{end+1} = 'zimmerman_11\zimmerman_11-949\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch03'; sex{end+1} = 'Female'; code{end+1} = '03-042';
fname{end+1} = 'zimmerman_11\zimmerman_11-988\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch03'; sex{end+1} = 'Female'; code{end+1} = '03-043';
fname{end+1} = 'zimmerman_11\zimmerman_11-989\imaging_request_2'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch03'; sex{end+1} = 'Female'; code{end+1} = '03-044';
fname{end+1} = 'zimmerman_11\zimmerman_11-990\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch03'; sex{end+1} = 'Female'; code{end+1} = '03-045';
fname{end+1} = 'zimmerman_11\zimmerman_11-998\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch03'; sex{end+1} = 'Female'; code{end+1} = '03-046';
fname{end+1} = 'zimmerman_11\zimmerman_11-999\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise'; batch{end+1} = 'Batch03'; sex{end+1} = 'Female'; code{end+1} = '03-047';

% %%%%% BATCH 7 %%%%%
% % Familiar/Retrieval/Male
fname{end+1} = 'zimmerman_16\zimmerman_16-065\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch07'; sex{end+1} = 'Male'; code{end+1} = '04-048';
fname{end+1} = 'zimmerman_16\zimmerman_16-066\imaging_request_2'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch07'; sex{end+1} = 'Male'; code{end+1} = '04-049';
fname{end+1} = 'zimmerman_16\zimmerman_16-068\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch07'; sex{end+1} = 'Male'; code{end+1} = '04-050';
fname{end+1} = 'zimmerman_16\zimmerman_16-069\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch07'; sex{end+1} = 'Male'; code{end+1} = '04-051';
fname{end+1} = 'zimmerman_16\zimmerman_16-070\imaging_request_2'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch07'; sex{end+1} = 'Male'; code{end+1} = '04-052';
fname{end+1} = 'zimmerman_16\zimmerman_16-071\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch07'; sex{end+1} = 'Male'; code{end+1} = '04-053';
% % Novel/Retrieval/Male
fname{end+1} = 'zimmerman_16\zimmerman_16-064\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch07'; sex{end+1} = 'Male'; code{end+1} = '04-054';
fname{end+1} = 'zimmerman_16\zimmerman_16-067\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch07'; sex{end+1} = 'Male'; code{end+1} = '04-055';
fname{end+1} = 'zimmerman_16\zimmerman_16-072\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch07'; sex{end+1} = 'Male'; code{end+1} = '04-056';
fname{end+1} = 'zimmerman_16\zimmerman_16-073\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch07'; sex{end+1} = 'Male'; code{end+1} = '04-057';
fname{end+1} = 'zimmerman_16\zimmerman_16-074\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch07'; sex{end+1} = 'Male'; code{end+1} = '04-058';
fname{end+1} = 'zimmerman_16\zimmerman_16-075\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch07'; sex{end+1} = 'Male'; code{end+1} = '04-059';

%%%%% BATCH 6 %%%%%
% Familiar/Retrieval/Female
fname{end+1} = 'zimmerman_15\zimmerman_15-564\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch06'; sex{end+1} = 'Female'; code{end+1} = '05-060';
fname{end+1} = 'zimmerman_15\zimmerman_15-565\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch06'; sex{end+1} = 'Female'; code{end+1} = '05-061';
fname{end+1} = 'zimmerman_15\zimmerman_15-566\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch06'; sex{end+1} = 'Female'; code{end+1} = '05-062';
fname{end+1} = 'zimmerman_15\zimmerman_15-570\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch06'; sex{end+1} = 'Female'; code{end+1} = '05-063';
fname{end+1} = 'zimmerman_15\zimmerman_15-571\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch06'; sex{end+1} = 'Female'; code{end+1} = '05-064';
fname{end+1} = 'zimmerman_15\zimmerman_15-572\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch06'; sex{end+1} = 'Female'; code{end+1} = '05-065';
% Novel/Retrieval/Female
fname{end+1} = 'zimmerman_15\zimmerman_15-567\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch06'; sex{end+1} = 'Female'; code{end+1} = '05-066';
fname{end+1} = 'zimmerman_15\zimmerman_15-568\imaging_request_2'; flavor{end+1} = 'Novel'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch06'; sex{end+1} = 'Female'; code{end+1} = '05-067';
fname{end+1} = 'zimmerman_15\zimmerman_15-569\imaging_request_2'; flavor{end+1} = 'Novel'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch06'; sex{end+1} = 'Female'; code{end+1} = '05-068';
fname{end+1} = 'zimmerman_15\zimmerman_15-573\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch06'; sex{end+1} = 'Female'; code{end+1} = '05-069';
fname{end+1} = 'zimmerman_15\zimmerman_15-574\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch06'; sex{end+1} = 'Female'; code{end+1} = '05-070';
fname{end+1} = 'zimmerman_15\zimmerman_15-575\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Retrieval'; batch{end+1} = 'Batch06'; sex{end+1} = 'Female'; code{end+1} = '05-071';

%%%% BATCH 4 %%%%%
% Familiar/CGRP/Male
fname{end+1} = 'zimmerman_12\zimmerman_12-233\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch04'; sex{end+1} = 'Male'; code{end+1} = '06-072';
fname{end+1} = 'zimmerman_12\zimmerman_12-234\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch04'; sex{end+1} = 'Male'; code{end+1} = '06-073';
fname{end+1} = 'zimmerman_12\zimmerman_12-235\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch04'; sex{end+1} = 'Male'; code{end+1} = '06-074';
% Familiar/CGRP/Female
fname{end+1} = 'zimmerman_12\zimmerman_12-246\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch04'; sex{end+1} = 'Female'; code{end+1} = '06-075';
fname{end+1} = 'zimmerman_12\zimmerman_12-247\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch04'; sex{end+1} = 'Female'; code{end+1} = '06-076';
fname{end+1} = 'zimmerman_12\zimmerman_12-248\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch04'; sex{end+1} = 'Female'; code{end+1} = '06-077';
% Novel/CGRP/Male
fname{end+1} = 'zimmerman_12\zimmerman_12-242\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch04'; sex{end+1} = 'Male'; code{end+1} = '06-078';
fname{end+1} = 'zimmerman_12\zimmerman_12-243\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch04'; sex{end+1} = 'Male'; code{end+1} = '06-079';
fname{end+1} = 'zimmerman_12\zimmerman_12-244\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch04'; sex{end+1} = 'Male'; code{end+1} = '06-080';
fname{end+1} = 'zimmerman_12\zimmerman_12-245\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch04'; sex{end+1} = 'Male'; code{end+1} = '06-081';
% Novel/CGRP/Female
fname{end+1} = 'zimmerman_12\zimmerman_12-240\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch04'; sex{end+1} = 'Female'; code{end+1} = '06-082';
fname{end+1} = 'zimmerman_12\zimmerman_12-241\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch04'; sex{end+1} = 'Female'; code{end+1} = '06-083';

%%%% BATCH 5 %%%%%
% Familiar/CGRP/Male
fname{end+1} = 'zimmerman_14\zimmerman_14-249\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch05'; sex{end+1} = 'Male'; code{end+1} = '07-084';
fname{end+1} = 'zimmerman_14\zimmerman_14-250\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch05'; sex{end+1} = 'Male'; code{end+1} = '07-085';
fname{end+1} = 'zimmerman_14\zimmerman_14-252\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch05'; sex{end+1} = 'Male'; code{end+1} = '07-086';
% Familiar/CGRP/Female
fname{end+1} = 'zimmerman_14\zimmerman_14-239\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch05'; sex{end+1} = 'Female'; code{end+1} = '07-087';
fname{end+1} = 'zimmerman_14\zimmerman_14-251\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch05'; sex{end+1} = 'Female'; code{end+1} = '07-088';
fname{end+1} = 'zimmerman_14\zimmerman_14-253\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch05'; sex{end+1} = 'Female'; code{end+1} = '07-089';
fname{end+1} = 'zimmerman_14\zimmerman_14-254\imaging_request_1'; flavor{end+1} = 'Familiar'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch05'; sex{end+1} = 'Female'; code{end+1} = '07-090';
% Novel/CGRP/Male
fname{end+1} = 'zimmerman_14\zimmerman_14-255\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch05'; sex{end+1} = 'Male'; code{end+1} = '07-091';
fname{end+1} = 'zimmerman_14\zimmerman_14-256\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch05'; sex{end+1} = 'Male'; code{end+1} = '07-092';
fname{end+1} = 'zimmerman_14\zimmerman_14-257\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch05'; sex{end+1} = 'Male'; code{end+1} = '07-093';
% Novel/CGRP/Female
fname{end+1} = 'zimmerman_14\zimmerman_14-258\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch05'; sex{end+1} = 'Female'; code{end+1} = '07-094';
fname{end+1} = 'zimmerman_14\zimmerman_14-259\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch05'; sex{end+1} = 'Female'; code{end+1} = '07-095';
fname{end+1} = 'zimmerman_14\zimmerman_14-260\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch05'; sex{end+1} = 'Female'; code{end+1} = '07-096';
fname{end+1} = 'zimmerman_14\zimmerman_14-261\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch05'; sex{end+1} = 'Female'; code{end+1} = '07-097';
fname{end+1} = 'zimmerman_14\zimmerman_14-262\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'CGRP stim'; batch{end+1} = 'Batch05'; sex{end+1} = 'Female'; code{end+1} = '07-098';

%%%%% BATCH 0 %%%%%
fname{end+1} = 'zimmerman_26\zimmerman_26-579\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch08'; sex{end+1} = 'Male';   code{end+1} = '08-099';
fname{end+1} = 'zimmerman_26\zimmerman_26-598\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch08'; sex{end+1} = 'Male';   code{end+1} = '08-100';
fname{end+1} = 'zimmerman_26\zimmerman_26-599\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch08'; sex{end+1} = 'Male';   code{end+1} = '08-101';
fname{end+1} = 'zimmerman_26\zimmerman_26-224\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch08'; sex{end+1} = 'Female'; code{end+1} = '08-102';
fname{end+1} = 'zimmerman_26\zimmerman_26-578\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch08'; sex{end+1} = 'Female'; code{end+1} = '08-103';
fname{end+1} = 'zimmerman_26\zimmerman_26-580\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch08'; sex{end+1} = 'Female'; code{end+1} = '08-104';

fname{end+1} = 'zimmerman_27\zimmerman_27-220\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch08'; sex{end+1} = 'Male';   code{end+1} = '08-105';
fname{end+1} = 'zimmerman_27\zimmerman_27-221\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch08'; sex{end+1} = 'Male';   code{end+1} = '08-106';
fname{end+1} = 'zimmerman_27\zimmerman_27-581\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch08'; sex{end+1} = 'Male';   code{end+1} = '08-107';
fname{end+1} = 'zimmerman_27\zimmerman_27-595\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch08'; sex{end+1} = 'Female'; code{end+1} = '08-108';
fname{end+1} = 'zimmerman_27\zimmerman_27-596\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch08'; sex{end+1} = 'Female'; code{end+1} = '08-109';
fname{end+1} = 'zimmerman_27\zimmerman_27-597\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch08'; sex{end+1} = 'Female'; code{end+1} = '08-110';

%%%%% BATCH 1 %%%%%
fname{end+1} = 'zimmerman_28\zimmerman_28-276\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch09'; sex{end+1} = 'Male';   code{end+1} = '09-111';
fname{end+1} = 'zimmerman_28\zimmerman_28-277\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch09'; sex{end+1} = 'Male';   code{end+1} = '09-112';
fname{end+1} = 'zimmerman_28\zimmerman_28-297\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch09'; sex{end+1} = 'Male';   code{end+1} = '09-113';
fname{end+1} = 'zimmerman_28\zimmerman_28-298\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch09'; sex{end+1} = 'Female'; code{end+1} = '09-114';
fname{end+1} = 'zimmerman_28\zimmerman_28-299\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch09'; sex{end+1} = 'Female'; code{end+1} = '09-115';
fname{end+1} = 'zimmerman_28\zimmerman_28-300\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch09'; sex{end+1} = 'Female'; code{end+1} = '09-116';

fname{end+1} = 'zimmerman_29\zimmerman_29-295\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch09'; sex{end+1} = 'Male';   code{end+1} = '09-117';
fname{end+1} = 'zimmerman_29\zimmerman_29-296\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch09'; sex{end+1} = 'Male';   code{end+1} = '09-118';
fname{end+1} = 'zimmerman_29\zimmerman_29-235\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch09'; sex{end+1} = 'Female'; code{end+1} = '09-119';
fname{end+1} = 'zimmerman_29\zimmerman_29-236\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch09'; sex{end+1} = 'Female'; code{end+1} = '09-120';
fname{end+1} = 'zimmerman_29\zimmerman_29-278\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch09'; sex{end+1} = 'Female'; code{end+1} = '09-121';
fname{end+1} = 'zimmerman_29\zimmerman_29-294\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch09'; sex{end+1} = 'Female'; code{end+1} = '09-122';

x = [];
for i = 1:length(fname)
    try
    x(:,i) = csvread([path.main,fname{i},'\rawdata\resolution_3.6x\cellfinder_stats.csv']);
    catch
        disp(fname{i})
    end
end
x = (x(2,:)-x(3,:))./x(2,:)*100;
disp(['mean: ',num2str(nanmean(x))])
disp(['s.e.m.: ',num2str(nanstd(x)./sqrt(length(x)))])
disp(['max: ',num2str(nanmax(x))])
disp(['min: ',num2str(nanmin(x))])
    
warning('off')
for i = 1:length(fname)
    try
        raw_data{i} = readtable([path.main,fname{i},path.file]);
        %disp([size(raw_data{i},2) sum(strcmp('CentralAmygdalarNucleus_VentralPart',raw_data{i}.Properties.VariableNames))])
    catch
        disp(['Missing file: ',path.main,fname{i},path.file]);
    end
end
warning('on')
counts = [];
for j = 1:length(fname)
    counts(j,:) = raw_data{j}{1,:};
    %counts(j,:) = raw_data{j}{2,:};
end

[~,idx] = sort(fname);
A = cellfun(@(x) x(1:12),fname,'UniformOutput',false)';
B = cellfun(@(x) x(27:29),fname,'UniformOutput',false)';
C = cellfun(@(x) x(end),fname,'UniformOutput',false)';
D = counts(:,find(cellfun(@(x) isequal("ParabrachialNucleus",x),RegionLibrary.full{:,2})));
E = counts(:,1) - counts(:,find(cellfun(@(x) isequal("Cerebellum",x),RegionLibrary.full{:,2})));
tbl = table([A(idx)],[B(idx)],[C(idx)],[flavor(idx)]',[phase(idx)]',[sex(idx)]',[batch(idx)]',D,E);
tbl.Properties.VariableNames = {'Cohort','Sample','Imaging session','Flavor','Phase','Sex','Batch code','PBN Fos','Total Fos'};
writetable(tbl,'sample_info.csv')
%% Save data for GLMM

% x = size(RegionLibrary.full,1);
% RegionLibrary.full(x+1,:) = RegionLibrary.full(find(cellfun(@(x) isequal("CentralAmygdalarNucleus_LateralPart",x),RegionLibrary.full{:,2})),:);
% RegionLibrary.full{x+1,1} = x;
% RegionLibrary.full{x+1,2} = "CentralAmygdalarNucleus_CombinedWithoutShell";
% RegionLibrary.full{x+1,3} = "Cerebral Nuclei";
% RegionLibrary.full{x+1,4} = true;
% RegionLibrary.full{x+1,5} = RegionLibrary.full{find(cellfun(@(x) isequal("CentralAmygdalarNucleus_MedialPart",x),RegionLibrary.full{:,2})),5}+RegionLibrary.full{find(cellfun(@(x) isequal("CentralAmygdalarNucleus_LateralPart",x),RegionLibrary.full{:,2})),5}+RegionLibrary.full{find(cellfun(@(x) isequal("CentralAmygdalarNucleus_CapsularPart",x),RegionLibrary.full{:,2})),5};
% RegionLibrary.reduced(end+1,:) = RegionLibrary.full(x+1,:);
% RegionLibrary.reduced{end,2} = "Central amygdalar nucleus, combined without shell";
% counts(:,x+1) = counts(:,find(cellfun(@(x) isequal("CentralAmygdalarNucleus_MedialPart",x),RegionLibrary.full{:,2})))+counts(:,find(cellfun(@(x) isequal("CentralAmygdalarNucleus_LateralPart",x),RegionLibrary.full{:,2})))+counts(:,find(cellfun(@(x) isequal("CentralAmygdalarNucleus_CapsularPart",x),RegionLibrary.full{:,2})));
%
% x = size(RegionLibrary.full,1);
% RegionLibrary.full(x+1,:) = RegionLibrary.full(find(cellfun(@(x) isequal("CentralAmygdalarNucleus_LateralPart",x),RegionLibrary.full{:,2})),:);
% RegionLibrary.full{x+1,1} = x;
% RegionLibrary.full{x+1,2} = "CentralAmygdalarNucleus_CombinedWithShell";
% RegionLibrary.full{x+1,3} = "Cerebral Nuclei";
% RegionLibrary.full{x+1,4} = true;
% RegionLibrary.full{x+1,5} = RegionLibrary.full{find(cellfun(@(x) isequal("CentralAmygdalarNucleus_MedialPart",x),RegionLibrary.full{:,2})),5}+RegionLibrary.full{find(cellfun(@(x) isequal("CentralAmygdalarNucleus_LateralPart",x),RegionLibrary.full{:,2})),5}+RegionLibrary.full{find(cellfun(@(x) isequal("CentralAmygdalarNucleus_CapsularPart",x),RegionLibrary.full{:,2})),5}+RegionLibrary.full{find(cellfun(@(x) isequal("CentralAmygdalarNucleus_DorsalPart",x),RegionLibrary.full{:,2})),5}+RegionLibrary.full{find(cellfun(@(x) isequal("CentralAmygdalarNucleus_VentralPart",x),RegionLibrary.full{:,2})),5};
% RegionLibrary.reduced(end+1,:) = RegionLibrary.full(x+1,:);
% RegionLibrary.reduced{end,2} = "Central amygdalar nucleus, combined with shell";
% counts(:,x+1) = counts(:,find(cellfun(@(x) isequal("CentralAmygdalarNucleus_MedialPart",x),RegionLibrary.full{:,2})))+counts(:,find(cellfun(@(x) isequal("CentralAmygdalarNucleus_LateralPart",x),RegionLibrary.full{:,2})))+counts(:,find(cellfun(@(x) isequal("CentralAmygdalarNucleus_CapsularPart",x),RegionLibrary.full{:,2})))+counts(:,find(cellfun(@(x) isequal("CentralAmygdalarNucleus_VentralPart",x),RegionLibrary.full{:,2})))+counts(:,find(cellfun(@(x) isequal("CentralAmygdalarNucleus_DorsalPart",x),RegionLibrary.full{:,2})));
%
% output_dir = GLMMdir;
% if exist(output_dir,'dir')
%     rmdir(output_dir,'s')
% end
% mkdir(output_dir)
% writematrix(RegionLibrary.reduced{:,2},[output_dir,'regions.csv']);
% writematrix(RegionLibrary.reduced{:,5},[output_dir,'region_sizes.csv']);
%
% y = cell(0,0);
% for i = 1:length(flavor)
%     y{i} = [phase{i},flavor{i}];
% end
%
% w0 = grp2idx(y);
% w = nan(size(w0));
% for i = 1:length(w0)
%     w(i) = 1./sum(w0==w0(i));
% end
% w = w./sum(w)*length(w);
%
% i = 0;
% c = counts;
% if contains(library,'pons')
%     offset = c(:,1) - c(:,find(cellfun(@(x) isequal("Cerebellum",x),RegionLibrary.full{:,2}))) - c(:,find(cellfun(@(x) isequal("Hindbrain",x),RegionLibrary.full{:,2})));
% elseif contains(library,'medulla')
%     offset = c(:,1) - c(:,find(cellfun(@(x) isequal("Cerebellum",x),RegionLibrary.full{:,2}))) - c(:,find(cellfun(@(x) isequal("Medulla",x),RegionLibrary.full{:,2})));
% elseif contains(library,'cerebellum')
%     offset = c(:,1) - c(:,find(cellfun(@(x) isequal("Cerebellum",x),RegionLibrary.full{:,2})));
% elseif contains(library,'all')
%     offset = c(:,1);
% else
%     disp('Problem disaster!')
%     return
% end
%
% PBN = c(:,find(cellfun(@(x) isequal("ParabrachialNucleus",x),RegionLibrary.full{:,2})));
% output_dir_i = [output_dir,num2str(i,'%04.f'),'/'];
% if exist(output_dir_i,'dir')
%     rmdir(output_dir_i,'s')
% end
% mkdir(output_dir_i)
%
% writecell(batch,[output_dir_i,'/batch.csv']);
% writecell(flavor,[output_dir_i,'/flavor.csv']);
% writecell(phase,[output_dir_i,'/phase.csv']);
% writecell(batch,[output_dir_i,'/batch.csv']);
% writecell(sex,[output_dir_i,'/sex.csv']);
% writematrix(offset,[output_dir_i,'/offset.csv']);
% writematrix(PBN,[output_dir_i,'/pbn.csv']);
% for j = 1:size(RegionLibrary.reduced,1)
%     writematrix(c(:,RegionLibrary.reduced.index(j)+1),[output_dir_i,'/counts_region_',num2str(j,'%04.f'),'.csv']);
% end
% writematrix(w,[output_dir_i,'/weights.csv']);
% return
%% Load GLMM results
clc
FDR = 0.1;
FWER = 0.05;

output = struct;
% output.regions_all = struct;
%output.regions.libraryname = library;
output.regions.name = RegionLibrary.reduced.region;
output.regions.parent = RegionLibrary.reduced.parent;
output.regions.size = RegionLibrary.reduced{:,5};
output.regions.index = RegionLibrary.reduced.index;
output.regions.counts = counts(:,RegionLibrary.reduced.index+1)';

output.GLMMinput.flavor = flavor;
output.GLMMinput.phase = phase;
output.GLMMinput.sex = sex;
output.GLMMinput.batch = batch;
if contains(library,'pons')
    output.GLMMinput.offset = counts(:,1) - counts(:,find(cellfun(@(x) isequal("Cerebellum",x),RegionLibrary.full{:,2}))) - counts(:,find(cellfun(@(x) isequal("Hindbrain",x),RegionLibrary.full{:,2})));
elseif contains(library,'medulla')
    output.GLMMinput.offset = counts(:,1) - counts(:,find(cellfun(@(x) isequal("Cerebellum",x),RegionLibrary.full{:,2}))) - counts(:,find(cellfun(@(x) isequal("Medulla",x),RegionLibrary.full{:,2})));
elseif contains(library,'cerebellum')
    output.GLMMinput.offset = counts(:,1) - counts(:,find(cellfun(@(x) isequal("Cerebellum",x),RegionLibrary.full{:,2})));
elseif contains(library,'all')
    output.GLMMinput.offset = counts(:,1);
else
    disp('Problem disaster!')
    return
end
output.GLMMinput.offset = counts(:,1);
output.GLMMinput.pbn = counts(:,find(cellfun(@(x) isequal("ParabrachialNucleus",x),RegionLibrary.full{:,2})));
output.GLMMinput.counts = counts(:,RegionLibrary.reduced.index+1);

idx = find(contains(RegionLibrary.full.region,'IntercalatedAmygdalarNucleus'));
RegionLibrary.full.index([idx-2 idx-1]) = -1; RegionLibrary.full.index(idx:idx+1) = RegionLibrary.full.index(idx:idx+1)-2;
idx = find(contains(RegionLibrary.reduced.region,'Intercalated amygdalar nucleus'));
RegionLibrary.reduced.index(idx:idx+1) = RegionLibrary.reduced.index(idx:idx+1)-2;
output.regions.name = RegionLibrary.reduced.region;
output.regions.parent = RegionLibrary.reduced.parent;
output.regions.size = RegionLibrary.reduced{:,5};
output.regions.index = RegionLibrary.reduced.index;

% %%%%%
% T = load(['regions_allen_no_cerebellum_NEW.mat'],'RegionLibrary');
% output.regions_all = struct;
% output.regions_all.name = T.RegionLibrary.full.region;
% output.regions_all.parent = T.RegionLibrary.full.parent;
% output.regions_all.size = T.RegionLibrary.full{:,5};
% output.regions_all.index = T.RegionLibrary.full.index;
% output.regions_all.counts = counts';
% %%%%%

%clear counts raw_data

results = struct;
results.Eq2.modelstats.Chi2stat = [];
results.Eq2.modelstats.pvalue_raw = [];

results.Eq2.coefficients.estimates = [];
results.Eq2.coefficients.std_errors = [];
results.Eq2.coefficients.Chi2stat = [];
results.Eq2.coefficients.pvalues_raw = [];
results.Eq2.coefficients.pvalues_corrected = [];

results.Eq2.flavor.estimates = [];
results.Eq2.flavor.std_errors = [];
results.Eq2.flavor.Zstat = [];
results.Eq2.flavor.pvalues_raw = [];
results.Eq2.flavor.pvalues_corrected = [];

results.Eq3.coefficients.estimates = [];
results.Eq3.coefficients.std_errors = [];
results.Eq3.coefficients.Zstat = [];
results.Eq3.coefficients.pvalues_raw = [];

results.Eq4.flavor.estimates = [];
results.Eq4.flavor.std_errors = [];
results.Eq4.flavor.Zstat = [];
results.Eq4.flavor.pvalues_raw = [];

k = 0;
m = 0;
for i = 1:size(RegionLibrary.reduced,1)
    file1 = [GLMMdir,num2str(k,'%04.f'),'/GLMM_output_region_',num2str(i,'%04.f'),'_coefficients.csv'];
    file2 = [GLMMdir,num2str(k,'%04.f'),'/GLMM_output_region_',num2str(i,'%04.f'),'_reducedmodels.csv'];
    file3 = [GLMMdir,num2str(k,'%04.f'),'/GLMM_output_region_',num2str(i,'%04.f'),'_marginaleffects.csv'];
    file4 = [GLMMdir,num2str(k,'%04.f'),'/GLMM_output_region_',num2str(i,'%04.f'),'_CGRP_coefficients.csv'];
    file5 = [GLMMdir,num2str(k,'%04.f'),'/GLMM_output_region_',num2str(i,'%04.f'),'_CGRP_marginaleffects.csv'];
    if exist(file1) & exist(file2) & exist(file3) & exist(file4) & exist(file5)
        
        temp = readmatrix(file1,'Delimiter',' ');
        t = [];
        t(1,:) = rmmissing(temp(1,:));
        t(2,:) = rmmissing(temp(2,:));
        t(3,:) = rmmissing(temp(3,:));
        t(4,:) = rmmissing(temp(4,:));
        t(5,:) = rmmissing(temp(5,:));
        t(6,:) = rmmissing(temp(6,:));
        t(7,:) = rmmissing(temp(7,:));
        try
            results.Eq2.coefficients.estimates(i,:) = t(:,1);
            results.Eq2.coefficients.std_errors(i,:) = t(:,2);
            
            temp = readmatrix(file2);
            results.Eq2.coefficients.Chi2stat(i,1:4) = temp([9,10,11,12],3);
            results.Eq2.coefficients.pvalues_raw(i,1:4) = temp([4,5,6,7],3);
             results.Eq2.coefficients.pvalues_corrected(i,1:4) = multicmp(results.Eq2.coefficients.pvalues_raw(i,:),'up',FWER);
%             results.Eq2.modelstats.BIC(i,1) = temp(1,3);
%             results.Eq2.modelstats.LL(i,1) = temp(2,3);
            results.Eq2.modelstats.Chi2stat(i,1) = temp(8,3);
            results.Eq2.modelstats.pvalue_raw(i,1) = temp(3,3);
            
            temp = readmatrix(file3,'Delimiter',' ','Range','A1:ZZZ10');
            results.Eq2.flavor.estimates(i,1:3) = rmmissing(temp(1,:));
            results.Eq2.flavor.std_errors(i,1:3) = rmmissing(temp(2,:));
            results.Eq2.flavor.pvalues_raw(i,rmmissing(temp(2,:))~=0) = rmmissing(temp(3,:));
            results.Eq2.flavor.pvalues_raw(i,rmmissing(temp(2,:))==0) = NaN;
            results.Eq2.flavor.Zstat(i,1:3) = results.Eq2.flavor.estimates(i,1:3)./results.Eq2.flavor.std_errors(i,1:3);
            results.Eq2.flavor.pvalues_corrected(i,1:3) = multicmp(results.Eq2.flavor.pvalues_raw(i,:),'up',FWER);
%             results.marginaleffects.phase.estimates(i,1:3) = rmmissing(temp(4,:));
%             results.marginaleffects.phase.std_errors(i,1:3) = rmmissing(temp(5,:));
%             results.marginaleffects.phase.pvalues_raw(i,1:3) = rmmissing(temp(6,:));
%             results.marginaleffects.phase.pvalues_corrected(i,1:3) = multicmp(results.marginaleffects.phase.pvalues_raw(i,:),'up',FWER);
            
            temp = readmatrix(file4,'Delimiter',' ');
            t = [];
            t(1,:) = rmmissing(temp(1,:));
            t(2,:) = rmmissing(temp(2,:));
            t(3,:) = rmmissing(temp(3,:));
            t(4,:) = rmmissing(temp(4,:));
            t(5,:) = rmmissing(temp(5,:));
            results.Eq3.coefficients.estimates(i,:) = t(:,1);
            results.Eq3.coefficients.std_errors(i,:) = t(:,2);
            results.Eq3.coefficients.Zstat(i,:) = t(:,3);
            results.Eq3.coefficients.pvalues_raw(i,:) = t(:,4);
            temp = readmatrix(file5,'Delimiter',' ','Range','A1:ZZZ10');
            results.Eq4.flavor.estimates(i,:) = rmmissing(temp(1,:));
            results.Eq4.flavor.std_errors(i,:) = rmmissing(temp(2,:));
            results.Eq4.flavor.pvalues_raw(i,rmmissing(temp(2,:))~=0) = rmmissing(temp(3,:));
            results.Eq4.flavor.pvalues_raw(i,rmmissing(temp(2,:))==0) = NaN;
            results.Eq4.flavor.Zstat(i,:) = results.Eq4.flavor.estimates(i,:)./results.Eq4.flavor.std_errors(i,:);
%             results.Eq4.flavor.pvalues_corrected(i,:) = multicmp(results.Eq4.flavor.pvalues_raw(i,:),'up',FWER);
%             results.Eq4.phase.estimates(i,:) = rmmissing(temp(4,:));
%             results.Eq4.phase.std_errors(i,:) = rmmissing(temp(5,:));
%             if size(temp,1)==6
%                 results.Eq4.phase.pvalues_raw(i,:) = rmmissing(temp(6,:));
%             else
%                 results.Eq4.phase.pvalues_raw(i,1:5) = rmmissing(temp(6,:));
%                 results.Eq4.phase.pvalues_raw(i,6) = rmmissing(temp(7,:));
%             end
%             results.Eq4.phase.pvalues_corrected(i,:) = multicmp(results.Eq4.phase.pvalues_raw(i,:),'up',FWER);
            
            disp(['Finished loading data for region #',num2str(i,'%03.f'),': ',output.regions.name{i,:}])
        catch
            temp = readmatrix(file2);
            results.Eq2.modelstats.pvalue_raw(i,1) = temp(3,3);
            if isnan(results.Eq2.modelstats.pvalue_raw(i,1))
                results.Eq2.modelstats.pvalue_raw(i,1) = 1;
            end
            disp(['<WARNING> Error loading GLMM output for region #',num2str(i,'%03.f'),': ',output.regions.name{i,:},' <WARNING>'])
        end
    else
        m = m + 1;
    end
end

h = fdr_bky(results.Eq2.modelstats.pvalue_raw,FDR);
% for i = 1:size(results.Eq2.coefficients.pvalues_corrected,1)
%     if h(i)
%         results.Eq2.modelstats.significant{i} = '*';
%     else
%         results.Eq2.modelstats.significant{i} = '';
%     end
%     for j = 1:size(results.Eq2.coefficients.pvalues_corrected,2)
%         if results.Eq2.coefficients.pvalues_corrected(i,j) <= 0.0001
%             results.Eq2.coefficients.significant{i,j} = '****';
%         elseif results.Eq2.coefficients.pvalues_corrected(i,j) <= 0.001
%             results.Eq2.coefficients.significant{i,j} = '***';
%         elseif results.Eq2.coefficients.pvalues_corrected(i,j) <= 0.01
%             results.Eq2.coefficients.significant{i,j} = '**';
%         elseif results.Eq2.coefficients.pvalues_corrected(i,j) <= 0.05
%             results.Eq2.coefficients.significant{i,j} = '*';
%         else
%             results.Eq2.coefficients.significant{i,j} = '';
%         end
%         if j < 4
%             if results.Eq2.flavor.pvalues_corrected(i,j) <= 0.0001
%                 results.Eq2.flavor.significant{i,j} = '****';
%             elseif results.Eq2.flavor.pvalues_corrected(i,j) <= 0.001
%                 results.Eq2.flavor.significant{i,j} = '***';
%             elseif results.Eq2.flavor.pvalues_corrected(i,j) <= 0.01
%                 results.Eq2.flavor.significant{i,j} = '**';
%             elseif results.Eq2.flavor.pvalues_corrected(i,j) <= 0.05
%                 results.Eq2.flavor.significant{i,j} = '*';
%             else
%                 results.Eq2.flavor.significant{i,j} = '';
%             end
%         end
% %         if j < 4
% %             if results.marginaleffects.phase.pvalues_corrected(i,j) <= 0.0001
% %                 results.marginaleffects.phase.significant{i,j} = '****';
% %             elseif results.marginaleffects.phase.pvalues_corrected(i,j) <= 0.001
% %                 results.marginaleffects.phase.significant{i,j} = '***';
% %             elseif results.marginaleffects.phase.pvalues_corrected(i,j) <= 0.01
% %                 results.marginaleffects.phase.significant{i,j} = '**';
% %             elseif results.marginaleffects.phase.pvalues_corrected(i,j) <= 0.05
% %                 results.marginaleffects.phase.significant{i,j} = '*';
% %             else
% %                 results.marginaleffects.phase.significant{i,j} = '';
% %             end
% %         end
%     end
% end
results.Eq2.modelstats.significant = h';

results.CEA.counts = counts(:,find(cellfun(@(x) isequal("CentralAmygdalarNucleus",x),RegionLibrary.full{:,2})));
% results.CEA.Eq2.modelstats.Chi2stat = [];
% results.CEA.Eq2.modelstats.pvalue_raw = [];
% results.CEA.Eq2.coefficients.estimates = [];
% results.CEA.Eq2.coefficients.std_errors = [];
% results.CEA.Eq2.coefficients.Chi2stat = [];
% results.CEA.Eq2.coefficients.pvalues_raw = [];
% results.CEA.Eq2.coefficients.pvalues_corrected = [];
results.CEA.Eq2.flavor.estimates = [];
results.CEA.Eq2.flavor.std_errors = [];
results.CEA.Eq2.flavor.Zstat = [];
results.CEA.Eq2.flavor.pvalues_raw = [];
results.CEA.Eq2.flavor.pvalues_corrected = [];
% results.CEA.Eq3.coefficients.estimates = [];
% results.CEA.Eq3.coefficients.std_errors = [];
% results.CEA.Eq3.coefficients.Zstat = [];
% results.CEA.Eq3.coefficients.pvalues_raw = [];
results.CEA.Eq4.flavor.estimates = [];
results.CEA.Eq4.flavor.std_errors = [];
results.CEA.Eq4.flavor.Zstat = [];
results.CEA.Eq4.flavor.pvalues_raw = [];

for i = 202
    k=0;
    file1 = [GLMMdir,num2str(k,'%04.f'),'/GLMM_output_region_',num2str(i,'%04.f'),'_coefficients.csv'];
    file2 = [GLMMdir,num2str(k,'%04.f'),'/GLMM_output_region_',num2str(i,'%04.f'),'_reducedmodels.csv'];
    file3 = [GLMMdir,num2str(k,'%04.f'),'/GLMM_output_region_',num2str(i,'%04.f'),'_marginaleffects.csv'];
    file4 = [GLMMdir,num2str(k,'%04.f'),'/GLMM_output_region_',num2str(i,'%04.f'),'_CGRP_coefficients.csv'];
    file5 = [GLMMdir,num2str(k,'%04.f'),'/GLMM_output_region_',num2str(i,'%04.f'),'_CGRP_marginaleffects.csv'];
    temp = readmatrix(file1,'Delimiter',' ');
    t = [];
    t(1,:) = rmmissing(temp(1,:));
    t(2,:) = rmmissing(temp(2,:));
    t(3,:) = rmmissing(temp(3,:));
    t(4,:) = rmmissing(temp(4,:));
    t(5,:) = rmmissing(temp(5,:));
    t(6,:) = rmmissing(temp(6,:));
    t(7,:) = rmmissing(temp(7,:));
%     results.CEA.Eq2.coefficients.estimates(1,:) = t(:,1);
%     results.CEA.Eq2.coefficients.std_errors(1,:) = t(:,2);
    
    temp = readmatrix(file2);
%     results.CEA.Eq2.modelstats.Chi2stat(1,1) = temp(8,3);
%     results.CEA.Eq2.modelstats.pvalue_raw(1,1) = temp(3,3);
%     results.CEA.Eq2.coefficients.Chi2stat(1,1:4) = temp([9,10,11,12],3);
%     results.CEA.Eq2.coefficients.pvalues_raw(1,1:4) = temp([4,5,6,7],3);
%     results.CEA.Eq2.coefficients.pvalues_corrected(1,1:4) = multicmp(results.CEA.Eq2.coefficients.pvalues_raw(1,:),'up',FWER);
    
    temp = readmatrix(file3,'Delimiter',' ','Range','A1:ZZZ10');
    results.CEA.Eq2.flavor.estimates(1,1:3) = rmmissing(temp(1,:));
    results.CEA.Eq2.flavor.std_errors(1,1:3) = rmmissing(temp(2,:));
    results.CEA.Eq2.flavor.pvalues_raw(1,rmmissing(temp(2,:))~=0) = rmmissing(temp(3,:));
    results.CEA.Eq2.flavor.pvalues_raw(1,rmmissing(temp(2,:))==0) = NaN;
    results.CEA.Eq2.flavor.Zstat(1,1:3) = results.CEA.Eq2.flavor.estimates(1,1:3)./results.CEA.Eq2.flavor.std_errors(1,1:3);
    results.CEA.Eq2.flavor.pvalues_corrected(1,1:3) = multicmp(results.CEA.Eq2.flavor.pvalues_raw(1,:),'up',FWER);
    
    temp = readmatrix(file4,'Delimiter',' ');
    t = [];
    t(1,:) = rmmissing(temp(1,:));
    t(2,:) = rmmissing(temp(2,:));
    t(3,:) = rmmissing(temp(3,:));
    t(4,:) = rmmissing(temp(4,:));
    t(5,:) = rmmissing(temp(5,:));
%     results.CEA.Eq3.coefficients.estimates(1,:) = t(:,1);
%     results.CEA.Eq3.coefficients.std_errors(1,:) = t(:,2);
%     results.CEA.Eq3.coefficients.Zstat(1,:) = t(:,3);
%     results.CEA.Eq3.coefficients.pvalues_raw(1,:) = t(:,4);
    
    temp = readmatrix(file5,'Delimiter',' ','Range','A1:ZZZ10');
    results.CEA.Eq4.flavor.estimates(1,:) = rmmissing(temp(1,:));
    results.CEA.Eq4.flavor.std_errors(1,:) = rmmissing(temp(2,:));
    results.CEA.Eq4.flavor.pvalues_raw(1,rmmissing(temp(2,:))~=0) = rmmissing(temp(3,:));
    results.CEA.Eq4.flavor.pvalues_raw(1,rmmissing(temp(2,:))==0) = NaN;
    results.CEA.Eq4.flavor.Zstat(1,:) = results.CEA.Eq4.flavor.estimates(1,:)./results.CEA.Eq4.flavor.std_errors(1,:);
end

output.GLMMoutput = results;
data = output;
save(['sharing/Fos-GLMM-statistics.mat'],'data')
mkdir(['plots-png/',path.file(63:end-4)])
mkdir(['plots-eps/',path.file(63:end-4)])
save(['plots-png/',path.file(63:end-4),'/GLMM-statistics.mat'],'output')
save(['plots-eps/',path.file(63:end-4),'/GLMM-statistics.mat'],'output')
%% Summary table
[~,idx1] = sort(output.GLMMoutput.modelstats.pvalue_raw(:,1));

tbl0 = table(output.regions.index(idx1),output.regions.name(idx1),output.regions.parent(idx1),'VariableNames',{'Index','Name','Parent'});
tbl1 = table(output.GLMMoutput.modelstats.pvalue_raw(idx1),output.GLMMoutput.modelstats.significant(idx1),'VariableNames',{'P-value','Significant'});
tbl2 = table(output.GLMMoutput.coefficients.significant(idx1,1),output.GLMMoutput.coefficients.significant(idx1,2),output.GLMMoutput.coefficients.significant(idx1,3),output.GLMMoutput.coefficients.significant(idx1,4),'VariableNames',{'Flavor','Phase','Flavor:Phase','Sex'});
tbl3 = table(output.GLMMoutput.flavorME.significant(idx1,1),output.GLMMoutput.flavorME.significant(idx1,2),output.GLMMoutput.flavorME.significant(idx1,3),'VariableNames',{'Drinking','LiCl','Retrieval'});
tbl3a = table(output.GLMMoutput.flavorME.pvalues_corrected(idx1,1),output.GLMMoutput.flavorME.pvalues_corrected(idx1,2),output.GLMMoutput.flavorME.pvalues_corrected(idx1,3),'VariableNames',{'Drinking','LiCl','Retrieval'});
tbl4 = table(output.GLMMoutput.marginaleffects.phase.significant(idx1,1),output.GLMMoutput.marginaleffects.phase.significant(idx1,2),output.GLMMoutput.marginaleffects.phase.significant(idx1,3),'VariableNames',{'Drinking/LiCl','Drinking/Retrieval','LiCl/Retrieval'});
tbl4a = table(output.GLMMoutput.marginaleffects.phase.pvalues_corrected(idx1,1),output.GLMMoutput.marginaleffects.phase.pvalues_corrected(idx1,2),output.GLMMoutput.marginaleffects.phase.pvalues_corrected(idx1,3),'VariableNames',{'Drinking/LiCl','Drinking/Retrieval','LiCl/Retrieval'});
tbl5 = table(output.GLMMoutput.flavorME.estimates(idx1,1),output.GLMMoutput.flavorME.estimates(idx1,2),output.GLMMoutput.flavorME.estimates(idx1,3),'VariableNames',{'Drinking','LiCl','Retrieval'});
tbl6 = table(output.GLMMoutput.flavorME.std_errors(idx1,1),output.GLMMoutput.flavorME.std_errors(idx1,2),output.GLMMoutput.flavorME.std_errors(idx1,3),'VariableNames',{'Drinking','LiCl','Retrieval'});
tbl7 = table(output.GLMMoutput.marginaleffects.phase.estimates(idx1,1)./output.regions.size(idx1),output.GLMMoutput.marginaleffects.phase.estimates(idx1,2)./output.regions.size(idx1),output.GLMMoutput.marginaleffects.phase.estimates(idx1,3)./output.regions.size(idx1),'VariableNames',{'Drinking','LiCl','Retrieval'});
tbl8 = table(output.GLMMoutput.marginaleffects.phase.std_errors(idx1,1)./output.regions.size(idx1),output.GLMMoutput.marginaleffects.phase.std_errors(idx1,2)./output.regions.size(idx1),output.GLMMoutput.marginaleffects.phase.std_errors(idx1,3)./output.regions.size(idx1),'VariableNames',{'Drinking','LiCl','Retrieval'});
tbl = table(tbl0,tbl1,tbl3,tbl3a,tbl4,tbl4a,tbl5,tbl6,tbl7,tbl8);
%tbl = table(tbl0,tbl1,tbl2,tbl3,tbl4,tbl1b,tbl2b,tbl3b,tbl4b,tbl5,tbl6,tbl7,tbl8);
tbl.Properties.VariableNames = {'Region','Full model','Novel-Familiar Significant','Novel-Familiar P-value','Phase Significant','Phase P-value','Novel-Familiar AME estimate','Novel-Familiar AME std. error','Phase MM estimate','Phase MM std. error'};

clc
%disp(tbl)
%disp(' ')
%disp(['Total significant regions (full regression model): ',num2str(sum(output.GLMMoutput.modelstats.significant=='*')),'/',num2str(length(output.GLMMoutput.modelstats.significant)),' (',num2str(size(output.regions.name,1)-length(output.GLMMoutput.modelstats.significant)),' missing)'])

idx1 = find(output.GLMMoutput.modelstats.significant=='*');

tbl0 = table(output.regions.index(idx1),output.regions.name(idx1),output.regions.parent(idx1),'VariableNames',{'Index','Name','Parent'});
tbl1 = table(output.GLMMoutput.modelstats.pvalue_raw(idx1),output.GLMMoutput.modelstats.significant(idx1),'VariableNames',{'P-value','Significant'});
tbl5 = table(output.GLMMoutput.flavorME.estimates(idx1,1)./output.GLMMoutput.flavorME.std_errors(idx1,1),output.GLMMoutput.flavorME.estimates(idx1,2)./output.GLMMoutput.flavorME.std_errors(idx1,2),output.GLMMoutput.flavorME.estimates(idx1,3)./output.GLMMoutput.flavorME.std_errors(idx1,3),output.GLMMoutput.cgrp_stim.flavorME.estimates(idx1,2)./output.GLMMoutput.cgrp_stim.flavorME.std_errors(idx1,2),'VariableNames',{'Drinking','LiCl','Retrieval','CGRP'});

tbl = table(tbl0,tbl1,tbl5);
tbl.Properties.VariableNames = {'Region','Full model','Novel-Familiar'};

clc
disp(tbl)
disp(['Total significant regions: ',num2str(sum(output.GLMMoutput.modelstats.significant=='*')),'/',num2str(length(output.GLMMoutput.modelstats.significant)),' (',num2str(size(output.regions.name,1)-length(output.GLMMoutput.modelstats.significant)),' missing)'])
writetable(splitvars(tbl),['GLMM-statistics.csv'])
%% Figure: Drinking timepoint

cmap2 = flipud(cbrewer('div','RdBu',1000,'spline')); cmap2(cmap2<0) = 0;
cmap1 = cbrewer('seq','Greys',1000,'spline');

significant = find(output.GLMMoutput.modelstats.significant=='*');
Aa = output.GLMMoutput.flavorME.estimates(significant,1)/mean(output.GLMMinput.offset(find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase))))./output.regions.size(significant)*100;
Ab = output.GLMMoutput.flavorME.std_errors(significant,1)/mean(output.GLMMinput.offset(find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase))))./output.regions.size(significant)*100;
B = contains(output.GLMMoutput.flavorME.significant(significant,1),'*');

idx_fam = find(contains(output.GLMMoutput.flavorME.significant(:,1),'*')&output.GLMMoutput.flavorME.estimates(:,1)<0);
idx_nov = find(contains(output.GLMMoutput.flavorME.significant(:,1),'*')&output.GLMMoutput.flavorME.estimates(:,1)>0);

Aa1 = Aa(find(B & Aa<0));
Aa2 = Aa(find(~B));
Aa3 = Aa(find(B & Aa>0));

Ab1 = Ab(find(B & Aa<0));
Ab2 = Ab(find(~B));
Ab3 = Ab(find(B & Aa>0));

[Aa1,i1] = sort(Aa1);
[Aa2,i2] = sort(Aa2);
[Aa3,i3] = sort(Aa3);

Ab1 = Ab1(i1);
Ab2 = Ab2(i2);
Ab3 = Ab3(i3);

figure('Position', get(0, 'Screensize'))

subplot(2,2,1)
hold on
axis square
A = [Aa1;Aa2;Aa3];
B = [Ab1;Ab2;Ab3];
for i = 1:length(A)
    scatter(i,A(i),64,cmap2(round((A(i)/B(i)+5)*100),:),'o','filled','MarkerEdgeColor','none')
    plot([i i],[A(i)-B(i) A(i)+B(i)],'Color',cmap2(round((A(i)/B(i)+5)*100),:),'LineWidth',1)
end
xlabel('Brain Regions')
ylabel('Novel–Familiar (Δ Fos^+ Cells, %/mm^3)')
ylim([-1 1])
yticks(-1:.5:1)
xlim([0 length(A)+1])
xticks(0:25:200)
for i = 1:length(idx_fam)
    text(0.05,0.95-.035*(i-1),RegionLibrary.reduced{idx_fam(i),2},'Units','Normalized','FontSize',8,'VerticalAlignment','top','Color','b')
end
for i = 1:length(idx_nov)
    text(0.95,0.05+.035*(i-1),RegionLibrary.reduced{idx_nov(i),2},'Units','Normalized','FontSize',8,'VerticalAlignment','bottom','Color','r','HorizontalAlignment','right')
end
title('Whole Brain')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

subplot(2,2,2)
hold on
axis square
for i = 1:length(A)
    scatter(i,A(i)/B(i),64,cmap2(round((A(i)/B(i)+5)*100),:),'o','filled','MarkerEdgeColor','k')
end
plot([1 length(idx_fam)],[4.5 4.5],'k','LineWidth',1)
plot([length(A)-length(idx_nov)+1 length(A)],[4.5 4.5],'k','LineWidth',1)
ylim([-5 5])
xlim([0 length(A)+1])
yticks(-5:2.5:5)
xticks(0:25:200)
xlabel('Brain Regions')
ylabel('Novel–Familiar (z)')
% for i = 1:length(idx_fam)
% text(0.05,0.95-.035*(i-1),RegionLibrary.reduced{idx_fam(i),2},'Units','Normalized','FontSize',8,'VerticalAlignment','top','Color','b')
% end
% for i = 1:length(idx_nov)
% text(0.95,0.05+.035*(i-1),RegionLibrary.reduced{idx_nov(i),2},'Units','Normalized','FontSize',8,'VerticalAlignment','bottom','Color','r','HorizontalAlignment','right')
% end
title('Whole Brain')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

phases = {'Drinking','LiCl','Retrieval'};
counts_norm = [output.GLMMinput.counts./output.GLMMinput.offset./output.regions.size']*100;

subplot(2,2,3)
p_label = [];
hold on
%axis square
FL = cell(0,0);
out = [];
for j = 1:length(idx_nov)
    PHASEidx = find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase)&cellfun(@(x) isequal('Familiar',x),output.GLMMinput.flavor));
    x = counts_norm(PHASEidx,idx_nov(j));
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05-.125,x,64,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','w')
    plot([-0.125 0.125]+j-.125,[xx(3) xx(3)],'k','LineWidth',2)
    plot([j j]-.125,[xx(2) xx(4)],'k','LineWidth',2)
    plot([j j]-.125,[xx(1) xx(5)],'k','LineWidth',2)
    
    out(1:length(x),j) = x;
        FL(1:length(x)) = {'Familiar'};
    
    PHASEidx = find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase)&cellfun(@(x) isequal('Novel',x),output.GLMMinput.flavor));
    x = counts_norm(PHASEidx,idx_nov(j));
    xx = prctile(x,[10 25 50 75 90]);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05+.125,x,64,'filled','MarkerFaceColor',[1 .75 .75],'MarkerEdgeColor','w')
    plot([-0.125 0.125]+j+.125,[xx(3) xx(3)],'r','LineWidth',2)
    plot([j j]+.125,[xx(2) xx(4)],'r','LineWidth',2)
    plot([j j]+.125,[xx(1) xx(5)],'r','LineWidth',2)

        out(1+length(x):2*length(x),j) = x;
        FL(1+length(x):2*length(x)) = {'Novel'};

        %ylim([0 1])
    p_label = [p_label,'P = ',num2str(output.GLMMoutput.flavorME.pvalues_corrected(idx_nov(j),1),2),char(10)];
end
text(0.1,0.95,p_label,'Units','Normalized','FontSize',8,'VerticalAlignment','top')
xticks(1:length(idx_nov))
xticklabels(output.regions.name(idx_nov,:))
% ylim([0 2000])
% yticks(0:500:2000)
xlim([0.5 length(idx_nov)+.5])
ylabel('Fos^+ Cells (% of Total) per mm^3')
title('Novel regions')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

tbl = table(FL',out(:,1),out(:,2),out(:,3),out(:,4),out(:,5),out(:,6),out(:,7),out(:,8),out(:,9),out(:,10));
writetable(tbl,['Z:\Chris\matlab\cz\cta-source-data\Fig-ED1a.csv'])

subplot(2,2,4)
hold on
%axis square
p_label = [];
FL = cell(0,0);
out = [];
for j = 1:length(idx_fam)
    PHASEidx = find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase)&cellfun(@(x) isequal('Familiar',x),output.GLMMinput.flavor));
    x = counts_norm(PHASEidx,idx_fam(j));
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05-.125,x,64,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','w')
    plot([-0.125 0.125]+j-.125,[xx(3) xx(3)],'k','LineWidth',2)
    plot([j j]-.125,[xx(2) xx(4)],'k','LineWidth',2)
    plot([j j]-.125,[xx(1) xx(5)],'k','LineWidth',2)
    
    out(1:length(x),j) = x;
        FL(1:length(x)) = {'Familiar'};

        PHASEidx = find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase)&cellfun(@(x) isequal('Novel',x),output.GLMMinput.flavor));
    x = counts_norm(PHASEidx,idx_fam(j));
    xx = prctile(x,[10 25 50 75 90]);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05+.125,x,64,'filled','MarkerFaceColor',[1 .75 .75],'MarkerEdgeColor','w')
    plot([-0.125 0.125]+j+.125,[xx(3) xx(3)],'r','LineWidth',2)
    plot([j j]+.125,[xx(2) xx(4)],'r','LineWidth',2)
    plot([j j]+.125,[xx(1) xx(5)],'r','LineWidth',2)

            out(1+length(x):2*length(x),j) = x;
        FL(1+length(x):2*length(x)) = {'Novel'};

        %ylim([0 1])
    p_label = [p_label,'P = ',num2str(output.GLMMoutput.flavorME.pvalues_corrected(idx_fam(j),1),2),char(10)];
end
text(0.1,0.95,p_label,'Units','Normalized','FontSize',8,'VerticalAlignment','top')
xticks(1:length(idx_fam))
xticklabels(output.regions.name(idx_fam,:))
% ylim([0 2000])
% yticks(0:500:2000)
xlim([0.5 length(idx_fam)+.5])
ylabel('Fos^+ Cells (% of Total) per mm^3')
title('Familiar regions')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

tbl = table(FL',out(:,1),out(:,2),out(:,3),out(:,4),out(:,5),out(:,6),out(:,7),out(:,8),out(:,9),out(:,10));
writetable(tbl,['Z:\Chris\matlab\cz\cta-source-data\Fig-ED1b.csv'])

saveas(gcf,['plots-png/',path.file(63:end-4),'/plot-1-drinking-timepoint'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/plot-1-drinking-timepoint'],'epsc')
%% Figure: Brainwide shift

bins = -5:.0001:5;
Z1 = results.Eq2.flavor.Zstat(:,1);
Z2 = results.Eq2.flavor.Zstat(:,2);
Z3 = results.Eq2.flavor.Zstat(:,3);
N = length(results.Eq2.flavor.Zstat);

figure('Position', get(0, 'Screensize'))

subplot(1,4,1)
idx = find(output.GLMMoutput.modelstats.significant=='*');
N = length(idx);
[~,p1] = kstest2(Z1(idx),Z2(idx));
[~,p2] = kstest2(Z1(idx),Z3(idx));
[~,p3] = kstest2(Z2(idx),Z3(idx));
p = multicmp([p1 p2 p3],'up',0.05);
hold on
axis square
x = struct;
x.Drinking = Z1(idx);
x.Malaise = Z2(idx);
x.Retrieval = Z3(idx);
y=violinplot(x);
y(1).ViolinColor = [0 0 0];
y(2).ViolinColor = [.5 0 0];
y(3).ViolinColor = [1 0 0];
ylim([-5 5])
yticks(-5:2.5:5)
xlim([0.5 3.5])
ylabel('Novel–Familiar (z)')
title('Whole Brain')
text(0.05,0.95,['Drinking vs. LiCl: ',num2str(p(1),2),char(10),'Drinking vs. Retrieval: ',num2str(p(2),2),char(10),'LiCl vs. Retrieval: ',num2str(p(3),2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

subplot(1,4,2)
idx = intersect(find(output.GLMMoutput.modelstats.significant=='*'),find(output.regions.parent=='Cerebral Cortex'));
N = length(idx)
[~,p1] = kstest2(Z1(idx),Z2(idx));
[~,p2] = kstest2(Z1(idx),Z3(idx));
[~,p3] = kstest2(Z2(idx),Z3(idx));
p = multicmp([p1 p2 p3],'up',0.05);
hold on
axis square
x = struct;
x.Drinking = Z1(idx);
x.Malaise = Z2(idx);
x.Retrieval = Z3(idx);
y=violinplot(x);
y(1).ViolinColor = [0 0 0];
y(2).ViolinColor = [.5 0 0];
y(3).ViolinColor = [1 0 0];
ylim([-5 5])
yticks(-5:2.5:5)
xlim([0.5 3.5])
ylabel('Novel–Familiar (z)')
text(0.05,0.95,['Drinking vs. LiCl: ',num2str(p(1),2),char(10),'Drinking vs. Retrieval: ',num2str(p(2),2),char(10),'LiCl vs. Retrieval: ',num2str(p(3),2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
title('Forebrain/Cortex')
hold off

subplot(1,4,3)
idx = intersect(find(output.GLMMoutput.modelstats.significant=='*'),find([output.regions.parent=='Cerebral Nuclei'|output.regions.parent=='Thalamus'|output.regions.parent=='Hypothalamus']));
N = length(idx)
[~,p1] = kstest2(Z1(idx),Z2(idx));
[~,p2] = kstest2(Z1(idx),Z3(idx));
[~,p3] = kstest2(Z2(idx),Z3(idx));
p = multicmp([p1 p2 p3],'up',0.05);
hold on
axis square
x = struct;
x.Drinking = Z1(idx);
x.Malaise = Z2(idx);
x.Retrieval = Z3(idx);
y=violinplot(x);
y(1).ViolinColor = [0 0 0];
y(2).ViolinColor = [.5 0 0];
y(3).ViolinColor = [1 0 0];
ylim([-5 5])
yticks(-5:2.5:5)
xlim([0.5 3.5])
ylabel('Novel–Familiar (z)')
text(0.05,0.95,['Drinking vs. LiCl: ',num2str(p(1),2),char(10),'Drinking vs. Retrieval: ',num2str(p(2),2),char(10),'LiCl vs. Retrieval: ',num2str(p(3),2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
title('Forebrain/Subcortex')
hold off

subplot(1,4,4)
idx = intersect(find(output.GLMMoutput.modelstats.significant=='*'),find([output.regions.parent=='Pons'|output.regions.parent=='Medulla'|output.regions.parent=='Midbrain']));
N = length(idx)
[~,p1] = kstest2(Z1(idx),Z2(idx));
[~,p2] = kstest2(Z1(idx),Z3(idx));
[~,p3] = kstest2(Z2(idx),Z3(idx));
p = multicmp([p1 p2 p3],'up',0.05);
hold on
axis square
x = struct;
x.Drinking = Z1(idx);
x.Malaise = Z2(idx);
x.Retrieval = Z3(idx);
y=violinplot(x);
y(1).ViolinColor = [0 0 0];
y(2).ViolinColor = [.5 0 0];
y(3).ViolinColor = [1 0 0];
ylim([-5 5])
yticks(-5:2.5:5)
xlim([0.5 3.5])
ylabel('Novel–Familiar (z)')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
title('Midbrain/Hindbrain')
text(0.05,0.95,['Drinking vs. LiCl: ',num2str(p(1),2),char(10),'Drinking vs. Retrieval: ',num2str(p(2),2),char(10),'LiCl vs. Retrieval: ',num2str(p(3),2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top')
hold off

saveas(gcf,['plots-png/',path.file(63:end-4),'/plot-2-brainwide-shift'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/plot-2-brainwide-shift'],'epsc')
%% Figure: Clustering

cmap2 = flipud(cbrewer('div','RdBu',1000,'spline')); cmap2(cmap2<0) = 0;
cmap1 = cbrewer('seq','Greys',1000,'spline');

significant = find(output.GLMMoutput.modelstats.significant=='*');
treedata = results.Eq2.flavor.Zstat(significant,:);
treenames = output.regions.name(significant);
tree = linkage(treedata,'ward','chebychev');

figure('Position', get(0,'Screensize'))

subplot(1,2,1)
D = pdist(treedata);
leafOrder = optimalleaforder(tree,D);
[hLines,T,outperm] = dendrogram(tree,length(significant)+1,'Labels',[],'ColorThreshold',4.7,'Reorder',leafOrder,'Orientation','left');
for i = 1:length(hLines)
    zz(i,:) = hLines(i).Color;
end
N_clusters = length(unique(zz,'rows'))-1;
idx = cluster(tree,'maxclust',N_clusters);
idx = idx(outperm);
list = unique(idx,'stable');
amygdala_cluster = find(list==idx(find(RegionLibrary.reduced{significant(leafOrder),2}=="Central amygdalar nucleus, medial part")));
amygdala_regions = significant(leafOrder(idx==list(amygdala_cluster)));
xtickangle(90)
yticklabels([])
set(hLines,'LineWidth',1)
xlabel('Distance')
ylim([0.5 length(treedata)+0.5])
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')

ax=subplot(1,2,2);
hold on
heatmap(treedata(leafOrder,1:3),[],[],[],'Colormap',cmap2,'ColorLevels',1000,'MaxColorValue',5,'MinColorValue',-5,'NaNColor',[1 1 1]);
idx = cluster(tree,'maxclust',N_clusters);
idx = idx(outperm);
list = unique(idx,'stable');
counter = [];
for i = 1:length(list)-1
    counter(i) = sum(idx==list(i));
    plot([0 3]+0.5,[sum(counter) sum(counter)]+0.5,'k','LineWidth',1)
end
yticks(1:length(leafOrder))
yticklabels(treenames(leafOrder))
xticks(1:3)
xticklabels({'Drinking','LiCl','Retrieval'})
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
ax.XAxis.FontSize=16;
h = colorbar('FontSize',8);
ylabel(h,'Novel–Familiar (z)');
h.Ticks = -5:2.5:5;
h.TickLength = 0;
h.Box = 'on';
hold off

saveas(gcf,['plots-png/',path.file(63:end-4),'/plot-3-clustering'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/plot-3-clustering'],'epsc')
%% Table with cluster info

idx1 = find(output.GLMMoutput.modelstats.significant=='*');

tbl0 = table(output.regions.index(idx1),output.regions.name(idx1),output.regions.parent(idx1),'VariableNames',{'Index','Name','Parent'});
tbl1 = table(output.GLMMoutput.modelstats.pvalue_raw(idx1),output.GLMMoutput.modelstats.significant(idx1),'VariableNames',{'P-value','Significant'});
tbl5 = table(output.GLMMoutput.flavorME.estimates(idx1,1)./output.GLMMoutput.flavorME.std_errors(idx1,1),output.GLMMoutput.flavorME.estimates(idx1,2)./output.GLMMoutput.flavorME.std_errors(idx1,2),output.GLMMoutput.flavorME.estimates(idx1,3)./output.GLMMoutput.flavorME.std_errors(idx1,3),output.GLMMoutput.cgrp_stim.flavorME.estimates(idx1,2)./output.GLMMoutput.cgrp_stim.flavorME.std_errors(idx1,2),'VariableNames',{'Drinking','LiCl','Retrieval','CGRP'});


D = pdist(treedata);
leafOrder = optimalleaforder(tree,D);
[hLines,T,outperm] = dendrogram(tree,length(significant)+1,'Labels',[],'ColorThreshold',4.7,'Reorder',leafOrder,'Orientation','left');
for i = 1:length(hLines)
    zz(i,:) = hLines(i).Color;
end
N_clusters = length(unique(zz,'rows'))-1;
idx = cluster(tree,'maxclust',N_clusters);
idx = idx(outperm);
list = unique(idx,'stable');
%idx = list(idx);
[~,loc] = ismember(cluster(tree,'maxclust',N_clusters),list);
tbl6 = table(11-loc,'VariableNames',{'Cluster'});

tbl = table(tbl0,tbl1,tbl6,tbl5);
tbl.Properties.VariableNames = {'Region','Full model','Cluster','Novel-Familiar'};
disp(tbl)
tbl_sd = splitvars(tbl);

writetable(splitvars(tbl_sd),['Z:\Chris\matlab\cz\cta-source-data\Fig1e.csv'])
%% Figure: LiCl/CGRP 1

idx1 = find(cellfun(@(x) isequal('LiCl',x),output.GLMMinput.phase));
idx2 = find(cellfun(@(x) isequal('PBNCGRP',x),output.GLMMinput.phase));
significant = find(output.GLMMoutput.modelstats.significant=='*');


counts_norm = [output.GLMMinput.counts./output.GLMMinput.offset./output.regions.size']*100;
figure
hold on
axis square
idx = find(output.regions.name=='Parabrachial nucleus');
phases = {'Drinking','LiCl','PBNCGRP','Retrieval'};
out = [];
PH = cell(0,0);
for j = 1:length(phases)
    PHASEidx = find(cellfun(@(x) isequal(phases{j},x),output.GLMMinput.phase));
    x = counts_norm(PHASEidx,idx);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05,x,64,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','w')
    plot([-0.125 0.125]+j,[xx(3) xx(3)],'k','LineWidth',2)
    plot([j j],[xx(2) xx(4)],'k','LineWidth',2)
    plot([j j],[xx(1) xx(5)],'k','LineWidth',2)
    
    out(end+1:end+length(x)) = x;
    PH(end+1:end+length(x)) = phases(j);
end
text(0.05,0.95,['Drinking vs. LiCl: ',    num2str(output.GLMMoutput.cgrp_stim.phase.pvalues_corrected(idx,1),2),char(10),...
    'Drinking vs. CGRP: ',     num2str(output.GLMMoutput.cgrp_stim.phase.pvalues_corrected(idx,2),2),char(10),...
    'Drinking vs. Retrieval: ',num2str(output.GLMMoutput.cgrp_stim.phase.pvalues_corrected(idx,3),2),char(10),...
    'LiCl vs. CGRP: ',         num2str(output.GLMMoutput.cgrp_stim.phase.pvalues_corrected(idx,4),2),char(10),...
    'LiCl vs. Retrieval: ',   num2str(output.GLMMoutput.cgrp_stim.phase.pvalues_corrected(idx,5),2),char(10),...
    'CGRP vs. Retrieval: ',    num2str(output.GLMMoutput.cgrp_stim.phase.pvalues_corrected(idx,6),2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top')
xticks(1:4)
xticklabels({'Drinking','LiCl','CGRP','Retrieval'})
% ylim([0 5])
% yticks(0:1:5)
xlim([.5 4.5])
ylabel('Fos^+ Cells (% of Total) per mm^3')
title('Parabrachial nucleus')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

tbl = table(PH',out','VariableNames',{'Timepoint','PB Fos'});
writetable(tbl,['Z:\Chris\matlab\cz\cta-source-data\Fig-ED6b.csv'])

saveas(gcf,['plots-png/',path.file(63:end-4),'/cgrp-pb-summary'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/cgrp-pb-summary'],'epsc')
%% Figure: LiCl/CGRP 2
figure('Position', get(0,'Screensize'))
subplot(2,2,3)
hold on
axis square
A = results.Eq3.coefficients.Zstat(setdiff(significant,[find(output.regions.name=='Parabrachial nucleus');amygdala_regions]),2);
B = results.Eq3.coefficients.Zstat(setdiff(significant,[find(output.regions.name=='Parabrachial nucleus');amygdala_regions]),3);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([floor(min([A;B])) ceil(max([A;B]))])
ylim([floor(min([A;B])) ceil(max([A;B]))])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[r_O,p_O]=corr(A,B);
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.9 .9 .9],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,64,'k','filled','MarkerEdgeColor','w')
xlim([floor(min([A;B])) ceil(max([A;B]))])
ylim([floor(min([A;B])) ceil(max([A;B]))])
xticks([floor(min([A;B])),0,ceil(max([A;B]))])
yticks([floor(min([A;B])),0,ceil(max([A;B]))])
xlabel('LiCl Average Fos (Z/GLMM)')
ylabel('CGRP Average Fos (Z/GLMM)')
text(0.05,0.95,['r = ',num2str(r_O,3),char(10),'p = ',num2str(p_O,2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top','Color','k')
title('Other regions')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

out = [A';B'];

subplot(2,2,1)
hold on
axis square
A = results.Eq3.coefficients.Zstat(amygdala_regions,2);
B = results.Eq3.coefficients.Zstat(amygdala_regions,3);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([floor(min(A)) ceil(max(A))])
ylim([floor(min(B)) ceil(max(B))])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[r_A,p_A]=corr(A,B);
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[1 .9 .9],'LineStyle','none'),
plot(x,y_fit,'r','LineWidth',1)
scatter(A,B,64,'r','filled','MarkerEdgeColor','w')
xlim([floor(min(A)) ceil(max(A))])
ylim([floor(min(B)) ceil(max(B))])
xticks([floor(min(A)),ceil(max(A))])
yticks([floor(min(B)),ceil(max(B))])
hold off
xlabel('LiCl Average Fos (Z/GLMM)')
ylabel('CGRP Average Fos (Z/GLMM)')
text(0.05,0.95,['r = ',num2str(r_A,3),char(10),'p = ',num2str(p_A,2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top','Color','r')
title('Fos cluster 1')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

out = [A,B;out'];
tbl = table(out(:,1),out(:,2),'VariableNames',{'Malaise','CGRP stim'});
writetable(tbl,'Z:\Chris\matlab\cz\cta-source-data\Fig-ED6c.csv')

subplot(2,2,4)
hold on
axis square
A = results.Eq4.flavor.Zstat(setdiff(significant,[find(output.regions.name=='Parabrachial nucleus');amygdala_regions]),1);
B = results.Eq4.flavor.Zstat(setdiff(significant,[find(output.regions.name=='Parabrachial nucleus');amygdala_regions]),2);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([floor(min(A)) ceil(max(A))])
ylim([floor(min(B)) ceil(max(B))])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[r_O,p_O]=corr(A,B);
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.9 .9 .9],'LineStyle','none'),
xlim([floor(min(A)) ceil(max(A))])
ylim([floor(min(B)) ceil(max(B))])
xticks([floor(min(A)),ceil(max(A))])
yticks([floor(min(B)),ceil(max(B))])
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,64,'k','filled','MarkerEdgeColor','w')
xlabel('LiCl Novel–Familiar Fos (Z/GLMM)')
ylabel('CGRP Novel–Familiar Fos (Z/GLMM)')
text(0.05,0.95,['r = ',num2str(r_O,3),char(10),'p = ',num2str(p_O,2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top','Color','k')
title('Other regions')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

out = [A';B'];

subplot(2,2,2)
hold on
axis square
A = results.Eq4.flavor.Zstat(amygdala_regions,1);
B = results.Eq4.flavor.Zstat(amygdala_regions,2);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xticks([floor(min(A)),ceil(max(A))])
yticks([floor(min(B)),ceil(max(B))])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[r_A,p_A]=corr(A,B);
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[1 .9 .9],'LineStyle','none'),
plot(x,y_fit,'r','LineWidth',1)
scatter(A,B,64,'r','filled','MarkerEdgeColor','w')
xlim([floor(min(A)) ceil(max(A))])
ylim([floor(min(B)) ceil(max(B))])
xticks([floor(min(A)),ceil(max(A))])
yticks([floor(min(B)),ceil(max(B))])
hold off
xlabel('LiCl Novel–Familiar Fos (Z/GLMM)')
ylabel('CGRP Novel–Familiar Fos (Z/GLMM)')
text(0.05,0.95,['r = ',num2str(r_A,3),char(10),'p = ',num2str(p_A,2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top','Color','r')
title('Fos cluster 1')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

out = [A,B;out'];
tbl = table(out(:,1),out(:,2),'VariableNames',{'Malaise','CGRP stim'});
writetable(tbl,'Z:\Chris\matlab\cz\cta-source-data\Fig-ED6d.csv')

saveas(gcf,['plots-png/',path.file(63:end-4),'/plot-5-cgrp'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/plot-5-cgrp'],'epsc')
%% Figure: LiCl/CGRP 3

idx_drinking = find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase));
idx_retrieval = find(cellfun(@(x) isequal('Retrieval',x),output.GLMMinput.phase));
idx_licl = find(cellfun(@(x) isequal('LiCl',x),output.GLMMinput.phase));
idx_cgrp = find(cellfun(@(x) isequal('PBNCGRP',x),output.GLMMinput.phase));
idx_nov = find(cellfun(@(x) isequal('Novel',x),output.GLMMinput.flavor));
idx_fam = find(cellfun(@(x) isequal('Familiar',x),output.GLMMinput.flavor));

counts_norm = [];
total = output.GLMMinput.offset./output.GLMMinput.pbn;

idx = sort([idx_nov,idx_fam]);
counts_norm(intersect(idx_drinking,idx),:) = output.GLMMinput.counts(intersect(idx_drinking,idx),:)./output.GLMMinput.pbn(intersect(idx_drinking,idx));
counts_norm(intersect(idx_drinking,idx),:) = counts_norm(intersect(idx_drinking,idx),:)./mean(total(intersect(idx_drinking,idx)));

idx = sort([idx_nov,idx_fam]);
counts_norm(intersect(idx_licl,idx),:) = output.GLMMinput.counts(intersect(idx_licl,idx),:)./output.GLMMinput.pbn(intersect(idx_licl,idx));
counts_norm(intersect(idx_licl,idx),:) = counts_norm(intersect(idx_licl,idx),:)./mean(total(intersect(idx_licl,idx)));

idx = sort([idx_nov,idx_fam]);
counts_norm(intersect(idx_retrieval,idx),:) = output.GLMMinput.counts(intersect(idx_retrieval,idx),:)./output.GLMMinput.pbn(intersect(idx_retrieval,idx));
counts_norm(intersect(idx_retrieval,idx),:) = counts_norm(intersect(idx_retrieval,idx),:)./mean(total(intersect(idx_retrieval,idx)));

idx = sort([idx_nov,idx_fam]);
counts_norm(intersect(idx_cgrp,idx),:) = output.GLMMinput.counts(intersect(idx_cgrp,idx),:)./output.GLMMinput.pbn(intersect(idx_cgrp,idx));
counts_norm(intersect(idx_cgrp,idx),:) = counts_norm(intersect(idx_cgrp,idx),:)./mean(total(intersect(idx_cgrp,idx)));

counts_norm = [counts_norm./output.regions.size']*100;

significant = find(output.GLMMoutput.modelstats.significant=='*');
idx_licl = find(cellfun(@(x) isequal('LiCl',x),output.GLMMinput.phase));
idx_cgrp = find(cellfun(@(x) isequal('PBNCGRP',x),output.GLMMinput.phase));
idx_nov = find(cellfun(@(x) isequal('Novel',x),output.GLMMinput.flavor));
idx_fam = find(cellfun(@(x) isequal('Familiar',x),output.GLMMinput.flavor));
idx_amygdala = amygdala_regions;
idx_other = setdiff(significant,[find(output.regions.name=='Parabrachial nucleus');amygdala_regions]);

figure('Position', get(0,'Screensize'))
subplot(2,2,3)
hold on
axis square
A = mean(counts_norm(idx_licl,idx_other));
B = mean(counts_norm(idx_cgrp,idx_other));
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 ceil(max(A)*10)/10])
ylim([0 ceil(max(B)*10)/10])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[r_O,p_O]=corr(A',B');
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.9 .9 .9],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,64,'k','filled','MarkerEdgeColor','w')
xlim([0 ceil(max(A)*10)/10])
ylim([0 ceil(max(B)*10)/10])
xticks([floor(min(A)*10)/10 ceil(max(A)*10)/10])
yticks([floor(min(B)*10)/10 ceil(max(B)*10)/10])
xlabel('LiCl Average Fos (% per mm^3)')
ylabel('CGRP Average Fos (% per mm^3)')
text(0.05,0.95,['r = ',num2str(r_O,3),char(10),'p = ',num2str(p_O,2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top','Color','k')
title('Other regions')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

out = [A;B];

subplot(2,2,1)
hold on
axis square
A = mean(counts_norm(idx_licl,idx_amygdala));
B = mean(counts_norm(idx_cgrp,idx_amygdala));
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1.2])
ylim([0 ceil(max(B)*10)/10])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[r_A,p_A]=corr(A',B');
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[1 .9 .9],'LineStyle','none'),
plot(x,y_fit,'r','LineWidth',1)
scatter(A,B,64,'r','filled','MarkerEdgeColor','w')
xlim([0 1.2])
ylim([0 ceil(max(B)*10)/10])
xticks([floor(min(A)*10)/10 ceil(max(A)*10)/10])
yticks([floor(min(B)*10)/10 ceil(max(B)*10)/10])
hold off
xlabel('LiCl Average Fos (% per mm^3)')
ylabel('CGRP Average Fos (% per mm^3)')
text(0.05,0.95,['r = ',num2str(r_A,3),char(10),'p = ',num2str(p_A,2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top','Color','r')
title('Fos cluster 1')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

out = [A',B';out'];
tbl = table(out(:,1),out(:,2),'VariableNames',{'Malaise','CGRP stim'});
writetable(tbl,'Z:\Chris\matlab\cz\cta-source-data\Fig-2i.csv')

subplot(2,2,4)
hold on
axis square
A = nanmean(counts_norm(intersect(idx_licl,idx_nov),idx_other))-nanmean(counts_norm(intersect(idx_licl,idx_fam),idx_other));
B = nanmean(counts_norm(intersect(idx_cgrp,idx_nov),idx_other))-nanmean(counts_norm(intersect(idx_cgrp,idx_fam),idx_other));
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([floor(min(A)*10)/10 ceil(max(A)*10)/10])
ylim([floor(min(B)*10)/10 0.5])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[r_O,p_O]=corr(A',B');
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.9 .9 .9],'LineStyle','none'),
xlim([floor(min(A)*10)/10 ceil(max(A)*10)/10])
ylim([floor(min(B)*10)/10 0.5])
xticks([floor(min(A)*10)/10 ceil(max(A)*10)/10])
yticks([floor(min(B)*10)/10 ceil(max(B)*10)/10])
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,64,'k','filled','MarkerEdgeColor','w')
xlabel('LiCl Novel–Familiar Fos (% per mm^3)')
ylabel('CGRP Novel–Familiar Fos (% per mm^3)')
text(0.05,0.95,['r = ',num2str(r_O,3),char(10),'p = ',num2str(p_O,2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top','Color','k')
title('Other regions')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

out = [A;B];

subplot(2,2,2)
hold on
axis square
A = mean(counts_norm(intersect(idx_licl,idx_nov),idx_amygdala))-mean(counts_norm(intersect(idx_licl,idx_fam),idx_amygdala));
B = mean(counts_norm(intersect(idx_cgrp,idx_nov),idx_amygdala))-mean(counts_norm(intersect(idx_cgrp,idx_fam),idx_amygdala));
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([floor(min(A)*10)/10 ceil(max(A)*10)/10])
ylim([-0.1 ceil(max(B)*10)/10])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[r_A,p_A]=corr(A',B');
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[1 .9 .9],'LineStyle','none'),
plot(x,y_fit,'r','LineWidth',1)
scatter(A,B,64,'r','filled','MarkerEdgeColor','w')
xlim([floor(min(A)*10)/10 ceil(max(A)*10)/10])
ylim([-0.1 ceil(max(B)*10)/10])
xticks([floor(min(A)*10)/10 ceil(max(A)*10)/10])
yticks([floor(min(B)*10)/10 ceil(max(B)*10)/10])
hold off
xlabel('LiCl Novel–Familiar Fos (% per mm^3)')
ylabel('CGRP Novel–Familiar Fos (% per mm^3)')
text(0.05,0.95,['r = ',num2str(r_A,3),char(10),'p = ',num2str(p_A,2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top','Color','r')
title('Fos cluster 1')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

out = [A',B';out'];
tbl = table(out(:,1),out(:,2),'VariableNames',{'Malaise','CGRP stim'});
writetable(tbl,'Z:\Chris\matlab\cz\cta-source-data\Fig-2j.csv')

saveas(gcf,['plots-png/',path.file(63:end-4),'/plot-5bb-cgrp'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/plot-5bb-cgrp'],'epsc')
%% Figure: Correlation matrix all regions

counts_norm = [output.GLMMinput.counts./output.GLMMinput.offset./output.regions.size']*100;
significant = find(output.GLMMoutput.modelstats.significant=='*');

figure('Position', get(0,'Screensize'))

subplot(1,5,1)
tree = linkage(treedata,'ward','chebychev');
treedata = results.Eq2.flavor.Zstat(significant,:);
treenames = output.regions.name(significant);
D = pdist(treedata);
leafOrder = optimalleaforder(tree,D);
[hLines,T,outperm] = dendrogram(tree,length(significant)+1,'Labels',[],'ColorThreshold',4.7,'Reorder',leafOrder,'Orientation','left');
for i = 1:length(hLines)
    zz(i,:) = hLines(i).Color;
end
N_clusters = length(unique(zz,'rows'))-1;
set(hLines,'LineWidth',1)
xtickangle(90)
yticklabels([])
xlabel('Distance')
ylim([0.5 length(treedata)+0.5])
title('Clustering')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')

ax=subplot(1,5,2);
hold on
heatmap(treedata(leafOrder,1:3),[],[],[],'Colormap',cmap2,'ColorLevels',1000,'MaxColorValue',5,'MinColorValue',-5,'NaNColor',[1 1 1]);
idx = cluster(tree,'maxclust',N_clusters);
idx = idx(outperm);
list = unique(idx,'stable');
counter = [];
for i = 1:length(list)-1
    counter(i) = sum(idx==list(i));
    plot([0 3]+0.5,[sum(counter) sum(counter)]+0.5,'k','LineWidth',1)
end
yticks(1:length(leafOrder))
yticklabels(treenames(leafOrder))
xticks(1:3)
xticklabels({'Drinking','LiCl','Retrieval'})
title('GLMM estimates')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
ax.XAxis.FontSize=16;
ax.Title.FontSize=16;
h = colorbar('FontSize',8);
ylabel(h,'Novel–Familiar (z)');
h.Ticks = -5:2.5:5;
h.TickLength = 0;
h.Box = 'on';
hold off

A = counts_norm(find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase)),significant(leafOrder));
ax=subplot(1,5,3);
hold on
heatmap(corr(A),[],[],[],'Colormap',cmap2,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',-1,'NaNColor',[1 1 1]);
idx = cluster(tree,'maxclust',N_clusters);
idx = idx(outperm);
list = unique(idx,'stable');
counter = [];
for i = 1:length(list)-1
    counter(i) = sum(idx==list(i));
    plot([0 size(A,2)]+0.5,[sum(counter) sum(counter)]+0.5,'k','LineWidth',1)
    plot([sum(counter) sum(counter)]+0.5,[0 size(A,2)]+0.5,'k','LineWidth',1)
end
% yticks(1:length(leafOrder))
% yticklabels(treenames(leafOrder))
% xticks(1:3)
% xticklabels({'Drinking','LiCl','Retrieval'})
title('Drinking')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
ax.XAxis.FontSize=16;
h = colorbar('FontSize',8);
ylabel(h,'Correlation coefficient');
h.Ticks = -1:.5:1;
h.TickLength = 0;
h.Box = 'on';
hold off

A = counts_norm(find(cellfun(@(x) isequal('LiCl',x),output.GLMMinput.phase)),significant(leafOrder));
ax=subplot(1,5,4);
hold on
heatmap(corr(A),[],[],[],'Colormap',cmap2,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',-1,'NaNColor',[1 1 1]);
idx = cluster(tree,'maxclust',N_clusters);
idx = idx(outperm);
list = unique(idx,'stable');
counter = [];
for i = 1:length(list)-1
    counter(i) = sum(idx==list(i));
    plot([0 size(A,2)]+0.5,[sum(counter) sum(counter)]+0.5,'k','LineWidth',1)
    plot([sum(counter) sum(counter)]+0.5,[0 size(A,2)]+0.5,'k','LineWidth',1)
end
% yticks(1:length(leafOrder))
% yticklabels(treenames(leafOrder))
% xticks(1:3)
% xticklabels({'Drinking','LiCl','Retrieval'})
title('Malaise')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
ax.XAxis.FontSize=16;
h = colorbar('FontSize',8);
ylabel(h,'Correlation coefficient');
h.Ticks = -1:.5:1;
h.TickLength = 0;
h.Box = 'on';
hold off

A = counts_norm(find(cellfun(@(x) isequal('Retrieval',x),output.GLMMinput.phase)),significant(leafOrder));
ax=subplot(1,5,5);
hold on
heatmap(corr(A),[],[],[],'Colormap',cmap2,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',-1,'NaNColor',[1 1 1]);
idx = cluster(tree,'maxclust',N_clusters);
idx = idx(outperm);
list = unique(idx,'stable');
counter = [];
for i = 1:length(list)-1
    counter(i) = sum(idx==list(i));
    plot([0 size(A,2)]+0.5,[sum(counter) sum(counter)]+0.5,'k','LineWidth',1)
    plot([sum(counter) sum(counter)]+0.5,[0 size(A,2)]+0.5,'k','LineWidth',1)
end
% yticks(1:length(leafOrder))
% yticklabels(treenames(leafOrder))
% xticks(1:3)
% xticklabels({'Drinking','LiCl','Retrieval'})
title('Retrieval')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
ax.XAxis.FontSize=16;
h = colorbar('FontSize',8);
ylabel(h,'Correlation coefficient');
h.Ticks = -1:.5:1;
h.TickLength = 0;
h.Box = 'on';
hold off
saveas(gcf,['plots-png/',path.file(63:end-4),'/plot-6-correlations'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/plot-6-correlations'],'epsc')
%% Figure: Correlation matrix by cluster
figure('Position', get(0,'Screensize'))

idx = cluster(tree,'maxclust',N_clusters);
idx = idx(outperm);
list = unique(idx,'stable');


ax=subplot(1,4,1);
axis square
pl = [];
for i = 1:length(list)
    pl(i,:) = mean(treedata(leafOrder(find(idx==list(i))),:));
end
AME = pl;
hold on
heatmap(pl,[],[],[],'Colormap',cmap2,'ColorLevels',1000,'MaxColorValue',5,'MinColorValue',-5,'NaNColor',[1 1 1]);
xticks(1:3)
xticklabels({'Drinking','LiCl','Retrieval'})
title('GLMM estimates')
ylabel('Clusters')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
ax.XAxis.FontSize=16;
ax.Title.FontSize=16;
h = colorbar('FontSize',8);
ylabel(h,'Novel–Familiar (z)');
h.Ticks = -5:2.5:5;
h.TickLength = 0;
h.Box = 'on';
hold off


A = counts_norm(find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase)),significant(leafOrder));
B = corr(A);
pl = [];
list = unique(idx,'stable');
for i = 1:length(list)
    for j = 1:length(list)
        pl(i,j) = mean(B(find(idx==list(i)),find(idx==list(j))),'all');
    end
end
AvCor(:,1) = diag(pl);
AmygCor(:,1) = pl(:,find(list==idx(find(RegionLibrary.reduced{significant(leafOrder),2}=="Central amygdalar nucleus, medial part"))));
ax=subplot(1,4,2);
hold on
axis square
heatmap(pl,[],[],[],'Colormap',cmap2,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
% yticks(1:length(leafOrder))
% yticklabels(treenames(leafOrder))
% xticks(1:3)
% xticklabels({'Drinking','LiCl','Retrieval'})
xlabel('Clusters')
ylabel('Clusters')
title('Drinking')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
ax.XAxis.FontSize=16;
h = colorbar('FontSize',8);
ylabel(h,'Correlation coefficient');
h.Ticks = -1:.25:1;
h.TickLength = 0;
h.Box = 'on';
hold off

A = counts_norm(find(cellfun(@(x) isequal('LiCl',x),output.GLMMinput.phase)),significant(leafOrder));
B = corr(A);
pl = [];
list = unique(idx,'stable');
for i = 1:length(list)
    for j = 1:length(list)
        pl(i,j) = mean(B(find(idx==list(i)),find(idx==list(j))),'all');
    end
end
AvCor(:,2) = diag(pl);
AmygCor(:,2) = pl(:,find(list==idx(find(RegionLibrary.reduced{significant(leafOrder),2}=="Central amygdalar nucleus, medial part"))));
ax=subplot(1,4,3);
axis square
hold on
heatmap(pl,[],[],[],'Colormap',cmap2,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
% yticks(1:length(leafOrder))
% yticklabels(treenames(leafOrder))
% xticks(1:3)
% xticklabels({'Drinking','LiCl','Retrieval'})
xlabel('Clusters')
ylabel('Clusters')
title('Malaise')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
ax.XAxis.FontSize=16;
h = colorbar('FontSize',8);
ylabel(h,'Correlation coefficient');
h.Ticks = -1:.25:1;
h.TickLength = 0;
h.Box = 'on';
hold off

A = counts_norm(find(cellfun(@(x) isequal('Retrieval',x),output.GLMMinput.phase)),significant(leafOrder));
B = corr(A);
pl = [];
list = unique(idx,'stable');
for i = 1:length(list)
    for j = 1:length(list)
        pl(i,j) = mean(B(find(idx==list(i)),find(idx==list(j))),'all');
    end
end
AvCor(:,3) = diag(pl);
AmygCor(:,3) = pl(:,find(list==idx(find(RegionLibrary.reduced{significant(leafOrder),2}=="Central amygdalar nucleus, medial part"))));

ax=subplot(1,4,4);
axis square
hold on
heatmap(pl,[],[],[],'Colormap',cmap2,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
% yticks(1:length(leafOrder))
% yticklabels(treenames(leafOrder))
% xticks(1:3)
% xticklabels({'Drinking','LiCl','Retrieval'})
xlabel('Clusters')
ylabel('Clusters')
title('Retrieval')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
ax.XAxis.FontSize=16;
h = colorbar('FontSize',8);
ylabel(h,'Correlation coefficient');
h.Ticks = -1:.25:1;
h.TickLength = 0;
h.Box = 'on';
hold off
saveas(gcf,['plots-png/',path.file(63:end-4),'/plot-7-correlations'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/plot-7-correlations'],'epsc')
%% Figure: Correlation by cluster scatter plot
idx = cluster(tree,'maxclust',N_clusters);
idx = idx(outperm);

idx = setdiff(1:N_clusters,find(list==idx(find(RegionLibrary.reduced{significant(leafOrder),2}=="Central amygdalar nucleus, medial part"))));
figure
hold on
axis square
A = AME(idx,:);
B = AmygCor(idx,:);
a = scatter(A(:),B(:),64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A(:)',B(:),1);
xlim([-4 4])
ylim([-.4 .4])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[r,p]=corr(A(:),B(:));
[y_fit] = polyval(p2,x,S);
fitresult = fit(A(:),B(:),'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.9 .9 .9],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)

out = [];
for i = length(A):-1:1
    scatter(A(i,1),B(i,1),128,cmap2(round((A(i,1)+5)*100),:),'o','filled','MarkerEdgeColor','k')
    scatter(A(i,2),B(i,2),128,cmap2(round((A(i,2)+5)*100),:),'^','filled','MarkerEdgeColor','k')
    scatter(A(i,3),B(i,3),128,cmap2(round((A(i,3)+5)*100),:),'square','filled','MarkerEdgeColor','k')
    
    out(end+1,:) = [11-i A(i,1) B(i,1)];
    out(end+1,:) = [11-i A(i,2) B(i,2)];
    out(end+1,:) = [11-i A(i,3) B(i,3)];
    
end
% scatter(A(:,1),B(:,1),128,'k','o','filled')
% scatter(A(:,2),B(:,2),128,'k','square','filled')
% scatter(A(:,3),B(:,3),128,'k','diamond','filled')
xlim([-4 4])
ylim([-.4 .4])
xlabel('Novel–Familiar (z)')
ylabel('Amygdala-cluster Fos correlation')
text(0.05,0.95,['r = ',num2str(r,3),char(10),'p = ',num2str(p,2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

tbl = table(out(:,1),out(:,2),out(:,3),'VariableNames',{'Cluster','X','Y'});
writetable(tbl,['Z:\Chris\matlab\cz\cta-source-data\Fig-ED4c.csv'])

saveas(gcf,['plots-png/',path.file(63:end-4),'/plot-7b-correlations'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/plot-7b-correlations'],'epsc')
%% Figure: Correlation matrix examples


sizes = readmatrix([GLMMdir,'region_sizes.csv']);
counts_norm_CEA = [readmatrix([GLMMdir,'0000/counts_region_0202.csv'])./output.GLMMinput.offset./sizes(202)]*100;
counts_norm = [output.GLMMinput.counts./output.GLMMinput.offset./output.regions.size']*100;

figure('Position', get(0,'Screensize'))

idx = cluster(tree,'maxclust',N_clusters);
idx = idx(outperm);

ax=subplot(2,3,1);
axis square
hold on
list = unique(idx,'stable');
data = struct;
a = [];
for i = 1:length(list)
    x = treedata(leafOrder(find(idx==list(length(list)-i+1))),1);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    xval = (i-1)+1;
    scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,16,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','w')
    plot([-0.5 0.5]+xval,[xx(3) xx(3)],'k','LineWidth',1)
    plot([xval xval],[xx(2) xx(4)],'k','LineWidth',1)
    a(i) = signrank(x);
    x1(i,1) = mean(x);
    if i == 1
        signrank(x)
        x
    end
end
xlim([0 N_clusters+1])
ylim([-5 5])
xticks([1:N_clusters])
yticks(-5:2.5:5)
xticklabels(num2str([1:N_clusters]'))
xtickangle(0)
xlabel('Cluster')
ylabel('Novel–Familiar (z)')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
a = multicmp(a,'up',FWER)

ax=subplot(2,3,2);
axis square
hold on
list = unique(idx,'stable');
data = struct;
a = [];
for i = 1:length(list)
    x = treedata(leafOrder(find(idx==list(length(list)-i+1))),2);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    xval = (i-1)+1;
    scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,16,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','w')
    plot([-0.5 0.5]+xval,[xx(3) xx(3)],'k','LineWidth',1)
    plot([xval xval],[xx(2) xx(4)],'k','LineWidth',1)
    a(i) = signrank(x);
    x1(i,1) = mean(x);
    if i == 1
        signrank(x)
        x
    end
end
xlim([0 N_clusters+1])
ylim([-5 5])
xticks([1:N_clusters])
yticks(-5:2.5:5)
xticklabels(num2str([1:N_clusters]'))
xtickangle(0)
xlabel('Cluster')
ylabel('Novel–Familiar (z)')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
a = multicmp(a,'up',FWER)

ax=subplot(2,3,3);
axis square
hold on
list = unique(idx,'stable');
data = struct;
a = [];
for i = 1:length(list)
    x = treedata(leafOrder(find(idx==list(length(list)-i+1))),3);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    xval = (i-1)+1;
    scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,16,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','w')
    plot([-0.5 0.5]+xval,[xx(3) xx(3)],'k','LineWidth',1)
    plot([xval xval],[xx(2) xx(4)],'k','LineWidth',1)
    a(i) = signrank(x);
    if i == 1
        signrank(x)
        x
    end
    x1(i,1) = mean(x);
end
xlim([0 N_clusters+1])
ylim([-5 5])
xticks([1:N_clusters])
yticks(-5:2.5:5)
xticklabels(num2str([1:N_clusters]'))
xtickangle(0)
xlabel('Cluster')
ylabel('Novel–Familiar (z)')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
a = multicmp(a,'up',FWER)

% ax=subplot(2,3,1);
% axis square
% hold on
% list = unique(idx,'stable');
% data = struct;
% plot([0 N_clusters*4],[0 0],'k','LineWidth',1)
% a = []; b = []; c = [];
% x1 = [];
% for i = 1:length(list)
%     x = treedata(leafOrder(find(idx==list(i))),1);
%     xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
%     xval = (i-1)*4+1;
%     scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,16,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','w')
%     plot([-0.5 0.5]+xval,[xx(3) xx(3)],'k','LineWidth',1)
%     plot([xval xval],[xx(2) xx(4)],'k','LineWidth',1)
%     a(i) = signrank(x);
%     x1(i,1) = mean(x);
%
%     x = treedata(leafOrder(find(idx==list(i))),2);
%     xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
%     xval = (i-1)*4+2;
%     scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,16,'filled','MarkerFaceColor',[.85 .75 .75],'MarkerEdgeColor','w')
%     plot([-0.5 0.5]+xval,[xx(3) xx(3)],'Color',[.5 0 0],'LineWidth',1)
%     plot([xval xval],[xx(2) xx(4)],'Color',[.5 0 0],'LineWidth',1)
%     b(i) = signrank(x);
%     x1(i,2) = mean(x);
%
%     x = treedata(leafOrder(find(idx==list(i))),3);
%     xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
%     xval = (i-1)*4+3;
%     scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,16,'filled','MarkerFaceColor',[1 .75 .75],'MarkerEdgeColor','w')
%     plot([-0.5 0.5]+xval,[xx(3) xx(3)],'Color',[1 0 0],'LineWidth',1)
%     plot([xval xval],[xx(2) xx(4)],'Color',[1 0 0],'LineWidth',1)
%     c(i) = signrank(x);
%     x1(i,3) = mean(x);
% end
% xlim([0 N_clusters*4])
% ylim([-5 5])
% xticks([2:4:N_clusters*4-2])
% yticks(-5:2.5:5)
% xticklabels(num2str([1:N_clusters]'))
% xtickangle(0)
% xlabel('Cluster')
% ylabel('Novel–Familiar (z)')
% set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
% a = multicmp(a,'up',FWER)
% b = multicmp(b,'up',FWER)
% c = multicmp(c,'up',FWER)
%
A = counts_norm(:,significant(leafOrder));
[~,F]=mode(idx); CorrStruct = nan(N_clusters,N_clusters,3,F);
for i = 1:length(list)
    in = find(idx==list(i));
    for j = 1:length(in)
        for k = 1:length(list)
            temp1 = [];
            temp2 = [];
            temp3 = [];
            test = setdiff(find(idx==list(k)),in(j));
            for l = 1:length(test)
                temp1(1) = corr(A(find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase)),in(j)),A(find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase)),test(l)));
                temp2(1) = corr(A(find(cellfun(@(x) isequal('LiCl',x),output.GLMMinput.phase)),in(j)),A(find(cellfun(@(x) isequal('LiCl',x),output.GLMMinput.phase)),test(l)));
                temp3(1) = corr(A(find(cellfun(@(x) isequal('Retrieval',x),output.GLMMinput.phase)),in(j)),A(find(cellfun(@(x) isequal('Retrieval',x),output.GLMMinput.phase)),test(l)));
            end
            CorrStruct(i,k,1,j) = mean(temp1);
            CorrStruct(i,k,2,j) = mean(temp2);
            CorrStruct(i,k,3,j) = mean(temp3);
        end
    end
end
%
% subplot(2,3,2)
% axis square
% hold on
% plot([0 N_clusters*4],[0 0],'k','LineWidth',1)
% y1 = [];
for i = 1:length(list)
    
    x = squeeze(CorrStruct(i,i,1,:)); x(isnan(x)) = [];
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    xval = (i-1)*4+1;
    %     scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,16,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','w')
    %     plot([-0.5 0.5]+xval,[xx(3) xx(3)],'k','LineWidth',1)
    %     plot([xval xval],[xx(2) xx(4)],'k','LineWidth',1)
    y1(i,1) = mean(x);
    if i == 10
        x
        signrank(x)
    end
    
    x = squeeze(CorrStruct(i,i,2,:)); x(isnan(x)) = [];
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    xval = (i-1)*4+2;
    %     scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,16,'filled','MarkerFaceColor',[.85 .75 .75],'MarkerEdgeColor','w')
    %     plot([-0.5 0.5]+xval,[xx(3) xx(3)],'Color',[.5 0 0],'LineWidth',1)
    %     plot([xval xval],[xx(2) xx(4)],'Color',[.5 0 0],'LineWidth',1)
    y1(i,2) = mean(x);
    if i == 10
        x
        signrank(x)
    end
    
    x = squeeze(CorrStruct(i,i,3,:)); x(isnan(x)) = [];
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    xval = (i-1)*4+3;
    %     scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,16,'filled','MarkerFaceColor',[1 .75 .75],'MarkerEdgeColor','w')
    %     plot([-0.5 0.5]+xval,[xx(3) xx(3)],'Color',[1 0 0],'LineWidth',1)
    %     plot([xval xval],[xx(2) xx(4)],'Color',[1 0 0],'LineWidth',1)
    y1(i,3) = mean(x);
    if i == 10
        x
        signrank(x)
    end
    
end
% xlim([0 N_clusters*4])
% ylim([-1 1])
% xticks([2:4:N_clusters*4-2])
% yticks(-1:.5:1)
% xticklabels(num2str([1:N_clusters]'))
% xtickangle(0)
% xlabel('Cluster')
% ylabel('Within-cluster correlation')
% set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
%
% amygdala_cluster = find(list==idx(find(RegionLibrary.reduced{significant(leafOrder),2}=="Central amygdalar nucleus, medial part")));
% subplot(2,3,3)
% axis square
% hold on
% plot([0 N_clusters*4],[0 0],'k','LineWidth',1)
% for i = setdiff(1:length(list),amygdala_cluster,'stable')
%
%     x = squeeze(CorrStruct(i,amygdala_cluster,1,:)); x(isnan(x)) = [];
%     xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
%     xval = (i-1)*4+1;
%     scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,16,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','w')
%     plot([-0.5 0.5]+xval,[xx(3) xx(3)],'k','LineWidth',1)
%     plot([xval xval],[xx(2) xx(4)],'k','LineWidth',1)
%
%     x = squeeze(CorrStruct(i,amygdala_cluster,2,:)); x(isnan(x)) = [];
%     xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
%     xval = (i-1)*4+2;
%     scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,16,'filled','MarkerFaceColor',[.85 .75 .75],'MarkerEdgeColor','w')
%     plot([-0.5 0.5]+xval,[xx(3) xx(3)],'Color',[.5 0 0],'LineWidth',1)
%     plot([xval xval],[xx(2) xx(4)],'Color',[.5 0 0],'LineWidth',1)
%
%     x = squeeze(CorrStruct(i,amygdala_cluster,3,:)); x(isnan(x)) = [];
%     xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
%     xval = (i-1)*4+3;
%     scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,16,'filled','MarkerFaceColor',[1 .75 .75],'MarkerEdgeColor','w')
%     plot([-0.5 0.5]+xval,[xx(3) xx(3)],'Color',[1 0 0],'LineWidth',1)
%     plot([xval xval],[xx(2) xx(4)],'Color',[1 0 0],'LineWidth',1)
%
% end
% xlim([0 N_clusters*4])
% ylim([-1 1])
% xticks([2:4:N_clusters*4-2])
% yticks(-1:.5:1)
% xticklabels(num2str([1:N_clusters]'))
% xtickangle(0)
% xlabel('Cluster')
% ylabel(['Across-cluster correlation (Seed: CEAm, cl',num2str(amygdala_cluster),')'])
% set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
%%
figure('Position', get(0,'Screensize'))
subplot(2,3,4)
regionA = 202;%find(output.regions.name=='Central amygdalar nucleus, medial part');
regionB = find(output.regions.name=='Agranular insular area, posterior part');
hold on
axis square
A = counts_norm_CEA(find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase)));
B = counts_norm(find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase)),regionB);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1])
ylim([0 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[r,p]=corr(A,B);
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.9 .9 .9],'LineStyle','none','FaceAlpha',0.5)
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,64,'k','filled','MarkerEdgeColor','w')

A = counts_norm_CEA(find(cellfun(@(x) isequal('LiCl',x),output.GLMMinput.phase)));
B = counts_norm(find(cellfun(@(x) isequal('LiCl',x),output.GLMMinput.phase)),regionB);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1])
ylim([0 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[rb,pb]=corr(A,B);
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[1 .9 .9],'LineStyle','none'),
plot(x,y_fit,'r','LineWidth',1)
scatter(A,B,64,'r','filled','MarkerEdgeColor','w')

A = counts_norm_CEA(find(cellfun(@(x) isequal('Retrieval',x),output.GLMMinput.phase)));
B = counts_norm(find(cellfun(@(x) isequal('Retrieval',x),output.GLMMinput.phase)),regionB);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1])
ylim([0 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[rc,pc]=corr(A,B);
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[0.9 0.9 1],'LineStyle','none','FaceAlpha',0.5)
plot(x,y_fit,'b','LineWidth',1)
scatter(A,B,64,'b','filled','MarkerEdgeColor','w')
hold off
xlabel('CEA Fos')
ylabel('AIp Fos')
text(0.05,0.95,['Drinking:  r = ',num2str(r,3),', p = ',num2str(p,2),char(10),...
    'Malaise:   r = ',num2str(rb,3),', p = ',num2str(pb,2),char(10),...
    'Retrieval: r = ',num2str(rc,3),', p = ',num2str(pc,2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top')
title('CEA vs. AIp')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

subplot(2,3,5)
regionA = 202;%find(output.regions.name=='Central amygdalar nucleus, medial part');
regionB = find(output.regions.name=='Basolateral amygdalar nucleus, posterior part');
hold on
axis square
A = counts_norm_CEA(find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase)));
B = counts_norm(find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase)),regionB);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1])
ylim([0 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[r,p]=corr(A,B);
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.9 .9 .9],'LineStyle','none','FaceAlpha',0.5),
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,64,'k','filled','MarkerEdgeColor','w')

A = counts_norm_CEA(find(cellfun(@(x) isequal('LiCl',x),output.GLMMinput.phase)));
B = counts_norm(find(cellfun(@(x) isequal('LiCl',x),output.GLMMinput.phase)),regionB);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1])
ylim([0 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[rb,pb]=corr(A,B);
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[1 .9 .9],'LineStyle','none','FaceAlpha',0.5),
plot(x,y_fit,'r','LineWidth',1)
scatter(A,B,64,'r','filled','MarkerEdgeColor','w')

A = counts_norm_CEA(find(cellfun(@(x) isequal('Retrieval',x),output.GLMMinput.phase)));
B = counts_norm(find(cellfun(@(x) isequal('Retrieval',x),output.GLMMinput.phase)),regionB);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1])
ylim([0 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[rc,pc]=corr(A,B);
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[0.9 0.9 1],'LineStyle','none','FaceAlpha',0.5)
plot(x,y_fit,'b','LineWidth',1)
scatter(A,B,64,'b','filled','MarkerEdgeColor','w')
hold off
xlabel('CEA Fos')
ylabel('BLAp Fos')
text(0.05,0.95,['Drinking:  r = ',num2str(r,3),', p = ',num2str(p,2),char(10),...
    'Malaise:   r = ',num2str(rb,3),', p = ',num2str(pb,2),char(10),...
    'Retrieval: r = ',num2str(rc,3),', p = ',num2str(pc,2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top')
title('CEA vs. BLAp')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

subplot(2,3,6)
regionA = 202;%find(output.regions.name=='Central amygdalar nucleus, medial part');
regionB = find(output.regions.name=='Bed nuclei of the stria terminalis');
hold on
axis square
A = counts_norm_CEA(find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase)));
B = counts_norm(find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase)),regionB);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1])
ylim([0 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[r,p]=corr(A,B);
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.9 .9 .9],'LineStyle','none','FaceAlpha',0.5)
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,64,'k','filled','MarkerEdgeColor','w')

A = counts_norm_CEA(find(cellfun(@(x) isequal('LiCl',x),output.GLMMinput.phase)));
B = counts_norm(find(cellfun(@(x) isequal('LiCl',x),output.GLMMinput.phase)),regionB);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1])
ylim([0 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[rb,pb]=corr(A,B);
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[1 .9 .9],'LineStyle','none','FaceAlpha',0.5)
plot(x,y_fit,'r','LineWidth',1)
scatter(A,B,64,'r','filled','MarkerEdgeColor','w')

A = counts_norm_CEA(find(cellfun(@(x) isequal('Retrieval',x),output.GLMMinput.phase)));
B = counts_norm(find(cellfun(@(x) isequal('Retrieval',x),output.GLMMinput.phase)),regionB);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1])
ylim([0 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[rc,pc]=corr(A,B);
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[0.9 0.9 1],'LineStyle','none','FaceAlpha',0.5)
plot(x,y_fit,'b','LineWidth',1)
scatter(A,B,64,'b','filled','MarkerEdgeColor','w')
hold off
xlabel('CEA Fos')
ylabel('BST Fos')
text(0.05,0.95,['Drinking:  r = ',num2str(r,3),', p = ',num2str(p,2),char(10),...
    'Malaise:   r = ',num2str(rb,3),', p = ',num2str(pb,2),char(10),...
    'Retrieval: r = ',num2str(rc,3),', p = ',num2str(pc,2)],'Units','Normalized','FontSize',8,'VerticalAlignment','top')
title('CEA vs. BST')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')


A = counts_norm_CEA;
B = counts_norm(:,find(output.regions.name=='Agranular insular area, posterior part'));
C = counts_norm(:,find(output.regions.name=='Bed nuclei of the stria terminalis'));

tbl = table(output.GLMMinput.phase',A,B,C,'VariableNames',{'TimePoint','CEA','AIp','BST'});
writetable(tbl,['Z:\Chris\matlab\cz\cta-source-data\Fig-ED4d.csv'])


saveas(gcf,['plots-png/',path.file(63:end-4),'/plot-8-correlations'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/plot-8-correlations'],'epsc')
%% Figure: CEA bar plots
for i = 202
    k=0;
    file1 = [GLMMdir,num2str(k,'%04.f'),'/GLMM_output_region_',num2str(i,'%04.f'),'_coefficients.csv'];
    file2 = [GLMMdir,num2str(k,'%04.f'),'/GLMM_output_region_',num2str(i,'%04.f'),'_reducedmodels.csv'];
    file3 = [GLMMdir,num2str(k,'%04.f'),'/GLMM_output_region_',num2str(i,'%04.f'),'_marginaleffects.csv'];
    file4 = [GLMMdir,num2str(k,'%04.f'),'/GLMM_output_region_',num2str(i,'%04.f'),'_CGRP_coefficients.csv'];
    file5 = [GLMMdir,num2str(k,'%04.f'),'/GLMM_output_region_',num2str(i,'%04.f'),'_CGRP_marginaleffects.csv'];
    temp = readmatrix(file1,'Delimiter',' ');
    t = [];
    t(1,:) = rmmissing(temp(1,:));
    t(2,:) = rmmissing(temp(2,:));
    t(3,:) = rmmissing(temp(3,:));
    t(4,:) = rmmissing(temp(4,:));
    t(5,:) = rmmissing(temp(5,:));
    t(6,:) = rmmissing(temp(6,:));
    t(7,:) = rmmissing(temp(7,:));
    results.Eq2.coefficients.estimates(i,:) = t(:,1);
    results.Eq2.coefficients.std_errors(i,:) = t(:,2);
    
    temp = readmatrix(file2);
%     results.CEA.Eq2.modelstats.BIC(i,1) = temp(1,3);
%     results.CEA.Eq2.modelstats.LL(i,1) = temp(2,3);
    results.CEA.Eq2.modelstats.Chi2stat(i,1) = temp(8,3);
    results.CEA.Eq2.modelstats.pvalue_raw(i,1) = temp(3,3);
    results.CEA.Eq2.coefficients.Chi2stat(i,1:4) = temp([9,10,11,12],3);
    results.CEA.Eq2.coefficients.pvalues_raw(i,1:4) = temp([4,5,6,7],3);
    results.CEA.Eq2.coefficients.pvalues_corrected(i,1:4) = multicmp(results.Eq2.coefficients.pvalues_raw(i,:),'up',FWER);
    
    temp = readmatrix(file3,'Delimiter',' ','Range','A1:ZZZ10');
    results.CEA.Eq2.flavor.estimates(i,1:3) = rmmissing(temp(1,:));
    results.CEA.Eq2.flavor.std_errors(i,1:3) = rmmissing(temp(2,:));
    results.CEA.Eq2.flavor.pvalues_raw(i,rmmissing(temp(2,:))~=0) = rmmissing(temp(3,:));
    results.CEA.Eq2.flavor.pvalues_raw(i,rmmissing(temp(2,:))==0) = NaN;
    results.CEA.Eq2.flavor.Zstat(i,1:3) = results.Eq2.flavor.estimates(i,1:3)./results.Eq2.flavor.std_errors(i,1:3);
    results.CEA.Eq2.flavor.pvalues_corrected(i,1:3) = multicmp(results.Eq2.flavor.pvalues_raw(i,:),'up',FWER);
%     results.marginaleffects.phase.estimates(i,1:3) = rmmissing(temp(4,:));
%     results.marginaleffects.phase.std_errors(i,1:3) = rmmissing(temp(5,:));
%     results.marginaleffects.phase.pvalues_raw(i,1:3) = rmmissing(temp(6,:));
%     results.marginaleffects.phase.pvalues_corrected(i,1:3) = multicmp(results.marginaleffects.phase.pvalues_raw(i,:),'up',FWER);
    
    temp = readmatrix(file4,'Delimiter',' ');
    t = [];
    t(1,:) = rmmissing(temp(1,:));
    t(2,:) = rmmissing(temp(2,:));
    t(3,:) = rmmissing(temp(3,:));
    t(4,:) = rmmissing(temp(4,:));
    t(5,:) = rmmissing(temp(5,:));
    results.CEA.Eq3.coefficients.estimates(i,:) = t(:,1);
    results.CEA.Eq3.coefficients.std_errors(i,:) = t(:,2);
    results.CEA.Eq3.coefficients.Zstat(i,:) = t(:,3);
    results.CEA.Eq3.coefficients.pvalues_raw(i,:) = t(:,4);
    temp = readmatrix(file5,'Delimiter',' ','Range','A1:ZZZ10');
    results.CEA.Eq4.flavor.estimates(i,:) = rmmissing(temp(1,:));
    results.CEA.Eq4.flavor.std_errors(i,:) = rmmissing(temp(2,:));
    results.CEA.Eq4.flavor.pvalues_raw(i,rmmissing(temp(2,:))~=0) = rmmissing(temp(3,:));
    results.CEA.Eq4.flavor.pvalues_raw(i,rmmissing(temp(2,:))==0) = NaN;
    results.CEA.Eq4.flavor.Zstat(i,:) = results.Eq4.flavor.estimates(i,:)./results.Eq4.flavor.std_errors(i,:);
%     results.Eq4.flavor.pvalues_corrected(i,:) = multicmp(results.Eq4.flavor.pvalues_raw(i,:),'up',FWER);
%     results.Eq4.phase.estimates(i,:) = rmmissing(temp(4,:));
%     results.Eq4.phase.std_errors(i,:) = rmmissing(temp(5,:));
%     if size(temp,1)==6
%         results.Eq4.phase.pvalues_raw(i,:) = rmmissing(temp(6,:));
%     else
%         results.Eq4.phase.pvalues_raw(i,1:5) = rmmissing(temp(6,:));
%         results.Eq4.phase.pvalues_raw(i,6) = rmmissing(temp(7,:));
%     end
%     results.Eq4.phase.pvalues_corrected(i,:) = multicmp(results.Eq4.phase.pvalues_raw(i,:),'up',FWER);
end
%%
figure('Position', get(0,'Screensize'))

sizes = readmatrix([GLMMdir,'region_sizes.csv']);

subplot(1,3,1)
axis square
p_label = [];
hold on
counts_norm = [readmatrix([GLMMdir,'0000/counts_region_0202.csv'])./output.GLMMinput.offset./sizes(202)]*100;
phaselist = {'Drinking','LiCl','Retrieval','PBNCGRP'};

PH = cell(0,0);
FL = cell(0,0);
out = [];
for j = 1:3%:length(phaselist)

    
    PHASEidx = find(cellfun(@(x) isequal(phaselist{j},x),output.GLMMinput.phase)&cellfun(@(x) isequal('Novel',x),output.GLMMinput.flavor));
    x = counts_norm(PHASEidx);
    xx = prctile(x,[10 25 50 75 90]);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05+.125,x,64,'filled','MarkerFaceColor',[1 .75 .75],'MarkerEdgeColor','w')
    plot([-0.125 0.125]+j+.125,[xx(3) xx(3)],'r','LineWidth',2)
    plot([j j]+.125,[xx(2) xx(4)],'r','LineWidth',2)
    plot([j j]+.125,[xx(1) xx(5)],'r','LineWidth',2)
    FL(end+1:end+length(x)) = {'Novel'};
    PH(end+1:end+length(x)) = phaselist(j);
    out(end+1:end+length(x)) = x;

    PHASEidx = find(cellfun(@(x) isequal(phaselist{j},x),output.GLMMinput.phase)&cellfun(@(x) isequal('Familiar',x),output.GLMMinput.flavor));
    x = counts_norm(PHASEidx);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05-.125,x,64,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','w')
    plot([-0.125 0.125]+j-.125,[xx(3) xx(3)],'k','LineWidth',2)
    plot([j j]-.125,[xx(2) xx(4)],'k','LineWidth',2)
    plot([j j]-.125,[xx(1) xx(5)],'k','LineWidth',2)
    FL(end+1:end+length(x)) = {'Familiar'};
    PH(end+1:end+length(x)) = phaselist(j);
    out(end+1:end+length(x)) = x;
    
    %ylim([0 1])
    if j<4
        p_label = [p_label,'P = ',num2str(results.Eq2.flavor.pvalues_corrected(202,j),2),char(10)];
    else
        p_label = [p_label,char(10)];
    end
end
text(0.1,0.95,p_label,'Units','Normalized','FontSize',8,'VerticalAlignment','top')
xticks(1:length(phaselist))
xticklabels(phaselist)
% ylim([0 2000])
% yticks(0:500:2000)
xlim([0.5 length(phaselist)+.5])
ylabel('Fos^+ Cells (% of Total) per mm^3')
title('With Shell, % of Total')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

tbl = table(FL',PH',out','VariableNames',{'Flavor','Timepoint','CEA Fos'});
writetable(tbl,['Z:\Chris\matlab\cz\cta-source-data\Fig-1i.csv'])

subplot(1,3,2)
p_label = [];
hold on

counts2 = readmatrix([GLMMdir,'0000/counts_region_0202.csv']);
idx_drinking = find(cellfun(@(x) isequal('Drinking',x),output.GLMMinput.phase));
idx_retrieval = find(cellfun(@(x) isequal('Retrieval',x),output.GLMMinput.phase));
idx_licl = find(cellfun(@(x) isequal('LiCl',x),output.GLMMinput.phase));
idx_cgrp = find(cellfun(@(x) isequal('PBNCGRP',x),output.GLMMinput.phase));
idx_nov = find(cellfun(@(x) isequal('Novel',x),output.GLMMinput.flavor));
idx_fam = find(cellfun(@(x) isequal('Familiar',x),output.GLMMinput.flavor));
counts_norm = nan(size(counts2));
total = output.GLMMinput.offset./output.GLMMinput.pbn;

idx = sort([idx_nov,idx_fam]);
counts_norm(intersect(idx_drinking,idx),:) = counts2(intersect(idx_drinking,idx),:)./output.GLMMinput.pbn(intersect(idx_drinking,idx));
counts_norm(intersect(idx_drinking,idx),:) = counts_norm(intersect(idx_drinking,idx),:)./mean(total(intersect(idx_drinking,idx)));

idx = sort([idx_nov,idx_fam]);
counts_norm(intersect(idx_licl,idx),:) = counts2(intersect(idx_licl,idx),:)./output.GLMMinput.pbn(intersect(idx_licl,idx));
counts_norm(intersect(idx_licl,idx),:) = counts_norm(intersect(idx_licl,idx),:)./mean(total(intersect(idx_licl,idx)));

idx = sort([idx_nov,idx_fam]);
counts_norm(intersect(idx_retrieval,idx),:) = counts2(intersect(idx_retrieval,idx),:)./output.GLMMinput.pbn(intersect(idx_retrieval,idx));
counts_norm(intersect(idx_retrieval,idx),:) = counts_norm(intersect(idx_retrieval,idx),:)./mean(total(intersect(idx_retrieval,idx)));

idx = sort([idx_nov,idx_fam]);
counts_norm(intersect(idx_cgrp,idx),:) = counts2(intersect(idx_cgrp,idx),:)./output.GLMMinput.pbn(intersect(idx_cgrp,idx));
counts_norm(intersect(idx_cgrp,idx),:) = counts_norm(intersect(idx_cgrp,idx),:)./mean(total(intersect(idx_cgrp,idx)));

counts_norm = [counts_norm./sizes(202)']*100;

axis square
PH = cell(0,0);
FL = cell(0,0);
out = [];
for j = 4%1:length(phaselist)
    PHASEidx = find(cellfun(@(x) isequal(phaselist{j},x),output.GLMMinput.phase)&cellfun(@(x) isequal('Familiar',x),output.GLMMinput.flavor));
    x = counts_norm(PHASEidx);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05-.125,x,64,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','w')
    plot([-0.125 0.125]+j-.125,[xx(3) xx(3)],'k','LineWidth',2)
    plot([j j]-.125,[xx(2) xx(4)],'k','LineWidth',2)
    plot([j j]-.125,[xx(1) xx(5)],'k','LineWidth',2)
    
        plot([j j]-.125,[xx(1) xx(5)],'k','LineWidth',2)
    FL(end+1:end+length(x)) = {'Familiar'};
    PH(end+1:end+length(x)) = phaselist(j);
    out(end+1:end+length(x)) = x;
    
    PHASEidx = find(cellfun(@(x) isequal(phaselist{j},x),output.GLMMinput.phase)&cellfun(@(x) isequal('Novel',x),output.GLMMinput.flavor));
    x = counts_norm(PHASEidx);
    xx = prctile(x,[10 25 50 75 90]);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05+.125,x,64,'filled','MarkerFaceColor',[1 .75 .75],'MarkerEdgeColor','w')
    plot([-0.125 0.125]+j+.125,[xx(3) xx(3)],'r','LineWidth',2)
    plot([j j]+.125,[xx(2) xx(4)],'r','LineWidth',2)
    plot([j j]+.125,[xx(1) xx(5)],'r','LineWidth',2)

        FL(end+1:end+length(x)) = {'Novel'};
    PH(end+1:end+length(x)) = phaselist(j);
    out(end+1:end+length(x)) = x;

    %ylim([0 1])
    if j<4
        p_label = [p_label,char(10)];
    else
        p_label = [p_label,'P = ',num2str(results.Eq4.flavor.pvalues_raw(202,2),2),char(10)];
    end
end
text(0.1,0.95,p_label,'Units','Normalized','FontSize',8,'VerticalAlignment','top')
xticks(1:length(phaselist))
xticklabels(phaselist)
% ylim([0 2000])
% yticks(0:500:2000)
xlim([0.5 length(phaselist)+.5])
ylabel('Fos^+ Cells (% of PB) per mm^3')
title('With Shell, % of PB')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

tbl = table(FL',PH',out','VariableNames',{'Flavor','Timepoint','CEA Fos'});
writetable(tbl,['Z:\Chris\matlab\cz\cta-source-data\Fig-2h.csv'])

subplot(1,3,3)
hold on
axis square
% plot([0 N_clusters*4],[0 0],'k','LineWidth',1)
y1 = [];
idx = cluster(tree,'maxclust',N_clusters);
idx = idx(outperm);
list = unique(idx,'stable');
amygdala_cluster = find(list==idx(find(RegionLibrary.reduced{significant(leafOrder),2}=="Central amygdalar nucleus, medial part")));
out = [];
for i = 1%:length(list)
    
    x = squeeze(CorrStruct(amygdala_cluster,amygdala_cluster,1,:)); x(isnan(x)) = [];
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    xval = (i-1)*4+1;
        scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,64,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','w')
        plot([-0.25 0.25]+xval,[xx(3) xx(3)],'k','LineWidth',1)
        plot([xval xval],[xx(2) xx(4)],'k','LineWidth',1)
    y1(i,1) = mean(x);
        [p,~,stats]=signrank(x)
        
        out(:,1) = x;
    
    x = squeeze(CorrStruct(amygdala_cluster,amygdala_cluster,2,:)); x(isnan(x)) = [];
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    xval = (i-1)*4+2;
        scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,64,'filled','MarkerFaceColor',[.85 .75 .75],'MarkerEdgeColor','w')
        plot([-0.25 0.25]+xval,[xx(3) xx(3)],'Color',[.5 0 0],'LineWidth',1)
        plot([xval xval],[xx(2) xx(4)],'Color',[.5 0 0],'LineWidth',1)
    y1(i,2) = mean(x);
        [p,~,stats]=signrank(x)
    
        out(:,2) = x;

        x = squeeze(CorrStruct(amygdala_cluster,amygdala_cluster,3,:)); x(isnan(x)) = [];
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    xval = (i-1)*4+3;
        scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,64,'filled','MarkerFaceColor',[1 .75 .75],'MarkerEdgeColor','w')
        plot([-0.25 0.25]+xval,[xx(3) xx(3)],'Color',[1 0 0],'LineWidth',1)
        plot([xval xval],[xx(2) xx(4)],'Color',[1 0 0],'LineWidth',1)
    y1(i,3) = mean(x);
        [p,~,stats]=signrank(x)

                out(:,3) = x;

end
xlim([0.5 3.4])
ylim([0 1])
xticks([1:3])
yticks(-1:.5:1)
xticklabels({'Meal','Malaise','Retrieval'})
xtickangle(90)
ylabel('Within-cluster correlation')
title('Fos cluster 1')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
%

tbl = table(out(:,1),out(:,2),out(:,3),'VariableNames',{'Consume','Malaise','Retrieval'});
writetable(tbl,['Z:\Chris\matlab\cz\cta-source-data\Fig-ED4b.csv'])

saveas(gcf,['plots-png/',path.file(63:end-4),'/plot-9-CEA'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/plot-9-CEA'],'epsc')
%% Figure: Phase summary by cluster

figure('Position', get(0,'Screensize'))

idx = cluster(tree,'maxclust',N_clusters);
idx = idx(outperm);

list = unique(idx,'stable');
data = struct;
a = [];
for i = 1:length(list)
ax=subplot(2,5,i);
axis square
hold on

    x = treedata(leafOrder(find(idx==list(length(list)-i+1))),1);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    xval = 1;
    scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,16,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','w')
    plot([-0.5 0.5]+xval,[xx(3) xx(3)],'k','LineWidth',1)
    plot([xval xval],[xx(2) xx(4)],'k','LineWidth',1)
    a(i,1) = signrank(x);

        x = treedata(leafOrder(find(idx==list(length(list)-i+1))),2);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    xval = 2;
    scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,16,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','w')
    plot([-0.5 0.5]+xval,[xx(3) xx(3)],'k','LineWidth',1)
    plot([xval xval],[xx(2) xx(4)],'k','LineWidth',1)
    a(i,2) = signrank(x);

        x = treedata(leafOrder(find(idx==list(length(list)-i+1))),3);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    xval = 3;
    scatter(ones(1,length(x))*xval+rand(1,length(x))*.1-.05,x,16,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','w')
    plot([-0.5 0.5]+xval,[xx(3) xx(3)],'k','LineWidth',1)
    plot([xval xval],[xx(2) xx(4)],'k','LineWidth',1)

    a(i,3) = signrank(x);
xlim([0 3+1])
ylim([-5 5])
xticks([1:3])
yticks(-5:2.5:5)
xticklabels(num2str([1:N_clusters]'))
xtickangle(0)
xlabel('Timepoint')
ylabel('Novel–Familiar (z)')
set(gca,'FontSize',8,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
a = multicmp(a(i,:),'up',FWER);
disp([i size(x,1) a])
end
saveas(gcf,['plots-png/',path.file(63:end-4),'/cluster-summary'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/cluster-summary'],'epsc')
%%
return
%% 3D brain renderings

idx_fam = find(contains(output.GLMMoutput.flavorME.significant(:,1),'*')&output.GLMMoutput.flavorME.estimates(:,1)<0);
idx_nov = find(contains(output.GLMMoutput.flavorME.significant(:,1),'*')&output.GLMMoutput.flavorME.estimates(:,1)>0);
amygdala_cluster = find(list==idx(find(RegionLibrary.reduced{significant(leafOrder),2}=="Central amygdalar nucleus, medial part")));
load('Z:\Chris\matlab\cz\neuropixels-utils\Allen_CCFv3_Fos.mat','atlas')

figure('Position', get(0, 'Screensize'))
hold on
amygdala_cluster = find(list==idx(find(RegionLibrary.reduced{significant(leafOrder),2}=="Central amygdalar nucleus, medial part")));
AMYG = significant(leafOrder(idx==list(amygdala_cluster)));
AMYG = sort(AMYG); AMYG = AMYG(1:end-1);
for i = 1:length(AMYG)
    t = atlas==(RegionLibrary.reduced.index(AMYG(i))+1);
    %t(1:size(t,1)/2,:,:) = 0;
    p = patch(isosurface(t));
    p.FaceColor = [229 45 38]/255;
    p.FaceAlpha = 0.4;
    p.LineStyle = 'none';
end
p = patch(isosurface(atlas>0));
p.FaceAlpha = 0.05;
p.FaceColor = [44 31 22]/255;
p.LineStyle = 'none';

t = atlas==(903+1);
p = patch(isosurface(t));
p.FaceColor = [44 31 22]/255;
p.FaceAlpha = 0.2;
p.LineStyle = 'none';
t = atlas==(950+1);
p = patch(isosurface(t));
p.FaceColor = [44 31 22]/255;
p.FaceAlpha = 0.2;
p.LineStyle = 'none';

h = get(gca,'DataAspectRatio');
set(gca,'DataAspectRatio',[1 1 h(3)],'CameraViewAngleMode','Manual');
set(gca,'Ydir','reverse','Zdir','reverse')
axis off
xlim('auto')
ylim('auto')
zlim('auto')
view(-20,20)
hold off
% set(gcf,'renderer','Painters')
exportgraphics(gcf,['plots-png/',path.file(63:end-4),'/plot-9-amygdala-3d.png'],'Resolution',300)

figure('Position', get(0, 'Screensize'))
hold on
for i = 1:length(idx_nov)
    t = atlas==(RegionLibrary.reduced.index(idx_nov(i))+1);
    %t(1:size(t,1)/2,:,:) = 0;
    p = patch(isosurface(t));
    p.FaceColor = [229 45 38]/255;
    p.FaceAlpha = 0.4;
    p.LineStyle = 'none';
end
p = patch(isosurface(atlas>0));
p.FaceAlpha = 0.05;
p.FaceColor = [44 31 22]/255;
p.LineStyle = 'none';

h = get(gca,'DataAspectRatio');
set(gca,'DataAspectRatio',[1 1 h(3)],'CameraViewAngleMode','Manual');
set(gca,'Ydir','reverse','Zdir','reverse')
axis off
xlim('auto')
ylim('auto')
zlim('auto')
view(-20,20)
hold off
% set(gcf,'renderer','Painters')
exportgraphics(gcf,['plots-png/',path.file(63:end-4),'/plot-10-novel-3d.png'],'Resolution',300)

figure('Position', get(0, 'Screensize'))
hold on
for i = 1:length(idx_nov)
    t = atlas==(RegionLibrary.reduced.index(idx_nov(i))+1);
    %t(1:size(t,1)/2,:,:) = 0;
    p = patch(isosurface(t));
    p.FaceColor = [229 45 38]/255;
    p.FaceAlpha = 0.4;
    p.LineStyle = 'none';
end
for i = 1:length(idx_fam)
    t = atlas==(RegionLibrary.reduced.index(idx_fam(i))+1);
    %t(1:size(t,1)/2,:,:) = 0;
    p = patch(isosurface(t));
    p.FaceColor = [55 136 193]/255;
    p.FaceAlpha = 0.4;
    p.LineStyle = 'none';
end
p = patch(isosurface(atlas>0));
p.FaceAlpha = 0.05;
p.FaceColor = [44 31 22]/255;
p.LineStyle = 'none';

h = get(gca,'DataAspectRatio');
set(gca,'DataAspectRatio',[1 1 h(3)],'CameraViewAngleMode','Manual');
set(gca,'Ydir','reverse','Zdir','reverse')
axis off
xlim('auto')
ylim('auto')
zlim('auto')
view(-20,20)
hold off
% set(gcf,'renderer','Painters')
exportgraphics(gcf,['plots-png/',path.file(63:end-4),'/plot-11-novel-familiar-3d.png'],'Resolution',300)

figure('Position', get(0, 'Screensize'))
hold on
for i = 1:length(idx_fam)
    t = atlas==(RegionLibrary.reduced.index(idx_fam(i))+1);
    %t(1:size(t,1)/2,:,:) = 0;
    p = patch(isosurface(t));
    p.FaceColor = [55 136 193]/255;
    p.FaceAlpha = 0.4;
    p.LineStyle = 'none';
end
p = patch(isosurface(atlas>0));
p.FaceAlpha = 0.05;
p.FaceColor = [44 31 22]/255;
p.LineStyle = 'none';

h = get(gca,'DataAspectRatio');
set(gca,'DataAspectRatio',[1 1 h(3)],'CameraViewAngleMode','Manual');
set(gca,'Ydir','reverse','Zdir','reverse')
axis off
xlim('auto')
ylim('auto')
zlim('auto')
view(-20,20)
hold off
% set(gcf,'renderer','Painters')
exportgraphics(gcf,['plots-png/',path.file(63:end-4),'/plot-11-familiar-3d.png'],'Resolution',300)
%%
load('Z:\Chris\matlab\cz\neuropixels-utils\Allen_CCFv3_Fos.mat','atlas')
idx = cluster(tree,'maxclust',N_clusters);
idx = idx(outperm);
list = unique(idx,'stable');

for i = 1:10
    figure('Position', get(0, 'Screensize'))
    hold on
    amygdala_cluster = i;
    AMYG = significant(leafOrder(idx==list(amygdala_cluster)));
    for j = 1:length(AMYG)
        t = atlas==(RegionLibrary.reduced.index(AMYG(j))+1);
        %t(1:size(t,1)/2,:,:) = 0;
        p = patch(isosurface(t));
        p.FaceColor = [44 31 22]/255;
        p.FaceAlpha = 0.4;
        p.LineStyle = 'none';
    end
    p = patch(isosurface(atlas>0));
    p.FaceAlpha = 0.05;
    p.FaceColor = [44 31 22]/255;
    p.LineStyle = 'none';
    
    h = get(gca,'DataAspectRatio');
    set(gca,'DataAspectRatio',[1 1 h(3)],'CameraViewAngleMode','Manual');
    set(gca,'Ydir','reverse','Zdir','reverse')
    axis off
    xlim('auto')
    ylim('auto')
    zlim('auto')
    view(-20,20)
    hold off
    % set(gcf,'renderer','Painters')
    exportgraphics(gcf,['plots-png/',path.file(63:end-4),'/cluster_',num2str(i),'-3d.png'],'Resolution',300)
end
%%
addpath(genpath('Z:\Chris\matlab\numpy-matlab\'));

fpath = 'Z:\Chris\data\clearmap2\kde-visualization\cta-final\';

fname = 'kde_20um_cta_drinking_novel_weighted_decr75pc.npy';
KDE.Drinking.Novel = readNPY([fpath,fname]);
fname = 'kde_20um_cta_drinking_familiar_weighted_decr75pc.npy';
KDE.Drinking.Familiar = readNPY([fpath,fname]);
fname = 'kde_20um_cta_licl_novel_weighted_decr75pc.npy';
KDE.LiCl.Novel = readNPY([fpath,fname]);
fname = 'kde_20um_cta_licl_familiar_weighted_decr75pc.npy';
KDE.LiCl.Familiar = readNPY([fpath,fname]);
fname = 'cgrp_pct_of_pbn\kde_20um_cta_cgrp_novel_weighted_decr75pc - FullWeight.npy';
KDE.CGRP.Novel = readNPY([fpath,fname]);
fname = 'cgrp_pct_of_pbn\kde_20um_cta_cgrp_familiar_weighted_decr75pc - FullWeight.npy';
KDE.CGRP.Familiar = readNPY([fpath,fname]);
fname = 'kde_20um_cta_retrieval_novel_weighted_decr75pc.npy';
KDE.Retrieval.Novel = readNPY([fpath,fname]);
fname = 'kde_20um_cta_retrieval_familiar_weighted_decr75pc.npy';
KDE.Retrieval.Familiar = readNPY([fpath,fname]);

fpath = 'Z:\Chris\data\clearmap2\utilities\allen-atlas-cz\';
fname = 'annotation_2017_25um_sagittal_16bit_hierarch_labels_fillmissing_cz_v2.tif';
filename = [fpath,fname];
tstack  = Tiff(filename);
[I,J] = size(tstack.read());
K = length(imfinfo(filename));
data = zeros(K,I,J);
data(1,:,:)  = tstack.read();
for n = 2:K
    tstack.nextDirectory()
    data(n,:,:) = tstack.read();
end
atlas = data;

atlasmask  = atlas>=1028 | atlas<=1;
KDE.Drinking.Novel(atlasmask) = NaN;
KDE.Drinking.Familiar(atlasmask) = NaN;
KDE.LiCl.Novel(atlasmask) = NaN;
KDE.LiCl.Familiar(atlasmask) = NaN;
KDE.CGRP.Novel(atlasmask) = NaN;
KDE.CGRP.Familiar(atlasmask) = NaN;
KDE.Retrieval.Novel(atlasmask) = NaN;
KDE.Retrieval.Familiar(atlasmask) = NaN;

method = 'diff';
cmap1 = cbrewer('seq','Greys',1000,'spline');
cmap2 = flipud(cbrewer('div','RdBu',1000,'spline')); cmap2(cmap2<0) = 0;
cmap3 = hsv(size(RegionLibrary.reduced,1)); cmap3 = cmap3(randperm(RandStream('mt19937ar','seed',12345),size(cmap3,1)),:)*0.5+.5;
lims1 = [0 1];
if isequal(method,'log')
    lims2 = [-2 2];
    lambda = 1E-2;
elseif isequal(method,'diff')
    lims2 = [-0.5 0.5];
else
    disp('Invalid comparison method.')
    return
end
%%
close all

bregma = 239;
planes = [-3.5:.5:6]/.025+bregma;

figure('Position', get(0,'Screensize'))
for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,10,pl)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    
    data = [KDE.Drinking.Novel(:,plane,:),KDE.Drinking.Familiar(:,plane,:)];
    if isequal(method,'log')
        data = log(squeeze([data(:,1,:)+lambda])./squeeze([data(:,2,:)+lambda]));
    elseif isequal(method,'diff')
        data = squeeze(data(:,1,:))-squeeze(data(:,2,:));
    end
    colormap(gcf,cmap2)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims2),'MinColorValue',min(lims2),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',0.25)
            end
    end
    for i = 1:length(B2a)
        fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',0.5)
    end
    for i = 1:length(B2b)
        fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B2d)
        fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B)
        plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',0.5)
    end
    title([num2str(plane)])
    set(gca,'FontSize',8)
    xlim([.5 228.5])
    axis off
    hold off

end
saveas(gcf,['plots-png/',path.file(63:end-4),'/drinking-heatmap-novelfamiliar'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/drinking-heatmap-novelfamiliar'],'epsc')

figure('Position', get(0,'Screensize'))
for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,10,pl)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    
    data = [KDE.LiCl.Novel(:,plane,:),KDE.LiCl.Familiar(:,plane,:)];
    if isequal(method,'log')
        data = log(squeeze([data(:,1,:)+lambda])./squeeze([data(:,2,:)+lambda]));
    elseif isequal(method,'diff')
        data = squeeze(data(:,1,:))-squeeze(data(:,2,:));
    end
    colormap(gcf,cmap2)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims2),'MinColorValue',min(lims2),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',0.25)
            end
    end
    for i = 1:length(B2a)
        fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',0.5)
    end
    for i = 1:length(B2b)
        fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B2d)
        fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B)
        plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',0.5)
    end
    title([num2str(plane)])
    set(gca,'FontSize',8)
    xlim([.5 228.5])
    axis off
    hold off

end
saveas(gcf,['plots-png/',path.file(63:end-4),'/malaise-heatmap-novelfamiliar'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/malaise-heatmap-novelfamiliar'],'epsc')

figure('Position', get(0,'Screensize'))
for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,10,pl)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    
    data = [KDE.Retrieval.Novel(:,plane,:),KDE.Retrieval.Familiar(:,plane,:)];
    if isequal(method,'log')
        data = log(squeeze([data(:,1,:)+lambda])./squeeze([data(:,2,:)+lambda]));
    elseif isequal(method,'diff')
        data = squeeze(data(:,1,:))-squeeze(data(:,2,:));
    end
    colormap(gcf,cmap2)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims2),'MinColorValue',min(lims2),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',0.25)
            end
    end
    for i = 1:length(B2a)
        fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',0.5)
    end
    for i = 1:length(B2b)
        fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B2d)
        fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B)
        plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',0.5)
    end
    title([num2str(plane)])
    set(gca,'FontSize',8)
    xlim([.5 228.5])
    axis off
    hold off

end
saveas(gcf,['plots-png/',path.file(63:end-4),'/retrieval-heatmap-novelfamiliar'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/retrieval-heatmap-novelfamiliar'],'epsc')
%%

bregma = 239;
planes = [-3.5:.5:6]/.025+bregma;

figure('Position', get(0,'Screensize'))
for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,10,pl)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    
    data = [KDE.Drinking.Novel(:,plane,:),KDE.Drinking.Familiar(:,plane,:)];
    data = squeeze(mean(data,2));
    colormap(gcf,cmap1)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims1),'MinColorValue',min(lims1),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',0.25)
            end
    end
    for i = 1:length(B2a)
        fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',0.5)
    end
    for i = 1:length(B2b)
        fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B2d)
        fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B)
        plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',0.5)
    end
    title([num2str(plane)])
    set(gca,'FontSize',8)
    xlim([.5 228.5])
    axis off
    hold off

end
saveas(gcf,['plots-png/',path.file(63:end-4),'/drinking-heatmap-average'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/drinking-heatmap-average'],'epsc')

figure('Position', get(0,'Screensize'))
for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,10,pl)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    
    data = [KDE.LiCl.Novel(:,plane,:),KDE.LiCl.Familiar(:,plane,:)];
    data = squeeze(mean(data,2));
    colormap(gcf,cmap1)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims1),'MinColorValue',min(lims1),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',0.25)
            end
    end
    for i = 1:length(B2a)
        fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',0.5)
    end
    for i = 1:length(B2b)
        fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B2d)
        fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B)
        plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',0.5)
    end
    title([num2str(plane)])
    set(gca,'FontSize',8)
    xlim([.5 228.5])
    axis off
    hold off

end
saveas(gcf,['plots-png/',path.file(63:end-4),'/malaise-heatmap-average'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/malaise-heatmap-average'],'epsc')

figure('Position', get(0,'Screensize'))
for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,10,pl)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    
    data = [KDE.Retrieval.Novel(:,plane,:),KDE.Retrieval.Familiar(:,plane,:)];
    data = squeeze(mean(data,2));
    colormap(gcf,cmap1)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims1),'MinColorValue',min(lims1),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',0.25)
            end
    end
    for i = 1:length(B2a)
        fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',0.5)
    end
    for i = 1:length(B2b)
        fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B2d)
        fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B)
        plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',0.5)
    end
    title([num2str(plane)])
    set(gca,'FontSize',8)
    xlim([.5 228.5])
    axis off
    hold off

end
saveas(gcf,['plots-png/',path.file(63:end-4),'/retrieval-heatmap-average'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/retrieval-heatmap-average'],'epsc')
%%
close all

bregma = 239;
planes = [-3.5:.5:6]/.025+bregma;

figure('Position', get(0,'Screensize'))
for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,10,pl)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    
    data = [KDE.CGRP.Novel(:,plane,:),KDE.CGRP.Familiar(:,plane,:)];
    if isequal(method,'log')
        data = log(squeeze([data(:,1,:)+lambda])./squeeze([data(:,2,:)+lambda]));
    elseif isequal(method,'diff')
        data = squeeze(data(:,1,:))-squeeze(data(:,2,:));
    end
    colormap(gcf,cmap2)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',100,'MinColorValue',-100,'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',0.25)
            end
    end
    for i = 1:length(B2a)
        fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',0.5)
    end
    for i = 1:length(B2b)
        fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B2d)
        fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B)
        plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',0.5)
    end
    title([num2str(plane)])
    set(gca,'FontSize',8)
    xlim([.5 228.5])
    axis off
    hold off

end
saveas(gcf,['plots-png/',path.file(63:end-4),'/cgrp-heatmap-novelfamiliar'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/cgrp-heatmap-novelfamiliar'],'epsc')

figure('Position', get(0,'Screensize'))
for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,10,pl)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    
    data = [KDE.CGRP.Novel(:,plane,:),KDE.CGRP.Familiar(:,plane,:)];
    data = squeeze(mean(data,2));
    colormap(gcf,cmap1)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',200,'MinColorValue',0,'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',0.25)
            end
    end
    for i = 1:length(B2a)
        fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',0.5)
    end
    for i = 1:length(B2b)
        fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B2d)
        fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B)
        plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',0.5)
    end
    title([num2str(plane)])
    set(gca,'FontSize',8)
    xlim([.5 228.5])
    axis off
    hold off

end
saveas(gcf,['plots-png/',path.file(63:end-4),'/cgrp-heatmap-average'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/cgrp-heatmap-average'],'epsc')
%%

bregma = 217;
planes = [-1.5:.5:5]/.025+bregma;

lims1 = [0 1.5];
lims2 = [-0.5 0.5];

figure('Position', get(0,'Screensize'))
for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,7,pl)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    
    data = [KDE.Drinking.Novel(:,plane,:),KDE.Drinking.Familiar(:,plane,:)];
    if isequal(method,'log')
        data = log(squeeze([data(:,1,:)+lambda])./squeeze([data(:,2,:)+lambda]));
    elseif isequal(method,'diff')
        data = squeeze(data(:,1,:))-squeeze(data(:,2,:));
    end
    colormap(gcf,cmap2)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims2),'MinColorValue',min(lims2),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',0.25)
            end
    end
    for i = 1:length(B2a)
        fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',0.5)
    end
    for i = 1:length(B2b)
        fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B2d)
        fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B)
        plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',0.5)
    end
    title([num2str(plane)])
    set(gca,'FontSize',8)
    %xlim([.5 228.5])
    axis off
    hold off

end
saveas(gcf,['plots-png/',path.file(63:end-4),'/drinking-heatmap-novelfamiliar-talk'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/drinking-heatmap-novelfamiliar-talk'],'epsc')

figure('Position', get(0,'Screensize'))
for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,7,pl)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    
    data = [KDE.Drinking.Novel(:,plane,:),KDE.Drinking.Familiar(:,plane,:)];
    data = squeeze(mean(data,2));
    colormap(gcf,cmap1)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims1),'MinColorValue',min(lims1),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',0.25)
            end
    end
    for i = 1:length(B2a)
        fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',0.5)
    end
    for i = 1:length(B2b)
        fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B2d)
        fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B)
        plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',0.5)
    end
    title([num2str(plane)])
    set(gca,'FontSize',8)
    axis off
    hold off

end
saveas(gcf,['plots-png/',path.file(63:end-4),'/drinking-heatmap-average-talk'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/drinking-heatmap-average-talk'],'epsc')
%%
close all
figure('Position', get(0,'Screensize'))

idx_fam = find(contains(output.GLMMoutput.flavorME.significant(:,1),'*')&output.GLMMoutput.flavorME.estimates(:,1)<0);
idx_nov = find(contains(output.GLMMoutput.flavorME.significant(:,1),'*')&output.GLMMoutput.flavorME.estimates(:,1)>0);
idx_amygdala = amygdala_regions;

planes_fam = [135,180,330];
planes_nov = [210,265,325];

for pl = 1:length(planes_fam)
    plane = planes_fam(pl);
    subplot(2,3,pl)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    
    data = [KDE.Drinking.Novel(:,plane,:),KDE.Drinking.Familiar(:,plane,:)];
    if isequal(method,'log')
        data = log(squeeze([data(:,1,:)+lambda])./squeeze([data(:,2,:)+lambda]));
    elseif isequal(method,'diff')
        data = squeeze(data(:,1,:))-squeeze(data(:,2,:));
    end
    colormap(gcf,cmap2)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims2),'MinColorValue',min(lims2),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
        if ~ismember(i,idx_fam)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',0.25)
            end
        end
    end
    for i = 1:size(RegionLibrary.reduced,1)
        if ismember(i,idx_fam)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[55 135 192]/255,'LineWidth',0.5)
            end
        end
    end
    for i = 1:length(B2a)
        fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',0.5)
    end
    for i = 1:length(B2b)
        fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B2d)
        fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B)
        plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',0.5)
    end
    title(['Drinking: Plane ',num2str(plane)])
    c2 = colorbar('Location','westoutside','FontSize',8);
    c2.Position = c2.Position+[-.03 0 0 0];
    c2.Ticks = [-.5 0 .5];
    c2.Label.String = '\DeltaFos (% per mm^3)';
    set(gca,'FontSize',8)
    xlim([.5 228.5])
    axis off
    hold off
end

for pl = 1:length(planes_nov)
    plane = planes_nov(pl);
    subplot(2,3,pl+3)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    
    data = [KDE.Drinking.Novel(:,plane,:),KDE.Drinking.Familiar(:,plane,:)];
    if isequal(method,'log')
        data = log(squeeze([data(:,1,:)+lambda])./squeeze([data(:,2,:)+lambda]));
    elseif isequal(method,'diff')
        data = squeeze(data(:,1,:))-squeeze(data(:,2,:));
    end
    colormap(gcf,cmap2)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims2),'MinColorValue',min(lims2),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
        if ~ismember(i,idx_fam)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',0.25)
            end
        end
    end
    for i = 1:size(RegionLibrary.reduced,1)
        if ismember(i,idx_nov)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[228 45 38]/255,'LineWidth',0.5)
            end
        end
    end
    for i = 1:length(B2a)
        fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',0.5)
    end
    for i = 1:length(B2b)
        fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B2d)
        fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B)
        plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',0.5)
    end
    title(['Drinking: Plane ',num2str(plane)])
    c2 = colorbar('Location','westoutside','FontSize',8);
    c2.Position = c2.Position+[-.03 0 0 0];
    c2.Ticks = [-.5 0 .5];
    c2.Label.String = '\DeltaFos (% per mm^3)';
    set(gca,'FontSize',8)
    xlim([.5 228.5])
    axis off
    hold off
end
saveas(gcf,['plots-png/',path.file(63:end-4),'/heatmap-drinking-clusters'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/heatmap-drinking-clusters'],'epsc')


figure('Position', get(0,'Screensize'))
for pl = 1:2
    plane = planes_nov(pl+1);
    subplot(2,3,1+(pl-1)*3)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    
    data = [KDE.LiCl.Novel(:,plane,:),KDE.LiCl.Familiar(:,plane,:)];
    if isequal(method,'log')
        data = log(squeeze([data(:,1,:)+lambda])./squeeze([data(:,2,:)+lambda]));
    elseif isequal(method,'diff')
        data = squeeze(data(:,1,:))-squeeze(data(:,2,:));
    end
    colormap(gcf,cmap2)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims2),'MinColorValue',min(lims2),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
        if ~ismember(i,idx_fam)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',0.25)
            end
        end
    end
    for i = 1:size(RegionLibrary.reduced,1)
        if ismember(i,idx_amygdala)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[228 45 38]/255,'LineWidth',0.5)
            end
        end
    end
    for i = 1:length(B2a)
        fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',0.5)
    end
    for i = 1:length(B2b)
        fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B2d)
        fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B)
        plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',0.5)
    end
    title(['Malaise: Plane ',num2str(plane)])
    c2 = colorbar('Location','westoutside','FontSize',8);
    c2.Position = c2.Position+[-.03 0 0 0];
    c2.Ticks = [-.5 0 .5];
    c2.Label.String = '\DeltaFos (% per mm^3)';
    set(gca,'FontSize',8)
    xlim([.5 228.5])
    axis off
    hold off

    subplot(2,3,2+(pl-1)*3)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    data = [KDE.Retrieval.Novel(:,plane,:),KDE.Retrieval.Familiar(:,plane,:)];
    if isequal(method,'log')
        data = log(squeeze([data(:,1,:)+lambda])./squeeze([data(:,2,:)+lambda]));
    elseif isequal(method,'diff')
        data = squeeze(data(:,1,:))-squeeze(data(:,2,:));
    end
    colormap(gcf,cmap2)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims2),'MinColorValue',min(lims2),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
        if ~ismember(i,idx_fam)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',0.25)
            end
        end
    end
    for i = 1:size(RegionLibrary.reduced,1)
        if ismember(i,idx_amygdala)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[228 45 38]/255,'LineWidth',0.5)
            end
        end
    end
    for i = 1:length(B2a)
        fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',0.5)
    end
    for i = 1:length(B2b)
        fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B2d)
        fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B)
        plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',0.5)
    end
    title(['Retrieval: Plane ',num2str(plane)])
    c2 = colorbar('Location','westoutside','FontSize',8);
    c2.Position = c2.Position+[-.03 0 0 0];
    c2.Ticks = [-.5 0 .5];
    c2.Label.String = '\DeltaFos (% per mm^3)';
    set(gca,'FontSize',8)
    xlim([.5 228.5])
    axis off
    hold off

    subplot(2,3,3+(pl-1)*3)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));    data = [KDE.CGRP.Novel(:,plane,:),KDE.CGRP.Familiar(:,plane,:)];
    if isequal(method,'log')
        data = log(squeeze([data(:,1,:)+lambda])./squeeze([data(:,2,:)+lambda]));
    elseif isequal(method,'diff')
        data = squeeze(data(:,1,:))-squeeze(data(:,2,:));
    end
    colormap(gcf,cmap2)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',100,'MinColorValue',-100,'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
        if ~ismember(i,idx_fam)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',0.25)
            end
        end
    end
    for i = 1:size(RegionLibrary.reduced,1)
        if ismember(i,idx_amygdala)
            idx = RegionLibrary.reduced{i,1}+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[228 45 38]/255,'LineWidth',0.5)
            end
        end
    end
    for i = 1:length(B2a)
        fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',0.5)
    end
    for i = 1:length(B2b)
        fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B2d)
        fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',0.5)
    end
    for i = 1:length(B)
        plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',0.5)
    end
    title(['CGRP: Plane ',num2str(plane)])
    c2 = colorbar('Location','westoutside','FontSize',8);
    c2.Position = c2.Position+[-.03 0 0 0];
    c2.Ticks = [-.5 0 .5];
    c2.Label.String = '\DeltaFos (% per mm^3)';
    set(gca,'FontSize',8)
    xlim([.5 228.5])
    axis off
    hold off

end
saveas(gcf,['plots-png/',path.file(63:end-4),'/heatmap-amygdala-cluster'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/',path.file(63:end-4),'/heatmap-amygdala-cluster'],'epsc')