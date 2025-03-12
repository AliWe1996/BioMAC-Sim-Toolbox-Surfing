%======================================================================
%> @file +IMUTracking3D/script_IMU_PE.m
%> @brief Script to run simulations with IMU tracking and PE angles
%>
%> @author Alexander Weiss
%> @date December, 2024
%======================================================================
clear all;
close all;


% path2repo
workDirectory = (what('BioMAC-Sim-Toolbox-Surfing').path); % If not already in the BioMAC-Sim-Toolbox folder
% cd(workDirectory);
% addpath(genpath('src'));
% savepath;
%% Settings
% Get path of this script
filePath = fileparts(mfilename('fullpath'));
% Path to your repository
path2repo = [filePath filesep '..' filesep '..' filesep];
workDirectory = path2repo;

% Specifiy participant numbers
participants = [1:7];
% No of camera
camView = 1;
%Tracking data weight
weights = "1";
% Angle of wave inclination
slope = 20;
% Flag for initial guess if == 1, standing without contact dynamics, if i
% == 0 initGuess from previous simulation
makeInitGuess = 0;

% Participants and camera view
noPar = numel(participants);
viewCount = 1;

if makeInitGuess == 1
    for iPar = participants
        % Weights and flags
        WAcc = '0';
        WGyro = '0';
        WReg = '1e-2';
        doStanding = 1;
        nNodes = 100;
        WAngles = '10';
        WTrans = '10';
        WEffort = 1;

        %InitGuess for the actual surfing simulation which excludes the contact
        % model and IMU tracking
        Surfing.run_motion_IMU_Surf(workDirectory, iPar, camView, ...
            WEffort, WTrans, WAcc, WGyro, WAngles, WReg, slope, doStanding, nNodes);
    end
end

disp("InitGuess Files created");

for iPar = participants
    % Weights and flags
    WAcc = weights;
    WGyro = weights;
    WReg = '1e-3';
    doStanding = 0;
    nNodes = 100;
    WAngles = weights;
    WTrans = weights;
    WEffort = 1e-1;

    % The actual surfing simulation which includes the contact
    % model and IMU tracking
    Surfing.run_motion_IMU_Surf(workDirectory, iPar, camView, ...
        WEffort, WTrans, WAcc, WGyro, WAngles, WReg, slope, doStanding, nNodes);
end

