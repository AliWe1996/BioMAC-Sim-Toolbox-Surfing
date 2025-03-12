%======================================================================
%> @file +Surfing/setup_motion_IMU_Surf.m
%> @brief Function to setup motion problem with IMU and PE tracking of Surfing
%>
%> @author Alexander Weiss
%> @date October, 2024
%======================================================================

% ======================================================================
%> @brief Function to setup motion problem with IMU and PE tracking
%>
%> @param   model               Gait3d: Model used for simulation
%> @param   dataFile_IMU        String: Filename containg data struct for IMU tracking
%> @param   dataFile_PE         String: Filename containg data struct for PE tracking
%> @param   initialFile         String: Result used for initial guess
%> @param   resultFile          String: Filename to log files
%> @param   W                   Struct: Weights for objective terms
%> @param   currRearFoot        String: Rear foot of the surfer
%> @param   center_board_foot   Double array: Center position of surf board relative to rear foot                       
%> @param   long_axis_board     Double: Long axis radius of Ellipse aka Board length in frontal plane
%> @param   short_axis_board    Double: Short axis radius of Ellipse aka Board length in sagittal plane
%> @param   doStanding          Int: Flag to specify if a standing simulation is needed or not
%> @param   nNodes              Int: Number of collocation nodes
%> @retval  problem             Collocation: Generated standing problem
% ======================================================================
function problem = setup_motion_IMU_Surf(model,dataFile_IMU, dataFile_PE,initialFile,resultFile,W, currRearFoot, center_board_foot, long_axis_board, short_axis_board, nNodes, doStanding)

%% Number of Nodes and portion of whole movement to simulate
N = nNodes; 
% Periodicity constraint, using full cycle works stable
portionSim = 1; % 1 = fullCircle, 2 = halfMov, 4 = quartMov

%% Tracking Data

%% Add translation to tracking data
load(dataFile_PE);

% Step 1: Identify rows containing 'translation' in their name
translation_rows = contains(dataStruct.variables.type, 'translation');

% Step 2: Remove those rows from the table
dataStruct.variables(translation_rows, :) = [];

%Check if 'translation_z' already exists
translation_z_exists = any(strcmp(dataStruct.variables.name, 'pelvis_tz'));

if ~translation_z_exists
    % Step 1: Extract length of the mean vector of pelvis_tilt
    pelvis_tilt_index = strcmp(dataStruct.variables.name, 'pelvis_tilt');
    mean_length = length(dataStruct.variables.mean{pelvis_tilt_index});
    % Step 2: Create the new mean curve for translation_z
    t = linspace(0, 1, mean_length);  % Time vector
    translation_z_mean = 3 * sin(2 * pi * t);  % Mean curve as specified
    % Step 3: Define the new row
    new_row = {
        'pelvis_tz', ...  % Name
        'translation',.... % Type
        translation_z_mean', ...  % Mean
        0.5 * ones(mean_length, 1), ...  % Var
        'm' ...  % Unit
        };
    % Step 4: Add the new row to the table
    dataStruct.variables = [dataStruct.variables; cell2table(new_row, 'VariableNames', dataStruct.variables.Properties.VariableNames)];
end
%dataStruct.variables.mean{4,1} = -dataStruct.variables.mean{4,1};
save(dataFile_PE, 'dataStruct');

%TrackingData IMU
trackingData_IMU = TrackingData.loadStruct(dataFile_IMU);
indicesStartEnd(1) = 1;
indicesStartEnd(2) = round(size(trackingData_IMU.variables.mean{1,1},1)/portionSim);
trackingData_IMU.trimData(indicesStartEnd(1), indicesStartEnd(2)-1);
targetframes = trackingData_IMU.nSamples; 
trackingData_IMU.resampleData(N);

%Tracking Data PE
trackingData_PE = TrackingData.loadStruct(dataFile_PE);
indicesStartEnd(1) = 1;
indicesStartEnd(2) = round(size(trackingData_PE.variables.mean{3,1},1)/portionSim);
trackingData_PE.trimData(indicesStartEnd(1), indicesStartEnd(2)-1);
trackingData_PE.resampleData(N);

% Set initial state of the surf cycle
for iInit = 2:size(trackingData_PE.variables.mean)
    initStateNames{iInit-1,1} = trackingData_PE.variables.name{iInit,1};
    initState(iInit-1,1) = trackingData_PE.variables.mean{iInit,1}(1,1);
    initVar(iInit-1,1) = trackingData_PE.variables.var{iInit,1}(1,1);
end

%% Problem formulation
% Create problem
Euler = 'BE';
plotLog = 1;
problem = Collocation(model,N,Euler,resultFile,plotLog);

% Add variables which are optimized
states_min = repmat(model.states.xmin,1,N+1);
states_max = repmat(model.states.xmax,1,N+1);
controls_min = repmat(model.controls.xmin, 1, N+1); 
controls_max = repmat(model.controls.xmax, 1, N+1);

idxCPxc = model.extractState('xc');
idxCPzc = model.extractState('zc');
p_global_x = [-5, 5];
p_global_z = [-5, 5];
states_min(idxCPxc, :) = p_global_x(1);
states_max(idxCPxc, :) = p_global_x(2);
states_min(idxCPzc, :) = p_global_z(1);
states_max(idxCPzc, :) = p_global_z(2);

% Bounds for global coordinates
idxPelvisX = model.extractState('q', 'pelvis_tx');
idxPelvisY = model.extractState('q', 'pelvis_ty');
idxPelvisZ = model.extractState('q', 'pelvis_tz');
idxdotPelvisZ = model.extractState('qdot', 'pelvis_tz');
idxPelvisTilt = model.extractState('q', 'pelvis_tilt');
idxPelvisRot = model.extractState('q', 'pelvis_rotation');
idxPelvisObliq = model.extractState('q', 'pelvis_obliquity');
idxHipFlexR = model.extractState('q', 'hip_flexion_r');
idxHipFlexL = model.extractState('q', 'hip_flexion_l');
idxKneeFlexR = model.extractState('q', 'knee_angle_r');
idxKneeFlexL = model.extractState('q', 'knee_angle_l');
idxAnkleFlexR = model.extractState('q', 'ankle_angle_r');
idxAnkleFlexL = model.extractState('q', 'ankle_angle_l');

% Bounds for pelvis orientation
states_min(idxPelvisTilt,:) = -pi;
states_max(idxPelvisTilt,:) = pi;
states_min(idxPelvisRot,:) = -2*pi;
states_max(idxPelvisRot,:) = 2*pi;
states_min(idxPelvisObliq,:) = -pi;
states_max(idxPelvisObliq,:) = pi;

% Bounds for ankle flexion
states_min(idxAnkleFlexR,:) = -20 /180*pi;
states_max(idxAnkleFlexR,:) = 50 /180*pi;
states_min(idxAnkleFlexL,:) = -20 /180*pi;
states_max(idxAnkleFlexL,:) = 50 /180*pi;

%Bounds for pelvis position
states_min(idxPelvisX,:) = -1.5;
states_max(idxPelvisX,:) = 1.5;
states_min(idxPelvisZ,:) = -5;
states_max(idxPelvisZ,:) = 5;
states_min(idxPelvisY,:) = 0;
states_max(idxPelvisY,:) = 2;

% Bounds for hip joint
states_min(model.extractState('q', 'hip_adduction_r'), :)     = -40 /180*pi; 
states_max(model.extractState('q', 'hip_adduction_r'), :)     = -15 /180*pi;
states_min(model.extractState('q', 'hip_rotation_r'), :)       = -40 /180*pi; 
states_max(model.extractState('q', 'hip_rotation_r'), :)       =  0 /180*pi;
states_min(model.extractState('q', 'hip_adduction_l'), :)     = -40 /180*pi; 
states_max(model.extractState('q', 'hip_adduction_l'), :)     = -15 /180*pi;
states_min(model.extractState('q', 'hip_rotation_l'), :)       = -40 /180*pi; 
states_max(model.extractState('q', 'hip_rotation_l'), :)       =  0 /180*pi;

% Add states and controls to problem
problem.addOptimVar('states',states_min,states_max);
problem.addOptimVar('controls',controls_min, controls_max);

% targetLength
h = 1/1600; % 1600 Hz of IMUs
targetdur =  h*(targetframes-1);

% Check if length of PE and IMU tracking data is equally long
load(dataFile_PE);
if abs(targetdur - dataStruct.variables.mean{1}/portionSim) > 0.3
    error('Propably wrong cycle allocation');
end

problem.addOptimVar('dur',targetdur,targetdur);

targetspeed = 0;
if portionSim == 1
    problem.addOptimVar('speed', targetspeed, targetspeed);
end
% Initialize the problem with an old result specified in initialFile
problem.makeinitialguess(initialFile);


% Add translation term
trackingDataTranslation = trackingData_PE.extractData('translation');
% Add tracking terms
% sensSignals = {'torso', 'femur_r', 'femur_l', 'tibia_r', 'tibia_l', 'toes_r', 'toes_l', 'ulna_r', 'ulna_l'};
trackingDataAcc = trackingData_IMU.extractData('acc');%, sensSignals);
if portionSim > 1
    trackingDataAcc.trimData(1, trackingDataAcc.nSamples-1); % remove last sample since it cannot be simulated wihout nNode+1 -1 for dyn FirstNode
end
trackingDataGyro = trackingData_IMU.extractData('gyro');%, sensSignals);

problem.addObjective(@trackTranslations,W.trackTrans,trackingDataTranslation);
if doStanding == 0
    problem.addObjective(@trackAcc,W.trackAcc,trackingDataAcc, 'minMaxSquared');
    problem.addObjective(@trackGyro,W.trackGyro,trackingDataGyro, 'minMaxSquared');
end


% Add tracking terms
%AngleSignals = {'pelvis_rotation', 'pelvis_tilt', 'hip_flexion_r', 'knee_angle_r', 'hip_flexion_l', 'knee_angle_l'};
trackingDataAngles = trackingData_PE.extractData('angle');%, AngleSignals);
trackingDataAngles.trimData(1, trackingDataAngles.nSamples); % remove last sample since it cannot be simulated wihout
problem.addObjective(@trackAngles,W.trackAngles,trackingDataAngles);
%problem.addObjective(@trackFirstNode,W.trackFirstNode,trackingDataAngles);


% Add effort terms
speedWeighting = 0;
problem.addObjective(@effortTermMuscles,W.effMuscles,'volumeweighted',3);%,speedWeighting);
problem.addObjective(@effortTermTorques,W.effTorques,2);%, speedWeighting);%volumeweighted
% Add regularization term
problem.addObjective(@regTerm,W.reg)

sym = 0;
% Update constraints to use forceFootEquilibrium for contact force
% calculation

% Add constraints
if portionSim > 1
    problem.addConstraint(@dynamicsFirstNodeConstraint,model.constraints.fmin,model.constraints.fmax)
    problem.addConstraint(@dynamicConstraints,repmat(model.constraints.fmin,1,N-1),repmat(model.constraints.fmax,1,N-1))

else
     problem.addConstraint(@periodicityConstraint, zeros(model.nStates+model.nControls,1),zeros(model.nStates+model.nControls,1),sym)
     problem.addConstraint(@dynamicConstraints,repmat(model.constraints.fmin,1,N),repmat(model.constraints.fmax,1,N))    
 end
    
problem.addConstraint(@maintainRelFootDistance, zeros(1,10),zeros(1,10), 1) % if orient = 1 zeros(1,10), zeros(1,1) else
problem.addConstraint(@maintainFootOnPlane, zeros(1,3*N),zeros(1,3*N), currRearFoot)
if doStanding == 0
    problem.addConstraint(@forceFootBoardEquilibrium, zeros(1,3*N),zeros(1,3*N), currRearFoot, center_board_foot, long_axis_board, short_axis_board)
end
end