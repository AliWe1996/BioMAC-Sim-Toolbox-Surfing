%======================================================================
%> @file Surfing/standing3D_pelvis213_surf.m
%> @brief Function to specify the optimization problem for 3D standing
%>         prior to surfing simulations
%>
%> @author Alexander Weiss
%> @date June, 2024
%======================================================================

% ======================================================================
%> @brief Function to specify the optimization problem for 3D standing
%>
%> @param   model               Gait3d: Model which should be used for the simulation
%> @param   resultfile          String: Name of the resultfile including path
%> @param   currRearFoot        String: Rear foot of the surfer
%> @param   center_board_foot   Double array: Center position of surf board
%>                              relative to rear foot
%> @param   long_axis_board     Double: Long axis radius of Ellipse aka Board length in frontal plane
%> @param   short_axis_board    Double: Short axis radius of Ellipse aka Board length in sagittal plane
%> @param   trackingData_PE     String: Name of tracking data from PE file
%>                                       to define initial state
%> @retval  problem             Collocation: Optimization problem for 3D standing
% ======================================================================
function problem = standing3D_pelvis213_surf(model, resultfile, curRearFoot, center_board_foot, long_axis_board, short_axis_board, trackingData_PE)

%% Fixed settings
% Number of collocation nodes is one since we want to simulate static standing 
% at one time point
nNodes = 1;   
% Discretization method is unimportant here since we do not have to compute deriatives. 
% However, most of the time we use backard euler which is encoded with 'BE'.
Euler = 'BE';
% We usually use the name of the resultfile for the name of the logfile
logfile = resultfile;
% We want to plot intermediate results during solving of the problem.
plotLog = 1;

%% Create collocation problem
problem = Collocation(model, nNodes, Euler, logfile, plotLog);

%% Add states and controls including their bounds and initial values to the problem
% Get upper and lower bounds for the states
xmin = model.states.xmin;
xmax = model.states.xmax;

trackingData_PE = load(trackingData_PE);
for iInit = 2:size(trackingData_PE.dataStruct.variables.mean)
    initStateNames{iInit-1,1} = trackingData_PE.dataStruct.variables.name{iInit,1};
    initState(iInit-1,1) = trackingData_PE.dataStruct.variables.mean{iInit,1}(1,1);
    initVar(iInit-1,1) = trackingData_PE.dataStruct.variables.var{iInit,1}(1,1);
end

% Bounds for global coordinates: Start at (x, z) = (0, 0) with rot = 0
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

% Initial states from Pose Estimation
xmin(idxPelvisTilt,1) = initState(strcmp(initStateNames, 'pelvis_tilt'), 1) - initVar(strcmp(initStateNames, 'pelvis_tilt'), 1);
xmax(idxPelvisTilt,1) = initState(strcmp(initStateNames, 'pelvis_tilt'), 1) + initVar(strcmp(initStateNames, 'pelvis_tilt'), 1);
xmin(idxPelvisRot,1) = initState(strcmp(initStateNames, 'pelvis_rotation'), 1) - initVar(strcmp(initStateNames, 'pelvis_rotation'), 1);
xmax(idxPelvisRot,1) = initState(strcmp(initStateNames, 'pelvis_rotation'), 1) + initVar(strcmp(initStateNames, 'pelvis_rotation'), 1);
xmin(idxPelvisObliq,1) = initState(strcmp(initStateNames, 'pelvis_obliquity'), 1) - initVar(strcmp(initStateNames, 'pelvis_obliquity'), 1);
xmax(idxPelvisObliq,1) = initState(strcmp(initStateNames, 'pelvis_obliquity'), 1) + initVar(strcmp(initStateNames, 'pelvis_obliquity'), 1);

xmin(idxHipFlexR,1) = initState(strcmp(initStateNames, 'hip_flexion_r'), 1) - initVar(strcmp(initStateNames, 'hip_flexion_r'), 1);
xmax(idxHipFlexR,1) = initState(strcmp(initStateNames, 'hip_flexion_r'), 1) + initVar(strcmp(initStateNames, 'hip_flexion_r'), 1);
xmin(idxHipFlexL,1) = initState(strcmp(initStateNames, 'hip_flexion_l'), 1) - initVar(strcmp(initStateNames, 'hip_flexion_l'), 1);
xmax(idxHipFlexL,1) = initState(strcmp(initStateNames, 'hip_flexion_l'), 1) + initVar(strcmp(initStateNames, 'hip_flexion_l'), 1);
xmin(idxKneeFlexR,1) = initState(strcmp(initStateNames, 'knee_angle_r'), 1) - initVar(strcmp(initStateNames, 'knee_angle_r'), 1);
xmax(idxKneeFlexR,1) = initState(strcmp(initStateNames, 'knee_angle_r'), 1) + initVar(strcmp(initStateNames, 'knee_angle_r'), 1);
xmin(idxKneeFlexL,1) = initState(strcmp(initStateNames, 'knee_angle_l'), 1) - initVar(strcmp(initStateNames, 'knee_angle_l'), 1);
xmax(idxKneeFlexL,1) = initState(strcmp(initStateNames, 'knee_angle_l'), 1) + initVar(strcmp(initStateNames, 'knee_angle_l'), 1);

% Before adding the states, we want to change the bounds of the DOFs to make
% it easier for the solver to find a proper solution.

xmin(model.extractState('q', 'pelvis_tx'))          =   0        ; % should be zero due to translation invariance in horizontal plane
xmax(model.extractState('q', 'pelvis_tx'))          =   0        ;
xmin(model.extractState('q', 'pelvis_ty'))          =   0     ; % pelvis should be about 1 m above ground during standing
xmax(model.extractState('q', 'pelvis_ty'))          =   0.8      ;
xmin(model.extractState('q', 'pelvis_tz'))          =   0        ; % should be zero due to translation invariance in horizontal plane
xmax(model.extractState('q', 'pelvis_tz'))          =   0        ;
xmin(model.extractState('q', {'hip_adduction_r', 'hip_adduction_l'}))     = -40 /180*pi; % some may be needed
xmax(model.extractState('q', {'hip_adduction_r', 'hip_adduction_l'}))     = -15 /180*pi;
xmin(model.extractState('q', {'hip_rotation_r', 'hip_rotation_l'}))       = -40 /180*pi; % should be small
xmax(model.extractState('q', {'hip_rotation_r', 'hip_rotation_l'}))       =  0 /180*pi;
xmin(model.extractState('q', {'ankle_angle_r', 'ankle_angle_l'}))         = 0 /180*pi; % some may be needed
xmax(model.extractState('q', {'ankle_angle_r', 'ankle_angle_l'}))         = 70 /180*pi;
xmin(model.extractState('q', {'subtalar_angle_r', 'subtalar_angle_l'}))   = 10 /180*pi; % some may be needed
xmax(model.extractState('q', {'subtalar_angle_r', 'subtalar_angle_l'}))   =  35 /180*pi; %!!!!
xmin(model.extractState('q', {'mtp_angle_r', 'mtp_angle_l'}))             = -10 /180*pi; % some may be needed
xmax(model.extractState('q', {'mtp_angle_r', 'mtp_angle_l'}))             =  10 /180*pi;
xmin(model.extractState('q', 'lumbar_extension'))   = -10 /180*pi; % some may be needed
xmax(model.extractState('q', 'lumbar_extension'))   =  10 /180*pi;
xmin(model.extractState('q', 'lumbar_bending'))     = -10 /180*pi; % should be zero due symmetry of model
xmax(model.extractState('q', 'lumbar_bending'))     =  20 /180*pi;
xmin(model.extractState('q', 'lumbar_rotation'))    = -25 /180*pi; % should be zero due symmetry of model
xmax(model.extractState('q', 'lumbar_rotation'))    =  -20 /180*pi;
xmin(model.extractState('q', {'arm_flex_r', 'arm_flex_l'}))        = -20 /180*pi; % some may be needed
xmax(model.extractState('q', {'arm_flex_r', 'arm_flex_l'}))        =  20 /180*pi;
xmin(model.extractState('q', {'arm_add_r', 'arm_add_l'}))          = -70 /180*pi; % some may be needed
xmax(model.extractState('q', {'arm_add_r', 'arm_add_l'}))          = -40 /180*pi;
xmin(model.extractState('q', {'arm_rot_r', 'arm_rot_l'}))          = -88 /180*pi; % some may be needed
xmax(model.extractState('q', {'arm_rot_r', 'arm_rot_l'}))          =  88 /180*pi;
xmin(model.extractState('q', {'elbow_flex_r', 'elbow_flex_l'}))    = -10 /180*pi; % some may be needed
xmax(model.extractState('q', {'elbow_flex_r', 'elbow_flex_l'}))    =  88 /180*pi;
xmin(model.extractState('q', {'pro_sup_r', 'pro_sup_l'}))          = -10 /180*pi; % some may be needed
xmax(model.extractState('q', {'pro_sup_r', 'pro_sup_l'}))          = 150 /180*pi;

% Get a random initial guess for the states (within [xmin, xmax])
xinit = xmin + (xmax-xmin) .* rand(size(xmin));

% Add states to the problem using the bounds and inital values which were
% adapted above
problem.addOptimVar('states', xmin, xmax, xinit);

% Add controls to the problem using the default bounds and neutral controls
problem.addOptimVar('controls', model.controls.xmin, model.controls.xmax, model.controls.xneutral); 
%problem.addOptimVar('dur',1,1,1)

%% Add objective term(s) to the problem
% For this simulation, we want to minimize muscle effort but we do not want
% to track data. 
% We have one term for muscle effort, for which we use a weight of 1.
Weff = 1; 
% In the effort term, we use volume dependent weighting to account for
% different muscle sizes.
weightsType = 'volumeweighted'; 
% In the effort term, we want to use the squared neural excitation of the muscles.
exponent = 2; 
% Now we can add the objective term by specifiying the function handle 'effortTermMuscles'.
problem.addObjective(@effortTermMuscles, Weff, weightsType, exponent); 

% Additionally, we have arms in our 3D model which are actuated using
% torques. We use here the same weighting and the same exponent as above.
problem.addObjective(@effortTermTorques, Weff, exponent); 

%% Add constraint(s) to the problem
%problem.addConstraint(@dynamicConstraints,repmat(model.constraints.fmin,1,nNodes-1),repmat(model.constraints.fmax,1,nNodes-1))

problem.addConstraint(@equilibriumConstraints, model.constraints.fmin, model.constraints.fmax); 
problem.addConstraint(@forceFootBoardEquilibrium, zeros(1,3*nNodes),zeros(1,3*nNodes), curRearFoot,center_board_foot, long_axis_board, short_axis_board);%
problem.addConstraint(@maintainRelFootDistance, zeros(10,1),zeros(10,1), 1)
problem.addConstraint(@maintainFootOnPlane, zeros(1,3*1),zeros(1,3*1), curRearFoot)


end