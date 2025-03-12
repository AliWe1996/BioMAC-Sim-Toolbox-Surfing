%======================================================================
%> @file +Surfing/run_motion_IMU_Surf.m
%> @brief Function to do simulation with IMU/PE tracking for Surfing
%>
%> @author Alexander Weiss
%> @date October, 2024
%======================================================================

%======================================================================
%> @brief Function to do surfing simulation with IMU and pose estimation angles tracking
%>
%> @details
%> This function tracks data to simulate a surfing motion.
%>
%> @param  workDirectory    String: Current work directory
%> @param  iPar             Int   : Number of participant
%> @param  iView            Int   : Number of camera view pose estimation
%> @param  WEffort          String: Weight of effort term in objective
%> @param  WTrans           String: Weight of translation tracking term in objective
%> @param  WAcc             String: Weight of acceleration tracking term in objective
%> @param  WGyro            String: Weight of gyroscope tracking term in objective
%> @param  WAng             String: Weight of angle tracking term in objective
%> @param  WReg             String: Weight of the regularization term in objective
%> @param  slope            Double: Inclination angle of wave
%> @param  doStanding       Int: Flag to specify if a standing simulation is needed or not
%> @param  nNodes           Int: Number of collocation nodes
%======================================================================
function run_motion_IMU_Surf(workDirectory, iPar, iView, WEffort, WTrans, WAcc, WGyro, WAngles, WReg, slope, doStanding, nNodes)

close all
% Flags
doShowResults = 0;
% number of max. iterations
run_iter = 10;

% Fixed settings
date                =  datetime('today');
dataFolder          = ['data' filesep 'SurfingData']; % Relative from the work directorygait3d_pelvis213.osim'MarkerTracking3D' filesep subject
dataFolderIMU       = 'AvgData_IMU';
dataFolderPE        = 'AvgData_Angle';
dataInfoFile_IMU    = sprintf('Subject_0%s_Surf_IMU_avgTurns.mat', int2str(iPar));                           % IMU tracking File
dataInfoFile_PE     = sprintf('Subject_0%s_C%s_Surf_PE_avgTurns.mat', int2str(iPar), int2str(iView));        % PE tracking file
modelFile           = 'gait3d_pelvis213.osim'; %'Subject.osim';

standingFile  = sprintf('Standing_Par_0%s_%s', int2str(iPar), date);

if doStanding == 1
    resultFolderStand  = ['results' filesep 'Surfing_initGuess'];
    resultFolder       = ['results' filesep 'Surfing_initGuess'];
    resultFileStand    = [workDirectory filesep resultFolderStand filesep standingFile];
    resultFileSurfing  = sprintf('Surfing_PE_IMU_WAcc_%s_WGyro_%s_WAngle_%s_Date_%s_Sub_%s_InitStep', WAcc, WGyro, WAngles, date, int2str(iPar));
    resultFileSurfing  = [workDirectory filesep resultFolder filesep resultFileSurfing];
else
    initFolder         = ['results' filesep 'Surfing_initGuess'];
    resultFolder       = ['results' filesep 'Surfing'];  % Relative from the path of the repository
    resultFileSurfing  = sprintf('Surfing_PE_IMU_WAcc_%s_WGyro_%s_WAngle_%s_Date_%s_Sub_%s_SecStep', WAcc, WGyro, WAngles, date, int2str(iPar));
    resultFileSurfing  = [workDirectory filesep resultFolder filesep resultFileSurfing];
end

% Get absolute file names
dataTrackFile_IMU  = [workDirectory filesep dataFolder filesep dataFolderIMU filesep dataInfoFile_IMU];
dataTrackFile_PE   = [workDirectory filesep dataFolder filesep dataFolderPE  filesep dataInfoFile_PE];

% Create resultfolder if it does not exist
if ~exist([workDirectory,filesep,resultFolder], 'dir')
    mkdir([workDirectory,filesep,resultFolder]);
end

% Specific participant parameters
load(dataTrackFile_IMU);
bodyheight          = dataStruct.subjectHeight;
bodymass            = dataStruct.subjectMass;
center_board        = dataStruct.center_board;
long_board_axis     = dataStruct.longBoardAxis;
short_board_axis    = dataStruct.shortBoardAxis;
curRearFoot         = dataStruct.rearFoot;

% Create an instane of our 3D model class_slope
model = Surf3d(modelFile, bodyheight, bodymass);
% Apply slope
model.slope = [slope 0];
[model.gravity(1), model.gravity(2), model.gravity(3)] = model.applySlope(model.gravity(1), model.gravity(2), model.gravity(3));

if doStanding == 1
        % Define standing problem    
        problemStanding = Surfing.standing3D_pelvis213_surf(model, resultFileStand, curRearFoot, center_board, long_board_axis, short_board_axis, dataTrackFile_PE);
        
        % Create an object of class solver. We use most of the time the IPOPT here.
        solver = IPOPT();
        
        % Change settings of the solver
        solver.setOptionField('tol', 0.001);
        solver.setOptionField('constr_viol_tol', 0.001);
        solver.setOptionField('max_iter', 50);
        
        % Solve the optimization problem
        resultStanding = solver.solve(problemStanding);
        
        % Save the result
        resultStanding.save(resultFileStand);

        resultStanding = load(resultFileStand);
        resultStanding = resultStanding.result;
        % To plot the result we have to extract the states x from the result vector X
        x = resultStanding.X(resultStanding.problem.idx.states);
        
        % Now, we can plot the stick figure visualizing the result
        figure();
        resultStanding.problem.model.showStick(x, [], 1, 1, 1, 1, 1, 125, 30, curRearFoot, center_board, long_board_axis, short_board_axis);
        title('3D Standing');
        result = resultFileStand;
end

% Specify the optimizaton problem
W.trackTrans = str2num(WTrans)/10;
W.trackAcc   = str2num(WAcc);   % Weight of acceleration tracking term in objective
W.trackGyro  = str2num(WGyro);  % Weight of gyroscope tracking term in objective5593
W.trackAngles= str2num(WAngles);% Weight of joint angles tracking term
W.effMuscles = WEffort;         % Weight of effort term for muscles in objective
W.effTorques = 1e-01;           % Weight of effort term for torques in objective
W.reg        = str2num(WReg);   % Weight of regularization term in objective180

if doStanding == 1
    % Initial guess from standing
    initialGuessPRE = resultFileStand;
else
    % Test to use a previous trial as initialGuess
    initFileEnding  = sprintf('_Sub_%s_InitStep.mat', int2str(iPar));
    initFolderPath  = [workDirectory filesep initFolder];
    initFile        = dir(fullfile(initFolderPath, ['*', initFileEnding]));
    if size(initFile,1) > 1
        warning("Folder should only contain a single InitGuess");
        keyboard
    end
    initialGuessPRE = [initFolderPath filesep initFile(1).name];
end

% Define Surfing problem
problemSurfing = Surfing.setup_motion_IMU_Surf(model,dataTrackFile_IMU,dataTrackFile_PE,initialGuessPRE,resultFileSurfing,W,curRearFoot,center_board,long_board_axis, short_board_axis, nNodes, doStanding);

% Solve
solver = IPOPT();
solver.setOptionField('max_iter', run_iter);
solver.setOptionField('tol', 0.0005);
resultSurfing = solver.solve(problemSurfing);
resultSurfing.save(resultFileSurfing);

if doShowResults
    figure();
    x = resultSurfing.X(resultSurfing.problem.idx.states);
    resultSurfing.problem.model.showStick(x, [], 1, 0, 1, 0, 1, 140, 35,curRearFoot, center_board, long_board_axis, short_board_axis); % plot feet and CP
    title('3D Surfing');
    
    settings.plotInitialGuess = 1;
    style.figureSize = [0 0 16 26];
    resultSurfing.problem.writeMovie(resultSurfing.X, resultSurfing.filename, [],[],{1, 1, 1, 0, 1, 140, 35, curRearFoot, center_board, long_board_axis, short_board_axis});
    resultSurfing.problem.getMetabolicCost(resultSurfing.X);
    % Plot convergence
    figure(); hold on;
    plot(resultSurfing.problem.objectiveTerms(1).weightedValueHist);
    plot(resultSurfing.problem.objectiveTerms(2).weightedValueHist);
    plot(resultSurfing.problem.objectiveTerms(3).weightedValueHist);
    plot(resultSurfing.problem.objectiveTerms(4).weightedValueHist);
    %plot(resultSurfing.problem.objectiveTerms(5).weightedValueHist);
    plot(sum([resultSurfing.problem.objectiveTerms.weightedValueHist], 2));
    legend({resultSurfing.problem.objectiveTerms.name, 'Sum'});
    
    IMUTrackingFile = dataTrackFile_IMU;
    
    % Adapt style to your needs
    style = struct();
    style.xLabelText  = 'Motion in \%';
    style.subFigSettings.nCol = 9;
    style.subFigSettings.width       = 5.5; % in cm
    style.subFigSettings.height      = 5;   % in cm
    style.subFigSettings.originUp    = 3;   % in cm
    style.subFigSettings.originRight = 1.5; % in cm
    style.subFigSettings.oneUp       = 1.0; % in cm
    style.subFigSettings.oneRight    = 1.5; % in cm
    
    dofNames = resultSurfing.problem.model.dofs.Properties.RowNames;
    idxTrans = ismember(dofNames, {'pelvis_tx', 'pelvis_ty', 'pelvis_tz'});
    settings.plotInitialGuess = 1;
    style.figureSize = [0 0 16 26];
    % When also giving a filename as input, the function will automatically
    % create a pdf summarizing all information and plots. Have a look at it!
    %resultSurfing.report(settings, style, resultFileSurfing);
    settings.translation = dofNames(idxTrans);
    settings.angle = dofNames(~idxTrans); % global pelvis orientation and joint angles
    settings.torque = dofNames(resultSurfing.problem.model.idxTorqueDof); % all arm DOFs
    settings.acc = resultSurfing.problem.objectiveTerms(strcmp({resultSurfing.problem.objectiveTerms.name}, 'trackAcc')).varargin{1}.variables;
    settings.gyro = resultSurfing.problem.objectiveTerms(strcmp({resultSurfing.problem.objectiveTerms.name}, 'trackGyro')).varargin{1}.variables;
    
    simVarTableIMU = resultSurfing.problem.extractData(resultSurfing.X, settings);%, simVarTableMarker);
    
    plotVarTable(simVarTableIMU, style);
    %resultSurfing.problem.writeMotionToOsim(resultSurfing.X, [workDirectory filesep dataFolder filesep resultFolder sprintf('Osim_Sol_', int2str(iPar)];
    
    [metCost, ~, ~, metRate, CoT] = resultSurfing.problem.getMetabolicCost(resultSurfing.X);
    fprintf('metCost: %4.3f\n', metCost);
    fprintf('metRate: %4.3f\n', metRate);
end

end