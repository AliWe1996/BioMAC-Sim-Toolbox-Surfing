%======================================================================
%> @file maximizeSpeed.m
%> @brief Collocation function to minimize difference between summed CP
%> forces and board force from intersection with water
%> @details
%> Details: Collocation::forceFootBoardEquilibrium()
%>
%> @author Alexander Weiss
%> @date October, 2023
%======================================================================

%======================================================================
%> @brief
%> Computes the squared difference of computed reaction forces and the ones
%> from feet
%> 
%>
%> @param   obj                    Collocation class object
%> @para    option                String parsing the demanded output: 'objval' or 'gradient'
%>                              (or 'init' for initialization)
%> @param   X                    Double array: State vector containing at least 'states'
%> @param   currRearFoot        String: Rear foot of the surfer
%> @param   center_board_foot   Double array: Center position of surf board relative to rear foot                       
%> @param   long_axis_board     Double: Long axis radius of Ellipse aka Board length in frontal plane
%> @param   short_axis_board    Double: Short axis radius of Ellipse aka Board length in sagittal plane
%> @retval  output               Objective values for input option 'objval' or vector
%>                              with gradient
%======================================================================
function output = forceFootBoardEquilibrium(obj, option, X, rearFoot, center_board_foot, long_axis_board, short_axis_board)

%% check input parameter
if ~isfield(obj.idx,'states') % check whether speed is stored in X
    error('Model states are not stored in state vector X.')
end

numNodes = obj.nNodes;
Ncontactpoints = obj.model.nCPs;

% Case standing
if isfield(obj.idx, 'dur')
    dur = X(obj.idx.dur);
else
    dur = 0.1;
end

for iNode = 1:numNodes
    % States current Node
    x = X(obj.idx.states(:,iNode));
    % Case standing
    if iNode > 1
        x_prev = X(obj.idx.states(:,iNode-1));
    else
        x_prev = X(obj.idx.states(:,numNodes));
    end

    for i = 1:Ncontactpoints
        %Indices
        CPfx(i) = obj.model.idxForward(i+1);
        CPfy(i) = obj.model.idxUpward(i+1);
        CPfz(i) = obj.model.idxSideward(i+1);
        statFx(i) = CPfx(i)-3;
        statFy(i) = CPfy(i)-3;
        statFz(i) = CPfz(i)-3;
    
        Fyc(i,iNode) = x(statFy(i));
        Fxc(i,iNode) = x(statFx(i));
        Fzc(i,iNode) = x(statFz(i));
    
        yc(i,iNode) = x(CPfy(i));
        xc(i,iNode) = x(CPfx(i));
        zc(i,iNode) = x(CPfz(i));
    
        Myc(i,iNode) = Fyc(i,iNode) * yc(i,iNode);
        Mxc(i,iNode) = Fxc(i,iNode) * xc(i,iNode);
        Mzc(i,iNode) = Fzc(i,iNode) * zc(i,iNode);
    end
    
    % Get acting forces
    [Fx, Fy, Fz, dFx_dq, dFy_dq, dFz_dq, dFx_ddur, dFy_ddur, dFz_ddur] = obj.model.getBoardForce(x, x_prev, numNodes, rearFoot, center_board_foot, long_axis_board, short_axis_board, dur);%,...
        %point_s(iNode,:), d_point_s_dq(iNode,:,:), center_s(iNode,:), d_center_dq(iNode,:,:), center_s_prev(iNode,:), d_center_dq_prev(iNode,:,:));
    
    Fpp_x(iNode) = Fx;
    Fpp_y(iNode) = Fy;
    Fpp_z(iNode) = Fz;

    dFpp_x_dq(iNode,:) = dFx_dq;
    dFpp_y_dq(iNode,:) = dFy_dq;
    dFpp_z_dq(iNode,:) = dFz_dq;

    dFpp_x_ddur(iNode,:) = dFx_ddur;
    dFpp_y_ddur(iNode,:) = dFy_ddur;
    dFpp_z_ddur(iNode,:) = dFz_ddur;


    %sum over all contactpoint forces
    FCP_x_all(iNode) = sum(Fxc(:,iNode));
    FCP_y_all(iNode) = sum(Fyc(:,iNode));
    FCP_z_all(iNode) = sum(Fzc(:,iNode));

    % Derivatives of themselves = 1
    dFCP_x_all_dFCP(iNode,:) = -ones(Ncontactpoints,1);
    dFCP_y_all_dFCP(iNode,:) = -ones(Ncontactpoints,1);
    dFCP_z_all_dFCP(iNode,:) = -ones(Ncontactpoints,1);


    % Sum over Moments of CPs
    MCP_x_all(iNode) = sum(Mxc(:,iNode));
    MCP_y_all(iNode) = sum(Myc(:,iNode));
    MCP_z_all(iNode) = sum(Mzc(:,iNode));

    dMCP_x_all_dCP(iNode) = sum(xc(:,iNode));
    dMCP_y_all_dCP(iNode) = sum(yc(:,iNode));
    dMCP_z_all_dCP(iNode) = sum(zc(:,iNode));

    dMCP_x_all_dCP(iNode) = FCP_x_all(iNode);
    dMCP_y_all_dCP(iNode) = FCP_y_all(iNode);
    dMCP_z_all_dCP(iNode) = FCP_z_all(iNode);
end

ix = obj.idx.states(obj.model.extractState('q'), 1:obj.nNodes);
ixCP_x = obj.idx.states(obj.model.extractState('Fx'), 1:obj.nNodes);
ixCP_y = obj.idx.states(obj.model.extractState('Fy'), 1:obj.nNodes);
ixCP_z = obj.idx.states(obj.model.extractState('Fz'), 1:obj.nNodes);

if strcmp(option,'confun') % objective value
    for iNode = 1:numNodes
        output(3*(iNode-1)+1,1) = Fpp_x(iNode) - FCP_x_all(iNode);
        output(3*(iNode-1)+2,1) = Fpp_y(iNode) - FCP_y_all(iNode);
        output(3*(iNode-1)+3,1) = Fpp_z(iNode) - FCP_z_all(iNode);

        %output(1,1) = 0;
    end

    %% TODO: Add moment equilibrium
    % output(4,1) = Mpp_x - MCP_x_all;
    % output(5,1) = Mpp_y - MCP_y_all;
    % output(6,1) = Mpp_z - MCP_z_all;

elseif strcmp(option,'jacobian') % gradient_all;

    output = zeros(3,size(X,1));
    for iNode = 1:numNodes
        output(3*(iNode-1)+1,ix(:,iNode)) = dFpp_x_dq(iNode,:)';
        output(3*(iNode-1)+2,ix(:,iNode)) = dFpp_y_dq(iNode,:)';
        output(3*(iNode-1)+3,ix(:,iNode)) = dFpp_z_dq(iNode,:)';

        output(3*(iNode-1)+1, ixCP_x(:,iNode)) = dFCP_x_all_dFCP(iNode,:)';
        output(3*(iNode-1)+2, ixCP_y(:,iNode)) = dFCP_y_all_dFCP(iNode,:)';
        output(3*(iNode-1)+3, ixCP_z(:,iNode)) = dFCP_z_all_dFCP(iNode,:)';

        if isfield(obj.idx, 'dur')
            ixDur = obj.idx.dur;
            output(3*(iNode-1)+1, ixDur) = dFpp_x_ddur(iNode,:);
            output(3*(iNode-1)+2, ixDur) = dFpp_y_ddur(iNode,:);
            output(3*(iNode-1)+3, ixDur) = dFpp_z_ddur(iNode,:);
        else
            %output(3*(iNode-1)+1, ixDur) = 0;%dFpp_x_ddur(iNode,:);
            %output(3*(iNode-1)+2, ixDur) = 0;%dFpp_y_ddur(iNode,:);
            %output(3*(iNode-1)+3, ixDur) = 0;%dFpp_z_ddur(iNode,:);
        end

        %output(4,ix(:,iNode)) = Mpp_x - MCP_x_all;
        %output(5,ix(:,iNode)) = Mpp_y - MCP_y_all;
        %output(6,ix(:,iNode)) = Mpp_z - MCP_z_all;
    end

else
    error('Unknown option.');
end

end