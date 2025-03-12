%======================================================================
%> @file maximizeSpeed.m
%> @brief Collocation function to ensure front foot is on plane of surf
%>        board
%> @details
%> Details: Collocation::maintainFootOnPlane()
%>
%> @author Alexander Weiss
%> @date October, 2023
%======================================================================

%======================================================================
%> @brief
%> Computes the perpendiculat vector from upward calcaneus to the one
%> between both calcaneus
%>
%> @todo This function is not yet working if there are multiple entries for
%> the speed (i.e. in the 3D model).
%>
%> @param  obj           Collocation class object
%> @param  option        String parsing the demanded output: 'objval' or 'gradient'
%>                       (or 'init' for initialization)
%> @param  X             Double array: State vector containing at least 'states'
%> @param  currRearFoot  String: Rear foot of the surfer
%> @retval output        Objective values for input option 'objval' or vector
%>                       with gradient
%======================================================================
function output = maintainFootOnPlane(obj, option, X, rearFoot)

%% check input parameter
if ~isfield(obj.idx,'states') % check whether speed is stored in X
    error('Model states are not stored in state vector X.')
end    
    
numNodes = obj.nNodes;
% Calceneus index
calcR_no = find(strcmp(obj.model.segments.Properties.RowNames, 'calcn_r'));
calcL_no = find(strcmp(obj.model.segments.Properties.RowNames, 'calcn_l'));

if strcmp(rearFoot, 'left')

    for iNode = 1:numNodes
        x = X(obj.idx.states(:,iNode));
        % Foot position and orientation
        [FK,dFKdq]= obj.model.getFkin(x(1:obj.model.nDofs));
        p_calc_r(iNode,:) = FK((calcR_no-2)*12 + (1:3));
        dp_calc_r_dq(iNode,:,:) = dFKdq((calcR_no-2)*12 +(1:3),:);
        p_calc_l(iNode,:) = FK((calcL_no-2)*12 + (1:3));
        dp_calc_l_dq(iNode,:,:) = dFKdq((calcL_no-2)*12 +(1:3),:);
        
        R_calc_r(iNode,:,:) = reshape(FK((calcR_no-2)*12 +(4:12)), 3, 3)';
        R_calc_l(iNode,:,:) = reshape(FK((calcL_no-2)*12 +(4:12)), 3, 3)';
        dR_calc_rdq(iNode,:,:) = dFKdq((calcR_no-2)*12 +(4:12),:);
        dR_calc_ldq(iNode,:,:) = dFKdq((calcL_no-2)*12 +(4:12),:);
    
        % Vector between feet
        dp_calc(iNode,:) = p_calc_r(iNode,:) - p_calc_l(iNode,:);

        vec_rl_x(iNode) = p_calc_r(iNode,1) - p_calc_l(iNode,1);
        vec_rl_y(iNode) = p_calc_r(iNode,2) - p_calc_l(iNode,2);
        vec_rl_z(iNode) = p_calc_r(iNode,3) - p_calc_l(iNode,3);
    
        vec_rl(iNode) = sqrt(vec_rl_x(iNode)^2 + vec_rl_y(iNode)^2 + vec_rl_z(iNode)^2);
    
        dvec_rl_dqx(iNode,:) = 2 * (vec_rl_x(iNode)) .* (dp_calc_r_dq(iNode,1,:) - dp_calc_l_dq(iNode,1,:));
        dvec_rl_dqy(iNode,:) = 2 * (vec_rl_y(iNode)) .* (dp_calc_r_dq(iNode,2,:) - dp_calc_l_dq(iNode,2,:));
        dvec_rl_dqz(iNode,:) = 2 * (vec_rl_z(iNode)) .* (dp_calc_r_dq(iNode,3,:) - dp_calc_l_dq(iNode,3,:));
    
        dvec_rl_dq(iNode,:) = (dvec_rl_dqx(iNode,:) + dvec_rl_dqy(iNode,:) + dvec_rl_dqz(iNode,:)) ./ (2 * vec_rl(iNode));
        
        %Derivatives
        dpp_calc_dq(iNode,:,:) = dp_calc_r_dq(iNode,:,:) - dp_calc_l_dq(iNode,:,:);
    
        % Normal vector 
        N_right_cur = reshape(R_calc_r(iNode,:,2),3,1);
        N_left_cur = reshape(R_calc_l(iNode,:,2),3,1);
       
        normN_right(iNode) = norm(N_right_cur);
        normN_left(iNode) = norm(N_left_cur);
        
        N_right(iNode,:) = N_right_cur;
        N_left(iNode,:) = N_left_cur;

        N_left_x(iNode,:) = reshape(R_calc_l(iNode,:,1),3,1);
        N_left_z(iNode,:) = reshape(R_calc_l(iNode,:,3),3,1);
       
    end
    
    for iNode = 2:numNodes
        dp_calc_sum(iNode) = (vec_rl(1) - vec_rl(iNode));      
    end
  
    ix = obj.idx.states(obj.model.extractState('q'), 1:obj.nNodes);  
    
    if strcmp(option,'confun') % objective value
        
        for iNode = 1:numNodes
            % Ensure parallelity in the both feets (calacaneus) y-axis
            output(3*(iNode-1)+1,1) = dot(reshape(N_right(iNode,:),3,1), reshape(N_left(iNode,:),3,1)) - 1;
            % Ensure that the left foots normal x-vector is perpendicular to the vector between the feet 
            output(3*(iNode-1)+2,1) = dot(reshape(N_left_x(iNode,:),3,1), dp_calc(iNode,:));
            % Ensure that the left foots normal z-vector is perpendicular to the vector between the feet 
            output(3*(iNode-1)+3,1) = dot(reshape(N_left(iNode,:),3,1), dp_calc(iNode,:));
        end
        
    elseif strcmp(option,'jacobian') % gradient
      
        output = zeros(1,size(X,1));

        for iNode = 1:numNodes
            
            output(3*(iNode-1)+1, ix(:,iNode)) = (reshape(dR_calc_rdq(iNode,2,:), 1, size(ix,1)) .* R_calc_l(iNode, 1,2) + R_calc_r(iNode, 1,2) .* reshape(dR_calc_ldq(iNode,2,:), 1, size(ix,1)) + ...
                                      reshape(dR_calc_rdq(iNode,5,:), 1, size(ix,1)) .* R_calc_l(iNode,2,2) + R_calc_r(iNode,2,2) .* reshape(dR_calc_ldq(iNode,5,:), 1, size(ix,1)) + ... 
                                      reshape(dR_calc_rdq(iNode,8,:), 1, size(ix,1)) .* R_calc_l(iNode,3,2) + R_calc_r(iNode, 3,2) .* reshape(dR_calc_ldq(iNode,8,:), 1, size(ix,1)))';
    
            output(3*(iNode-1)+2, ix(:,iNode)) = (reshape(dR_calc_ldq(iNode,1,:), 1, size(ix,1)) .* dp_calc(iNode, 1) + R_calc_l(iNode, 1,1) .* reshape(dpp_calc_dq(iNode, 1, :), 1, size(ix,1)) + ...
                                      reshape(dR_calc_ldq(iNode,4,:), 1, size(ix,1)) .* dp_calc(iNode, 2) + R_calc_l(iNode, 2,1) .* reshape(dpp_calc_dq(iNode, 2, :), 1, size(ix,1))  + ... 
                                      reshape(dR_calc_ldq(iNode,7,:), 1, size(ix,1)) .* dp_calc(iNode, 3) + R_calc_l(iNode, 3,1) .* reshape(dpp_calc_dq(iNode, 3, :), 1, size(ix,1)))' ;
    
            output(3*(iNode-1)+3, ix(:,iNode)) = (reshape(dR_calc_ldq(iNode,2,:), 1, size(ix,1)) .* dp_calc(iNode, 1) + R_calc_l(iNode, 1,2) .* reshape(dpp_calc_dq(iNode, 1, :), 1, size(ix,1)) + ...
                                      reshape(dR_calc_ldq(iNode,5,:), 1, size(ix,1)) .* dp_calc(iNode, 2) + R_calc_l(iNode, 2,2) .* reshape(dpp_calc_dq(iNode, 2, :), 1, size(ix,1))  + ... 
                                      reshape(dR_calc_ldq(iNode,8,:), 1, size(ix,1)) .* dp_calc(iNode, 3) + R_calc_l(iNode, 3,2) .* reshape(dpp_calc_dq(iNode, 3, :), 1, size(ix,1)))' ;
        end
   
    else
        error('Unknown option');
    end
    
else %RIGHT rear foot
    for iNode = 1:numNodes
        x = X(obj.idx.states(:,iNode));
        [FK,dFKdq]= obj.model.getFkin(x(1:obj.model.nDofs));
        p_calc_r(iNode,:) = FK((calcR_no-2)*12 + (1:3));
        dp_calc_r_dq(iNode,:,:) = dFKdq((calcR_no-2)*12 +(1:3),:);
        p_calc_l(iNode,:) = FK((calcL_no-2)*12 + (1:3));
        dp_calc_l_dq(iNode,:,:) = dFKdq((calcL_no-2)*12 +(1:3),:);
        
        R_calc_r(iNode,:,:) = reshape(FK((calcR_no-2)*12 +(4:12)), 3, 3)';
        R_calc_l(iNode,:,:) = reshape(FK((calcL_no-2)*12 +(4:12)), 3, 3)';
        dR_calc_rdq(iNode,:,:) = dFKdq((calcR_no-2)*12 +(4:12),:);
        dR_calc_ldq(iNode,:,:) = dFKdq((calcL_no-2)*12 +(4:12),:);
    
        % Vector between feet
        dp_calc(iNode,:) = p_calc_l(iNode,:) - p_calc_r(iNode,:);
        
    
        vec_rl_x(iNode) = p_calc_l(iNode,1) - p_calc_r(iNode,1);
        vec_rl_y(iNode) = p_calc_l(iNode,2) - p_calc_r(iNode,2);
        vec_rl_z(iNode) = p_calc_l(iNode,3) - p_calc_r(iNode,3);
    
        vec_rl(iNode) = sqrt(vec_rl_x(iNode)^2 + vec_rl_y(iNode)^2 + vec_rl_z(iNode)^2);
    
        dvec_rl_dqx(iNode,:) = 2 * (vec_rl_x(iNode)) .* (dp_calc_l_dq(iNode,1,:) - dp_calc_r_dq(iNode,1,:));
        dvec_rl_dqy(iNode,:) = 2 * (vec_rl_y(iNode)) .* (dp_calc_l_dq(iNode,2,:) - dp_calc_r_dq(iNode,2,:));
        dvec_rl_dqz(iNode,:) = 2 * (vec_rl_z(iNode)) .* (dp_calc_l_dq(iNode,3,:) - dp_calc_r_dq(iNode,3,:));
    
        dvec_rl_dq(iNode,:) = (dvec_rl_dqx(iNode,:) + dvec_rl_dqy(iNode,:) + dvec_rl_dqz(iNode,:)) ./ (2 * vec_rl(iNode));

        %Derivatives
        dpp_calc_dq(iNode,:,:) = dp_calc_l_dq(iNode,:,:) - dp_calc_r_dq(iNode,:,:);
    
        N_right_cur = reshape(R_calc_r(iNode,:,2),3,1);
        N_left_cur = reshape(R_calc_l(iNode,:,2),3,1);
       
    
        normN_right(iNode) = norm(N_right_cur);
        normN_left(iNode) = norm(N_left_cur);
        
        N_right(iNode,:) = N_right_cur;%/normN_right(iNode);
        N_right_x(iNode,:) = reshape(R_calc_r(iNode,:,1),3,1);
        N_right_z(iNode,:) = reshape(R_calc_r(iNode,:,3),3,1);
       
        N_left(iNode,:) = N_left_cur;%/normN_left(iNode);
        N_left_x(iNode,:) = reshape(R_calc_l(iNode,:,1),3,1);
        N_left_z(iNode,:) = reshape(R_calc_l(iNode,:,3),3,1);
       
    end
    
    for iNode = 2:numNodes
        dp_calc_sum(iNode) = (vec_rl(1) - vec_rl(iNode));      
    end
    
    ix = obj.idx.states(obj.model.extractState('q'), 1:obj.nNodes);
    
    
    if strcmp(option,'confun') % objective value
        
        for iNode = 1:numNodes
            % Ensure parallelity in the both feets (calacaneus) y-axis
            output(3*(iNode-1)+1,1) = dot(reshape(N_left(iNode,:),3,1), reshape(N_right(iNode,:),3,1)) - 1;
            % Ensure that the left foots normal x-vector is perpendicular to the vector between the feet 
            output(3*(iNode-1)+2,1) = dot(reshape(N_right_x(iNode,:),3,1), dp_calc(iNode,:));
            % Ensure that the left foots normal z-vector is perpendicular to the vector between the feet
            output(3*(iNode-1)+3,1) = dot(reshape(N_right(iNode,:),3,1), dp_calc(iNode,:));
        end

    elseif strcmp(option,'jacobian') % gradient

        output = zeros(3,size(X,1));

        for iNode = 1:numNodes

            output(3*(iNode-1)+1, ix(:,iNode)) = (reshape(dR_calc_ldq(iNode,2,:), 1, size(ix,1)) .* R_calc_r(iNode, 1,2) + R_calc_l(iNode, 1,2) .* reshape(dR_calc_rdq(iNode,2,:), 1, size(ix,1)) + ...
                reshape(dR_calc_ldq(iNode,5,:), 1, size(ix,1)) .* R_calc_r(iNode,2,2) + R_calc_l(iNode,2,2) .* reshape(dR_calc_rdq(iNode,5,:), 1, size(ix,1)) + ...
                reshape(dR_calc_ldq(iNode,8,:), 1, size(ix,1)) .* R_calc_r(iNode,3,2) + R_calc_l(iNode, 3,2) .* reshape(dR_calc_rdq(iNode,8,:), 1, size(ix,1)))';

            output(3*(iNode-1)+2, ix(:,iNode)) = (reshape(dR_calc_rdq(iNode,1,:), 1, size(ix,1)) .* dp_calc(iNode, 1) + R_calc_r(iNode, 1,1) .* reshape(dpp_calc_dq(iNode, 1, :), 1, size(ix,1)) + ...
                reshape(dR_calc_rdq(iNode,4,:), 1, size(ix,1)) .* dp_calc(iNode, 2) + R_calc_r(iNode, 2,1) .* reshape(dpp_calc_dq(iNode, 2, :), 1, size(ix,1))  + ...
                reshape(dR_calc_rdq(iNode,7,:), 1, size(ix,1)) .* dp_calc(iNode, 3) + R_calc_r(iNode, 3,1) .* reshape(dpp_calc_dq(iNode, 3, :), 1, size(ix,1)))' ;

            output(3*(iNode-1)+3, ix(:,iNode)) = (reshape(dR_calc_rdq(iNode,2,:), 1, size(ix,1)) .* dp_calc(iNode, 1) + R_calc_r(iNode, 1,2) .* reshape(dpp_calc_dq(iNode, 1, :), 1, size(ix,1)) + ...
                reshape(dR_calc_rdq(iNode,5,:), 1, size(ix,1)) .* dp_calc(iNode, 2) + R_calc_r(iNode, 2,2) .* reshape(dpp_calc_dq(iNode, 2, :), 1, size(ix,1))  + ...
                reshape(dR_calc_rdq(iNode,8,:), 1, size(ix,1)) .* dp_calc(iNode, 3) + R_calc_r(iNode, 3,2) .* reshape(dpp_calc_dq(iNode, 3, :), 1, size(ix,1)))' ;
        end
    else
        error('Unknown option');
    end
end
end