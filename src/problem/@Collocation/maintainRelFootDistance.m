%======================================================================
%> @file maximizeSpeed.m
%> @brief Collocation function to maintain a constant relative foot
%>  distance and optionally orientation
%> @details
%> Details: Collocation::maintainRelFootDistance()
%>
%> @author Alexander Weiss
%> @date October,2023
%======================================================================

%======================================================================
%> @brief
%> Computes the squared summed difference between first node foot distance and
%> (optioanlly orientation) and all nodes
%>
%>
%> @param  obj           Collocation class object
%> @param  option        String parsing the demanded output: 'objval' or 'gradient'
%>                       (or 'init' for initialization)
%> @param  X             Double array: State vector containing at least 'states'
%> @param  orient        Bool: If relative foot orientation should be
%>                       constant (default: 0)
%> @retval output        Objective values for input option 'objval' or vector
%>                       with gradient
%======================================================================
function output = maintainRelFootDistance(obj, option, X, orient)

%% check input parameter
if ~isfield(obj.idx,'states') % check whether speed is stored in X
    error('Model states are not stored in state vector X.')
end

numNodes = obj.nNodes;

calcR_no = find(strcmp(obj.model.segments.Properties.RowNames, 'calcn_r'));
calcL_no = find(strcmp(obj.model.segments.Properties.RowNames, 'calcn_l'));

for iNode = 1:numNodes
    x = X(obj.idx.states(:,iNode));
    [FK,dFKdq]= obj.model.getFkin(x(1:obj.model.nDofs));
    p_calc_r(iNode,:) = FK((calcR_no-2)*12 + (1:3));
    dp_calc_r_dq(iNode,:,:) = dFKdq((calcR_no-2)*12 +(1:3),:);
    p_calc_l(iNode,:) = FK((calcL_no-2)*12 + (1:3));
    dp_calc_l_dq(iNode,:,:) = dFKdq((calcL_no-2)*12 +(1:3),:);

    R_calc_r = reshape(FK((calcR_no-2)*12 +(4:12)), 3, 3)';
    R_calc_l = reshape(FK((calcL_no-2)*12 +(4:12)), 3, 3)';
    dR_calc_rdq = dFKdq((calcR_no-2)*12 +(4:12),:);
    dR_calc_ldq = dFKdq((calcL_no-2)*12 +(4:12),:);

    dp_calc(iNode,:) = p_calc_r(iNode,:) - p_calc_l(iNode,:);
    dR_calc(iNode,:,:) = R_calc_r(:,:) * R_calc_l(:,:)';

    % Vector between both calcaneus
    vec_rl_x(iNode) = p_calc_r(iNode,1) - p_calc_l(iNode,1);
    vec_rl_y(iNode) = p_calc_r(iNode,2) - p_calc_l(iNode,2);
    vec_rl_z(iNode) = p_calc_r(iNode,3) - p_calc_l(iNode,3);

    vec_rl(iNode) = sqrt(vec_rl_x(iNode)^2 + vec_rl_y(iNode)^2 + vec_rl_z(iNode)^2);

    % Derivative
    dvec_rl_dqx(iNode,:) = 2 * (vec_rl_x(iNode)) .* (dp_calc_r_dq(iNode,1,:) - dp_calc_l_dq(iNode,1,:));
    dvec_rl_dqy(iNode,:) = 2 * (vec_rl_y(iNode)) .* (dp_calc_r_dq(iNode,2,:) - dp_calc_l_dq(iNode,2,:));
    dvec_rl_dqz(iNode,:) = 2 * (vec_rl_z(iNode)) .* (dp_calc_r_dq(iNode,3,:) - dp_calc_l_dq(iNode,3,:));

    dvec_rl_dq(iNode,:) = (dvec_rl_dqx(iNode,:) + dvec_rl_dqy(iNode,:) + dvec_rl_dqz(iNode,:)) ./ (2 * vec_rl(iNode));


    %Derivatives
    for iq = 1:obj.model.nDofs
        dR_calc_rdq_reshape = reshape(dR_calc_rdq(:,iq),3, 3)';
        dR_calc_ldq_reshape = reshape(dR_calc_ldq(:,iq),3, 3);

        %TODO: Make nicer
        dR_calc_dq(iq,1,1) = dR_calc_rdq_reshape(1,1) * R_calc_l(1,1) + dR_calc_rdq_reshape(1,2) * R_calc_l(1,2) + dR_calc_rdq_reshape(1,3) * R_calc_l(1,3) + R_calc_r(1,1) * dR_calc_ldq_reshape(1,1) + R_calc_r(1,2) * dR_calc_ldq_reshape(2,1) + R_calc_r(1,3) * dR_calc_ldq_reshape(3,1);
        dR_calc_dq(iq,1,2) = dR_calc_rdq_reshape(1,1) * R_calc_l(2,1) + dR_calc_rdq_reshape(1,2) * R_calc_l(2,2) + dR_calc_rdq_reshape(1,3) * R_calc_l(2,3) + R_calc_r(1,1) * dR_calc_ldq_reshape(1,2) + R_calc_r(1,2) * dR_calc_ldq_reshape(2,2) + R_calc_r(1,3) * dR_calc_ldq_reshape(3,2);
        dR_calc_dq(iq,1,3) = dR_calc_rdq_reshape(1,1) * R_calc_l(3,1) + dR_calc_rdq_reshape(1,2) * R_calc_l(3,2) + dR_calc_rdq_reshape(1,3) * R_calc_l(3,3) + R_calc_r(1,1) * dR_calc_ldq_reshape(1,3) + R_calc_r(1,2) * dR_calc_ldq_reshape(2,3) + R_calc_r(1,3) * dR_calc_ldq_reshape(3,3);

        dR_calc_dq(iq,2,1) = dR_calc_rdq_reshape(2,1) * R_calc_l(1,1) + dR_calc_rdq_reshape(2,2) * R_calc_l(1,2) + dR_calc_rdq_reshape(2,3) * R_calc_l(1,3) + R_calc_r(2,1) * dR_calc_ldq_reshape(1,1) + R_calc_r(2,2) * dR_calc_ldq_reshape(2,1) + R_calc_r(2,3) * dR_calc_ldq_reshape(3,1);
        dR_calc_dq(iq,2,2) = dR_calc_rdq_reshape(2,1) * R_calc_l(2,1) + dR_calc_rdq_reshape(2,2) * R_calc_l(2,2) + dR_calc_rdq_reshape(2,3) * R_calc_l(2,3) + R_calc_r(2,1) * dR_calc_ldq_reshape(1,2) + R_calc_r(2,2) * dR_calc_ldq_reshape(2,2) + R_calc_r(2,3) * dR_calc_ldq_reshape(3,2);
        dR_calc_dq(iq,2,3) = dR_calc_rdq_reshape(2,1) * R_calc_l(3,1) + dR_calc_rdq_reshape(2,2) * R_calc_l(3,2) + dR_calc_rdq_reshape(2,3) * R_calc_l(3,3) + R_calc_r(2,1) * dR_calc_ldq_reshape(1,3) + R_calc_r(2,2) * dR_calc_ldq_reshape(2,3) + R_calc_r(2,3) * dR_calc_ldq_reshape(3,3);

        dR_calc_dq(iq,3,1) = dR_calc_rdq_reshape(3,1) * R_calc_l(1,1) + dR_calc_rdq_reshape(3,2) * R_calc_l(1,2) + dR_calc_rdq_reshape(3,3) * R_calc_l(1,3) + R_calc_r(3,1) * dR_calc_ldq_reshape(1,1) + R_calc_r(3,2) * dR_calc_ldq_reshape(2,1) + R_calc_r(3,3) * dR_calc_ldq_reshape(3,1);
        dR_calc_dq(iq,3,2) = dR_calc_rdq_reshape(3,1) * R_calc_l(2,1) + dR_calc_rdq_reshape(3,2) * R_calc_l(2,2) + dR_calc_rdq_reshape(3,3) * R_calc_l(2,3) + R_calc_r(3,1) * dR_calc_ldq_reshape(1,2) + R_calc_r(3,2) * dR_calc_ldq_reshape(2,2) + R_calc_r(3,3) * dR_calc_ldq_reshape(3,2);
        dR_calc_dq(iq,3,3) = dR_calc_rdq_reshape(3,1) * R_calc_l(3,1) + dR_calc_rdq_reshape(3,2) * R_calc_l(3,2) + dR_calc_rdq_reshape(3,3) * R_calc_l(3,3) + R_calc_r(3,1) * dR_calc_ldq_reshape(1,3) + R_calc_r(3,2) * dR_calc_ldq_reshape(2,3) + R_calc_r(3,3) * dR_calc_ldq_reshape(3,3);
    end
    dR_calc_dq_all(iNode,:,:) = reshape(dR_calc_dq(:,:,:), 33,9)';
    ddp_calc_dq(iNode,:,:) = dp_calc_r_dq(iNode,:,:) - dp_calc_l_dq(iNode,:,:);

end

% Constant relative distance and orientation depending on first node
for iNode = 1:numNodes
    dp_calc_sum(iNode) = (vec_rl(1) - vec_rl(iNode));
    dR_calc_sum(iNode,:,:) = (dR_calc(1,:,:) - dR_calc(iNode,:,:));
end

ix = obj.idx.states(obj.model.extractState('q'), 1:obj.nNodes);

if strcmp(option,'confun') % objective value

    output(1,1) = sum(dp_calc_sum(1,:).^2)/obj.nNodes;

    if orient == 1
        output(2,1) = (sum(dR_calc_sum(:,1,1).^2))/obj.nNodes;
        output(3,1) = (sum(dR_calc_sum(:,2,1).^2))/obj.nNodes;
        output(4,1) = (sum(dR_calc_sum(:,3,1).^2))/obj.nNodes;
        output(5,1) = (sum(dR_calc_sum(:,1,2).^2))/obj.nNodes;
        output(6,1) = (sum(dR_calc_sum(:,2,2).^2))/obj.nNodes;
        output(7,1) = (sum(dR_calc_sum(:,3,2).^2))/obj.nNodes;
        output(8,1) = (sum(dR_calc_sum(:,1,3).^2))/obj.nNodes;
        output(9,1) = (sum(dR_calc_sum(:,2,3).^2))/obj.nNodes;
        output(10,1) = (sum(dR_calc_sum(:,3,3).^2))/obj.nNodes;
    end


elseif strcmp(option,'jacobian') % gradient

    % Constant relative foot orientation
    if orient == 1
        output = zeros(10,size(X,1));
        output(2, ix(:,1)) =  (2*((obj.nNodes-1) * dR_calc(1,1,1) - sum(dR_calc(2:obj.nNodes,1,1))) .*  dR_calc_dq_all(1,1,:)) / obj.nNodes;
        output(3, ix(:,1)) =  (2*((obj.nNodes-1) * dR_calc(1,2,1) - sum(dR_calc(2:obj.nNodes,2,1))) .*  dR_calc_dq_all(1,2,:)) / obj.nNodes;
        output(4, ix(:,1)) =  (2*((obj.nNodes-1) * dR_calc(1,3,1) - sum(dR_calc(2:obj.nNodes,3,1))) .*  dR_calc_dq_all(1,3,:)) / obj.nNodes;
        output(5, ix(:,1)) =  (2*((obj.nNodes-1) * dR_calc(1,1,2) - sum(dR_calc(2:obj.nNodes,1,2))) .*  dR_calc_dq_all(1,4,:)) / obj.nNodes;
        output(6, ix(:,1)) =  (2*((obj.nNodes-1) * dR_calc(1,2,2) - sum(dR_calc(2:obj.nNodes,2,2))) .*  dR_calc_dq_all(1,5,:)) / obj.nNodes;
        output(7, ix(:,1)) =  (2*((obj.nNodes-1) * dR_calc(1,3,2) - sum(dR_calc(2:obj.nNodes,3,2))) .*  dR_calc_dq_all(1,6,:)) / obj.nNodes;
        output(8, ix(:,1)) =  (2*((obj.nNodes-1) * dR_calc(1,1,3) - sum(dR_calc(2:obj.nNodes,1,3))) .*  dR_calc_dq_all(1,7,:)) / obj.nNodes;
        output(9, ix(:,1)) =  (2*((obj.nNodes-1) * dR_calc(1,2,3) - sum(dR_calc(2:obj.nNodes,2,3))) .*  dR_calc_dq_all(1,8,:)) / obj.nNodes;
        output(10, ix(:,1)) = (2*((obj.nNodes-1) * dR_calc(1,3,3) - sum(dR_calc(2:obj.nNodes,3,3))) .*  dR_calc_dq_all(1,9,:)) / obj.nNodes;

        for iNode = 2:numNodes
            output(2,ix(:,iNode)) = (2*(- dR_calc_sum(iNode,1,1)) .*  dR_calc_dq_all(iNode,1,:)) / obj.nNodes;
            output(3,ix(:,iNode)) = (2*(- dR_calc_sum(iNode,2,1)) .*  dR_calc_dq_all(iNode,2,:)) / obj.nNodes;
            output(4,ix(:,iNode)) = (2*(- dR_calc_sum(iNode,3,1)) .*  dR_calc_dq_all(iNode,3,:)) / obj.nNodes;
            output(5,ix(:,iNode)) = (2*(- dR_calc_sum(iNode,1,2)) .*  dR_calc_dq_all(iNode,4,:)) / obj.nNodes;
            output(6,ix(:,iNode)) = (2*(- dR_calc_sum(iNode,2,2)) .*  dR_calc_dq_all(iNode,5,:)) / obj.nNodes;
            output(7,ix(:,iNode)) = (2*(- dR_calc_sum(iNode,3,2)) .*  dR_calc_dq_all(iNode,6,:)) / obj.nNodes;
            output(8,ix(:,iNode)) = (2*(- dR_calc_sum(iNode,1,3)) .*  dR_calc_dq_all(iNode,7,:)) / obj.nNodes;
            output(9,ix(:,iNode)) = (2*(- dR_calc_sum(iNode,2,3)) .*  dR_calc_dq_all(iNode,8,:)) / obj.nNodes;
            output(10,ix(:,iNode)) = (2*(- dR_calc_sum(iNode,3,3)) .*  dR_calc_dq_all(iNode,9,:)) / obj.nNodes;
        end
    else
        output = zeros(1,size(X,1));
    end

    output(1,ix(:,1)) =  (2*((obj.nNodes-1) * vec_rl(1) - sum(vec_rl(2:obj.nNodes))) .*  dvec_rl_dq(1,:)) / obj.nNodes;%
    for iNode = 2:numNodes
        output(1,ix(:,iNode)) = (2*(- dp_calc_sum(iNode)) .*  dvec_rl_dq(iNode,:)) / obj.nNodes;
    end

else
    error('Unknown option');
end


end