% ======================================================================
%> @file @Surf3d.m
%> @brief Matlab class describing the Surf3d model
%>
%> @author Alexander, Ton, Eva, Anne, Marlies
%> @date October 2024
% ======================================================================

%======================================================================
%> @brief The class describes the Surf3d (for Surfing)
%>
%> @details This class tilts the ground by a defined degree. Ground contact
%> is modelled through lift and drag forces acting on the surfboard that is simulated through an ellipsoidal trajectory around the rear foot.
%>
%> - The model is defined by the .osim file. Scaling must be done
%>   previously on the .osim file.
%> - Currently only the following models can be used
%>    - '3D Gait Model with Simple Arms'
%>    - '3D Gait Model with Simple Arms and Pelvis Rotation-Obliquity-Tilt Sequence'
%>    - 'gait24dof22musc'
%======================================================================
%======================================================================
classdef Surf3d < Gait3d % the class is a subclass of Gait3d

    properties(SetAccess = protected)
        scalebodyheight
        scalebodymass
    end

    properties (SetObservable, AbortSet)
        slope
        %> Double: Slope of ground in deg (positive is uphill)
    end

    methods
        %======================================================================
        %> @brief Constructor setting default Surf3d object
        %>
        %> @details
        %> Initializes the model and builds and initializes the mex function
        %> by calling the constructor Gait3d.Surf3d().
        %>
        %> @param   opensimfile         String: Opensim file with path, name and extension
        %> @param   scalebodyheight     Double: Bodyheight for scaling
        %> @param   scalebodymass       Double: Bodymass for scaling
        %> @retval  obj                 Surf3d class object
        %======================================================================
        function [obj] = Surf3d(opensimfile, scalebodyheight, scalebodymass)  

            %call Gait3d constructor opensimfile,
            obj = obj@Gait3d(opensimfile);
         
            if nargin > 1
                obj.bodyheight = scalebodyheight;
                if nargin > 2
                    obj.bodymass = scalebodymass;
                end
            end
        end
    end

    methods(Access = public)

        %======================================================================
        %> @brief Function to compute implicit differential equation for 3D musculoskeletal model
        %>
        %> @details
        %> This function calls the mex file of Gait3d.c:
        %> [f, dfdx, dfdxdot, dfdumus, dfdMextra]  = Gait3d('Dynamics',x,xdot,umus,Mextra);
        %>
        %> with the neural excitation umus and the extra torques Mextra.
        %>
        %> The dynamic residuals will be between fmin and fmax when inputs
        %> satisfy system dynamics: fmin <= f(x,dx/dt,umus,Mextra) <= fmax
        %>
        %> The last four outputs are optional and some computation time is saved if you do
        %> not request all of them.
        %>
        %> Compared to Gait3d.getDynamics we overwrite the dynamics of the
        %> ground contact using the equations of Dieter Heinrich
        %> Contact equations are removed and later added again through
        %> forceFootBoardEquilibrium constraint.
        %>
        %>
        %> @param   obj     Surf3d class object
        %> @param   x       Double array: State of the model (Surf3d.nStates x 1)
        %> @param   xdot    Double array: State derivatives (Surf3d.nStates x 1)
        %> @param   u       Double array: Controls of the model (Surf3d.nControls x 1)
        %>
        %> @retval  f       Double array: Dynamic residuals (Surf3d.nConstraints x 1)
        %> @retval	dfdx	(optional) Double matrix: Transpose of Jacobian matrix df/dx 		(Surf3d.nStates x Surf3d.nConstraints)
        %> @retval	dfdxdot	(optional) Double matrix: Transpose of Jacobian matrix df/dxdot 	(Surf3d.nStates x Surf3d.nConstraints)
        %> @retval	dfdu	(optional) Double matrix: Transpose of Jacobian matrix df/du        (Surf3d.nControls x Surf3d.nConstraints)

        %======================================================================
        function [f, dfdx, dfdxdot, dfdu] = getDynamics(obj,x,xdot,u)

            %Get Dynamics to overwrite
            [f, dfdx, dfdxdot, dfdu] = getDynamics@Gait3d(obj,x,xdot,u);
            
            isCP = strcmp(obj.constraints.type, 'CP');
            idxEqu1 = find(isCP & strcmp(obj.constraints.equation, 'equality 1')); % Fx
            idxEqu2 = find(isCP & strcmp(obj.constraints.equation, 'equality 2')); % Fy
            idxEqu3 = find(isCP & strcmp(obj.constraints.equation, 'equality 3')); % Fz
            idxEqu4 = find(isCP & strcmp(obj.constraints.equation, 'equality 4')); % xc
            idxEqu5 = find(isCP & strcmp(obj.constraints.equation, 'equality 5')); % yc
            idxEqu6 = find(isCP & strcmp(obj.constraints.equation, 'equality 6')); % zc

            % Initialize a logical mask with zeros values for all indices.
            filtered_f = f;
            filtered_dfdx = dfdx;
            filtered_dfdxdot = dfdxdot;
            filtered_dfdu = dfdu;

            % Update the mask for each set of indices to remove.
            filtered_f(idxEqu4) = 0;
            filtered_f(idxEqu5) = 0;
            filtered_f(idxEqu6) = 0;

            filtered_dfdx(:,idxEqu4) = 0;
            filtered_dfdx(:,idxEqu5) = 0;
            filtered_dfdx(:,idxEqu6) = 0;

            filtered_dfdu(:,idxEqu4) = 0;
            filtered_dfdu(:,idxEqu5) = 0;
            filtered_dfdu(:,idxEqu6) = 0;

            filtered_dfdxdot(:,idxEqu4) = 0;
            filtered_dfdxdot(:,idxEqu5) = 0;
            filtered_dfdxdot(:,idxEqu6) = 0;

            f = filtered_f;
            dfdx = filtered_dfdx;
            dfdu = filtered_dfdu;
            dfdxdot = filtered_dfdxdot;
        end


        %======================================================================
        %> @brief Function to compute deepest point on surf board and
        %> center of board depending on the horizontal and vertical axis
        %> an derivative by q from both
        %>
        %> @details
        %>
        %> Contact dynamics of surfing model depends on penetration depth
        %> of surf board with ground
        %>
        %>
        %>
        %> @param   obj               Surf3d class object
        %> @param   x                 Double array: State of the model (Surf3d.nStates x 1)
        %> @param                     rearFoot String: Rear foot of surfer that defines the
        %>                            foot around which the elliptical surf board is defined by
        %> @param   center_board_foot Double array: Center position of surf board relative to rear foot                       
        %> @param   long_axis_board   Double: Long axis radius of Ellipse aka Board length in frontal plane
        %> @param   short_axis_board  Double: Short axis radius of Ellipse aka Board length in sagittal plane
        %>
        %>
        %>
        %> @retval  point_s          Double array: Global location of deepest point on board           (Surf3d.nNodes x 3)
        %> @retval	dpoint_s_dq	     Double matrix: Derivative for Jacobian matrix dpoint_s/dq         (Surf3d.nDoFs x 3)
        %> @retval	center_s	     Double array: Global location of center of elliptical surf board  (Surf3d.nNodes x 3)
        %> @retval	dcenter_dq	     Double matrix: Derivative for Jacobian matrix dcenter_s/dq        (Surf3d.nDoFs x 3)
        %> @retval	sag_Board_axis	 Double array: Global location of center of elliptical surf board  (Surf3d.nNodes x 3)
        %> @retval	front_Board_axis Double matrix: Derivative for Jacobian matrix dcenter_s/dq        (Surf3d.nDoFs x 3)

        %======================================================================
        function [point_s, dpoint_s_dq, point_s_max, dpoint_s_max_dq, center_s, dcenter_dq, short_axis_board, long_axis_board] = getPointOnBoard(obj, x, rearFoot, center_board_foot, long_axis_board, short_axis_board)

            % MAKE PARAMS
            % as input argument
            calcR_no = find(strcmp(obj.segments.Properties.RowNames, 'calcn_r'));
            calcL_no = find(strcmp(obj.segments.Properties.RowNames, 'calcn_l'));
           
            % Calceneus position
            [FK,dFKdq]= obj.getFkin(x(1:obj.nDofs));
            p_calc_r(:) = FK((calcR_no-2)*12 + (1:3));
            dp_calc_r_dq(:,:) = dFKdq((calcR_no-2)*12 +(1:3),:);
            p_calc_l(:) = FK((calcL_no-2)*12 + (1:3));
            dp_calc_l_dq(:,:) = dFKdq((calcL_no-2)*12 +(1:3),:);
            % Calcaneus orientation
            R_calc_r = reshape(FK((calcR_no-2)*12 +(4:12)), 3, 3)';
            R_calc_l = reshape(FK((calcL_no-2)*12 +(4:12)), 3, 3)';
            dR_calc_rdq = dFKdq((calcR_no-2)*12 +(4:12),:);
            dR_calc_ldq = dFKdq((calcL_no-2)*12 +(4:12),:);
            
            % Differentiate between right and left rear foot
            if strcmp(rearFoot, 'left')                
                center_board_foot = [1, 1, -1] .* center_board_foot;
                R_calc = R_calc_l;
                dR_calc_dq = dR_calc_ldq;
                dp_calc_dq = dp_calc_l_dq;
                center_s = p_calc_l' + R_calc * center_board_foot';
            elseif strcmp(rearFoot, 'right')
                R_calc = R_calc_r;
                dp_calc_dq = dp_calc_r_dq;
                dR_calc_dq = dR_calc_rdq;

                center_s = p_calc_r' + R_calc * center_board_foot';
            else
                Error('Choose either right or left foot');
            end

            % Rotation matrix entries
            x_q = R_calc(2,3);
            y_q = R_calc(2,1);
            dx_dq = dR_calc_dq(6,:);
            dy_dq = dR_calc_dq(4,:);

            % Sign of all angles
            theta_x = atan2(R_calc(3,2), R_calc(3,3));
            theta_y = atan2(-R_calc(3,1), sqrt(R_calc(3,2)^2+R_calc(3,3)^2));
            theta_z = atan2(R_calc(2,1),R_calc(1,1));

            % Angle of deepest point around ellipse (= surf board)
            phi = atan((long_axis_board/short_axis_board) * (x_q / y_q)) + pi * ((theta_z + abs(theta_z)) / (2*theta_z));% or sigmoid function
            phi_max = phi + pi;
         
            sign_all = [sign(theta_x);-sign(theta_y);1];

            % Ellipse cooridnates
            base_ell = [short_axis_board*cos(phi); 0; long_axis_board*sin(phi)];
            base_ell_max = [short_axis_board*cos(phi_max); 0; long_axis_board*sin(phi_max)];
            
            % Deepest penetration point
            point_s = center_s + (R_calc * base_ell);
            % Opposite highest point
            point_s_max = center_s + (R_calc * base_ell_max);
          
            
            %Derivatives
            drcp_dq_x(:) = (long_axis_board^2 * x_q * (y_q * dx_dq - x_q*dy_dq)) / (short_axis_board*y_q^3 * (((long_axis_board^2*x_q^2)/(short_axis_board^2*y_q^2)) + 1)^(3/2));
            drcp_dq_z(:) = -(short_axis_board*long_axis_board^2*(y_q*dx_dq - x_q*dy_dq))/(sqrt(((long_axis_board^2*x_q^2)/(short_axis_board^2*y_q^2))+1)*(short_axis_board^2*y_q^2+long_axis_board^2*x_q^2));

            drcp_dq_x_max(:) = (short_axis_board*long_axis_board^2 * x_q * (x_q * dy_dq - y_q*dx_dq)) / (sqrt(((long_axis_board^2*x_q^2)/(short_axis_board^2*y_q^2))+1) * (short_axis_board^2*y_q^3 + long_axis_board^2*x_q^2*y_q));
            drcp_dq_z_max(:) = (short_axis_board*long_axis_board^2 * (y_q * dx_dq - x_q*dy_dq)) / (sqrt(((long_axis_board^2*x_q^2)/(short_axis_board^2*y_q^2))+1) * (short_axis_board^2*y_q^2 + long_axis_board^2*x_q^2));

            d_center_x_dq = dp_calc_dq(1,:) + dR_calc_dq(1,:).* center_board_foot(1) +  dR_calc_dq(2,:).* center_board_foot(2) +  dR_calc_dq(3,:).* center_board_foot(3);
            d_center_y_dq = dp_calc_dq(2,:) + dR_calc_dq(4,:).* center_board_foot(1) +  dR_calc_dq(5,:).* center_board_foot(2) +  dR_calc_dq(6,:).* center_board_foot(3);
            d_center_z_dq = dp_calc_dq(3,:) + dR_calc_dq(7,:).* center_board_foot(1) +  dR_calc_dq(8,:).* center_board_foot(2) +  dR_calc_dq(9,:).* center_board_foot(3);

            dcenter_dq = [d_center_x_dq; d_center_y_dq; d_center_z_dq];

            d_point_s_x_dq = d_center_x_dq + dR_calc_dq(1,:) * (short_axis_board * cos(phi)) + R_calc(1,1) * drcp_dq_x + dR_calc_dq(3,:) * (long_axis_board * sin(phi)) + R_calc(1,3) * drcp_dq_z;
            d_point_s_y_dq = d_center_y_dq + dR_calc_dq(4,:) * (short_axis_board * cos(phi)) + R_calc(2,1) * drcp_dq_x + dR_calc_dq(6,:) * (long_axis_board * sin(phi)) + R_calc(2,3) * drcp_dq_z;
            d_point_s_z_dq = d_center_z_dq + dR_calc_dq(7,:) * (short_axis_board * cos(phi)) + R_calc(3,1) * drcp_dq_x + dR_calc_dq(9,:) * (long_axis_board * sin(phi)) + R_calc(3,3) * drcp_dq_z;

            d_point_s_max_x_dq = d_center_x_dq + dR_calc_dq(1,:) * (short_axis_board * cos(phi_max)) + R_calc(1,1) * drcp_dq_x_max + dR_calc_dq(3,:) * (long_axis_board * sin(phi_max)) + R_calc(1,3) * drcp_dq_z_max;
            d_point_s_max_y_dq = d_center_y_dq + dR_calc_dq(4,:) * (short_axis_board * cos(phi_max)) + R_calc(2,1) * drcp_dq_x_max + dR_calc_dq(6,:) * (long_axis_board * sin(phi_max)) + R_calc(2,3) * drcp_dq_z_max;
            d_point_s_max_z_dq = d_center_z_dq + dR_calc_dq(7,:) * (short_axis_board * cos(phi_max)) + R_calc(3,1) * drcp_dq_x_max + dR_calc_dq(9,:) * (long_axis_board * sin(phi_max)) + R_calc(3,3) * drcp_dq_z_max;

            point_s = [point_s(1), point_s(2), point_s(3)]';
            dpoint_s_dq = [d_point_s_x_dq; d_point_s_y_dq; d_point_s_z_dq];
            dpoint_s_max_dq = [d_point_s_max_x_dq; d_point_s_max_y_dq; d_point_s_max_z_dq];
        end


         %======================================================================
        %> @brief Function to compute the forces acting on the surf board
        %>  depending on the orientation and penetration depth of the surf
        %>  board. The forces are lift and drag forces. 
        %>
        %> @details
        %>
        %> Board water contact forces depends on various factors and are computed from lift and drag forces. 
        %> Constants are defined in the beginning of the function
        %>
        %>
        %>
        %> @param   obj      Surf3d class object
        %> @param   x        Double array: State of the model (Surf3d.nStates x 1)
        %> @param   rearFoot String: Rear foot of surfer that defines the
        %>          foot around which the elliptical surf board is defined by
        %>
        %>
        %>
        %> @retval  point_s          Double array: Global location of deepest point on board           (Surf3d.nNodes x 3)
        %> @retval	dpoint_s_dq	     Double matrix: Derivative for Jacobian matrix dpoint_s/dq         (Surf3d.nDoFs x 3)
        %> @retval	center_s	     Double array: Global location of center of elliptical surf board  (Surf3d.nNodes x 3)
        %> @retval	dcenter_dq	     Double matrix: Derivative for Jacobian matrix dcenter_s/dq        (Surf3d.nDoFs x 3)
        %> @retval	sag_Board_axis	 Double array: Global location of center of elliptical surf board  (Surf3d.nNodes x 3)
        %> @retval	front_Board_axis Double matrix: Derivative for Jacobian matrix dcenter_s/dq        (Surf3d.nDoFs x 3)

        %======================================================================
        function [Fpp_x, Fpp_y, Fpp_z, dFpp_x_dq, dFpp_y_dq, dFpp_z_dq, dFpp_x_ddur, dFpp_y_ddur, dFpp_z_ddur] = getBoardForce(obj, x, x_prev, numNodes, rearFoot, center_board_foot, b, a, dur)
                    
            %% Drag Parameters
            % Floe rate of wave
            Q = -8.2;
            % Width of wave
            b_wave = 8;
            %Depth of wave
            d_wave = 0.48;
            % Drag coefficients
            cdx_R = 0.04;
            cdz_R = 0.04;
            cdy_R = 0.8;
            % Density of water
            rho = 1 * 1e+3;

            %inidces of calc segments
            calcR_no = find(strcmp(obj.segments.Properties.RowNames, 'calcn_r'));
            calcL_no = find(strcmp(obj.segments.Properties.RowNames, 'calcn_l'));

            % Bodyweight for normalizing forces
            bodyweight = obj.bodymass * -obj.gravity(2);
    
            % Deepest penetration point
            [point_s(:), d_point_s_dq(:,:), ~,~, center_s(:), d_center_dq(:,:) ,~, ~] = obj.getPointOnBoard(x, rearFoot, center_board_foot, b, a);
            [~,~,~,~, center_s_prev(:), d_center_dq_prev(:,:) ,~, ~] = obj.getPointOnBoard(x_prev, rearFoot, center_board_foot, b, a);

            % Find point of deepest penetration for standing in global
            % coordinates
            if numNodes == 1
                center_s_prev(:) = zeros(size(center_s(:),1), size(center_s(:),2));
                d_center_dq_prev(:,:) = zeros(size(d_center_dq(:,:),1), size(d_center_dq(:,:),2), size(d_center_dq(:,:),3));
            end

            % Foot orientation
            [FK,dFKdq]= obj.getFkin(x(1:obj.nDofs));
            p_calc_r(:) = FK((calcR_no-2)*12 + (1:3));
            dp_calc_r_dq(:,:) = dFKdq((calcR_no-2)*12 +(1:3),:);
            p_calc_l(:) = FK((calcL_no-2)*12 + (1:3));
            dp_calc_l_dq(:,:) = dFKdq((calcL_no-2)*12 +(1:3),:);

            R_calc_r = reshape(FK((calcR_no-2)*12 +(4:12)), 3, 3)';
            R_calc_l = reshape(FK((calcL_no-2)*12 +(4:12)), 3, 3)';
            dR_calc_rdq = dFKdq((calcR_no-2)*12 +(4:12),:);
            dR_calc_ldq = dFKdq((calcL_no-2)*12 +(4:12),:);

            % Orientation and position of rear foot
            if strcmp(rearFoot, 'left')
                center_board_foot = [1, 1, -1] .* center_board_foot;
                R_calc = R_calc_l;
                dR_calc_dq = dR_calc_ldq;
                dp_calc_dq = dp_calc_l_dq;
            elseif strcmp(rearFoot, 'right')
                R_calc = R_calc_r;
                dp_calc_dq = dp_calc_r_dq;
                dR_calc_dq = dR_calc_rdq;
            else
                Error('Choose either right or left foot');
            end

            % cd_R =  [cdx; cdy; cdz];
            % cdx_R(iNode) = cd_R(1);
            % cdy_R(iNode) = cd_R(2);
            % cdz_R(iNode) = cd_R(3);
            cdx_R_dq(:) = 0;%(dR_calc_dq(1,:) * cdx + dR_calc_dq(2,:) * cdy + dR_calc_dq(3, :) * cdz);% * sign(cd_R(1));
            cdy_R_dq(:) = 0;%(dR_calc_dq(4,:) * cdx + dR_calc_dq(5,:) * cdy + dR_calc_dq(6, :) * cdz);% * sign(cd_R(1));
            cdz_R_dq(:) = 0;%(dR_calc_dq(7,:) * cdx + dR_calc_dq(8,:) * cdy + dR_calc_dq(9, :) * cdz);% * sign(cd_R(3));

            %% Factors of forces
            % Get velocity of surfer
            if numNodes < 2
                v_surf = 0;
            else
                v_surf = (center_s(1, 3) - center_s_prev(1, 3)) / (dur / numNodes);
            end

            % Derivatives of veloctiy
            dv_surf_dq = d_center_dq(3,:) / (dur / numNodes);
            dv_surf_ddur = -(center_s(1,3) - center_s_prev(1,3)) * numNodes / (dur)^2;

            % Force in global x-direction
            Sw = 0.5 .* cdx_R  .* rho .* (Q / (b_wave * d_wave))^2;
            dSw_dq = 0.5 .* cdx_R_dq(:) .* rho .* (Q / (b_wave * d_wave))^2;

            % Force in global z-direction (in movement direction of surfer)
            Sz = 0.5 * cdz_R  * rho * v_surf^2;
            dSz_dq = 0.5 * cdz_R_dq(:) * rho *  v_surf^2 + 0.5 * cdz_R  * rho * 2 * v_surf * dv_surf_dq(:);
            dSz_ddur = 0.5 * cdz_R  * rho *  2 * v_surf * dv_surf_ddur;

            %% Geometrical parameters A and L
            % Translation matrix rotation indices

            Ryx = R_calc(2,1);
            Ryz = R_calc(2,3);
            dRyx_dq = dR_calc_dq(4,:);
            dRyz_dq = dR_calc_dq(6,:);

            % Angle of intersection points on elliptical board
            theta1 = -atan((Ryx * a) / (Ryz * b));
            theta2 = -atan((Ryx * a) / (Ryz * b)) + pi;
            % Derivative of theta1 with respect to each q
            dtheta1_dq = (a*b*(Ryx*dRyz_dq-Ryz*dRyx_dq)/((Ryx * a)^2 + (Ryz * b)^2));
            dtheta2_dq = (a*b*(Ryx*dRyz_dq-Ryz*dRyx_dq)/((Ryx * a)^2 + (Ryz * b)^2));

            %phi(iNode) = atan((b/a) * (Rzy / Rzz));

            %Normal of the board
            normal = [0, 1, 0] * R_calc;
            normNormal = norm(normal);
            normal = normal/normNormal;
            %normal_plot(iNode) = normal(1,2);

            %Angle between normal and y-axis - center_s(2) and derivative
            theta_normal = acos(normal(1,2));
            d_theta_normal_dq(:) = -(dR_calc_dq(5,:)/normNormal) .* (1/sqrt(1-normal(1,2)^2));

            % Over water level point_s equals zero
            point_s_exp(1,2) = - 0.5 * (sign(point_s(1,2)) - 1) * point_s(1,2);
            d_point_s_exp_dq(2,:) = -0.5 * (sign(point_s(1,2)) - 1) * d_point_s_dq(2,:);

            % Intersection points
            ell_intersec_1 =  R_calc * [a * cos(theta1); 0; b * sin(theta1)];
            ell_intersec_2 =  R_calc * [a * cos(theta2); 0; b * sin(theta2)];
            vert_Length_depth  = -point_s_exp(1,2) / sin(theta_normal);
            d_vert_Length_depth_dq(:,:) = csc(theta_normal) .* -(d_point_s_exp_dq(2,:) - point_s_exp(1,2) .* d_theta_normal_dq .* cot(theta_normal)); % 1xn vector

            % Total length calculation
            L_pen = sqrt((ell_intersec_1(1)  - ell_intersec_2(1))^2 + ...
                (ell_intersec_1(2) - ell_intersec_2(2))^2 + ...
                (ell_intersec_1(3) - ell_intersec_2(3))^2);

            % Initialize the derivative of total length with respect to q
            dL_pen_dq = zeros(1, obj.nDofs);

            % Derivative of intersection points with respect to each q
            for j = 1:obj.nDofs
                dR_calc_dq_j = reshape(dR_calc_dq(:, j), 3, 3)';

                dEllIntersec1_dq_j =  dR_calc_dq_j * [a * cos(theta1); 0; b * sin(theta1)] + ...
                    R_calc * [-a * sin(theta1) * dtheta1_dq(j); 0; b * cos(theta1) * dtheta1_dq(j)];

                % Derivative of ell_intersec_2 with respect to q(j)
                dEllIntersec2_dq_j =  dR_calc_dq_j * [a * cos(theta2); 0; b * sin(theta2)] + ...
                    R_calc * [-a * sin(theta2) * dtheta2_dq(j); 0; b * cos(theta2) * dtheta2_dq(j)];

                % Compute derivative of total_length with respect to q(j)
                diff_vector = ell_intersec_1 - ell_intersec_2;
                dL_pen_dq(j) = (1 / L_pen) * ...
                    (diff_vector' * (dEllIntersec1_dq_j - dEllIntersec2_dq_j));
            end


            %% Forces
            % Penetration force
            Fpp_p = 0.5 * rho * cdy_R  * vert_Length_depth  * L_pen * (Q / (b_wave * d_wave))^2;
            Fpp_s = Sw * vert_Length_depth  * L_pen;
            Fpp_l = Sz * vert_Length_depth  *  L_pen;

            % Derivatives
            dFpp_p_dq = 0.5 * rho * (Q / (b_wave * d_wave))^2 * (cdy_R_dq(:)' * L_pen * vert_Length_depth ...
                + cdy_R  * dL_pen_dq * vert_Length_depth  ...
                + cdy_R  * L_pen *  d_vert_Length_depth_dq(:,:));
            dFpp_s_dq =   (Sw * d_vert_Length_depth_dq * L_pen + dL_pen_dq * Sw * vert_Length_depth  + dSw_dq' * vert_Length_depth  * L_pen);
            dFpp_l_dq =   (Sz * d_vert_Length_depth_dq * L_pen + dL_pen_dq * Sz * vert_Length_depth  + dSz_dq' * vert_Length_depth  * L_pen);
            dFpp_l_ddur = (dSz_ddur * vert_Length_depth  * L_pen);

            % Forces dependet on rotation of rear foot calcaneus
            rot_Forces = R_calc * [Fpp_s; Fpp_p; Fpp_l]; % Flow of wave in negative x-direction

            %% Rotated foreces
            Fpp_y  = rot_Forces(2) / bodyweight;
            dFpp_y_dq(1,:) = (R_calc(2,1) * dFpp_s_dq(:)' + dR_calc_dq(4,:) * Fpp_s + R_calc(2,2) * dFpp_p_dq(:)' + dR_calc_dq(5,:) * Fpp_p + R_calc(2,3) * dFpp_l_dq(:)' + dR_calc_dq(6,:) * Fpp_l) / bodyweight;
            dFpp_y_ddur(1,:) = ( R_calc(2,3) * dFpp_l_ddur(:)') / bodyweight;

            Fpp_z  = rot_Forces(3) / bodyweight;
            dFpp_z_dq(1,:) = (R_calc(3,1) * dFpp_s_dq(:)' + dR_calc_dq(7,:) * Fpp_s + R_calc(3,2) * dFpp_p_dq(:)' + dR_calc_dq(8,:) * Fpp_p + R_calc(3,3) * dFpp_l_dq(:)' + dR_calc_dq(9,:) * Fpp_l) / bodyweight;
            dFpp_z_ddur(1,:) = (R_calc(3,3) * dFpp_l_ddur(:)') / bodyweight;

            Fpp_x  = rot_Forces(1)  / bodyweight;
            dFpp_x_dq(1,:) = (R_calc(1,1) * dFpp_s_dq(:)' + dR_calc_dq(1,:) * Fpp_s + R_calc(1,2) * dFpp_p_dq(:)' + dR_calc_dq(2,:) * Fpp_p + R_calc(1,3) * dFpp_l_dq(:)' + dR_calc_dq(3,:) * Fpp_l) / bodyweight;
            dFpp_x_ddur(1,:) = (R_calc(1,3) * dFpp_l_ddur(:)') / bodyweight;
        end
    end

    methods

        %======================================================================
        %> @brief Function to show model as stick figure
        %>
        %> @param   obj              Surf3d class object
        %> @param   x                Double matrice: State vector of model for n time points (Surf3d.nStates x n)
        %> @param   range            (optional) Double matrice: Defining the range of the figure
        %>                           with [xmin, xmax; ymin, ymax, zmin, zmax]. (3 x 2)
        %>                           (default: [pelvisX-1, pelvisX+1; -0.2, 2, pelvisZ-1, pelvisZ+1])
        %> @param   plotFeet         (optional) Integer: Defines if and how the feet are plotted (default: 0)
        %>                           - 0: Toe segments are not plotted,
        %>                           - 1: Toe segments are plotted using double the center of mass in x direction,
        %>                           - 2: Feet segments are plotted using the predefined position of the CPs. In future,
        %>                                we could think of using the simulated positions saved in the states.
        %> @param   plotGRF          (optional) Bool: If true, it plots arrows for the GRFs (default: 0)
        %> @param   plotCPs          (optional) Bool: If true, it plots spheres for the CPs (default: 0)
        %> @param   plotJointCOSYs   (optional) Bool: If true, it plots the coordinate systems of the joints (default: 0)
        %> @param   plotSurfBoard    (optional) Bool: If true, it plots the surfboard attachesd to rear calcaneus (default: 0)
        %> @param   az               (optional) Double: Azimuth for the view of the figure (view(az, el)). (default: 90)
        %> @param   el               (optional) Double: Elevation for the view of the figure (view(az, el)). (default: 0)
        %> @param   rearFoot          String: Rear foot of surfer that defines the
        %>                            foot around which the elliptical surf board is defined by
        %> @param   center_board_foot Double array: Center position of surf board relative to rear foot                       
        %> @param   long_axis_board   Double: Long axis radius of Ellipse aka Board length in frontal plane
        %> @param   short_axis_board  Double: Short axis radius of Ellipse aka Board length in sagittal plane
        %======================================================================
        function showStick(obj,x, range, plotFeet, plotGRF, plotCPs, plotJointCOSYs, plotSurfBoard, az, el, rearFoot, center_board_foot, long_axis_board, short_axis_board)
            if nargin < 9
                az = 125;
            end
            if nargin < 10
                el = 30;
            end
            if nargin < 11
                plotSurfBoard = 1;
                rearFoot = 'right';
                center_board_foot = [0.125, -0.05, -0.6];%[0.1, 0, 0.3]; % as INPUT argument
                long_axis_board = 0.81; %Long axis radius of Ellipse aka Board length in frontal plane
                short_axis_board = 0.245; %Short axis radius of Ellipse aka Board length in sagittal plane
            end

            % make coordinates for a 1 cm radius sphere
            [xs,ys,zs] = sphere(10);
            xs = xs*0.01;
            ys = ys*0.01;
            zs = zs*0.01;

            % define lines or polygons to draw the segments
            if ~strcmp(obj.osim.name, 'gait14dof22musc') && ~strcmp(obj.osim.name, 'gait14dof22musc and Pelvis Rotation-Obliquity-Tilt Sequence')
                % all our usual 3D models
                polygons = {   ...
                    [obj.joints.location(12,:);	obj.joints.location(7,:); obj.joints.location(2,:); obj.joints.location(12,:)]; ...	% pelvis
                    [0 0 0; obj.joints.t1_coefs(3,5) obj.joints.t2_coefs(3,5) 0.0]; ...	% right femur
                    [0 0 0;	obj.joints.location(4,:)]; ...										% right shank
                    []; ...																		% right talus
                    [0 0 0;	obj.joints.location(6,:)]; ...										% right calcaneus
                    []; ...																		% right toes
                    [0 0 0; obj.joints.t1_coefs(8,5) obj.joints.t2_coefs(8,5) 0.0]; ...	% left femur
                    [0 0 0;	obj.joints.location(9,:)]; ...										% left shank
                    []; ...																		% left talus
                    [0 0 0;	obj.joints.location(11,:)]; ...										% left calcaneus
                    []; ...																		% left toes
                    [0 0 0; obj.joints.location(13,:); obj.joints.location(17,:); 0 0 0]; ...	% torso
                    [0 0 0;	obj.joints.location(14,:)]; ...										% right humerus
                    [0 0 0;	obj.joints.location(15,:)]; ...										% right ulna
                    [0 0 0;	obj.joints.location(16,:)]; ...										% right radius
                    []; ...																		% right hand
                    [0 0 0;	obj.joints.location(18,:)]; ...										% left humerus
                    [0 0 0;	obj.joints.location(19,:)]; ...										% left ulna
                    [0 0 0;	obj.joints.location(20,:)]; ...										% left radius
                    []; ...																		% left hand
                    };
                colors_ploygons = {'k', 'r', 'r', 'r', 'r', 'r', 'b', 'b', 'b', 'b', 'b', 'k', 'r', 'r', 'r', 'r', 'b', 'b', 'b', 'b'};
            else
                % model used for NeuIPS challenge 2019
                polygons = {   ...
                    [obj.joints.location(12,:);	obj.joints.location(7,:); obj.joints.location(2,:); obj.joints.location(12,:)]; ...	% pelvis
                    [0 0 0; obj.joints.t1_coefs(3,5) obj.joints.t2_coefs(3,5) 0.0]; ...	        % right femur
                    [0 0 0;	obj.joints.location(4,:)]; ...										% right shank
                    [0 0 0;	obj.joints.location(5,:)]; ...		     							% right talus
                    [0 0 0;	obj.joints.location(6,:)]; ...										% right calcaneus
                    []; ...																		% right toes
                    [0 0 0; obj.joints.t1_coefs(8,5) obj.joints.t2_coefs(8,5) 0.0]; ...         % left femur
                    [0 0 0;	obj.joints.location(9,:)]; ...										% left shank
                    [0 0 0;	obj.joints.location(10,:)]; ...			     						% left talus
                    [0 0 0;	obj.joints.location(11,:)]; ...										% left calcaneus
                    []; ...																		% left toes
                    [0 0 0; obj.joints.location(13,:)]; ...	                                    % torso
                    [] ...																		% head
                    };
                colors_ploygons = {'k', 'r', 'r', 'r', 'r', 'r', 'b', 'b', 'b', 'b', 'b', 'k', 'k'};

            end

            % get range of figures
            if nargin > 2 && size(range, 1) == 3 && size(range, 2) == 2
                xrange = range(1, :);
                yrange = range(2, :);
                zrange = range(3, :);
            else
                idxPelvisX = obj.extractState('q', 'pelvis_tx');
                idxPelvisZ = obj.extractState('q', 'pelvis_tz');
                xrange = [min(x(idxPelvisX,:))-1  , max(x(idxPelvisX,:))+1];
                yrange = [-0.2, 2];
                zrange = [min(x(idxPelvisZ,:))-1  , max(x(idxPelvisZ,:))+1];
            end

            % plot the ground
            %             [xrange, groundY, zrange] = obj.applySlope(xrange, [0 , 0], zrange);
            %
            %             yrange = yrange + sort(groundY);
            %
            %             fill3([zrange(1), zrange(2), zrange(2), zrange(1)], ... % sidewards
            %                   [xrange(1), xrange(1), xrange(2), xrange(2)], ... % forwards
            %                   [groundY(1), groundY(1), groundY(2), groundY(2)], ... % upwards
            %                   [0.5, 0.5, 0.5], 'EdgeColor', [0.5, 0.5, 0.5]);



            % plot the ground
            fill3([zrange(1), zrange(2), zrange(2), zrange(1)], ... % sidewards
                [xrange(1), xrange(1), xrange(2), xrange(2)], ... % forwards
                [0, 0, 0, 0], ...                                 % upwards
                [0.5, 0.5, 0.5], 'EdgeColor', [0.5, 0.5, 0.5]);

            % plot the stick figure
            hold on;
            for iTime = 1 : size(x,2)
                % run the forward kinematics
                fk = obj.getFkin(x(1:obj.nDofs, iTime));
                fk = reshape(fk, 12, obj.nSegments-1);			% forward kinematics output from MEX function does not have ground
                % fk = obj.applySlope(fk);

                for i=1:obj.nSegments-1

                    O = fk(1:3,i);							% position of origin
                    R = reshape(fk(4:12,i)',3,3)';			% rotation matrix

                    % [O_x, O_y, O_z] = (O(1), O(2), O(3));obj.applySlope
                    % [R_x(:), R_y(:), R_z(:)] = obj.applySlope(R(1,1:3)', R(2,1:3)', R(3,1:3)');

                    % O = [O_x; O_y; O_z];
                    % R = [R_x(:), R_y(:), R_z(:)]';


                    % draw the XYZ axes of the 18 segments in red, green, blue
                    if nargin > 6 && plotJointCOSYs
                        axislength = 0.1;
                        X = O + axislength*R(:,1);				% end of X axis
                        Y = O + axislength*R(:,2);				% end of Y axis
                        Z = O + axislength*R(:,3);				% end of Z axis
                        % lot ZXY as XYZ so Y will be up in the Matlab Window
                        plot3([O(3) X(3)],[O(1) X(1)],[O(2) X(2)],'r','LineWidth',2)	% show X axis
                        plot3([O(3) Y(3)],[O(1) Y(1)],[O(2) Y(2)],'g','LineWidth',2)	% show Y axis
                        plot3([O(3) Z(3)],[O(1) Z(1)],[O(2) Z(2)],'b','LineWidth',2)	% show Z axis
                    end

                    % Get polygon coordinates
                    p = polygons{i}';

                    % Adapt polygon if feet should be plotted
                    if nargin > 3 && plotFeet == 1
                        % Toe segments are plotted using double the center of mass in x direction
                        if ismember(obj.segments.Properties.RowNames{i+1}, {'toes_r', 'toes_l'})
                            p = [0 0 0;	2*obj.segments.mass_center(i+1, 1),0,0]';
                        end

                    elseif nargin > 3 && plotFeet == 2
                        % Feet segments are plotted using the CPs
                        if ismember(obj.segments.Properties.RowNames{i+1}, {'calcn_r', 'toes_r', 'calcn_l', 'toes_l'})
                            p = obj.CPs.position(obj.CPs.segmentindex==i+1, :)';
                            p = p(:, [1, 3, 4, 2]); %> @todo: Order of CPs should not be hard coded!
                        end
                    end

                    % Draw each segment a black line or polygon
                    % ->transform polygon coordinates to global, and plot it
                    np = size(p,2);
                    if (np>0)
                        pg = repmat(O,1,np) + R*p;
                        fill3(pg(3,:),pg(1,:),pg(2,:),'k', 'edgecolor', colors_ploygons{i});
                    end

                end

                % get the ground reaction forces and plot them
                if nargin > 4 && plotGRF
                    grf = obj.getGRF(x(:, iTime));
                    [CoP_r, CoP_l] = obj.getCoP(grf);
                    Gait3d.plotgrf(grf(1:3), CoP_r, 'r');  % right side GRF vector
                    Gait3d.plotgrf(grf(7:9), CoP_l, 'b');  % left side GRF vector
                end

                % plot the contact points as spheres
                if nargin > 5 && plotCPs
                    for i = 1 : obj.nCPs
                        curName = obj.CPs.Properties.RowNames{i};
                        xc = xs + x(obj.extractState('xc', curName), iTime);
                        yc = ys + x(obj.extractState('yc', curName), iTime);
                        zc = zs + x(obj.extractState('zc', curName), iTime);
                        %[xc, yc, zc] = obj.applySlope(xc, yc, zc);
                        if endsWith(curName, '_r') % right leg
                            cpColor = [1, 0, 0]; % see colors of polygons
                        elseif endsWith(curName, '_l') % left leg
                            cpColor = [0, 0, 1]; % see colors of polygons
                        else % unknown segment
                            cpColor = [0, 0, 0]; % see colors of polygons for other segment
                        end
                        surf(zc,xc,yc,'FaceColor',cpColor, 'EdgeColor', 'none');
                    end
                end

                if nargin > 6 && plotSurfBoard
                    Ncontactpoints = obj.nCPs;

                    for iTime = 1 : size(x,2)

                        % Get deepest penetration point on board
                        [point_s, d_point_s_dq, point_s_max, dpoint_s_max_dq,center_s, d_center_dq, short_axis_board,long_axis_board] = getPointOnBoard(obj,  x(1:obj.nDofs, iTime), rearFoot, center_board_foot, long_axis_board, short_axis_board);

                        surf(zs+point_s(3), xs+point_s(1), ys+point_s(2), 'FaceColor', [0, 1, 0], 'EdgeColor','none');
                        surf(zs+center_s(3), xs+center_s(1), ys+center_s(2), 'FaceColor', [0, 0, 0], 'EdgeColor','none');

                        [FK]= obj.getFkin(x(1:obj.nDofs, iTime));

                        % Calceneus segment indices
                        calcR_no = find(strcmp(obj.segments.Properties.RowNames, 'calcn_r'));
                        calcL_no = find(strcmp(obj.segments.Properties.RowNames, 'calcn_l'));

                        p_calc_r(:) = FK((calcR_no-2)*12 + (1:3));
                        p_calc_l(:) = FK((calcL_no-2)*12 + (1:3));
                        R_calc_r = reshape(FK((calcR_no-2)*12 +(4:12)), 3, 3)';
                        R_calc_l = reshape(FK((calcL_no-2)*12 +(4:12)), 3, 3)';

                        U = R_calc_l(3,2);
                        V = R_calc_l(1,2);
                        W = R_calc_l(2,2);

                        quiver3(p_calc_l(1,3), p_calc_l(1,1), p_calc_l(1,2) ,U, V, W);

                        U = R_calc_r(3,2);
                        V = R_calc_r(1,2);
                        W = R_calc_r(2,2);

                        % Set center cepending on rear foot
                        if strcmp(rearFoot, 'left')
                            center_board_foot = [1, 1, -1] .* center_board_foot;
                            R_calc = R_calc_l;
                        elseif strcmp(rearFoot, 'right')
                            R_calc = R_calc_r;
                        else
                            Error('Choose either right or left foot');
                        end

                        t_board = 0:0.01:2*pi;
                        board_ell = [short_axis_board*cos(t_board);0*ones(size(t_board));long_axis_board*sin(t_board)];

                        % rotate the board
                        board_ell = R_calc*board_ell;
                        x_board = center_s(1)+board_ell(1,:);
                        y_board = center_s(2)+board_ell(2,:);
                        z_board = center_s(3)+board_ell(3,:);

                        board = [x_board; y_board; z_board];
                    end
                end

            end

            % finish the plot
            grid on; box on;
            xlabel('Z');
            ylabel('X');
            zlabel('Y');

            axis equal
            ylim(xrange);
            xlim(zrange);
            zlim(yrange);
            view(az, el);
            hold off;

        end
    end

    methods (Access = public)

        %> @cond DO_NOT_DOCUMENT
        %======================================================================
        %> @brief Helper function to apply the slopes to the coordinates for showStick()
        %>
        %> @param   obj    Surf3d class object
        %> @param   x      Double (vector/matrix): x coordinate before applying the slope
        %> @param   y      Double (vector/matrix): y coordinate before applying the slope
        %> @param   z      Double (vector/matrix): z coordinate before applying the slope
        %> @retval  xRot   Double (vector/matrix): x coordinate after applying the slope
        %> @retval  yRot   Double (vector/matrix): y coordinate after applying the slope
        %> @retval  zRot   Double (vector/matrix): z coordinate after applying the slope

        %======================================================================
        function [xRot, yRot, zRot] = applySlope(obj, x, y, z) %
            if isempty(obj.slope)
                xRot = x;
                yRot = y;
                zRot = z;
            else
                xRot = cosd(obj.slope(1))*x - sind(obj.slope(1))*y;
                yRot = sind(obj.slope(1))*x + cosd(obj.slope(1))*y;
                zRot = z;
            end
        end
   
    end

    methods(Static, Access = protected)
        % % % Declared function to read Opensim file
        model = readOsim(osimfile)

        %======================================================================
        %> @brief Function used in showStick() to plot GRFs
        %>
        %> @details
        %> Plots a ground reaction force vector represented by 3D force and
        %> CoP.
        %>
        %> @param   F      Double array: Force (3 x 1)
        %> @param   CoP    Double array: Center of pressure in x, y, z coordinates (3 x 1)
        %> @param   color  String: Color (e.g. color = 'r')
        %======================================================================
        function plotgrf(F,CoP,color)

            scale = 1.0;   % scale of the force vector visualization (meters per BW)
            x = [CoP(1) CoP(1)+scale*F(1)];
            y = [CoP(2) CoP(2)+scale*F(2)];
            z = [CoP(3) CoP(3)+scale*F(3)];
            plot3(z,x,y,color,'LineWidth',2);
        end
    end
end
