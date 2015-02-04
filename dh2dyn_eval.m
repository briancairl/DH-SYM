function [D_eval, C_eval, V_eval] = dh2dyn_eval(D, C, V, q, q_dot, m, L, g_vec)
% Take Robot Manipulator symbolic matrices D, C, and V generated from the
% DH-SYM library (dh2dyn function) with numerical parameters q, q_dot, m, l, and g.
% Returns the D, C, and V numeric matrices evaluated from the symbolic structures.
%   Outputs:
%		
%		D_eval	N x N mass matrix, numerically evaluated
%
%		C_eval	N x N coriolis matrix, numerically evaluated
%
%		V_eval	N x 1 gravity vector, numerically evaluated
%
%	Inputs:
%
%		D   	N x N symbolic mass matrix
%
%		C   	N x N symbolic coriolis matrix
%
%		V   	N x 1 symbolic gravity vector
%
%       q       N x 1 numeric vector of joint angles
%
%       q_dot   N x 1 numeric vector of joint velocities
%
%       m       N x 1 numeric vector of link masses
%
%       L       N x 1 numeric vector of link lengths
%
%       g_vec   3 x 1 numeric vector describing the direction and magnitude
%               of gravity with respect to the robot arm. Typically, 
%               g*[gx,gy,gz] where g is the magnitude of gravity and the
%               elements are a unit vector in the direction of the gravity
%               vector.
%
%
% Example:
% 
%   % To see how to configure dh2dyn correctly, see the help file for 
%   % dh2dyn
%
% % Get manipulator matrices from dh2dyn
% [D, C, V, J, H, A] = dh2dyn(table,config,gravity);
%
% % Set parameter vectors for unit two-link arm at rest(origin)
% m = [1,1]';
% L = [1,1]';
% g = 9.8*[0,1,0]'; % Normal gravity, pointing "down"
% q = [0,0]';
% q_dot = [0,0]';
%
% [D_eval, C_eval, V_eval] = dh2dyn_eval(D, C, V, q, q_dot, m, L, g);
%
% Griswald Brooks
% griswald.brooks@gmail.com

% Check that a column vector has been passed in
if  size(q,1) < size(q,2)           || ...
    size(q_dot,1) < size(q_dot,2)   || ...
    size(m,1) < size(m,2)           || ...
    size(L,1) < size(L,2)           || ...
    size(g_vec,1) < size(g_vec,2)
    
    error('myApp:argChk', 'Parameter vectors must be column vectors.')
    
% Check that the q and q_dot are of the same dimension
elseif size(q,1) ~= size(q_dot,1)
    error('myApp:argChk', 'Vector q and q_dot must be of the same dimension.')
% Make sure gravity looks right
elseif size(g_vec,1) ~= 3
    error('myApp:argChk', 'Gravity vector must have three elements.')
else
    % Convert Symbolic Matrix to String
    C_char = char(C);

    % Replace differentials with symbols
    for i = 1:size(q,1)
        % Replace strings 'diff(q1(t), t)' with 'dq1', etc.
        C_char = strrep(C_char, ['diff(q',num2str(i),'(t), t)'], ['dq',num2str(i)]);
    end

    % Convert back to symbolic matrix
    C_eval = sym(C_char);

    %%% Evaluate matrices %%%

    % Get gravity terms
    g = norm(g_vec);
    g_norm = g_vec/g;
    gx = g_norm(1);
    gy = g_norm(2);
    gz = g_norm(3);
    
    % Create new variables q1, q2, etc. and set them to q(1), q(2), etc.
    % dq1, dq2, etc. and set them to q_dot(1), q_dot(2), etc.
    % m1, m2, etc. and set them to m(1), m(2), etc.
    % L1, L2, etc. and set them to l(1), l(2), etc.
    for i = 1:size(q,1)
        eval([genvarname(['q',num2str(i)]),     '= q(i);']);
        eval([genvarname(['dq',num2str(i)]),    '= q_dot(i);']);
        eval([genvarname(['m',num2str(i)]),     '= m(i);']);
        eval([genvarname(['L',num2str(i)]),     '= L(i);']);
    end
    
    % Numerically evaluate matrices
    D_eval = double(subs(D));
    C_eval = double(subs(C_eval));
    V_eval = double(subs(V));
end