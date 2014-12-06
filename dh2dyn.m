function [D, C, V, J, H, A] = dh2dyn(dhtable,cstr,gravity)
%	DH2DYN
%
%	Generates the dynamics of an N-link kinematic chain in Brunovksy form, i.e.
%
%		D(q)*ddq + C(q,dq)*dq + V(q) = tau
%
%	Kimatic chain is modeled with inertia-less links ("floating" point masses).
%
%	Outputs:
%		
%		D 		N x N mass matrix
%
%		C		N x N correolis matrix
%
%		C		N x 1 gravity vector
%
%		J 		N x 1 cell of Jacobian matrices. J{end} contains the jacobian of 
% 				the full kinematic chain.
%
%		H 		(N+1) x 1 cell of homogeneous transformations to each i-th link
% 				from the base of the kinematic chain.
%
%		A 		N x 1 cell of homogeneous transformations from each link. Used
% 				to form each successive matrix in "H."
%
%	Inputs:
%
%		DHTABLE	
%				A symbolic table of Davenit-Hartenberg parameters describing the
% 				kinematic chain. The table should be percribed as a symbolic arraw
% 				like so:
%
%				dhtable = [ sym('[ a1, alpha2, d1, theta1 ]') ;
%							sym('[ a2, alpha2, d2, theta2 ]') ;
%							...
%							sym('[ aN, alphaN, dN, thetaN ]') ]
%
%				See usage example.
%
%
%		CSTR 	
%				A configuration string used to specify the type of joint between
% 				each link. Accepted configuration specifiers are as follows:
%
%			P	: specifies a prismatic joint
%				* 	requires that the i-th row of DHTABLE specifies a 'z-offset' parameter
%					dependent of time, i.e. 'd(t)'
%
%			R	: specifies a revolute joint
%				* 	requires that the i-th row of DHTABLE specifies a 'theta' parameter
%					dependent of time, i.e. 'q(t)'
%
%			F	: specifies a fixed joint
%				* 	no special requirements
%
%		GRAVITY
%				3x1 unit-vector whhich specifies the direction of gravity relative to the 
%				manipulator base-frame. This parameter can be symbolic for more generality.
%
%
%	Example:
%
%		confstr = 'RRRR';
%
%		table   = [   
%					sym('[L1,pi/2,  0,  q1(t)]'); ...
%					sym('[L2,   0,  0,  q2(t)]'); ...
%					sym('[L3,   0,  0,  q3(t)]'); ...
%					sym('[L4,   0,  0,  q4(t)]')  ...
%				];			
%
%		gravity = sym('[gx;gy;zy]');
%
%		[D, C, V, J, H, A] = dh2dyn(table,confstr,gravity);
%
%		disp(D)
%		disp(C)
%		disp(V)
%
%
% 	Brian Cairl
%	briancairl@gmail.com

    % Get n-Frames
    nframes     = numel(cstr);
    A           = cell(nframes,1);
    [dhy,dhx]   = size(dhtable);
    
    % Table Width Check
    if dhx ~= 4
       error( 'DHTABLE format is incorrect. Must have DHTABLE(idx,:) = [a,alpha,d,theta].'); 
    end
    
    % Configuration/Table Length Check
    if dhy ~= nframes
        error( 'Configuration string and n-frames do not match.'); 
    end
    
    
    % Configure Joint Vars
    qnot  = sym(zeros(1,2));
    for i = 2:nframes+1
        if cstr(i-1) == 'R'
            H{i,2}      = 'ROTATIONAL';
            q(i-1)      = dhtable(i-1,4);
        elseif cstr(i-1) == 'P'
            H{i-1,2}    = 'PRISMATIC';
            q(i-1)      = dhtable(i-1,3);
        end
        charize = char(q(i-1));
        if strcmpi(charize(end-2:end),'(t)')
            qnot(i-1) = sym(charize(1:end-3));
        else
            error('Joint variables must be a function of time. i.e, q(t).')
        end
    end
    
    
    % Generate A-Matrices
    for i = 1:nframes
        
        % Configuration Check
        if      (cstr(i) == 'R')&&(isnumeric(dhtable(i,4))||(dhtable(i,4)==0))
            error( ['Rotational Joint-', num2str(i), ' must be a non-zero sym-variable.'] ); 
        elseif  (cstr(i) == 'P')&&(isnumeric(dhtable(i,3))||(dhtable(i,3)==0))
            error( ['Prismatic Joint-', num2str(i), ' must be a non-zero sym-variable.'] ); 
        elseif (cstr(i) ~= 'R')&&(cstr(i) ~= 'P')
            error( ['Joints may be prismatic (P) or rotational (R). Remove erroneous configuration specifications at cstr(' num2str(i) ')'] )
        end
        
        % A Mat
        if cstr(i) == 'R'
            A{i}    = dh2sym(char(qnot(i)),char(dhtable(i,2)),char(dhtable(i,1)),char(dhtable(i,3)));
        else
            A{i}    = dh2sym(char(dhtable(i,4)),char(dhtable(i,2)),char(dhtable(i,1)),char(qnot(i)));
        end
           
    end
    
    
    % Generate H-Matrices
    wait_idx = 0;
    wait_len = nframes^2;
    wb = waitbar(0,'Generating Homogeneous Transformations...');
    
    H     = cell(nframes+1,1);
    H{1,1}= sym(eye(4));
    H{1,2}= 'FIXED';
    for i = 2:nframes+1
        waitbar(wait_idx/wait_len);
        
        if cstr(i-1) == 'R'
            H{i,2}      = 'ROTATIONAL';
            q(i-1)      = dhtable(i-1,4);
        elseif cstr(i-1) == 'P'
            H{i-1,2}    = 'PRISMATIC';
            q(i-1)      = dhtable(i-1,3);
        end
        charize = char(q(i-1));
        if strcmpi(charize(end-2:end),'(t)')
            qnot(i-1) = sym(charize(1:end-3));
        else
            error('Joint variables must be a function of time. i.e, f(t).')
        end
        H{i,1} = simplify(H{i-1}*A{i-1});
    end
    delete(wb);
    
    
    
    % Mass Values
    M = sym;
    for i = 1:nframes
        M(i) = ['m', num2str(i)];
    end
    
    
    
    
    % Generate the Jacobian
    J = cell(nframes,1);
    D = sym(zeros(nframes));
    
    wait_idx = 0;
    wait_len = nframes^2;
    wb = waitbar(0,'Generating Mass Matrix...');
    for j = 1:nframes
        waitbar(wait_idx/wait_len);
        wait_idx = wait_idx + 1;
        
        EF  = simplify((H{j+1,1}(1:3,4)-H{j,1}(1:3,4))/2 + H{j,1}(1:3,4));
        
        for i = 1:j
            Z    = H{i}(1:3,3);
            if cstr(i) == 'R'
                J{j}(1:3,i) = simplify(cross(Z,(EF-H{i,1}(1:3,4))));
                J{j}(4:6,i) = Z;
            else
                J{j}(1:3,i) = Z;
                J{j}(4:6,i) = sym('0');
            end
        end
        
        % Append Zeros for Dimension Consistency
        J{j} = [ J{j}, sym(zeros(6,nframes-i)) ];
        D    = D + M(j)*transpose(J{j})*J{j};
    end
    D = simplify(D);
    delete(wb)
    

    
    
    %Generate Correolis Symbols
    C    = sym(zeros(nframes,nframes));
    wait_idx = 0;
    wait_len = nframes^3;
    wb = waitbar(0,'Generating Correolis Matrix...');
    for k = 1:nframes
        for j = 1:nframes
            for i = 1:nframes
                waitbar(wait_idx/wait_len);
                wait_idx = wait_idx + 1;
                C(k,j)  = C(k,j)  + (diff(D(k,j),qnot(i)) + diff(D(k,i),qnot(j)) - diff(D(i,j),qnot(k)))*diff(q(i),sym('t'));
            end
            C(k,j) = simplify(C(k,j));
        end
    end
    delete(wb);
    


                
    %Generate Generate Symbols
    V    = sym(zeros(nframes,1));
    wait_idx = 0;
    wait_len = nframes;
    wb = waitbar(0,'Generating Gravity Vector...');    
    for i = 1:nframes
        wait_idx = wait_idx + 1;
        waitbar(wait_idx/wait_len);
        V(i,1) = M(i)*sym('g')*diff(abs(gravity)'*H{i+1,1}(1:3,4),qnot(i));
    end
    delete(wb)
    
    
    
end