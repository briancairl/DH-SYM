function HMATOUT = dhsym(tht,alp,a,d)
% 	DH2SYM
%
% 	Creates a homegeneous transformation matrix given standard
% 	Devenin-Havenart Parameters for kinematic chains.
%
%
% 	Inputs:
%
%   	THT     offset about parent z-axis
%   	ALP     offest about parent x-axis
%   	A       offset along x-axis
%   	D       offset along z-axis
%
%
% 	Outputs:
%
%   	HMATOUT homogenous transformation matrix
%
%
% 	Brian Cairl
%	briancairl@gmail.com

  
    if isnumeric(tht)
        ROTz = [ cos(tht), -sin(tht), 0, 0 ;
                 sin(tht),  cos(tht), 0, 0 ;
                 0,         0,        1, 0 ;
                 0,         0,        0, 1 ];
    elseif ischar(tht)
        Ct = sym([ 'cos(', tht, ')' ]);
        St = sym([ 'sin(', tht, ')' ]);
        
        ROTz = [ Ct, -St,  0,  0 ;
                 St,  Ct,  0,  0 ;
                 0,    0,  1,  0 ;
                 0,    0,  0,  1 ];
    else
        ROTz = 0;
    end
    %////////////////////////////////////////////
    if isnumeric(alp)
        ROTx = [ 1,        0,         0, 0 ;
                 0, cos(alp), -sin(alp), 0 ;
                 0, sin(alp),  cos(alp), 0 ;
                 0,         0,        0, 1 ];
    elseif ischar(alp)
        Ca = sym([ 'cos(', alp, ')' ]);
        Sa = sym([ 'sin(', alp, ')' ]);
        
        ROTx = [ 1,  0,   0, 0 ;
                 0, Ca, -Sa, 0 ;
                 0, Sa,  Ca, 0 ;
                 0,  0,   0, 1 ];
    else
        ROTx = 0;
    end
    %////////////////////////////////////////////
    if isnumeric(a)
        TRNa = [ 1,  0,  0,  a ;
                 0,  1,  0,  0 ;
                 0,  0,  1,  0 ;
                 0,  0,  0,  1 ];
    elseif ischar(a)
        A = sym( a );
 
        TRNa = [ 1,  0,  0,  A ;
                 0,  1,  0,  0 ;
                 0,  0,  1,  0 ;
                 0,  0,  0,  1 ];
    else
        TRNa = 0;
    end
    %////////////////////////////////////////////
    if isnumeric(d)
        TRNd = [ 1,  0,  0,  0 ;
                 0,  1,  0,  0 ;
                 0,  0,  1,  d ;
                 0,  0,  0,  1 ];
    elseif ischar(d)
        D = sym( d );
 
        TRNd = [ 1,  0,  0,  0 ;
                 0,  1,  0,  0 ;
                 0,  0,  1,  D ;
                 0,  0,  0,  1 ];
    else
        TRNd = 0;
    end
    %////////////////////////////////////////////
    if( isnumeric(TRNd)&&...
        isnumeric(TRNd)&&...
        isnumeric(TRNd)&&...
        isnumeric(TRNd) )
        HMATOUT = TRNd*ROTz*TRNa*ROTx;
    else
        HMATOUT = simplify(TRNd*ROTz*TRNa*ROTx,50);
    end
end