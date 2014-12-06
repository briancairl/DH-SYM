% Unit-test for dh2dyn function

clear; clc

config  = 'RRRR';

table   = [   
            sym('[L1,pi/2,  0,  q1(t)]'); ...
            sym('[L2,   0,  0,  q2(t)]'); ...
            sym('[L3,   0,  0,  q3(t)]'); ...
            sym('[L4,   0,  0,  q4(t)]')  ...
        ];
    
size(table)

gravity = sym('[gx;gy;zy]');

[D, C, V, J, H, A] = dh2dyn(table,config,gravity);
