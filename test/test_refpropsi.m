clear all ;
clc ;
addpath('../refpropsi') ;

d = 1.807825145974120e+01 ;
u = 2.889369831786341e+05 ;
h = 3.000000000000000e+05 ;
p = 2.000000000000000e+05 ;
q = 5.504434242162227e-01 ;
s = 1.381340097152742e+03 ;
t = 2.630737275397669e+02 ;

dv = 7.681974423050868e+00 ;
uv = 3.682496177194109e+01 ;
hv = 4.500000000000000e+05 ;
pv = 2.000000000000000e+05 ;
qv = 9.980000000000000e+02 ;
sv = 1.927824051729751e+03 ;
tv = 3.288562284801340e+02 ;




fprintf('\n\n========== General properties ==========')

fprintf('\n\n=== DU ===\n')
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'D', d, 'U', u, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'D', d, 'U', u, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;


fprintf('\n=== DH ===\n')
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'D', d, 'H', h, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'D', d, 'H', h, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;


fprintf('\n=== DS ===\n')
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'D', d, 'S', s, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'D', d, 'S', s, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;


% fprintf('\n=== US ===\n')
% [D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'U', u, 'S', s, 'R134a') ;
% fprintf('refpropsi: %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;
% 
% [D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'S', s, 'U', u, 'R134a') ;
% fprintf('refpropm:  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;


fprintf('\n=== HS ===\n')
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'H', h, 'S', s, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'H', h, 'S', s, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;


fprintf('\n=== PD ===\n')
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'P', p, 'D', d, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'P', p/1000, 'D', d, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;


fprintf('\n=== PU ===\n')
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'P', p, 'U', u, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'P', p/1000, 'U', u, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;


fprintf('\n=== PH ===\n')
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'P', p, 'H', h, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'P', p/1000, 'H', h, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;


fprintf('\n=== PQ ===\n')
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'P', p, 'Q', q, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'P', p/1000, 'Q', q, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;


fprintf('\n=== PS ===\n')
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'P', p, 'S', s, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'P', p/1000, 'S', s, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;


fprintf('\n=== TD ===\n')
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'T', t, 'D', d, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'T', t, 'D', d, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;


fprintf('\n=== TU ===\n')
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'T', t, 'U', u, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'T', t, 'U', u, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;


fprintf('\n=== TH ===\n')
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'T', t, 'H', h, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'T', t, 'H', h, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;


fprintf('\n=== TP ===\n')
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'T', t, 'P', p, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'T', t, 'P', p/1000, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;


fprintf('\n=== TQ ===\n')
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'T', t, 'Q', q, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'T', t, 'Q', q, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;


fprintf('\n=== TS ===\n')
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'T', t, 'S', s, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'T', t, 'S', s, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;



fprintf('\n\n========== Transport properties ==========')

fprintf('\n=== TS ===\n')
[V,L,Pr] = refpropsi('VL^', 'T', tv, 'S', sv, 'R134a') ;
fprintf('refpropsi (VL^): %g %g %g\n', V, L, Pr) ;

[V,L,Pr] = refpropm('VL^', 'T', tv, 'S', sv, 'R134a') ;
fprintf('refpropm (VL^):  %g %g %g\n', V, L, Pr) ;

fprintf('\n=== PH ===\n')
[V,L,Pr] = refpropsi('VL^', 'P', pv, 'H', hv, 'R134a') ;
fprintf('refpropsi (VL^): %g %g %g\n', V, L, Pr) ;

[V,L,Pr] = refpropm('VL^', 'P', pv/1000, 'H', hv, 'R134a') ;
fprintf('refpropm (VL^):  %g %g %g\n', V, L, Pr) ;



fprintf('\n\n========== Critical point ==========\n')

[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'C', 0, ' ', 0, 'R134a') ;
fprintf('refpropsi (DUHPQST): %g %g %g %g %g %g %g\n', D, U, H, P, Q, S, T) ;

[D,U,H,P,Q,S,T] = refpropm('DUHPQST', 'C', 0, ' ', 0, 'R134a') ;
fprintf('refpropm (DUHPQST):  %g %g %g %g %g %g %g\n', D, U, H, P*1000, Q, S, T) ;



fprintf('\n\n========== Heat capacities ==========')

fprintf('\n=== PH ===\n')
[C,O,K] = refpropsi('COK', 'P', pv, 'H', hv, 'R134a') ;
fprintf('refpropsi (COK): %g %g %g\n', C,O,K) ;
[C,O,K] = refpropm('COK', 'P', pv/1000, 'H', hv, 'R134a') ;
fprintf('refpropm (COK):  %g %g %g\n', C,O,K) ;



fprintf('\n\n========== Surface tension ==========')
fprintf('\n=== PH ===\n')
[I] = refpropsi('I', 'P', p, 'H', h, 'R134a') ;
fprintf('refpropsi (I): %g\n', I) ;

[I] = refpropm('I', 'P', p/1000, 'H', h, 'R134a') ;
fprintf('refpropm (I):  %g\n', I) ;



fprintf('\n\n========== Molar mass ==========')
fprintf('\n=== PH ===\n')
[M] = refpropsi('M', 'P', p, 'H', h, 'R134a') ;
fprintf('refpropsi (M): %g\n', M) ;

[M] = refpropm('M', 'P', p/1000, 'H', h, 'R134a') ;
fprintf('refpropm (M):  %g\n', M/1000) ;



fprintf('\n\n========== Compressibilty factor ==========')
fprintf('\n=== PH ===\n')
[Z] = refpropsi('Z', 'P', p, 'H', h, 'R134a') ;
fprintf('refpropsi (Z): %g\n', Z) ;

[Z] = refpropm('Z', 'P', p/1000, 'H', h, 'R134a') ;
fprintf('refpropm (Z):  %g\n', Z) ;



fprintf('\n\n========== Mix fluid ==========')

fprintf('\n=== TD ===\n')
[P,H,S] = refpropsi('PHS','T',t, 'D', d, {'R134a', 'R152a'}, [0.95, 0.05]) ;
fprintf('refpropsi (PHS): %g %g %g\n', P, H, S) ;

[P,H,S] = refpropm('PHS','T',t, 'D', d, 'R134a', 'R152a', [0.95, 0.05]) ;
fprintf('refpropm (PHS):  %g %g %g\n', P*1000, H, S) ;


fprintf('\n=== TD ===\n')
[P,H,S,V,L] = refpropsi('PHSVL','T', tv, 'D', dv, {'R134a', 'R152a'}, [0.95, 0.05]) ;
fprintf('refpropsi (PHSVL): %g %g %g %g %g\n', P, H, S, V, L) ;

[P,H,S,V,L] = refpropm('PHSVL','T',tv, 'D', dv, 'R134a', 'R152a', [0.95, 0.05]) ;
fprintf('refpropm (PHSVL):  %g %g %g %g %g\n', P*1000, H, S, V, L) ;

fprintf('\n=== PH ===\n')
[T,D,S,V,L] = refpropsi('TDSVL','P', pv, 'H', hv, {'R134a', 'R152a'}, [0.95, 0.05]) ;
fprintf('refpropsi (TDSVL): %g %g %g %g %g\n', T, D, S, V, L) ;

[T,D,S,V,L] = refpropm('TDSVL','P', pv/1000, 'H', hv, 'R134a', 'R152a', [0.95, 0.05]) ;
fprintf('refpropm (TDSVL):  %g %g %g %g %g\n', T, D, S, V, L) ;



fprintf('\n\n========== Speed test ==========')

fprintf('\n=== PH ===\n')
tic
for k = 1:1000
    [~,~,~,~,~,~,~] = refpropsi('DUHPQST', 'P', p, 'H', h, 'R134a') ;
end
t = toc ;
fprintf('refpropsi time: %f\n', t) ;

tic
for k = 1:1000
    [~,~,~,~,~,~,~] = refpropm('DUHPQST', 'P', p/1000, 'H', h, 'R134a') ;
end
t = toc ;
fprintf('refpropm time:  %f\n', t) ;

fprintf('\n=== HS ===\n')
tic
for k = 1:1000
    [~,~,~,~,~,~,~] = refpropsi('DUHPQST', 'H', h, 'S', s, 'R134a') ;
end
t = toc ;
fprintf('refpropsi time: %f\n', t) ;

tic
for k = 1:1000
    [~,~,~,~,~,~,~] = refpropm('DUHPQST', 'H',h, 'S', s, 'R134a') ;
end
t = toc ;
fprintf('refpropm time:  %f\n', t) ;




fprintf('\n\n========== Output type ==========\n')

fprintf('Default (variables):\n') ;
[D,U,H,P,Q,S,T] = refpropsi('DUHPQST', 'P', p, 'H', h, 'R134a') ;
fprintf('D=%g U=%g H=%g P=%g Q=%g S=%g T=%g \n', D,U,H,P,Q,S,T) ;
fprintf('\n') ;

fprintf('Array:\n') ;
A = refpropsi('DUHPQST', 'P', p, 'H', h, 'R134a', 1, 'array') ;
disp(A)

fprintf('Cell:\n') ;
C = refpropsi('DUHPQST', 'P', p, 'H', h, 'R134a', 1, 'cell') ;
disp(C)

fprintf('Struct:\n') ;
S = refpropsi('DUHPQST', 'P', p, 'H', h, 'R134a', 1, 'struct') ;
disp(S)
