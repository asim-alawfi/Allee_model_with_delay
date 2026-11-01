%% allee_make_symbolic

clear;
close all;

base = [pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
addpath([base,filesep,'ddebiftool'], ...
        [base,filesep,'ddebiftool_extra_psol'], ...
        [base,filesep,'ddebiftool_utilities'], ...
        [base,filesep,'ddebiftool_extra_rotsym'], ...
        [base,filesep,'ddebiftool_extra_nmfm'], ...
        [base,filesep,'ddebiftool_extra_symbolic'], ...
        [base,filesep,'ddebiftool_coco'], ...
        pwd);

format compact
format short g

%% Parameters
parnames = {'a','b','lambda','tau'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:}); %#ok<NASGU>

syms a b lambda tau real
par = [a,b,lambda,tau];

%% One state and one discrete delay
% x(1,1) = u(t), x(1,2) = u(t-tau)
x = sym('x',[1,2],'real');

%% Right-hand side
rhs = a*x(1,1) - x(1,1)^3 + lambda - b*x(1,2);
%% Generate directional derivatives as an independent fallback
[fstr,derivs]=dde_sym2funcs(rhs,x,par, ...
    'filename','sym_allee_model', ...
    'maxorder',5, ...
    'directional_derivative',true);


