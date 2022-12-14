function results=fit_threshold_function(directory)
% Description: this function calculates the constrained fit of the density
% of the sponge as a function of the threshold value

% Input:  
% 1. directory (string): string that specifies the name of the folder with
% all input files, has to be inside of a subfolder in '2_Inputfiles'
% Output: 
% 1. results (struct): struct that contains all fit results, i.e. the
% density and the threshold values from the files, the fit function, the
% fitted paramters and the used tolerance for the constraints  

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

%% determine all files in the given folder
files=dir(sprintf('%s/*.mat',directory));
data=load(sprintf('%s/%s',directory,files.name));
density_vec = data.results(:,2);
threshold_vec = data.results(:,1);

%% sort values
[threshold_vec,inds_sort]=sort(threshold_vec);
density_vec=density_vec(inds_sort);

%% tolerance for inequality constraint
tol=1e-5;

%% define fit model 
zfunc_model=@(p,x) p(1)*tanh(p(2)*(x-p(3)))+...
                    p(4)*tanh(p(5)*(x-p(6)));

%% perform minimization
options = optimoptions('fmincon','MaxFunctionEvaluations',1e5,'MaxIterations',1e5,...
    'ConstraintTolerance',1e-12,'OptimalityTolerance',1e-12,'StepTolerance',1e-12);
nonlcon = @constraint_eq;
p_sol=fmincon(@(p)Fit_Data(density_vec,threshold_vec,zfunc_model,p),ones(1,6),...
    [],[],[],[],[],[],@(p)nonlcon(density_vec,zfunc_model,tol,p),options);

%% check results
zfunc_model(p_sol,0)
zfunc_model(p_sol,1)
min(zfunc_model(p_sol,density_vec))
max(zfunc_model(p_sol,density_vec))

%% store results in a single struct
results.density_vec=density_vec;
results.threshold_vec=threshold_vec;
results.zfunc_model=zfunc_model;
results.p_sol=p_sol;
results.tol=tol;

%% plot results
density_plot=linspace(0,1,1000);
figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
hold on
plot(density_plot,zfunc_model(p_sol,density_plot),'-b','LineWidth',1.5)
plot(density_vec,threshold_vec,'.r','MarkerSize',10)
xlabel('\rho_A')
ylabel('\phi_{T,1}')
title('Fit threshold function')
axis([0,1,-1.1 1.1])
grid on
drawnow
end

%% objective function
function Z=Fit_Data(x,y_soll,f_handle,p)
Z=sum((f_handle(p,x)-y_soll).^2);
end

%% constraints
function [c, ceq] = constraint_eq(x,f_handle,tol,p)
f1=f_handle(p,0);
f2=f_handle(p,1);

f_all=f_handle(p,x(2:end-1));

% function should be between -1+tol and 1-tol
c = abs(f_all)-(1-tol); 

% constraint f(0)=1 and f(1)=-1
ceq = [f1-1;f2+1];
end
