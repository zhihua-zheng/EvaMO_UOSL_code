function Q_inq = q_in_SFWB(fh,axm,r2o,xi_inq)
%
% q_in_SFWB
%==========================================================================
%
% USAGE:
%  Q_inq = q_in_SFWB(r2o,xi_inq,fh,axm)
%
% DESCRIPTION:
%  Compute a proxy of turbulent kinetic energy (q^3) in Second Moment 
%   Closure including surface TKE flux from surface wave breaking 
%   (Craig and Banner 1994 + buoyancy flux) by solving the differential
%   equation for TKE.
%
% INPUT:
%
%  r2o - ratio of roughness length to Obukhov length, z_0/L
%  xi  - Normalized vertical distance, |z|/z_0
%  fh  - Figure handle
%  axm - Axis handle
%
% OUTPUT:
%
%  Q_inq - q^3 normalized by u_*^3
%
% AUTHOR:
%  June 23 2020, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

%% Single constants

A1 = 0.92;
A2 = 0.74;
B1 = 16.6;
B2 = 10.1;
C1 = 0.08;
C2 = 0.7;
C3 = 0.2;

Sq = 0.2;
Sm = 0.39;

alpha_b = 100; % breaking energy factor
kappa   = 0.4;
mxi = 60;

%% Combined constants

a = 3/B1/Sq;
b = 3/Sq*r2o;
c = 1/A1/Sq;

R2 = 1/3 - 2*A1/B1;
R1 = R2 - C1;

Z2 = (6*A1 + B2*(1-C3))*r2o;
Z1 = (6*A1 + 3*A2*(1-C2))*r2o;

bef = 3*alpha_b/Sq/2; % breaking energy flux

%% Computation

Y     = -log(mxi)/kappa;
ymesh = linspace(0,Y,101);
solinit = bvpinit(ymesh,@guess);

% Q at lower boundary
Zeta = r2o*mxi;
log_Qbal = @(x) a*x + b*mxi*( 1 - ...
                x^(2/3)/(Z2*mxi-R2*x)/(Z1*mxi-R1*x) ) + ...
                c*x^(2/3)/(Z1*mxi-R1*x);            
[~,pM]    = get_emp_phi(Zeta,'Kansas');
QYlog_emp = B1*(pM - Zeta);
QYlog     = fzero(log_Qbal,QYlog_emp);

% dQdy at lower boundary
log_dQdybal = @(x) (a*x - b*kappa*mxi)*(Z1*mxi-R1*QYlog)*(Z2*mxi-R2*QYlog) + ...
                   (a*QYlog + b*mxi)*(-kappa*Z1*mxi-R1*x)*(Z2*mxi-R2*QYlog) + ...
                   (a*QYlog + b*mxi)*(Z1*mxi-R1*QYlog)*(-kappa*Z2*mxi-R2*x) - ...
                   kappa*(c*Z2-b)*mxi*c + ...
                   2/3*(c*Z2-b)*mxi*x*QYlog^(-1/3) - ...
                   5/3*c*R2*x*QYlog^(2/3);
dQdylog_emp = B1*Zeta*kappa*( 1 - 15/4*(1-15*Zeta) );
dQdyYlog = fzero(log_dQdybal,dQdylog_emp);


% numerical solution
sol = bvp4c(@bvpfcn,@bcfcn,solinit);
xi_mesh = exp(-kappa*sol.x); % |z|/z_0

y_inq = -log(xi_inq)/kappa;
if ~isempty(y_inq)
    if xi_inq <= mxi
        Q_inq = deval(sol,y_inq,1);
    else % extend solution to new interval
        Ynew    = floor(y_inq);
        solinit = bvpinit(sol,[0 Ynew]); 
        sol     = bvp4c(@bvpfcn,@bcfcn,solinit);
        Q_inq   = deval(sol,y_inq,1);
    end
end

%% Visualization of solutions

if ~isempty(fh)
    
j = 0;
Qlog = nan(size(xi_mesh));
for xi = xi_mesh
    j = j+1;
    log_Qbal = @(x) a*x + b*xi*( 1 - ...
                    x^(2/3)/(Z2*xi-R2*x)/(Z1*xi-R1*x) ) + ...
                    c*x^(2/3)/(Z1*xi-R1*x);
    [~,pM]   = get_emp_phi(xi*r2o,'Kansas');
    Qlog_emp = B1*(pM - xi*r2o);
    Qlog(j)  = fzero(log_Qbal,Qlog_emp);
end
apr_sol = Qlog + alpha_b/2*sqrt(3*B1/Sq)*exp(sqrt(a)*sol.x);
C96_sol = (B1/Sm)^(3/4) + alpha_b/2*sqrt(3*B1/Sq)*exp(sqrt(a)*sol.x);

plot(axm,(sol.y(1,:)).^(1/3),xi_mesh,'linewidth',2)
hold(axm,'on'); grid(axm,'on')
plot(axm,apr_sol.^(1/3),xi_mesh,'color',rgb('amber'),'linewidth',2);
plot(axm,C96_sol.^(1/3),xi_mesh,'linestyle','--','color',rgb('blood'),...
         'linewidth',2)
plot(axm,Qlog.^(1/3),xi_mesh,'linestyle','--','color',[.5 .5 .5],...
         'linewidth',2)
set(axm,'ydir','reverse','ylim',[1 10],'xlim',[2 10],...
        'TickDir','out','FontSize',14)
text(axm,0.92,0.02,'(b)','Units','Normalized','FontSize',22,...
     'HorizontalAlignment','right','VerticalAlignment','bottom')

lgd = legend(axm,'numerical sol.','approximate sol.',...
                 'neutral sol.','no wave breaking','fontsize',16,...
                 'Position',[0.562 0.33 0.18 0.17],'AutoUpdate','off');
set(lgd.BoxFace,'ColorType','truecoloralpha',...
    'ColorData',uint8(255*[1; 1; 1; .9]))

% add axis for comparison of predicted profiles
% axmP   = axm.Position;
% axProf = axes('Parent',fh,'Position',axmP);
% linkaxes([axm axProf],'y')

% adjust label position
off_r = 0.02;
ylh1 = ylabel(axm,'|z|/z_0','fontsize',18);
xlh1 = xlabel(axm,'q/u_*','fontsize',18);
xlh1.Units = 'Normalized';
ylh1.Units = 'Normalized';
axm.Position = axm.Position + off_r*[1 1 0 -2];
set(xlh1,'Position',xlh1.Position + [0 -off_r 0]);

end

%% Subfunctions for BVP

% ODE
function Qode = bvpfcn(y,Q)
Qode = [Q(2)
        a*Q(1) + b*exp(-kappa*y)*( 1 - ...
        Q(1)^(2/3)/(Z2*exp(-kappa*y)-R2*Q(1))/(Z1*exp(-kappa*y)-R1*Q(1)) ) + ...
        c*Q(1)^(2/3)/(Z1*exp(-kappa*y)-R1*Q(1))];
end

% Boundary conditions
function res = bcfcn(Q0,QY)
res = [Q0(2) - bef
       QY(2) - dQdyYlog];
end

% Initial guess
function g = guess(y)
g = [(c/a/R1)^(3/4)
     0];
end

end