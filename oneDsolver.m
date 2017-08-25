%% 1d Scalar Transport
% by Eray inanc
% define 2 profiles (sin/top-hat) and transport phi
% 4 numerical scheme for convective fluxes are considered, easily extendable
% Default settings: CFD=0.5, D=1e-5, dx=2e-2;
clear;clc;clf;

%% parameters run
u           = 1;            % axial velocity
CFL         = 0.5;          % timestep condition
tmax        = 1.0;          % max val of t
D           = 1e-5;         % viscosity

%% parameters net 
Ima         = 100;          % axial points
dx          = 2e-2;         % spacing between points
nG          = 1;            % ghost cells (depending on the scheme)

%% case 
profile     = 2;            % 1: sin / 2: top-hat 

%% parameters scheme
scheme      = 3;            % 0:cds / 1:uds / 2:dds / 3:udcds 
fac         = 0.5;          % udcds weighting (0 to cds - 1 to ud)

%% boundary conditions
periBC      = true;         % use periodic BC?

%% initialise and create class of Phi elements
x = linspace(0,2*pi,Ima+2*nG);
switch(profile)
    case 1
        phi = 0.5.*sin(x)+0.5;
    case 2
        phi = x.*0;
        phi(Ima/2-round(Ima/10,1):Ima/2+round(Ima/10,1)) = 1;
        exact = phi;
end
% predicted scalar
phiP = phi;

%% naming convention
% Ifim - Ifi - Ifip - .......- Ilam - Ila - Ilap (Ifim/Ilap ghost cells)
Ifim=1; Ifi=nG+1; Ifip=Ifi+1; Ilam=Ima; Ila=Ima+nG; Ilap=Ila+nG;

%% time integration, euler explicit
% stop button to end the loop
h = uicontrol('Style', 'PushButton', 'String', 'Stop', ...
              'Callback', 'delete(gcbo)');       

t=0;dt=CFL*dx/max(u); 
for n=1:round(tmax/dt);  
    %calculate the fluxes with scheme
    switch(scheme)
        case 0 % CDS
            fluxCon = 0.5*u*(phi(Ifim:Ila)+phi(Ifi:Ilap));
        case 1 % UDS
            fluxCon = u*phi(Ifim:Ila);
        case 2 % DDS
            fluxCon = u*phi(Ifi:Ilap);
        case 3 % UDCDS
            ww = fac*(u*dt/dx+1); % weight factor
            fluxCon = u*(ww.*phi(Ifim:Ila) + (1.0-ww).*phi(Ifi:Ilap));
    end
    fluxDif = D/dx.*(phi(Ifi:Ilap) - phi(Ifim:Ila));

    phiP(Ifi:Ila) = phiP(Ifi:Ila) - dt/dx*( (fluxCon(Ifi:Ila)-fluxCon(Ifim:Ilam)) ...
                                           - fluxDif(Ifi:Ila)-fluxDif(Ifim:Ilam));
                    
    % peridoic boundaries
    if periBC;
        phiP(Ifim) = phiP(Ila);
        phiP(Ilap) = phiP(Ifi);
    end
    
    % update t and phi
    t = t + dt;
    phi = phiP;
    
    % post-process
    plot(x,phiP,'-r','LineWidth',2); % numerical solution
    hold on 
    if profile==1; exact = 0.5.*sin(x-u*t)+0.5; end; % exact solution
    plot(x,exact,':k','LineWidth',2); 
    hold off;    
    xlabel('x'); ylabel('\phi'); title({['\itt = ',num2str(n*dt),'ms']});
    legend('computed','exact','Location','southwest'); legend('boxoff');
    axis([0 2*pi 0 1.05]); axis square;
    
    % Redraw
    shg;
    
    % Figure position
    bx = gcf; bx.Color = [1 1 1]; bx.Resize = 'off'; bx.ToolBar = 'none'; bx.MenuBar = 'none';

    % Pause in b/w iterations to see the evolvement (for fast cpus you need this)
    pause(dt/6400);
        
    % check the stop button
    drawnow
    if ~ishandle(h)
        break;
    end
    
    % printout some stuff
    display(['nt: ',num2str(n),' | t: ',num2str(t)]);
end;
