% Run script for fine scale saturation equation

% Initialize function handle for right hand side

gg = @(x,varargin)ones(size(x,1),1).*(x(:,2)>(1-2*h)).*(x(:,1)<2*h);
S = zeros(size(Mesh.Elements,1)/4,1);
coarse2fine;

% Assembling the advection matrix and load vector

Adv = assemMat_Adv_Sa3(Mesh,U_aux,Dofs,h);
%Lg = assemLoad_Vol_scalar(Mesh,h,QuadRule_2D,gg);
Lg = L;

% Time stepping (implicit
note=1;
clear sat
%==============================
    Pc = 0; Tt = 0;  pc = 32; Pc1 = 0;
    %==============================
for i = 1:Nt
    conv=0; IT=0; S00=S;
    
    while conv==0;
        ddt = dt/2^IT; % timestep
        I=0;
        while I<2^IT; % loop over sub?timesteps
            S0=S; dsn=1; it=0; I=I+1;
            while dsn>1e-3 && it<10
                K_w = S.^2/mu_w; %mobility
                K_o = (1-S).^2/mu_o;
                dK_w = 2*S/mu_w; %derivatives
                dK_o = 2*(S-1)/mu_o;
                df = dK_w./(K_w+K_o) - K_w./(K_w+K_o).^2.*(dK_w+dK_o); % df w/ds
                dR = speye(n^2/4) + (ddt/(2*h)^2)*Adv*spdiags(df,0,n^2/4,n^2/4); % R�(S)
                f = K_w./(K_w+K_o); % fractional flow
%               R = S - S0 + (ddt/(2*h)^2)*(Adv*f-Lg); % R(S)
                R = S - S0 + (ddt/(2*h)^2)*(Adv*f - min(Lg,0).*f - max(Lg,0));
                ds = -dR\R; % increment ds
                S = S+ds; % update S
                dsn = norm(ds); % norm of increment
                it = it+1; % number of N?R iterations
            end
            
            if dsn>1e-3
                I=2^IT; S=S00;
            end % check for convergence
        end
        
        if dsn<1e-3
            conv=1; % check for convergence
        else
            IT=IT+1;
        end % if not converged, decrease
    end
    Tt =[Tt, (1/500)*(i-1)*dt+I*ddt];
    Pc =[Pc; f(pc)];
    imagesc(reshape(S,32,32)); colorbar; drawnow;
    if mod(i,50)==0
        S2 = c2f*S;
        fine_update;
        Adv = assemMat_Adv_Sa3(Mesh,U_aux,Dofs,h);
    end
    
end

% Project the solution from the staggered DG space back to usual FEM space

S_fine = nodal_value_scalar(Mesh,S,h);

% Remove the saturation at the sink (not in the domain of consideration)

S_fine(1/h*(1/h-1)+1)=0;
S_fine(1/h*(1/h-2)+1)=0;
S_fine(1/h*(1/h-1)+2)=0;
S_fine(1/h*(1/h-2)+2)=0;

