%% Graphene Hamiltonian at K point, Kubo formula. 
clear all

%ref 

% Select approiate phenomenological broadening
% and you should adjust k mesh points to get convergent results

a=1.42*1e-10;  % a is lattic constant (nearest atom distance)
% t=2.8;  %% hopping parameter (eV)
% define k-line
vf=1e6;  % Fermi velocity m/sec
e = 1.602177e-19; % electron charge (C)
hbar = 1.054572e-34; % reduced Planck constant ( J . s )
sigma0=e^2/hbar; % normalized to e^2/hbar

spin=2; % Spin degeneracy 2
valey=2; % Graphene at K points has valley degeneracy of 2;
% Here the integration is only over K points, so valley=2 is necessary. 
scale=spin*valey/(2*pi)^2;  % here I figured out the unit problem. 

points=150; % k mesh points 150 is OK ofr eta=40e-3*e.

Ef=0.2*e;% Fermi level of graphene unit J
E=[0.01:0.02:1.4].*e; % Photo Energy unit J
omega=E./hbar; % omega unit s^-1
eta=40e-3*e;  % phenomenological broadening

alpha=1:1:2;  % how many bands? Hamiltonian of graphene is 2*2 matrix

Kx=2*pi/(a*3); % kx of K points
Kx_array=[-Kx/5:2*Kx/5/points:Kx/5]; % define Kx array
dkx=Kx_array(2)-Kx_array(1);
Ky_array=[-Kx/5:2*Kx/5/points:Kx/5]; % define Ky array
dky=Ky_array(2)-Ky_array(1);

% scan
for omega_count=1:1:length(omega)
    
    for alpha_count=1:1:length(alpha)
        alpha_loop=alpha(alpha_count);
        beta=1:1:2;
        beta(alpha_loop)=[];
        
        for beta_count=1:1:length(beta)
            beta_loop=beta(beta_count);
            
            for kx_count=1:1:length(Kx_array)
                kx=Kx_array(kx_count);
                
                for ky_count=1:1:length(Ky_array)
                    ky=Ky_array(ky_count);
                    
                    %% core part of the code
                    [Vector,E_eig]=eig(Hamil_K(kx,ky));
                    Band=diag(E_eig);
                    vx=(Hamil_K(kx+Kx/10000,ky)-Hamil_K(kx,ky))./(Kx/10000)./hbar;
                    
                    Mean_alpha_beta=Vector(:,alpha_loop)'*vx*Vector(:,beta_loop);
                    Mean_beta_alpha=Vector(:,beta_loop)'*vx*Vector(:,alpha_loop);
                    % S not imputed for simplicity
                    cond_ky(ky_count)=e^2*hbar/i*(fermi(Band(alpha_loop),Ef)-fermi(Band(beta_loop),Ef))...
                        /(Band(alpha_loop)-Band(beta_loop))...
                        *(Mean_alpha_beta*Mean_beta_alpha)...
                        ./(Band(alpha_loop)-Band(beta_loop)+hbar*omega(omega_count)+i*eta);
                    %%
                end
                cond_kx_ky(kx_count)=sum(cond_ky).*dky;
            end
            cond_beta(beta_count)=sum(cond_kx_ky).*dkx;
        end
        cond_alpha_beta(alpha_count)=sum(cond_beta);
    end
    cond_omega(omega_count)=sum(cond_alpha_beta);
    
    omega_count./length(omega)*100
end

% % %% Plot band structure
subplot(2,1,1)
plot(E./e,real(cond_omega)./sigma0*scale,'color','k','linewidth',2)
hold on 
subplot(2,1,2)
plot(E./e,imag(cond_omega)./sigma0*scale,'color','k','linewidth',2)
hold on 

