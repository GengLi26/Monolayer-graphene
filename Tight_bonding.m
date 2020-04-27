%% Tight bonding method calculation: band structure of bilayer graphene
% Ref: Rep. Prog. Phys. 76 (2013) 056503 (28pp)
clear all


a=1.42;  % a is nearest atom distance ( Angstrom)
t=2.8;  %% hopping parameter (eV)

% define k-line
Kx=2*pi/(a*3); %% x axis of K point.
points=1000; % Mesh points

% From K point to -K point
k_K_nK_x=[Kx:-Kx/points:-Kx];
k_K_nK_y=k_K_nK_x./sqrt(3);
% From -K point to -M point
k_nK_nM_x=-Kx.*ones(1,points+1);
k_nK_nM_y=[-Kx/sqrt(3):Kx/sqrt(3)/points:0];
% From -M point to M point
k_nM_M_x=[-Kx:Kx/points:Kx];
k_nM_M_y=k_nM_M_x.*0;
% From M point to K point
k_M_K_y=[0:Kx/sqrt(3)/points:Kx/sqrt(3)];
k_M_K_x=Kx.*ones(1,length(k_M_K_y));
% Generate kx and ky array 
Kx_array=[k_K_nK_x k_nK_nM_x  k_nM_M_x k_M_K_x];
Ky_array=[k_K_nK_y k_nK_nM_y k_nM_M_y k_M_K_y];
% plot kx and ky for checking the scan trace of k vector. 
%plot(Kx_array,Ky_array)

% x axis for ploting band structure
Kx_plot=[linspace(0,4,length(k_K_nK_x)) ...
         linspace(4,5,length(k_nK_nM_x))...
         linspace(5,5+2*sqrt(3),length(k_nM_M_x))...
         linspace(5+2*sqrt(3),6+2*sqrt(3),length(k_M_K_y))];


% Scan through k-line, solve eigenvalues for each k.
for count=1:1:length(Kx_array)
    
  kx=Kx_array(count); 
  ky=Ky_array(count);
  
  % Construct hamiltonian 
  fk=-t*exp(-i*kx*a)*(1+2*exp(i*3*kx*a/2)*cos(sqrt(3)/2*ky*a));
  Hamiltonian=[0 fk; fk' 0];
  
  % get eigen values of Hamiltonian
  Band(count,:)=sort(real(eig(Hamiltonian)),'ascend');

end 

%% Plot band structure

plot(Kx_plot,Band,'color','k','linewidth',2)
hold on 
set(gca,'fontsize',28)
set(gca,'XTick',[])
ylabel(['E (eV)'],'FontSize',28)
set(gcf,'Position',[500 300 800 600])
xlim([0,max(Kx_plot)])
% grid on 
% grid minor
y_l=10;
ylim([-y_l,y_l])

% Gamma points
x=2.*ones(1,points+1);
y=[-y_l:2*y_l./points:y_l];
plot(x,y,'LineWidth',1.5,'color',[0.2 0.2 0.2],'linestyle','--')
hold on 
% -K points
x=4.*ones(1,points+1);
y=[-y_l:2*y_l./points:y_l];
plot(x,y,'LineWidth',1.5,'color',[0.2 0.2 0.2],'linestyle','--')
hold on 
hold on 
% -M points
x=5.*ones(1,points+1);
y=[-y_l:2*y_l./points:y_l];
plot(x,y,'LineWidth',1.5,'color',[0.2 0.2 0.2],'linestyle','--')
hold on 
% Gamma points
x=(5+1*sqrt(3)).*ones(1,points+1);
y=[-y_l:2*y_l./points:y_l];
plot(x,y,'LineWidth',1.5,'color',[0.2 0.2 0.2],'linestyle','--')
hold on 
% M points
x=(5+2*sqrt(3)).*ones(1,points+1);
y=[-y_l:2*y_l./points:y_l];
plot(x,y,'LineWidth',1.5,'color',[0.2 0.2 0.2],'linestyle','--')
hold on 


% Create textbox
annotation('textbox',...
    [0.12 0.06 0.03 0.04],'String','K',...
    'LineWidth',4,...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation('textbox',...
    [0.28 0.06 0.03 0.04],'String','\Gamma',...
    'LineWidth',4,...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');


annotation('textbox',...
    [0.43 0.06 0.03 0.04],'String','-K',...
    'LineWidth',4,...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');


annotation('textbox',...
    [0.51 0.06 0.03 0.04],'String','-M',...
    'LineWidth',4,...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation('textbox',...
    [0.67 0.06 0.03 0.04],'String','\Gamma',...
    'LineWidth',4,...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation('textbox',...
    [0.8 0.06 0.03 0.04],'String','M',...
    'LineWidth',4,...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation('textbox',...
    [0.89 0.06 0.03 0.04],'String','K',...
    'LineWidth',4,...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% save figure as 'pdf'
set(gcf,'PaperOrientation','landscape')
print(gcf, 'PR 10band delta1_mag.pdf', '-dpdf','-r0','-bestfit')


