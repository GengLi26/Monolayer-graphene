%% DOS of monolayer graphene 

clc
clear all

a=1.42e-10;  % Anstrom  the equation in energy band is atom distance (m)
gamma0=2.8; % eV  Energy dispersion fitting parameter
sigma=0.01;  % Gausian broadening 


spin=2; % Spin degeneracy 2
valey=1; % Valley degeneracy: when intergration is over the entire
% first brolluin zone, valley=1.
scale=spin*valey/(2*pi)^2;  % here I figured out the unit problem. 


Kx=2*pi/(a*3); %% kx of K point
Ky=2*pi/(a*3)/sqrt(3); %% ky of K point
points=500; %% mesh point

% define ky array
ky=linspace(-Kx/(sqrt(3)/2),Kx/(sqrt(3)/2),points); 
dky=abs(ky(2)-ky(1)); 

%% Photo energy
E_boundary=10;
E=linspace(0,E_boundary,points);


% Generate k points in FBZ 
for count=1:1:length(E)
    for count_ky=1:1:length(ky)
        
        % limit kx in the first Brolluin zone 
        if ky(count_ky)<-Ky|ky(count_ky)>Ky
            Kx_BZ=abs(abs(ky(count_ky))-Kx/(sqrt(3)/2))*sqrt(3);
        else
            Kx_BZ=Kx;
        end
        % define kx array. 
        kx=linspace(-Kx_BZ,Kx_BZ,points);
        dkx=abs(kx(2)-kx(1));
        %scatter(kx(count_ky,:),count_ky.*ones(1,length(kx(count_ky,:))));
        %hold on
        
        % Graphene energy dispersion. 
        Ek=sqrt(gamma0^2*(1+4*((cos(ky(count_ky)*sqrt(3)*a/2)).^2)+...
            4*cos(ky(count_ky)*sqrt(3)*a/2).*cos(kx*3*a/2)));
        
        % delta function replaced by Gaussian function. 
        delta=exp(-(E(count)-Ek).^2./(2*sigma^2))./(sqrt(2*pi*sigma^2));
        
        % DOS
        DOS_ky(count_ky)=sum(delta)*dkx;

    end
    DOS(count)=sum(DOS_ky).*dky./(points^2);
end

DOS_final= scale.*DOS./(1e9)^2; % DOS unit eV-1.nm-2. 

%% plot 
plot(DOS_final,E,'color','k','linewidth',2)
hold on 
plot(DOS_final,-E,'color','k','linewidth',2)

set(gca,'fontsize',28)
xlabel(['DOS (a.u.)'],'FontSize',28)
ylabel(['E (eV)'],'FontSize',28)
set(gcf,'Position',[500 300 800 600])
%xlim([-7.9,7.9])


set(gcf,'PaperOrientation','landscape')
print(gcf, 'DOS_monolayer graphene.pdf', '-dpdf','-r0','-bestfit')




