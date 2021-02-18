
%% Initialize workspace
clear all; close all; clc; format long;

%% Fundamental physical constants
% Absolute vacuum permittivity
epsilon_0 = 8.854*1e-12;
% Absulute vacuum permeability
mu_0 = 4*pi*1e-7;
% Light speed in vacuum
c = 1/sqrt(epsilon_0*mu_0);
%% Main parameters
% Calculation area length per x and y axes
L = [4e-6, 4e-6];
delta=10e-9;
% Uniform grid points for x and y axes
nx = L(1)/delta; ny = L(1)/delta;

% End time [s]
t_end = 56619.5e-19;
% Excitation source amplitude and frequency [Hz]
E0 = 1.0;   
lambda=405e-9;
f = c/lambda;   
% Number of materials (include vacuum backrgound)
number_of_materials = 1;
% Background relative permittivity, relative permeability and absolute
% conductivity
eps_back = 1.0; mu_back = 1.0; sig_back = 0.0;
% Width of alinea between total field area and calculation area border
% (scattered field interface width) [m]
len_tf = L(1)/4;

% Grid calculation
X = linspace(0,L(1),nx)-L(1)/2;
Y = linspace(0,L(2),ny)-L(2)/2;
dx = X(nx)-X(nx-1);
dy = Y(ny)-Y(ny-1);
% Time step from CFL condition
dt = ( 1/c/sqrt( 1/(dx^2) + 1/(dy^2) ) )*0.99;
number_of_iterations = round(t_end/dt);

%% Geometry matrix 
Index = zeros(nx,ny);
IndexX = zeros(nx,ny);
IndexY = zeros(nx,ny);

%% Materials matrix 
Material = zeros(number_of_materials,3);
Material(1,1) = eps_back;
Material(1,2) = mu_back;
Material(1,3) = sig_back;

% %% Simple geometry (here metal round cylinder 1/2)
% % Diameter [m]
% d0 = 1.3;
% % Cylinder materials  
% Material(2,1) = 1;       % relative permittivity
% Material(2,2) = 1;       % relative permeability
% Material(2,3) = 1.0e+70; % absolute conductivity
% % Fill geometry matrix for 1/2 cylinder (no vectorization here), but for TE
% % mode we need three different Index for Ex, Ey and Hz due to leapfrog
% for I = 1:nx
%     for J = 1:ny
%         if ((J-nx/2)^2 + (I-ny/2)^2 <= 0.25*(d0/dx)^2 && (L(1)-I*dx<=J*dy))
%             Index(J,I) = 1;
%         end
%         if ((J-nx/2)^2 + (I+0.5-ny/2)^2 <= 0.25*(d0/dx)^2 && (L(1)-(I+0.5)*dx<=J*dy))
%             IndexX(J,I) = 1;
%         end
%         if ((J+0.5-nx/2)^2 + (I-ny/2)^2 <= 0.25*(d0/dx)^2 && (L(1)-I*dx<=(J+0.5)*dy))
%             IndexY(J,I) = 1;
%         end
%     end
% end

%% Calculate size of total field area in TF/SF formalism
%nx_a = round(len_tf/dx);
%nx_b = round((L(1)-len_tf)/dx);
nx_a=0;
ax=1;
nx_b=L(1)-0;
bx=ceil(nx_b/dx);
%ny_a = round(len_tf/dy);
%ny_b = round((L(2)-len_tf)/dy);
ny_a=0;
ay=1;
ny_b=L(2)-0;
by=ceil(ny_b/dy);
%% Allocate 2D arrays
% TE physical and auxiliary fields
% Fz = zeros(nx,ny); Tz = zeros(nx,ny);   Gx = zeros(nx,ny-1); Gy = zeros(nx-1,ny);
% Ez = zeros(nx,ny); Hx = zeros(nx,ny-1); Hy = zeros(nx-1,ny);
Ez=zeros(nx,ny); Hx=zeros(nx,ny); Hy=zeros(nx,ny);
m_hx=-1*dt/(Material(1,2)*mu_0)/dy*ones(nx,ny);
m_hy=-1*dt/(Material(1,2)*mu_0)/dx*ones(nx,ny);
m_ezhy=1*dt/(Material(1,1)*epsilon_0)/dx*ones(nx,ny);
m_ezhx=-1*dt/(Material(1,1)*epsilon_0)/dy*ones(nx,ny);
% % TM physical and auxiliary fields
% Wz = zeros(nx,ny); Mx = zeros(nx,ny-1); My = zeros(nx-1,ny);
% Hz = zeros(nx,ny); Ex = zeros(nx,ny-1); Ey = zeros(nx-1,ny);
Ex=zeros(nx,ny); Ey=zeros(nx,ny); Hz=zeros(nx,ny);
m_ex=1*dt/(Material(1,1)*epsilon_0)/dy*ones(nx,ny);
m_ey=-1*dt/(Material(1,1)*epsilon_0)/dx*ones(nx,ny);
m_hzey=-1*dt/(Material(1,2)*mu_0)/dx*ones(nx,ny);
m_hzex=1*dt/(Material(1,2)*mu_0)/dy*ones(nx,ny);
%%
% %try to plot this model
% rgb=zeros(nx,ny);
% rgb(nx_a:nx_b,ny_a:ny_b)=2;
% %heatmap(rgb)
% hmh = heatmap(rgb,'xlabel','\mu','ylabel','\nu');
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
%%update equation
movieVector(1:number_of_iterations) = struct('cdata',[],'colormap',[]);
f=figure('Color','w','Position',[100 100 800 800]);
%zbuffer or painter are recommended as renderer
set(gcf,'Renderer','zbuffer');


%  %Gaussian Source
%     f(n)= (-2*(n*dt-t0)*dt/(tw^2))*exp(-(n*dt-t0)^2/(tw^2))/dy;
%     Ez(srcx,srcy) = Ez(srcx,srcy) + f(n);
tw     = 16*dt;
source=zeros(number_of_iterations);
%%
for i=1:1:number_of_iterations
   %TE mode
   Hx(ax:bx,ay:by-1)=Hx(ax:bx,ay:by-1)+m_hx(ax:bx,ay:by-1).*(Ez(ax:bx,ay+1:by)-Ez(ax:bx,ay:by-1));
   Hx(ax:bx,by)=Hx(ax:bx,by)+m_hx(ax:bx,by).*(Ez(ax:bx,ay)-Ez(ax:bx,by));
   
   Hy(ax:bx-1,ay:by)=Hy(ax:bx-1,ay:by)+m_hy(ax:bx-1,ay:by).*(-Ez(ax+1:bx,ay:by)+Ez(ax:bx-1,ay:by));
   Hy(bx,ay:by)=Hy(bx,ay:by)+m_hy(bx,ay:by).*(-Ez(1,ay:by)+Ez(bx,ay:by));
   
   Ez(ax+1:bx,ay+1:by)=Ez(ax+1:bx,ay+1:by)+m_ezhy(ax+1:bx,ay+1:by).*(Hy(ax+1:bx,ay+1:by)-Hy(ax:bx-1,ay+1:by))+m_ezhx(ax+1:bx,ay+1:by).*(Hx(ax+1:bx,ay+1:by)-Hx(ax+1:bx,ay:by-1));
%   Ez(ax,ay+1:by)=Ez(ax,ay+1:by)+m_ezhy(ax,ay+1:by).*(hy(ax,ay+1:by)-hy(bx,ay+1:by))+m_ezhx(ax,ay+1:by).*(hx(ax,ay+1:by)-hx(ax,ay:by-1));
   Ez(ax,ay+1:by)=Ez(ax,ay+1:by)+m_ezhy(ax,ay+1:by).*(Hy(ax,ay+1:by)-Hy(bx,ay+1:by))+m_ezhx(ax,ay+1:by).*(Hx(ax,ay+1:by)-Hx(ax,ay:by-1));
   Ez(ax+1:bx,ay)=Ez(ax+1:bx,ay)+m_ezhy(ax+1:bx,ay).*(Hy(ax+1:bx,ay)-Hy(ax:bx-1,ay))+m_ezhx(ax+1:bx,ay).*(Hx(ax+1:bx,ay)-Hx(ax+1:bx,by));
   Ez(ax,ay)=Ez(ax,ay)+m_ezhy(ax,ay).*(Hy(ax,ay)-Hy(bx,ay))+m_ezhx(ax,ay).*(Hx(ax,ay)-Hx(ax,by));

   
   %TM mode
%    Hz(ax:bx-1,ay:by-1)=Hz(ax:bx-1,ay:by-1)+m_hzey(ax:bx-1,ay:by-1).*(Ey(ax+1:bx,ay:by-1)-Ey(ax:bx-1,ay:by-1))+m_hzex(ax:bx-1,ay:by-1).*(Ex(ax:bx-1,ay+1:by)-Ex(ax:bx-1,ay:by-1));
%    Hz(bx,ay:by-1)=Hz(bx,ay:by-1)+m_hzey(bx,ay:by-1).*(Ey(ax,ay:by-1)-Ey(bx,ay:by-1))+m_hzex(bx,ay:by-1).*(Ex(bx,ay+1:by)-Ex(bx,ay:by-1));
%    Hz(ax:bx-1,by)=Hz(ax:bx-1,by)+m_hzey(ax:bx-1,by).*(Ey(ax+1:bx,by)-Ey(ax:bx-1,by))+m_hzex(ax:bx-1,by).*(Ex(ax:bx-1,ay)-Ex(ax:bx-1,by));
%    Hz(bx,by)=Hz(bx,by)+m_hzey(bx,by).*(Ey(ax,by)-Ey(bx,by))+m_hzex(bx,by).*(Ex(bx,ay)-Ey(bx,by));
%    
%    Ex(ax:bx,ay+1:by)=Ex(ax:bx,ay+1:by)+m_ex(ax:bx,ay+1:by).*(Hz(ax:bx,ay+1:by)-Hz(ax:bx,ay:by-1));
%    Ex(ax:bx,ay)=Ex(ax:bx,ay)+m_ex(ax:bx,ay).*(Hz(ax:bx,ay)-Hz(ax:bx,by));
%    
%    Ey(ax+1:bx,ay:by)=Ey(ax+1:bx,ay:by)+m_ey(ax+1:bx,ay:by).*(Hz(ax+1:bx,ay:by)-Hz(ax:bx-1,ay:by));
%    Ey(ax,ay:by)=Ey(ax,ay:by)+m_ey(ax,ay:by).*(Hz(ax,ay:by)-Hz(bx,ay:by));
   Hz(ax+1:bx,ay+1:by)=Hz(ax+1:bx,ay+1:by)+m_hzey(ax+1:bx,ay+1:by).*(Ey(ax+1:bx,ay+1:by)-Ey(ax:bx-1,ay+1:by))+m_hzex(ax+1:bx,ay+1:by).*(Ex(ax+1:bx,ay+1:by)-Ex(ax+1:bx,ay:by-1));
   Hz(ax,ay+1:by)=Hz(ax,ay+1:by)+m_hzey(ax,ay+1:by).*(Ey(ax,ay+1:by)-Ey(bx,ay+1:by))+m_hzex(ax,ay+1:by).*(Ex(ax,ay+1:by)-Ex(ax,ay:by-1));
   Hz(ax+1:bx,ay)=Hz(ax+1:bx,ay)+m_hzey(ax+1:bx,ay).*(Ey(ax+1:bx,ay)-Ey(ax:bx-1,ay))+m_hzex(ax+1:bx,ay).*(Ex(ax+1:bx,ay)-Ex(ax+1:bx,by));
   Hz(ax,ay)=Hz(ax,ay)+m_hzey(ax,ay).*(Ey(ax,ay)-Ey(bx,ay))+m_hzex(ax,ay).*(Ex(ax,ay)-Ex(ax,by));
   
   Ex(ax:bx,ay:by-1)=Ex(ax:bx,ay:by-1)+m_ex(ax:bx,ay:by-1).*(Hz(ax:bx,ay+1:by)-Hz(ax:bx,ay:by-1));
   Ex(ax:bx,by)=Ex(ax:bx,by)+m_ex(ax:bx,by).*(Hz(ax:bx,ay)-Hz(ax:bx,by));
   
   Ey(ax:bx-1,ay:by)=Ey(ax:bx-1,ay:by)+m_ey(ax:bx-1,ay:by).*(Hz(ax+1:bx,ay:by)-Hz(ax:bx-1,ay:by));
   Ey(bx,ay:by)=Ey(bx,ay:by)+m_ey(bx,ay:by).*(Hz(ax,ay:by)-Hz(bx,ay:by));
   
   
    %Gaussian Source
    %f(n)= (-2*(n*dt-t0)*dt/(tw^2))*exp(-(n*dt-t0)^2/(tw^2))/dy;
 source(i)= (-2*(i*dt-0)*dt/(tw^2))*exp(-(i*dt-0)^2/(tw^2))/dy;
 Ez(100,100)=Ez(100,100)+source(i);
 %Hz(100,100)=Hz(100,100)+source(i);
%heatmapplot
% heatmap(X,Y,Ez,'GridVisible','off');
%TE mode
mesh(X,Y,Ez,'linewidth',2);
%mesh(X,Y,Hz,'linewidth',2);
%TM mode
%mesh(X,Y,Ex,'linewidth',2);
%mesh(X,Y,Ex,'linewidth',2);
xlabel('X \rightarrow');
ylabel('\leftarrow Y');
zlabel('E_z \rightarrow');
titlestring=['\fontsize{20}PBC 2D FDTD  at time  =',num2str(i*dt*1e18),'as'];
title(titlestring,'color','k');
%TE mode
axis([-L(1)/2 L(1)/2 -L(2)/2 L(2)/2 -2e5 2e5]);
caxis([-1e5 1e5]);
%TM mode

goodplot()
    
% hmh = heatmap(rgb,'xlabel','\mu','ylabel','\nu');
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
%drawnow
movieVector(i)=getframe(f);


end

%Create a VideoWriter object and set properties
myWriter = VideoWriter('2dfdtdwithpbc1');            %create an .avi file
% myWriter = VideoWriter('curve','MPEG-4');   %create an .mp4 file
myWriter.FrameRate = 20;
%Open the VideoWriter object, write the movie, and close the file
open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);
%if finished normally, display this
disp('DONE!')

save('fdtd2dwithPBCboudary.mat','Ez','Hx','Hy');



%%%%%%%%inverse%%%%%%%%%%%
%%
clf;
movieVector(1:number_of_iterations) = struct('cdata',[],'colormap',[]);
f=figure('Color','w','Position',[100 100 800 800]);
%zbuffer or painter are recommended as renderer
set(gcf,'Renderer','zbuffer');
for i=number_of_iterations:-1:1
   %TE mode
   Hx(ax:bx,ay:by-1)=Hx(ax:bx,ay:by-1)-m_hx(ax:bx,ay:by-1).*(Ez(ax:bx,ay+1:by)-Ez(ax:bx,ay:by-1));
   Hx(ax:bx,by)=Hx(ax:bx,by)-m_hx(ax:bx,by).*(Ez(ax:bx,ay)-Ez(ax:bx,by));
   
   Hy(ax:bx-1,ay:by)=Hy(ax:bx-1,ay:by)-m_hy(ax:bx-1,ay:by).*(-Ez(ax+1:bx,ay:by)+Ez(ax:bx-1,ay:by));
   Hy(bx,ay:by)=Hy(bx,ay:by)-m_hy(bx,ay:by).*(-Ez(1,ay:by)+Ez(bx,ay:by));
   
   Ez(ax+1:bx,ay+1:by)=Ez(ax+1:bx,ay+1:by)-m_ezhy(ax+1:bx,ay+1:by).*(Hy(ax+1:bx,ay+1:by)-Hy(ax:bx-1,ay+1:by))-m_ezhx(ax+1:bx,ay+1:by).*(Hx(ax+1:bx,ay+1:by)-Hx(ax+1:bx,ay:by-1));
%   Ez(ax,ay+1:by)=Ez(ax,ay+1:by)+m_ezhy(ax,ay+1:by).*(hy(ax,ay+1:by)-hy(bx,ay+1:by))+m_ezhx(ax,ay+1:by).*(hx(ax,ay+1:by)-hx(ax,ay:by-1));
   Ez(ax,ay+1:by)=Ez(ax,ay+1:by)-m_ezhy(ax,ay+1:by).*(Hy(ax,ay+1:by)-Hy(bx,ay+1:by))-m_ezhx(ax,ay+1:by).*(Hx(ax,ay+1:by)-Hx(ax,ay:by-1));
   Ez(ax+1:bx,ay)=Ez(ax+1:bx,ay)-m_ezhy(ax+1:bx,ay).*(Hy(ax+1:bx,ay)-Hy(ax:bx-1,ay))-m_ezhx(ax+1:bx,ay).*(Hx(ax+1:bx,ay)-Hx(ax+1:bx,by));
   Ez(ax,ay)=Ez(ax,ay)-m_ezhy(ax,ay).*(Hy(ax,ay)-Hy(bx,ay))-m_ezhx(ax,ay).*(Hx(ax,ay)-Hx(ax,by));

   
   %TM mode
%    Hz(ax:bx-1,ay:by-1)=Hz(ax:bx-1,ay:by-1)+m_hzey(ax:bx-1,ay:by-1).*(Ey(ax+1:bx,ay:by-1)-Ey(ax:bx-1,ay:by-1))+m_hzex(ax:bx-1,ay:by-1).*(Ex(ax:bx-1,ay+1:by)-Ex(ax:bx-1,ay:by-1));
%    Hz(bx,ay:by-1)=Hz(bx,ay:by-1)+m_hzey(bx,ay:by-1).*(Ey(ax,ay:by-1)-Ey(bx,ay:by-1))+m_hzex(bx,ay:by-1).*(Ex(bx,ay+1:by)-Ex(bx,ay:by-1));
%    Hz(ax:bx-1,by)=Hz(ax:bx-1,by)+m_hzey(ax:bx-1,by).*(Ey(ax+1:bx,by)-Ey(ax:bx-1,by))+m_hzex(ax:bx-1,by).*(Ex(ax:bx-1,ay)-Ex(ax:bx-1,by));
%    Hz(bx,by)=Hz(bx,by)+m_hzey(bx,by).*(Ey(ax,by)-Ey(bx,by))+m_hzex(bx,by).*(Ex(bx,ay)-Ey(bx,by));
%    
%    Ex(ax:bx,ay+1:by)=Ex(ax:bx,ay+1:by)+m_ex(ax:bx,ay+1:by).*(Hz(ax:bx,ay+1:by)-Hz(ax:bx,ay:by-1));
%    Ex(ax:bx,ay)=Ex(ax:bx,ay)+m_ex(ax:bx,ay).*(Hz(ax:bx,ay)-Hz(ax:bx,by));
%    
%    Ey(ax+1:bx,ay:by)=Ey(ax+1:bx,ay:by)+m_ey(ax+1:bx,ay:by).*(Hz(ax+1:bx,ay:by)-Hz(ax:bx-1,ay:by));
%    Ey(ax,ay:by)=Ey(ax,ay:by)+m_ey(ax,ay:by).*(Hz(ax,ay:by)-Hz(bx,ay:by));
   Hz(ax+1:bx,ay+1:by)=Hz(ax+1:bx,ay+1:by)-m_hzey(ax+1:bx,ay+1:by).*(Ey(ax+1:bx,ay+1:by)-Ey(ax:bx-1,ay+1:by))-m_hzex(ax+1:bx,ay+1:by).*(Ex(ax+1:bx,ay+1:by)-Ex(ax+1:bx,ay:by-1));
   Hz(ax,ay+1:by)=Hz(ax,ay+1:by)-m_hzey(ax,ay+1:by).*(Ey(ax,ay+1:by)-Ey(bx,ay+1:by))-m_hzex(ax,ay+1:by).*(Ex(ax,ay+1:by)-Ex(ax,ay:by-1));
   Hz(ax+1:bx,ay)=Hz(ax+1:bx,ay)-m_hzey(ax+1:bx,ay).*(Ey(ax+1:bx,ay)-Ey(ax:bx-1,ay))-m_hzex(ax+1:bx,ay).*(Ex(ax+1:bx,ay)-Ex(ax+1:bx,by));
   Hz(ax,ay)=Hz(ax,ay)-m_hzey(ax,ay).*(Ey(ax,ay)-Ey(bx,ay))-m_hzex(ax,ay).*(Ex(ax,ay)-Ex(ax,by));
   
   Ex(ax:bx,ay:by-1)=Ex(ax:bx,ay:by-1)-m_ex(ax:bx,ay:by-1).*(Hz(ax:bx,ay+1:by)-Hz(ax:bx,ay:by-1));
   Ex(ax:bx,by)=Ex(ax:bx,by)-m_ex(ax:bx,by).*(Hz(ax:bx,ay)-Hz(ax:bx,by));
   
   Ey(ax:bx-1,ay:by)=Ey(ax:bx-1,ay:by)-m_ey(ax:bx-1,ay:by).*(Hz(ax+1:bx,ay:by)-Hz(ax:bx-1,ay:by));
   Ey(bx,ay:by)=Ey(bx,ay:by)-m_ey(bx,ay:by).*(Hz(ax,ay:by)-Hz(bx,ay:by));
   
   
    %Gaussian Source
    %f(n)= (-2*(n*dt-t0)*dt/(tw^2))*exp(-(n*dt-t0)^2/(tw^2))/dy;
 %source(i)= (-2*(i*dt-0)*dt/(tw^2))*exp(-(i*dt-0)^2/(tw^2))/dy;
 %Ez(100,100)=Ez(100,100)+source(i);
 %Hz(100,100)=Hz(100,100)+source(i);
%heatmapplot
% heatmap(X,Y,Ez,'GridVisible','off');
%TE mode
mesh(X,Y,Ez,'linewidth',2);
%mesh(X,Y,Hz,'linewidth',2);
%TM mode
%mesh(X,Y,Ex,'linewidth',2);
%mesh(X,Y,Ex,'linewidth',2);
xlabel('X \rightarrow');
ylabel('\leftarrow Y');
zlabel('E_z \rightarrow');
titlestring=['\fontsize{20}Inverse 2d fdtd at time =',num2str(i*dt*1e18),'as'];
title(titlestring,'color','k');
%TE mode
axis([-L(1)/2 L(1)/2 -L(2)/2 L(2)/2 -4e5 4e5]);
caxis([-1e5 1e5]);
%TM mode

goodplot()
    
% hmh = heatmap(rgb,'xlabel','\mu','ylabel','\nu');
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
%drawnow
movieVector(number_of_iterations-i+1)=getframe(f);


end

%Create a VideoWriter object and set properties
myWriter = VideoWriter('2dfdtdinversewithpbc');            %create an .avi file
% myWriter = VideoWriter('curve','MPEG-4');   %create an .mp4 file
myWriter.FrameRate = 20;
%Open the VideoWriter object, write the movie, and close the file
open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);
%if finished normally, display this
disp('DONE2!')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function goodplot()
% function which produces a nice-looking plot
% and sets up the page for nice printing
set(get(gca,'xlabel'),'FontSize', 18, 'FontWeight', 'Bold');
set(get(gca,'ylabel'),'FontSize', 18, 'FontWeight', 'Bold');
set(get(gca,'zlabel'),'FontSize', 18, 'FontWeight', 'Bold');
set(get(gca,'title'),'FontSize', 18, 'FontWeight', 'Bold');
colorbar;
% box off; axis square;
set(gca,'LineWidth',2);
set(gca,'FontSize',14);
set(gca,'FontWeight','Bold');
set(gcf,'color','w');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [12 12]);
set(gcf,'PaperPosition',[0.5 0.5 7 7]);
set(gcf,'PaperPositionMode','Manual');
end