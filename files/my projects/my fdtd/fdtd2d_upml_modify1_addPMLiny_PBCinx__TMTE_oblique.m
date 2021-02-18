% fdtd2d for TM and TE mode
% supporting oblique incidence and PML boundary.

% INITIALIZE MATLAB
close all;
clc;
clear all;
TE=0;   % set to 1 for TE mode
TM=1;   % set to 1 for TM mode
plot=1; % set to 1 to show plot

angle=45; %set incident angle
angle=angle*pi/180;
color_range=[-1 1];
%Units
meters =1;
centimeters=1e-2*meters;
millimeters=1e-3*meters;
nm=1e-9*meters;
as=1e-18*meters;
seconds=1;
hertz=1/seconds;
gigahertz=1e9*hertz;

% CONSTANTS
e0=8.85418782e-12 * 1/meters;
u0=1.25663706e-6 * 1/meters;
N0=sqrt(u0/e0);
c0=299792458 * meters/seconds;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Source parameters
wavelength=488*nm;
fmax=c0/wavelength;
incidentangle=50;
incidentangle=incidentangle/180*pi;
wavevector=2*pi/wavelength;
%Grid Parameters
NRES=10; %resolution
nmax=1;
Nx=800;
Ny=500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate optimized grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute grid resolution

 xlength=4000*nm;
 ylength=2500*nm;
 dx=xlength/Nx;
 dy=ylength/Ny;


%compute grid axes (for graphics)
xa=[0:Nx-1]*dx;
ya=[0:Ny-1]*dy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build device on grid by modification of parameters of epsilon and mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize materials to free space, TE mode uses parameters with comments
ERzz=ones(Nx,Ny);   %epsilon
ERxx=ones(Nx,Ny);
ERyy=ones(Nx,Ny);
URxx=ones(Nx,Ny);   %mu in x
URyy=ones(Nx,Ny);   %mu in y
URzz=ones(Nx,Ny);
ERzz_ref=ones(Nx,Ny);   %epsilon
ERxx_ref=ones(Nx,Ny);
ERyy_ref=ones(Nx,Ny);
URxx_ref=ones(Nx,Ny);   %mu in x
URyy_ref=ones(Nx,Ny);   %mu in y
URzz_ref=ones(Nx,Ny);

%compute PML parameters on 2x grid
Nx2=2*Nx;
Ny2=2*Ny;
sigx=zeros(Nx2,Ny2);
sigy=zeros(Nx2,Ny2);
sigx_ref=zeros(Nx2,Ny2);
sigy_ref=zeros(Nx2,Ny2);
%add a silicon device
material=2; %number of materials
mx1=1;
mx2=Nx;
my1=1;
my2=200;
%add silicon
URxx(mx1:mx2,my1:my2)=1;
URyy(mx1:mx2,my1:my2)=1;
URzz(mx1:mx2,my1:my2)=1;
ERzz(mx1:mx2,my1:my2)=5.4376^2;
ERxx(mx1:mx2,my1:my2)=5.4376^2;
ERyy(mx1:mx2,my1:my2)=5.4376^2;

% add material conductivity
%add silicon
conductivity=1000;
sigx(mx1*2:mx2*2,my1*2:my2*2)=conductivity;
sigy(mx1*2:mx2*2,my1*2:my2*2)=conductivity;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%air W200D260
width=200*nm;
depth=350*nm;
mairx1=Nx/2-width/dx/2;
mairx2=Nx/2+width/dx/2;
mairy1=200-depth/dx;
mairy2=200;
%add air
URxx(mairx1:mairx2,mairy1:mairy2)=1;
URyy(mairx1:mairx2,mairy1:mairy2)=1;
URzz(mairx1:mairx2,mairy1:mairy2)=1;
ERzz(mairx1:mairx2,mairy1:mairy2)=1;
ERxx(mairx1:mairx2,mairy1:mairy2)=1;
ERyy(mairx1:mairx2,mairy1:mairy2)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute time step
% dmin=min([dx dy]);
% dt=dmin/(2*c0);
 t_end=10000*as;
 dt = ( 1/c0/sqrt( 1/((dx)^2) + 1/((dy)^2) ) )*0.5;
%dt=1.6e-10;
% calculate gaussian source parameters
tau = 0.5/fmax;%tau (pulse duration) of the gaussian
%tau=3.3e-9;
t0 = 3*tau ;%offset of the gaussian, delay 5 tau time

%courant number
S=c0*dt/dx;
%FDTD HARD SOURCE AND SOFT SOURCE REVIEWS
%AND MODIFICATIONS
alpha=1/S*1/(5*S-40*S^3+126*S^5-160*S^7+70*S^9);
% calculate number of iterations
STEPS = ceil(t_end/dt);
%STEPS=500;

%calculate gaussian source
t = [0:STEPS-1]*dt;
%g = exp(-((t - t0)/tau).^2);
nx_src=1+floor(0.25*Nx);         %the position of the source
ny_src=400;
N_lambda=wavelength/dx;
k=2*pi/N_lambda;
w=k*c0/dx;
%calculate complex CW
%h = sin(2*pi*fmax*t);
%h=S*alpha*exp(1i*(angularfreq*t+pi/2));

%h = cos(2*pi*fmax*t);
% plot(t,g);
% axis tight;
% return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate update coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NPML=[0 0 50 50]; % the thickness of PML boundaries
%% SOURCE
xcoordinate=(1+floor(0*Nx):1:floor(1*Nx));
nx_srcline=(1+NPML(1)+floor(0*Nx):1:-NPML(2)+floor(1*Nx));
nx_srclineshift=nx_srcline-min(nx_srcline)+1;
g=zeros(length(xcoordinate),STEPS);
%% NPML CONTINUE
for T = [1:STEPS]
    g(nx_srcline,T) = exp(((w)*T*dt*(nx_srcline./nx_srcline)-k*sin(angle)*(nx_srclineshift)+pi/2)*1i);
end


%left x (0 to NPML1)
for nx = 1 : 2*NPML(1)
   nx1 = 2*NPML(1) - nx + 1;
   sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(1))^3;
end
%right x (NPML2 to Nx)
for nx = 1 : 2*NPML(2)
   nx1 = Nx2-2*NPML(2) + nx;
   sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(2))^3;
end
%upper y (0 to NPML3)
for ny = 1 : 2*NPML(3)
   ny1 = 2*NPML(3) - ny + 1;
   sigy(:,ny1) = conductivity+(0.5*e0/dt)*(ny/2/NPML(3))^3;  %material is set at PML
end
%lower y (npml4 to Ny)
for ny = 1 : 2*NPML(4)
   ny1 = Ny2-2*NPML(4) + ny;
   sigy(:,ny1) =  (0.5*e0/dt)*(ny/2/NPML(4))^3; 
end

%left x (0 to NPML1)
for nx = 1 : 2*NPML(1)
   nx1 = 2*NPML(1) - nx + 1;
   sigx_ref(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(1))^3;
end
%right x (NPML2 to Nx)
for nx = 1 : 2*NPML(2)
   nx1 = Nx2-2*NPML(2) + nx;
   sigx_ref(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(2))^3;
end
%upper y (0 to NPML3)
for ny = 1 : 2*NPML(3)
   ny1 = 2*NPML(3) - ny + 1;
   sigy_ref(:,ny1) = conductivity+(0.5*e0/dt)*(ny/2/NPML(3))^3;  %material is set at PML
end
%lower y (npml4 to Ny)
for ny = 1 : 2*NPML(4)
   ny1 = Ny2-2*NPML(4) + ny;
   sigy_ref(:,ny1) =  (0.5*e0/dt)*(ny/2/NPML(4))^3; 
end






%add air
sigx(mairx1*2:mairx2*2,mairy1*2:mairy2*2)=0;
sigy(mairx1*2:mairx2*2,mairy1*2:mairy2*2)=0;
% imagesc(sigy.');
% axis equal tight;
% colorbar;
% return;
% 
%% Perform FDTD
%% TE
if TE==1
    %%%%%%%%%%%%%for soft source%%%%%%%%%%%%%%%%%%%%%%
    disp("TE mode Simulation");
    % 1D FDTD Initialization of field vectors%%%%%%%%%%%%%%%%%%%%%
    Ez_1D=zeros(Nx,1);
    %Dz_1D=zeros(Nx,1);
    H_1D=zeros(Nx,1);
    % Initialize curl arrays
    CE_1D=zeros(Nx,1);
    CH_1D=zeros(Nx,1);
    
    % Initialization of permittivity and permeability vectors
    ERzz_1D=ones(Nx,1);
    UR_1D=ones(Nx,1);
    %update coefficient
    mH_1D =  c0*dt./UR_1D;
    mEz_1D=c0*dt ./ERzz_1D;
    Ez_1D_receiver=zeros(Nx,1);
    %%%%%%%%%%%%%%%%1D end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2D parameters
    %update for mHx
    sigHx = sigx(1:2:Nx2,2:2:Ny2);
    sigHy = sigy(1:2:Nx2,2:2:Ny2);
    mHx0 = (1/dt) + sigHy/(2*e0);
    mHx1 = ((1/dt) - sigHy/(2*e0)) ./ mHx0;
    mHx2 = -c0 ./ URxx./mHx0;
    mHx3 = (-c0*dt/e0)*sigHx./URxx./mHx0;
    
    %update for mHy
    sigHx = sigx(2:2:Nx2,1:2:Ny2);
    sigHy = sigy(2:2:Nx2,1:2:Ny2);
    mHy0 = (1/dt) + sigHx/(2*e0);
    mHy1 = ((1/dt) - sigHx/(2*e0)) ./ mHy0;
    mHy2 = -c0 ./ URyy./mHy0;
    mHy3 = (-c0*dt/e0)*sigHy./URyy./mHy0;
    
    %update for Dz
    sigDx = sigx(1:2:Nx2,1:2:Ny2);
    sigDy = sigy(1:2:Nx2,1:2:Ny2);
    mDz0 = (1/dt) + (sigDx+sigDy)/(2*e0) + sigDx.*sigDy*(dt/4/e0^2);
    mDz1 = (1/dt) - (sigDx+sigDy)/(2*e0) - sigDx.*sigDy*(dt/4/e0^2);
    mDz1 = mDz1./mDz0;
    mDz2 = c0 ./mDz0;
    mDz4 = - (dt/e0^2)*sigDx.*sigDy./mDz0;
    
    mEz = 1./ERzz;
    mEz1=mEz;
    
    % Initialize fields
    Hx=zeros(Nx,Ny);
    Hy=zeros(Nx,Ny);
    Dz=zeros(Nx,Ny);
    Ez=zeros(Nx,Ny);
    Ezdx=zeros(Nx,Ny);
    Ezdy=zeros(Nx,Ny);
    % Initialize curl arrays
    CEx=zeros(Nx,Ny);
    CEy=zeros(Nx,Ny);
    CHz=zeros(Nx,Ny);
    
    % Initialize Intergration Terms
    ICEx=zeros(Nx,Ny);
    ICEy=zeros(Nx,Ny);
    IDz =zeros(Nx,Ny);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% intialize 2D reference field
    %update for mHx_ref
    sigHx_ref = sigx_ref(1:2:Nx2,2:2:Ny2);
    sigHy_ref = sigy_ref(1:2:Nx2,2:2:Ny2);
    mHx0_ref = (1/dt) + sigHy_ref/(2*e0);
    mHx1_ref = ((1/dt) - sigHy_ref/(2*e0)) ./ mHx0_ref;
    mHx2_ref = -c0 ./ URxx_ref./mHx0_ref;
    mHx3_ref = (-c0*dt/e0)*sigHx_ref./URxx_ref./mHx0_ref;
    
    %update for mHy_ref
    sigHx_ref = sigx_ref(2:2:Nx2,1:2:Ny2);
    sigHy_ref = sigy_ref(2:2:Nx2,1:2:Ny2);
    mHy0_ref = (1/dt) + sigHx_ref/(2*e0);
    mHy1_ref = ((1/dt) - sigHx_ref/(2*e0)) ./ mHy0_ref;
    mHy2_ref = -c0 ./ URyy_ref./mHy0_ref;
    mHy3_ref = (-c0*dt/e0)*sigHy_ref./URyy_ref./mHy0_ref;
    
    %update for Dz_ref
    sigDx_ref = sigx_ref(1:2:Nx2,1:2:Ny2);
    sigDy_ref = sigy_ref(1:2:Nx2,1:2:Ny2);
    mDz0_ref = (1/dt) + (sigDx_ref+sigDy_ref)/(2*e0) + sigDx_ref.*sigDy_ref*(dt/4/e0^2);
    mDz1_ref = (1/dt) - (sigDx_ref+sigDy_ref)/(2*e0) - sigDx_ref.*sigDy_ref*(dt/4/e0^2);
    mDz1_ref = mDz1_ref./mDz0_ref;
    mDz2_ref = c0 ./mDz0_ref;
    mDz4_ref = - (dt/e0^2)*sigDx_ref.*sigDy_ref./mDz0_ref;
    
    mEz_ref = 1./ERzz_ref;
    mEz1_ref=mEz_ref;
    % Initialize fields
    Hx_ref=zeros(Nx,Ny);
    Hy_ref=zeros(Nx,Ny);
    Dz_ref=zeros(Nx,Ny);
    Ez_ref=zeros(Nx,Ny);
    Ezdx_ref=zeros(Nx,Ny);
    Ezdy_ref=zeros(Nx,Ny);
    % Initialize curl arrays
    CEx_ref=zeros(Nx,Ny);
    CEy_ref=zeros(Nx,Ny);
    CHz_ref=zeros(Nx,Ny);
    
    % Initialize Intergration Terms
    ICEx_ref=zeros(Nx,Ny);
    ICEy_ref=zeros(Nx,Ny);
    IDz_ref =zeros(Nx,Ny);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Perform 2d fdtd
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %
    % MAIN LOOP
    %
    %Electricfieldsaver=zeros(Nx,Ny,STEPS);
    %movieVector(1:STEPS) = struct('cdata',[],'colormap',[]);
    f=figure('color','w','units','normalized','outerposition',[0 0 1 1]);
    %f=figure('Color','w','Position',[100 100 800 800],'visible','off');
    %S=0.5;
    %alpha=1/S*1/(5*S-40*S^2+126*S^5-160*S^7+70*S^9);
    set(gcf, 'Position', get(0, 'Screensize'));
    %entropy=zeros(STEPS,1);
    for T=1: STEPS
        %tic;
        %[Ez,Hx,Hy,ICEx,ICEy,IDz,Dz,Nx,Ny,dx,dy,mHx1,mHx2,mHx3,mDz1,mDz2,mDz4,mHy1,mHy2,mHy3,mEz]=fdtd(Ez,Hx,Hy,ICEx,ICEy,IDz,Dz,Nx,Ny,dx,dy,mHx1,mHx2,mHx3,mDz1,mDz2,mDz4,mHy1,mHy2,mHy3,mEz);
        %G=parfeval(@fdtd,21,Ez,Hx,Hy,ICEx,ICEy,IDz,Dz,Nx,Ny,dx,dy,mHx1,mHx2,mHx3,mDz1,mDz2,mDz4,mHy1,mHy2,mHy3,mEz);
        %[Ez,Hx,Hy,ICEx,ICEy,IDz,Dz,Nx,Ny,dx,dy,mHx1,mHx2,mHx3,mDz1,mDz2,mDz4,mHy1,mHy2,mHy3,mEz]=fetchOutputs(G);
        %batch('fdtd',21,{Ez,Hx,Hy,ICEx,ICEy,IDz,Dz,Nx,Ny,dx,dy,mHx1,mHx2,mHx3,mDz1,mDz2,mDz4,mHy1,mHy2,mHy3,mEz});
        %% fdtd 2D with device%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ezdy(1:Nx,1:Ny-1)=Ez(1:Nx,2:Ny+1-1);
        
        Ezdy(1:Nx,Ny)=Ez(1:Nx,1);
        
        CEx=(Ezdy-Ez)./dy;

        %Compute CEy, which is related to dHy/dt, and is related to
        %electric field dEz/dx
        %Update Integration terms
        ICEx=ICEx + CEx;
        
        Ezdx(1:Nx-1,1:Ny)=Ez(2:Nx+1-1,1:Ny);
        Ezdx(Nx,1:Ny)    =Ez(1,1:Ny);
        
        %Update Integration terms
        CEy = -(Ezdx-Ez)./dx;
        ICEy=ICEy + CEy;
        
        %Update Hx and Hy
        
        Hx = mHx1.*Hx + mHx2.*CEx + mHx3.*ICEx;
        
        Hy = mHy1.*Hy + mHy2.*CEy + mHy3.*ICEy;
        
        
        
        CHz(1,1)=(Hy(1,1)-Hy(Nx,1))/dx-(Hx(1,1)-Hy(1,Ny))/dy;
        
        CHz(2:Nx,1)=(Hy(2:Nx,1)-Hy(1:Nx-1,1))/dx ...
            -(Hx(2:Nx,1)-Hx(2:Nx,Ny))/dy;
        
        CHz(1,2:Ny)=(Hy(1,2:Ny)-Hy(Nx,2:Ny))/dx ...            %%PBC is applied Hy(Nx,ny) term
            -(Hx(1,2:Ny)-Hx(1,1:Ny-1))/dy;
        
        CHz(2:Nx,2:Ny)=(Hy(2:Nx,2:Ny)-Hy(1:Nx-1,2:Ny))/dx ...
            -(Hx(2:Nx,2:Ny)-Hx(2:Nx,1:Ny-1))/dy;
        %Update integration term
        IDz = IDz + Dz;
        
        
        %update Dz
        
        Dz = mDz1.* Dz + mDz2 .* CHz + mDz4.*IDz;
        
        
        %update Ez
        Ez=mEz.*Dz;
        %% fdtd 2D with device set source to eliminate horizontally propagation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for point0=1:Nx
            Ez_1D(:)=Ez_1D(:)+g(T)/Nx;
            Ez_1D(point0)=Ez_1D(point0)-g(T)/Nx;
            % Update H from E (Dirichlet Boundary Conditions)
            %for nx = 1 : Nx-1
            H_1D(1 : Nx-1) = H_1D(1 : Nx-1) + mH_1D(1 : Nx-1).*(Ez_1D(2 : Nx) - Ez_1D(1 : Nx-1))./dx;
            %end
            H_1D(Nx) = H_1D(Nx) + mH_1D(Nx)*(Ez_1D(1) - Ez_1D(Nx))/dx;
            % Update E from H (Dirichlet Boundary Conditions)
            Ez_1D(1) = Ez_1D(1) + mEz_1D(1)*(H_1D(1) - H_1D(Nx))/dx;
            %for nx = 2 : Nx
            Ez_1D( 2 : Nx) = Ez_1D( 2 : Nx)  + mEz_1D( 2 : Nx).*(H_1D( 2 : Nx) - H_1D( 1 : Nx-1))./dx;
            % [Ez_1D,H_1D,mH_1D,mEz_1D,dx,Nx]=fdtd_1D(Ez_1D,H_1D,mH_1D,mEz_1D,dx,Nx);
            %end
            Ez_1D_receiver(point0,1)=Ez_1D(point0);
        end
        
        %% fdtd 2D with device inject source%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if angle==0
            Ez(:,ny_src)= Ez(:,ny_src)+g(T)-Ez_1D_receiver(:,1);
        else
            Ez(:,ny_src)= g(:,T);
        end    
        %inject source, soft source
        %Ez(1:Nx,ny_src)= Ez(1:Nx,ny_src)+h(T);
        %inject source, hard source
        %Ez(1:Nx,ny_src)= h(T);
        
        %% fdtd 2D plot in figure1 with boundaries and devices%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %To check the update engine is fine, use an if plot
        if mod(T,2)==0 && plot==1
                % Draw field
                subplot(1,3,1)
                imagesc(xa,ya,real(Ez).');
        
                set(gca,'Ydir','normal');
                title(['STEP' num2str(T) ' of ' num2str(STEPS)]);
                %colorbar;
                caxis(color_range);
        
            hold on;
            if NPML(1)
                x1 = xa(1);
                x2 = xa(NPML(1));
                y1 = ya(1);
                y2 = ya(Ny);
                x  = [ x1 x2 x2 x1 x1 ];
                y  = [ y1 y1 y2 y2 y1 ];
                c  = 0.5*[1 1 1];
                fill (x,y,c,'FaceAlpha',0.5);
        
            end
            if NPML(2)
                x1 = xa(Nx-NPML(2));
                x2 = xa(Nx);
                y1 = ya(1);
                y2 = ya(Ny);
                x  = [ x1 x2 x2 x1 x1 ];
                y  = [ y1 y1 y2 y2 y1 ];
                c  = 0.5*[1 1 1];
                fill (x,y,c,'FaceAlpha',0.5);
        
            end
             if NPML(3)
                x1 = xa(1);
                x2 = xa(Nx);
                y1 = ya(1);
                y2 = ya(NPML(3));
        
                x  = [ x1 x2 x2 x1 x1 ];
                y  = [ y1 y1 y2 y2 y1 ];
                c  = 0.5*[1 1 1];
                fill (x,y,c,'FaceAlpha',0.5);
        
             end
             if NPML(4)
                 x1 = xa(1);
                 x2 = xa(Nx);
                 y1 = ya(Ny - NPML(4));
                 y2 = ya(Ny);
        
                 x  = [ x1 x2 x2 x1 x1 ];
                 y  = [ y1 y1 y2 y2 y1 ];
                 c  = 0.5*[1 1 1];
                 fill (x,y,c,'FaceAlpha',0.5);
        
             end
             %plot the device
             if material==1
                 x1 = xa(mx1);
                 x2 = xa(mx2);
                 y1 = ya(my1);
                 y2 = ya(my2);
        
                 x  = [ x1 x2 x2 x1 x1 ];
                 y  = [ y1 y1 y2 y2 y1 ];
                 c  = 0.8*[1 1 1];
                 fill (x,y,c,'FaceAlpha',0.5);
        
             end
             if material==2
                 %left
                 x1 = xa(mx1);
                 x2 = xa(mairx1);
                 y1 = ya(my1);
                 y2 = ya(my2);
                 x  = [ x1 x2 x2 x1 x1 ];
                 y  = [ y1 y1 y2 y2 y1 ];
                 c  = 0.8*[1 1 1];
                 fill (x,y,c,'FaceAlpha',0.5);
                 %middle
                 x1 = xa(mairx1);
                 x2 = xa(mairx2);
                 y1 = ya(my1);
                 y2 = ya(mairy1);
                 x  = [ x1 x2 x2 x1 x1 ];
                 y  = [ y1 y1 y2 y2 y1 ];
                 c  = 0.8*[1 1 1];
                 fill (x,y,c,'FaceAlpha',0.5);
                  %right
                 x1 = xa(mairx2);
                 x2 = xa(mx2);
                 y1 = ya(my1);
                 y2 = ya(my2);
                 x  = [ x1 x2 x2 x1 x1 ];
                 y  = [ y1 y1 y2 y2 y1 ];
                 c  = 0.8*[1 1 1];
                 fill (x,y,c,'FaceAlpha',0.5);
             end
        end
            hold off;
            goodplot();
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% fdtd 2D without device%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %use functions will be slower
        %F = parfeval(@fdtd,21,Ez_ref,Hx_ref,Hy_ref,ICEx_ref,ICEy_ref,IDz_ref,Dz_ref,Nx,Ny,dx,dy,mHx1_ref,mHx2_ref,mHx3_ref,mDz1_ref,mDz2_ref,mDz4_ref,mHy1_ref,mHy2_ref,mHy3_ref,mEz_ref);
        %[Ez_ref,Hx_ref,Hy_ref,ICEx_ref,ICEy_ref,IDz_ref,Dz_ref,Nx,Ny,dx,dy,mHx1_ref,mHx2_ref,mHx3_ref,mDz1_ref,mDz2_ref,mDz4_ref,mHy1_ref,mHy2_ref,mHy3_ref,mEz_ref]=fdtd(Ez_ref,Hx_ref,Hy_ref,ICEx_ref,ICEy_ref,IDz_ref,Dz_ref,Nx,Ny,dx,dy,mHx1_ref,mHx2_ref,mHx3_ref,mDz1_ref,mDz2_ref,mDz4_ref,mHy1_ref,mHy2_ref,mHy3_ref,mEz_ref);
        %[Ez_ref,Hx_ref,Hy_ref,ICEx_ref,ICEy_ref,IDz_ref,Dz_ref,Nx,Ny,dx,dy,mHx1_ref,mHx2_ref,mHx3_ref,mDz1_ref,mDz2_ref,mDz4_ref,mHy1_ref,mHy2_ref,mHy3_ref,mEz_ref] = fetchOutputs(F);
        %batch('fdtd_ref',21,{Ez,Hx,Hy,ICEx,ICEy,IDz,Dz,Nx,Ny,dx,dy,mHx1,mHx2,mHx3,mDz1,mDz2,mDz4,mHy1,mHy2,mHy3,mEz});
        %Ez_ref(1:Nx,ny_src)= h(T);
        %Ez_ref(1:Nx,ny_src)=Ez_ref(1:Nx,ny_src)+h(T);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fdtd%%%%%%%%%%%%
        Ezdy_ref(1:Nx,1:Ny-1)=Ez_ref(1:Nx,2:Ny+1-1);
        
        Ezdy_ref(1:Nx,Ny)=Ez_ref(1:Nx,1);
        %CEx(1:Nx,Ny) = (Ezdy(1:Nx,Ny) - Ez(1:Nx,Ny))/dy;
        CEx_ref=(Ezdy_ref-Ez_ref)./dy;
        
        
        %Compute CEy, which is related to dHy/dt, and is related to
        %electric field dEz/dx
        %Update Integration terms
        ICEx_ref=ICEx_ref + CEx_ref;
        
        Ezdx_ref(1:Nx-1,1:Ny)=Ez_ref(2:Nx+1-1,1:Ny);
        Ezdx_ref(Nx,1:Ny)    =Ez_ref(1,1:Ny);
        
        %Update Integration terms
        CEy_ref = -(Ezdx_ref-Ez_ref)./dx;
        ICEy_ref=ICEy_ref + CEy_ref;
        
        %Update Hx and Hy
        
        
        Hx_ref = mHx1_ref.*Hx_ref + mHx2_ref.*CEx_ref + mHx3_ref.*ICEx_ref;
        
        Hy_ref = mHy1_ref.*Hy_ref + mHy2_ref.*CEy_ref + mHy3_ref.*ICEy_ref;
        
        
        %Calculate CHz
        
        CHz_ref(1,1)=(Hy_ref(1,1)-Hy_ref(Nx,1))/dx-(Hx_ref(1,1)-Hy_ref(1,Ny))/dy;
        
        CHz_ref(2:Nx,1)=(Hy_ref(2:Nx,1)-Hy_ref(1:Nx-1,1))/dx ...
            -(Hx_ref(2:Nx,1)-Hx_ref(2:Nx,Ny))/dy;
        
        CHz_ref(1,2:Ny)=(Hy_ref(1,2:Ny)-Hy_ref(Nx,2:Ny))/dx ...            %%PBC is applied Hy(Nx,ny) term
            -(Hx_ref(1,2:Ny)-Hx_ref(1,1:Ny-1))/dy;
        
        CHz_ref(2:Nx,2:Ny)=(Hy_ref(2:Nx,2:Ny)-Hy_ref(1:Nx-1,2:Ny))/dx ...
            -(Hx_ref(2:Nx,2:Ny)-Hx_ref(2:Nx,1:Ny-1))/dy;
        %Update integration term
        IDz_ref = IDz_ref + Dz_ref;
        
        
        %update Dz
        
        Dz_ref = mDz1_ref.* Dz_ref + mDz2_ref .* CHz_ref + mDz4_ref.*IDz_ref;
        
        
        %update Ez
        Ez_ref=mEz_ref.*Dz_ref;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% inject source. Since it is same source, just inject.
        % inject source
        if angle==0
            Ez_ref(:,ny_src)= Ez_ref(:,ny_src)+g(T)-Ez_1D_receiver(:,1);
        else
            Ez_ref(:,ny_src)= g(:,T);
        end    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% plot boundaries without devices
        %plot the second figure of reference area
        if mod(T,2)==0 && plot==1
            subplot(1,3,2)
            imagesc(xa,ya,real(Ez_ref).');
            set(gca,'Ydir','normal');
            title(['STEP' num2str(T) ' of ' num2str(STEPS)]);
            hold on;
            if NPML(1)
                x1 = xa(1);
                x2 = xa(NPML(1));
                y1 = ya(1);
                y2 = ya(Ny);
                x  = [ x1 x2 x2 x1 x1 ];
                y  = [ y1 y1 y2 y2 y1 ];
                c  = 0.5*[1 1 1];
                fill (x,y,c,'FaceAlpha',0.5);
                
            end
            if NPML(2)
                x1 = xa(Nx-NPML(2));
                x2 = xa(Nx);
                y1 = ya(1);
                y2 = ya(Ny);
                x  = [ x1 x2 x2 x1 x1 ];
                y  = [ y1 y1 y2 y2 y1 ];
                c  = 0.5*[1 1 1];
                fill (x,y,c,'FaceAlpha',0.5);
                
            end
            if NPML(3)
                x1 = xa(1);
                x2 = xa(Nx);
                y1 = ya(1);
                y2 = ya(NPML(3));
                
                x  = [ x1 x2 x2 x1 x1 ];
                y  = [ y1 y1 y2 y2 y1 ];
                c  = 0.5*[1 1 1];
                fill (x,y,c,'FaceAlpha',0.5);
                
            end
            if NPML(4)
                x1 = xa(1);
                x2 = xa(Nx);
                y1 = ya(Ny - NPML(4));
                y2 = ya(Ny);
                
                x  = [ x1 x2 x2 x1 x1 ];
                y  = [ y1 y1 y2 y2 y1 ];
                c  = 0.5*[1 1 1];
                fill (x,y,c,'FaceAlpha',0.5);
                
            end
            hold off;
            %colorbar;
            caxis(color_range);
            goodplot();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %plot the second figure
            subplot(1,3,3)
            imagesc(xa,ya,real(Ez-Ez_ref).');
            set(gca,'Ydir','normal');
            title(['STEP' num2str(T) ' of ' num2str(STEPS)]);
            caxis(color_range);
            %Electricfieldsaver(:,:,T)=((Ez-Ez_ref));
            goodplot();
            %force graphics to show
            %drawnow;
        end
        movieVector(T)=getframe(f);
        
        %entropy(T,1)=entropy(T,1)+sum(Ez(1:Nx,1:Ny).^2,"all")^2/sum(Ez(1:Nx,1:Ny).^4,"all");
        %toc;
        %break;
    end
end
if TM==1
    disp("TM mode Simulation");
    
        % 1D FDTD Initialization of field vectors
        Ex_1D=zeros(Nx,1);
        Dx_1D=zeros(Nx,1);
        H_1D=zeros(Nx,1);
        % Initialize curl arrays
        CE_1D=zeros(Nx,1);
        CH_1D=zeros(Nx,1);
    
        % Initialization of permittivity and permeability vectors
        ERzz_1D=ones(Nx,1);
        UR_1D=ones(Nx,1);
        %update coefficient
        mH_1D =  c0*dt./UR_1D;
        mEz_1D=c0*dt ./ERzz_1D;
        Ex_1D_receiver=zeros(Nx,1);
    %%%%%%%%%%%1D end%%%%%%%%%%%%%%%%%%%%%%%
    %% initialize 2D parameters
    
    %update for mHz
    sigHx = sigx(1:2:Nx2,1:2:Ny2);
    sigHy = sigy(1:2:Nx2,1:2:Ny2);
    mHz0 = (1 / dt) + (sigHx + sigHy) / (2 * e0) + sigHx .* sigHy .* dt / (4*e0*e0);
    mHz1 = ((1 / dt) - (sigHx + sigHy) / (2 * e0) - sigHx .* sigHy .* dt / (4*e0*e0)) ./ mHz0;
    mHz2 = -c0 ./ URzz./mHz0;
    mHz3 = 0;
    mHz4 = - dt / (e0*e0) .* sigHx .* sigHy ./ mHz0;
    
    %update for mDx
    sigDx = sigx(1:2:Nx2,2:2:Ny2);
    sigDy = sigy(1:2:Nx2,2:2:Ny2);
    mDx0 = (1/dt) + (sigDy)/(2*e0) ;
    mDx1 = (1/dt) - (sigDy)/(2*e0) ;
    mDx1 = mDx1./mDx0;
    mDx2 = c0 ./mDx0;
    mDx3 = c0 * dt .* sigDx / e0./mDx0;
    mDx4 = 0;
    mEx = 1./ERxx;
    
    
    %update for mDy
    sigDx = sigx(2:2:Nx2,1:2:Ny2);
    sigDy = sigy(2:2:Nx2,1:2:Ny2);
    mDy0 = (1/dt) + (sigDx)/(2*e0) ;
    mDy1 = (1/dt) - (sigDx)/(2*e0) ;
    mDy1 = mDy1./mDy0;
    mDy2 = c0 ./mDy0;
    mDy3 = c0 * dt .* sigDy / e0./mDy0;
    mDy4 = 0;
    mEy = 1./ERyy;
  
    
    % TM mode fileds
    Hz=zeros(Nx,Ny);
    Dx=zeros(Nx,Ny);
    Dy=zeros(Nx,Ny);
    Ex=zeros(Nx,Ny);
    Ey=zeros(Nx,Ny);
    
    % Initialize curl arrays
    
    % TM mode curl arrays
    CHx=zeros(Nx,Ny);
    CHy=zeros(Nx,Ny);
    CEz=zeros(Nx,Ny);
    %ICEz=zeros(Nx,Ny);
    IHz =zeros(Nx,Ny);
    %IDx =zeros(Nx,Ny);
    ICHx=zeros(Nx,Ny);
    %IDy =zeros(Nx,Ny);
    ICHy=zeros(Nx,Ny);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% intialize reference field
    %update for mHz
    sigHx_ref = sigx_ref(1:2:Nx2,1:2:Ny2);
    sigHy_ref = sigy_ref(1:2:Nx2,1:2:Ny2);
    mHz0_ref = (1 / dt) + (sigHx_ref + sigHy_ref) / (2 * e0) + sigHx_ref .* sigHy_ref .* dt / (4*e0*e0);
    mHz1_ref = ((1 / dt) - (sigHx_ref + sigHy_ref) / (2 * e0) - sigHx_ref .* sigHy_ref .* dt / (4*e0*e0)) ./ mHz0_ref;
    mHz2_ref = -c0 ./ URzz_ref./mHz0_ref;
    mHz3_ref = 0;
    mHz4_ref = - dt / (e0*e0) .* sigHx_ref .* sigHy_ref ./ mHz0_ref;
    
    %update for mDx
    sigDx_ref = sigx_ref(1:2:Nx2,2:2:Ny2);
    sigDy_ref = sigy_ref(1:2:Nx2,2:2:Ny2);
    mDx0_ref = (1/dt) + (sigDy_ref)/(2*e0) ;
    mDx1_ref = (1/dt) - (sigDy_ref)/(2*e0) ;
    mDx1_ref = mDx1_ref./mDx0_ref;
    mDx2_ref = c0 ./mDx0_ref;
    mDx3_ref = c0 * dt .* sigDx_ref / e0./mDx0_ref;
    mDx4_ref = 0;
    mEx_ref = 1./ERxx_ref;
    
    %update for mDy
    sigDx_ref = sigx_ref(2:2:Nx2,1:2:Ny2);
    sigDy_ref = sigy_ref(2:2:Nx2,1:2:Ny2);
    mDy0_ref = (1/dt) + (sigDx_ref)/(2*e0) ;
    mDy1_ref = (1/dt) - (sigDx_ref)/(2*e0) ;
    mDy1_ref = mDy1_ref./mDy0_ref;
    mDy2_ref = c0 ./mDy0_ref;
    mDy3_ref = c0 * dt .* sigDy_ref / e0./mDy0_ref;
    mDy4_ref = 0;
    mEy_ref = 1./ERyy_ref;
    
    % TM mode fileds
    Hz_ref =zeros(Nx,Ny);
    Dx_ref =zeros(Nx,Ny);
    Dy_ref =zeros(Nx,Ny);
    Ex_ref =zeros(Nx,Ny);
    Ey_ref =zeros(Nx,Ny);
    
    % Initialize curl arrays
    
    % TM mode curl arrays
    CHx_ref =zeros(Nx,Ny);
    CHy_ref =zeros(Nx,Ny);
    CEz_ref =zeros(Nx,Ny);
    %ICEz=zeros(Nx,Ny);
    IHz_ref  =zeros(Nx,Ny);
    %IDx =zeros(Nx,Ny);
    ICHx_ref =zeros(Nx,Ny);
    %IDy =zeros(Nx,Ny);
    ICHy_ref =zeros(Nx,Ny);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Perform 2d fdtd
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %movieVector(1:STEPS) = struct('cdata',[],'colormap',[]);
    f=figure('color','w','units','normalized','outerposition',[0 0 1 1]);
    for T=1: STEPS
        
        %compute CEz
        CEz(1:Nx-1,1:Ny-1) =  (Ey(2:Nx-1+1,1:Ny-1)-Ey(1:Nx-1,1:Ny-1))/dx - (Ex(1:Nx-1,2:Ny-1+1)-Ex(1:Nx-1,1:Ny-1))/dy;
        CEz(Nx,1:Ny-1)     =  (Ey(1,1:Ny-1)-Ey(Nx,1:Ny-1))/dx - (Ex(Nx,2:Ny-1+1)-Ex(Nx,1:Ny-1))/dy;
        CEz(1:Nx-1,Ny) =  (Ey(2:Nx-1+1,Ny)-Ey(1:Nx-1,Ny))/dx - (Ex(1:Nx-1,1)-Ex(1:Nx-1,Ny))/dy;
        CEz(Nx,Ny)     =  (Ey(1,Ny)-Ey(Nx,Ny))/dx - (Ex(Nx,1)-Ex(Nx,Ny))/dy;
        
        %Update Integration terms
        IHz=IHz + Hz;
        
        %Update Hz
        Hz= mHz1.*Hz + mHz2.*CEz + mHz4.*IHz;
        
        %calculate CHx
        CHx(1:Nx,2:Ny)=(Hz(1:Nx,2:Ny)-Hz(1:Nx,1:Ny-1))/dy;
        CHx(1:Nx,1)   =(Hz(1:Nx,1)-Hz(1:Nx,Ny))/dy;
        
        %Update Integration terms
        %IDx=IDx + Dx;
        ICHx=ICHx + CHx;
        
        %update Dx
        Dx=mDx1.*Dx + mDx2.*CHx + mDx3.*ICHx;
        
        %calculate CHy
        CHy(2:Nx,1:Ny)=-(Hz(2:Nx,1:Ny)-Hz(1:Nx-1,1:Ny))/dx;
        CHy(1,1:Ny)=-(Hz(1,1:Ny)-Hz(Nx,1:Ny))/dx;
        
        %Update Integration terms
        ICHy=ICHy + CHy;
        
        %update Dy
        Dy=mDy1.*Dy + mDy2.*CHy + mDy3.*ICHy;
                
        %update Ex,Ey
        Ex=mEx.*Dx;
        Ey=mEy.*Dy;
        
        %inject source
        %Hz(nx_src,ny_src)= Hz(nx_src,ny_src)+g(T);
        %Hz(nx_src,ny_src)= +g(T);
        
        %To check the update engine is fine, use an if plot
        %% 1D source
        for point0=1:Nx
            Ex_1D(:)=Ex_1D(:)+g(T)/Nx;
            Ex_1D(point0)=Ex_1D(point0)-g(T)/Nx;
            % Update H from E (Dirichlet Boundary Conditions)
            %for nx = 1 : Nx-1
            H_1D(1 : Nx-1) = H_1D(1 : Nx-1) + mH_1D(1 : Nx-1).*(Ex_1D(2 : Nx) - Ex_1D(1 : Nx-1))./dx;
            %end
            H_1D(Nx) = H_1D(Nx) + mH_1D(Nx)*(Ex_1D(1) - Ex_1D(Nx))/dx;
            % Update E from H (Dirichlet Boundary Conditions)
            Ex_1D(1) = Ex_1D(1) + mEz_1D(1)*(H_1D(1) - H_1D(Nx))/dx;
            %for nx = 2 : Nx
            Ex_1D( 2 : Nx) = Ex_1D( 2 : Nx)  + mEz_1D( 2 : Nx).*(H_1D( 2 : Nx) - H_1D( 1 : Nx-1))./dx;
            % [Ez_1D,H_1D,mH_1D,mEz_1D,dx,Nx]=fdtd_1D(Ez_1D,H_1D,mH_1D,mEz_1D,dx,Nx);
            %end
            Ex_1D_receiver(point0,1)=Ex_1D(point0);
        end
        
        %% fdtd 2D with device inject source%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ex(:,ny_src)= Ex(:,ny_src)+g(T)-Ex_1D_receiver(:,1);
        if angle==0
            Ex(:,ny_src)= Ex(:,ny_src)+g(T)-Ez_1D_receiver(:,1);
        else
            Ex(:,ny_src)= g(:,T);
        end    
        %inject source, soft source
        %Ez(1:Nx,ny_src)= Ez(1:Nx,ny_src)+h(T);
        %inject source, hard source
        %Ez(1:Nx,ny_src)= h(T);
        
        %% fdtd 2D plot in figure1 with boundaries and devices%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %To check the update engine is fine, use an if plot
        if mod(T,1)==0 && plot==1
                % Draw field
                subplot(1,3,1)
                imagesc(xa,ya,real(Ex).');
        
                set(gca,'Ydir','normal');
                title(['STEP' num2str(T) ' of ' num2str(STEPS)]);
                %colorbar;
                caxis(color_range);
        
            hold on;
            if NPML(1)
                x1 = xa(1);
                x2 = xa(NPML(1));
                y1 = ya(1);
                y2 = ya(Ny);
                x  = [ x1 x2 x2 x1 x1 ];
                y  = [ y1 y1 y2 y2 y1 ];
                c  = 0.5*[1 1 1];
                fill (x,y,c,'FaceAlpha',0.5);
        
            end
            if NPML(2)
                x1 = xa(Nx-NPML(2));
                x2 = xa(Nx);
                y1 = ya(1);
                y2 = ya(Ny);
                x  = [ x1 x2 x2 x1 x1 ];
                y  = [ y1 y1 y2 y2 y1 ];
                c  = 0.5*[1 1 1];
                fill (x,y,c,'FaceAlpha',0.5);
        
            end
             if NPML(3)
                x1 = xa(1);
                x2 = xa(Nx);
                y1 = ya(1);
                y2 = ya(NPML(3));
        
                x  = [ x1 x2 x2 x1 x1 ];
                y  = [ y1 y1 y2 y2 y1 ];
                c  = 0.5*[1 1 1];
                fill (x,y,c,'FaceAlpha',0.5);
        
             end
             if NPML(4)
                 x1 = xa(1);
                 x2 = xa(Nx);
                 y1 = ya(Ny - NPML(4));
                 y2 = ya(Ny);
        
                 x  = [ x1 x2 x2 x1 x1 ];
                 y  = [ y1 y1 y2 y2 y1 ];
                 c  = 0.5*[1 1 1];
                 fill (x,y,c,'FaceAlpha',0.5);
        
             end
             %plot the device
             if material==1
                 x1 = xa(mx1);
                 x2 = xa(mx2);
                 y1 = ya(my1);
                 y2 = ya(my2);
        
                 x  = [ x1 x2 x2 x1 x1 ];
                 y  = [ y1 y1 y2 y2 y1 ];
                 c  = 0.8*[1 1 1];
                 fill (x,y,c,'FaceAlpha',0.5);
        
             end
             if material==2
                 %left
                 x1 = xa(mx1);
                 x2 = xa(mairx1);
                 y1 = ya(my1);
                 y2 = ya(my2);
                 x  = [ x1 x2 x2 x1 x1 ];
                 y  = [ y1 y1 y2 y2 y1 ];
                 c  = 0.8*[1 1 1];
                 fill (x,y,c,'FaceAlpha',0.5);
                 %middle
                 x1 = xa(mairx1);
                 x2 = xa(mairx2);
                 y1 = ya(my1);
                 y2 = ya(mairy1);
                 x  = [ x1 x2 x2 x1 x1 ];
                 y  = [ y1 y1 y2 y2 y1 ];
                 c  = 0.8*[1 1 1];
                 fill (x,y,c,'FaceAlpha',0.5);
                  %right
                 x1 = xa(mairx2);
                 x2 = xa(mx2);
                 y1 = ya(my1);
                 y2 = ya(my2);
                 x  = [ x1 x2 x2 x1 x1 ];
                 y  = [ y1 y1 y2 y2 y1 ];
                 c  = 0.8*[1 1 1];
                 fill (x,y,c,'FaceAlpha',0.5);
             end
        
            hold off;
            goodplot();
        end
        %% fdtd 2D without device%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
              %compute CEz
        CEz_ref(1:Nx-1,1:Ny-1) =  (Ey_ref(2:Nx-1+1,1:Ny-1)-Ey_ref(1:Nx-1,1:Ny-1))/dx - (Ex_ref(1:Nx-1,2:Ny-1+1)-Ex_ref(1:Nx-1,1:Ny-1))/dy;
        CEz_ref(Nx,1:Ny-1)     =  (Ey_ref(1,1:Ny-1)-Ey_ref(Nx,1:Ny-1))/dx - (Ex_ref(Nx,2:Ny-1+1)-Ex_ref(Nx,1:Ny-1))/dy;
        CEz_ref(1:Nx-1,Ny) =  (Ey_ref(2:Nx-1+1,Ny)-Ey_ref(1:Nx-1,Ny))/dx - (Ex_ref(1:Nx-1,1)-Ex_ref(1:Nx-1,Ny))/dy;
        CEz_ref(Nx,Ny)     =  (Ey_ref(1,Ny)-Ey_ref(Nx,Ny))/dx - (Ex_ref(Nx,1)-Ex_ref(Nx,Ny))/dy;
        
        %Update Integration terms
        IHz_ref=IHz_ref + Hz_ref;
        
        %Update Hz
        Hz_ref= mHz1_ref.*Hz_ref + mHz2_ref.*CEz_ref + mHz4_ref.*IHz_ref;
        
        %calculate CHx
        CHx_ref(1:Nx,2:Ny)=(Hz_ref(1:Nx,2:Ny)-Hz_ref(1:Nx,1:Ny-1))/dy;
        CHx_ref(1:Nx,1)   =(Hz_ref(1:Nx,1)-Hz_ref(1:Nx,Ny))/dy;
        
        %Update Integration terms
        %IDx=IDx + Dx;
        ICHx_ref=ICHx_ref + CHx_ref;
        
        %update Dx
        Dx_ref=mDx1_ref.*Dx_ref + mDx2_ref.*CHx_ref + mDx3_ref.*ICHx_ref;
        
        %calculate CHy
        CHy_ref(2:Nx,1:Ny)=-(Hz_ref(2:Nx,1:Ny)-Hz_ref(1:Nx-1,1:Ny))/dx;
        CHy_ref(1,1:Ny)=-(Hz_ref(1,1:Ny)-Hz_ref(Nx,1:Ny))/dx;
        
        %Update Integration terms
        ICHy_ref=ICHy_ref + CHy_ref;
        
        %update Dy
        Dy_ref=mDy1_ref.*Dy_ref + mDy2_ref.*CHy_ref + mDy3_ref.*ICHy_ref;
                
        %update Ex,Ey
        Ex_ref=mEx_ref.*Dx_ref;
        Ey_ref=mEy_ref.*Dy_ref;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% inject source. Since it is same source, just inject.
        % inject source
        %Ex_ref(:,ny_src)= Ex_ref(:,ny_src)+g(T)-Ex_1D_receiver(:,1);
        if angle==0
            Ex_ref(:,ny_src)= Ex_ref(:,ny_src)+g(T)-Ez_1D_receiver(:,1);
        else
            Ex_ref(:,ny_src)= g(:,T);
        end 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if mod(T,1)==0 && plot==1
            %plot the second figure of reference area
            subplot(1,3,2)
            %% plot boundaries without devices
            imagesc(xa,ya,real(Ex_ref).');
            set(gca,'Ydir','normal');
            title(['STEP' num2str(T) ' of ' num2str(STEPS)]);
            hold on;
            if NPML(1)
                x1 = xa(1);
                x2 = xa(NPML(1));
                y1 = ya(1);
                y2 = ya(Ny);
                x  = [ x1 x2 x2 x1 x1 ];
                y  = [ y1 y1 y2 y2 y1 ];
                c  = 0.5*[1 1 1];
                fill (x,y,c,'FaceAlpha',0.5);
                
            end
            if NPML(2)
                x1 = xa(Nx-NPML(2));
                x2 = xa(Nx);
                y1 = ya(1);
                y2 = ya(Ny);
                x  = [ x1 x2 x2 x1 x1 ];
                y  = [ y1 y1 y2 y2 y1 ];
                c  = 0.5*[1 1 1];
                fill (x,y,c,'FaceAlpha',0.5);
                
            end
            if NPML(3)
                x1 = xa(1);
                x2 = xa(Nx);
                y1 = ya(1);
                y2 = ya(NPML(3));
                
                x  = [ x1 x2 x2 x1 x1 ];
                y  = [ y1 y1 y2 y2 y1 ];
                c  = 0.5*[1 1 1];
                fill (x,y,c,'FaceAlpha',0.5);
                
            end
            if NPML(4)
                x1 = xa(1);
                x2 = xa(Nx);
                y1 = ya(Ny - NPML(4));
                y2 = ya(Ny);
                
                x  = [ x1 x2 x2 x1 x1 ];
                y  = [ y1 y1 y2 y2 y1 ];
                c  = 0.5*[1 1 1];
                fill (x,y,c,'FaceAlpha',0.5);
                
            end
            hold off;
            %colorbar;
            caxis(color_range);
            goodplot();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %plot the second figure
            subplot(1,3,3)
            imagesc(xa,ya,real(Ex-Ex_ref).');
            set(gca,'Ydir','normal');
            title(['STEP' num2str(T) ' of ' num2str(STEPS)]);
            caxis(color_range);
            %Electricfieldsaver(:,:,T)=((Ez-Ez_ref));
            goodplot();
            %force graphics to show
            %drawnow;
            movieVector(T)=getframe(f);
        end
        %entropy(T,1)=entropy(T,1)+sum(Ez(1:Nx,1:Ny).^2,"all")^2/sum(Ez(1:Nx,1:Ny).^4,"all");
        %toc;
        %break;
    end

    %
    %% MAIN LOOP
    %
    %Electricfieldsaver=zeros(Nx,Ny,STEPS);
    %movieVector(1:STEPS) = struct('cdata',[],'colormap',[]);
  
    %f=figure('Color','w','Position',[100 100 800 800],'visible','off');
    %S=0.5;
    %alpha=1/S*1/(5*S-40*S^2+126*S^5-160*S^7+70*S^9);
    set(gcf, 'Position', get(0, 'Screensize'));
    entropy=zeros(STEPS,1);
end


% %if finished normally, display this
disp('DONE2!');
if TE==1
    Enew=Ez-Ez_ref;
    Hxnew=Hx-Hx_ref;
    Hynew=Hy-Hy_ref;
    save('electricfield.mat','Enew');
    save('electricmagnetic.mat','Enew','Hxnew','Hynew');
    save('electricfieldnosubstract.mat','Ez','Hx','Hy');
end
if TM==1
    Exnew=Ex-Ex_ref;
    Hznew=Hz-Hz_ref;
    Eynew=Ey-Ey_ref;
    save('electricfield.mat','Exnew');
    save('electricmagnetic.mat','Exnew','Hznew','Eynew');
    save('electricfieldnosubstract.mat','Ex','Hz','Ey');
end


%% 
% Out put Video
% Create a VideoWriter object and set properties
myWriter = VideoWriter('2dfdtdwithdevice','MPEG-4');            %create an .avi file
% myWriter = VideoWriter('curve','MPEG-4');   %create an .mp4 file
myWriter.FrameRate = 20;
%Open the VideoWriter object, write the movie, and close the file
open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ez,Hx,Hy,ICEx,ICEy,IDz,Dz,Nx,Ny,dx,dy,mHx1,mHx2,mHx3,mDz1,mDz2,mDz4,mHy1,mHy2,mHy3,mEz]=fdtd_TE(Ez,Hx,Hy,ICEx,ICEy,IDz,Dz,Nx,Ny,dx,dy,mHx1,mHx2,mHx3,mDz1,mDz2,mDz4,mHy1,mHy2,mHy3,mEz)
    %Compute CEx, which is related to dHx/dt, and is related to
    %electric field dEz/dy+

    
    Ezdy(1:Nx,1:Ny-1)=Ez(1:Nx,2:Ny+1-1);
   
    Ezdy(1:Nx,Ny)=Ez(1:Nx,1);
    
    CEx=(Ezdy-Ez)./dy;
    
    
    %Compute CEy, which is related to dHy/dt, and is related to
    %electric field dEz/dx
    %Update Integration terms
     ICEx=ICEx + CEx;

    Ezdx(1:Nx-1,1:Ny)=Ez(2:Nx+1-1,1:Ny);
    Ezdx(Nx,1:Ny)    =Ez(1,1:Ny);

    %Update Integration terms
    CEy = -(Ezdx-Ez)./dx;
    ICEy=ICEy + CEy;
    
    %Update Hx and Hy

    
    Hx = mHx1.*Hx + mHx2.*CEx + mHx3.*ICEx;
    
    Hy = mHy1.*Hy + mHy2.*CEy + mHy3.*ICEy;
    
    
    %Calculate CHz
   
    CHz(1,1)=(Hy(1,1)-Hy(Nx,1))/dx-(Hx(1,1)-Hy(1,Ny))/dy;  
    
      CHz(2:Nx,1)=(Hy(2:Nx,1)-Hy(1:Nx-1,1))/dx ...
            -(Hx(2:Nx,1)-Hx(2:Nx,Ny))/dy;

      CHz(1,2:Ny)=(Hy(1,2:Ny)-Hy(Nx,2:Ny))/dx ...            %%PBC is applied Hy(Nx,ny) term
            -(Hx(1,2:Ny)-Hx(1,1:Ny-1))/dy;
        
      CHz(2:Nx,2:Ny)=(Hy(2:Nx,2:Ny)-Hy(1:Nx-1,2:Ny))/dx ...
                -(Hx(2:Nx,2:Ny)-Hx(2:Nx,1:Ny-1))/dy;
    %Update integration term
    IDz = IDz + Dz;
    
    
    %update Dz
    Dz = mDz1.* Dz + mDz2 .* CHz + mDz4.*IDz;
    
    
    %update Ez
    Ez=mEz.*Dz;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ez_1D,H_1D,mH_1D,mEz_1D,dx,Nx]=fdtd_1D(Ez_1D,H_1D,mH_1D,mEz_1D,dx,Nx)
     
        %for nx = 1 : Nx-1
            H_1D(1 : Nx-1) = H_1D(1 : Nx-1) + mH_1D(1 : Nx-1).*(Ez_1D(2 : Nx) - Ez_1D(1 : Nx-1))./dx;
        %end
        H_1D(Nx) = H_1D(Nx) + mH_1D(Nx)*(Ez_1D(1) - Ez_1D(Nx))/dx;
        % Update E from H (Dirichlet Boundary Conditions)
        Ez_1D(1) = Ez_1D(1) + mEz_1D(1)*(H_1D(1) - H_1D(Nx))/dx;
        %for nx = 2 : Nx
            Ez_1D( 2 : Nx) = Ez_1D( 2 : Nx)  + mEz_1D( 2 : Nx).*(H_1D( 2 : Nx) - H_1D( 1 : Nx-1))./dx;
            
       
end
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
