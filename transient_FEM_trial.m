function data=transient_FEM_trial(Lx,Ly,nelx,nely,t_target)
%TRANSIENT_FEM_TRIAL A Simple demo of Transient 2D Plate Conduction Solver using FEM & Implict Finite Difference Method. Adopted from top88.m by Ole Sigmund et al.
%   Input: Lx: length of plate in x-direction
%   Input: Ly: length of plate in y-direction
%   Input: nelx: FEM discretization in x-direction
%   Input: nely: FEM discretization in y-direction
%   Input: t_target: target time to reach
%   WARNING: Please MAKE SURE that Lx/Ly==nelx/nely, or the rendering of results may go wrong!

a=Lx/nelx;
b=Ly/nely;

%% MATERIAL PROPERTIES OF THE BASIC MATERIAL(SS316)
Kalpa=16.3;%thermal conductivity of basic material(isotropic, unit: W/(m.K))
c=0.45;%density of basic material(unit: J/(kg.K))
rho=7980;%density of basic material(unit: kg/(m^3))

%% PREPARE FINITE ELEMENT ANALYSIS
% KE = [ 2/3 -1/6 -1/3 -1/6
%        -1/6 2/3 -1/6 -1/3
%        -1/3 -1/6 2/3 -1/6
%        -1/6 -1/3 -1/6 2/3];
KE1=[b/(3*a) -b/(3*a) -b/(6*a) b/(6*a)
    -b/(3*a) b/(3*a) b/(6*a) -b/(6*a)
    -b/(6*a) b/(6*a) b/(3*a) -b/(3*a)
    b/(6*a) -b/(6*a) -b/(3*a) b/(3*a)];
KE2=[a/(3*b) a/(6*b) -a/(6*b) -a/(3*b)
    a/(6*b) a/(3*b) -a/(3*b) -a/(6*b)
    -a/(6*b) -a/(3*b) a/(3*b) a/(6*b)
    -a/(3*b) -a/(6*b) a/(6*b) a/(3*b)];
KE=KE1+KE2;
ME=0.25*a*b*[4/9 2/9 1/9 2/9
        2/9 4/9 2/9 1/9
        1/9 2/9 4/9 2/9
        2/9 1/9 2/9 4/9];

nodenrs=reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec=reshape(nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat=repmat(edofVec,1,4)+repmat([0,nely+[1,0],-1],nelx*nely,1);
iK=reshape(kron(edofMat,ones(4,1))',16*nelx*nely,1);
jK=reshape(kron(edofMat,ones(1,4))',16*nelx*nely,1);
iM=iK;
jM=jK;

%fixeddofs = [nely/2+1-(nely/10):nely/2+1+(nely/10)];
%fixeddofs=[];
%fixeddofs1 = [(nely+1)*nelx+nely/2+1-(nely/20):(nely+1)*nelx+nely/2+1+(nely/20)];
fixeddofs1 = [1:(nely+1)];
%fixeddofs1=[];
%fixeddofs2 = [(nely+1)*nelx+nely/2+1-(nely/20):(nely+1)*nelx+nely/2+1+(nely/20)];
fixeddofs2 = [(nely+1)*nelx+1:(nely+1)*(nelx+1)];
%fixeddofs = [(nely+1)*nelx+1:(nely+1)*(nelx+1)];
alldofs = [1:(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,[fixeddofs1,fixeddofs2]);

T = zeros((nelx+1)*(nely+1),1);
Q = sparse((nely+1)*(nelx+1),1); %External heat flux input, unit: K/(m^2)
Q(:,1) = 0;
%Q((nely/2)*(nelx+1)+nely/2+1,1) = 1;
%Q(nely/2+1-(nely/5):nely/2+1+(nely/5),1) = 1;
%% INITIALIZE ITERATION
xPhys=ones(nely,nelx);

sK=reshape(KE(:)*(xPhys(:)'*Kalpa),16*nelx*nely,1);
K=sparse(iK,jK,sK);
K=(K+K')/2;

sM=reshape(ME(:)*(xPhys(:)'*rho*c),16*nelx*nely,1);
M=sparse(iM,jM,sM);
M=(M+M')/2;

t_current=0;
T_lasttimestep = 273.15*ones((nelx+1)*(nely+1),1);% Initial temperature distribution, unit: K
T_lasttimestep(fixeddofs1,:)=373.15;
T_lasttimestep(fixeddofs2,:)=273.15;
dt=0.1;%In Finite Different methods, dt should always BE SMALLER THAN 1, or it will be unstable
beta=1;
data=[T_lasttimestep];

T_figure=reshape(T_lasttimestep,(1+nely),(1+nelx));
xdots=repmat([0:(nelx)],(nely+1),1);
ydots=repmat([nely:-1:0]',1,(nelx+1));
contourf(xdots,ydots,T_figure,50,'LineStyle','none');
shading interp;
colorbar;
colormap("jet");
axis equal;
axis off;
title('Temperature at time=',t_current);
pause(dt);

while t_current<t_target
    %% FE-ANALYSIS
    Kt=M+beta*dt*K;
    R=dt*Q+(M-(1-beta)*dt*K)*T_lasttimestep;
    T(fixeddofs1,:)= 373.15;
    T(fixeddofs2,:)= 273.15;
    T(freedofs,:) = Kt(freedofs,freedofs) \ (R(freedofs,:)-Kt(freedofs,fixeddofs1)*T(fixeddofs1,:)-Kt(freedofs,fixeddofs2)*T(fixeddofs2,:));
    
    %% PLOT RESULTS
    T_figure=reshape(T,(1+nely),(1+nelx));
    xdots=repmat([0:(nelx)],(nely+1),1);
    ydots=repmat([nely:-1:0]',1,(nelx+1));
    contourf(xdots,ydots,T_figure,50,'LineStyle','none');
    shading interp;
    colorbar;
    colormap("jet");
    axis equal;
    axis off;
    title('Temperature at time=',t_current+dt);
    pause(dt);

    t_current=t_current+dt;
    data=[data,T];
    T_lasttimestep=T;
end
end

