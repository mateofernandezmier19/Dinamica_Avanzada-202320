% Implementation of a Kinematic Model for a Passive Walker with
% hemispherical feet
% In this part, the evolution of the mass properties of the passive walker
% is considered. This include:
% * Evolution of the center of mass of the different bodies
% * Evolution of the inertia scalars of the bodies with respect to fixed
%   points in the bodies
% * Evolution of the inertia scalars of the bodies with respect to fixed
%   points

clear, clc, close all
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Definition of geometrical parameters
R=0.48939; %Radious of the sphere [m]
Lp=0.4418; %Length of each leg [m]
d3=35/1000;%Distance
drl3=30.236/1000; %Distance
drl1=316.887/1000; %Distance

%% Definition of the mass properties of the bodies
mD=0.444; %Mass of body D [kg]
IDD=[9252.59,0,0;0,7.990,0;0,0,9252.059]*(1e-6); %Inertia Matrix of body D (with respect to its center of mass) [kg*m^2]
mE=0.200; %Mass of body E [kg]
IEE=[4489.520,-0.011,-0.238;-0.011,4647.233,-361.076;-0.238,-361.076,532.433]*(1e-6); %Inertia Matrix of body E (with respect to its center of mass) [kg*m^2]
mF=0.200; %Mass of body F [kg]
IFF=[4489.520,-0.011,-0.238;-0.011,4647.233,-361.076;-0.238,-361.076,532.433]*(1e-6); %Inertia Matrix of body F (with respect to its center of mass) [kg*m^2]

%% Definition of the generalized coordinates and their derivatives
dt=0.001;
t=0:dt:10;
q1=zeros(1,length(t));
q2=zeros(1,length(t));
q3=zeros(1,length(t));
q4=zeros(1,length(t));
q5=zeros(1,length(t));
q6=zeros(1,length(t));
q7=zeros(1,length(t));
q8=zeros(1,length(t));
qx=zeros(1,length(t));
qy=zeros(1,length(t));

q1_d=zeros(1,length(t));
q2_d=zeros(1,length(t));
q3_d=zeros(1,length(t));
q4_d=zeros(1,length(t));
q5_d=zeros(1,length(t));
q6_d=zeros(1,length(t));
q7_d=zeros(1,length(t));
q8_d=zeros(1,length(t));
qx_d=zeros(1,length(t));
qy_d=zeros(1,length(t));

q1_dd=zeros(1,length(t));
q2_dd=zeros(1,length(t));
q3_dd=zeros(1,length(t));
q4_dd=zeros(1,length(t));
q5_dd=zeros(1,length(t));
q6_dd=zeros(1,length(t));
q7_dd=zeros(1,length(t));
q8_dd=zeros(1,length(t));
qx_dd=zeros(1,length(t));
qy_dd=zeros(1,length(t));

%% Initial conditions
q1(1)=0;
q2(1)=0;
q3(1)=0;
% q3(1)=Lp;
% q3(2)=Lp;
q4(1)=0;
q5(1)=0;
q6(1)=0;
q7(1)=0;
q8(1)=0;
qx(1)=0;
qy(1)=0;

q1_d(1)=0;
q2_d(1)=0;
q3_d(1)=0;
q4_d(1)=0;
q5_d(1)=0;
q6_d(1)=0;
q7_d(1)=0;
q8_d(1)=0;
qx_d(1)=0;
qy_d(1)=0;

q1_dd(1)=0;
q2_dd(1)=0;
q3_dd(1)=0;
q4_dd(1)=0;
q5_dd(1)=0;
q6_dd(1)=0;
q7_dd(1)=0;
q8_dd(1)=0;
qx_dd(1)=0;
qy_dd(1)=0;


%% Algorithm for the kinematically driven calculations

%Mass center of the bodies
rDstarO=cell(1,length(t));
rEstarO=cell(1,length(t));
rFstarO=cell(1,length(t));
rSO=cell(1,length(t));

%Inertia scalars with respect to fixed points in the bodies
DID1=cell(1,length(t));
EID1=cell(1,length(t));
FID2=cell(1,length(t));

%Inertia scalars with respect to fixed points
DIO=cell(1,length(t));
EIO=cell(1,length(t));
FIO=cell(1,length(t));

%Angular velocities
AwB=cell(1,length(t));
BwC=cell(1,length(t));
CwD=cell(1,length(t));
DwE=cell(1,length(t));
DwF=cell(1,length(t));

%Angular velocities in the frame A
AwC=cell(1,length(t));
AwD=cell(1,length(t));
AwE=cell(1,length(t));
AwF=cell(1,length(t));

%Angular aceleration
AalphaB=cell(1,length(t));
BalphaC=cell(1,length(t));
CalphaD=cell(1,length(t));
DalphaE=cell(1,length(t));
DalphaF=cell(1,length(t));

%Angular velocities in the frame A
AalphaC=cell(1,length(t));
AalphaD=cell(1,length(t));
AalphaE=cell(1,length(t));
AalphaF=cell(1,length(t));

%Velocity and acceleration of point Dstar
AvDstar=cell(1,length(t));
AaDstar=cell(1,length(t));

%Velocity and acceleration of point D1
AvD1=cell(1,length(t));
AaD1=cell(1,length(t));

%Velocity and acceleration of point D2
AvD2=cell(1,length(t));
AaD2=cell(1,length(t));

%Velocity and acceleration of point Estar
AvEstar=cell(1,length(t));
AaEstar=cell(1,length(t));

%Velocity and acceleration of point Fstar
AvFstar=cell(1,length(t));
AaFstar=cell(1,length(t));

for i=2:length(t)-1         
    q4_d(i)=0;
    q5_d(i)=0;
    q6_d(i)=0;
    q7_d(i)=0;
    q8_d(i)=0;
    
    q4_dd(i)=0;
    q5_dd(i)=0;
    q6_dd(i)=0;
    q7_dd(i)=0;
    q8_dd(i)=0;

    %Update of the generalized coordinates
    q4(i+1)=q4(i)+dt*q4_d(i);
    q5(i+1)=q5(i)+dt*q5_d(i);
    q6(i+1)=q6(i)+dt*q6_d(i);
    q7(i+1)=q7(i)+dt*q7_d(i);
    q8(i+1)=q8(i)+dt*q8_d(i);

    if q5(i)>=0
        %Udpate of the position and velocity of point Dstar when point Ehat is the
        %contact point
%         q1_d(i)=(-q4_d(i)*cos(q4(i))*cos(q5(i))+q5_d(i)*sin(q5(i))*sin(q4(i)))*d3+Lp*((q6_d(i)+q7_d(i))*cos(q4(i))*(cos(q5(i))^2)*cos(q6(i)+q7(i))-(q6_d(i)+q7_d(i))*sin(q5(i))*(sin(q4(i))*sin(q6(i)+q7(i))-cos(q4(i))*sin(q5(i))*cos(q6(i)+q7(i))));
%         q2_d(i)=(-q4_d(i)*sin(q4(i))*cos(q5(i))-q5_d(i)*sin(q5(i))*cos(q4(i)))*d3+Lp*(-(q6_d(i)+q7_d(i))*sin(q4(i))*(cos(q5(i))^2)*cos(q6(i)+q7(i))+(q6_d(i)+q7_d(i))*sin(q5(i))*(cos(q4(i))*sin(q6(i)+q7(i))+sin(q4(i))*sin(q5(i))*cos(q6(i)+q7(i))));
%         q3_d(i)=q5_d(i)*cos(q5(i))*d3+Lp*(-(q6_d(i)+q7_d(i))*sin(q4(i))*cos(q5(i))*(sin(q4(i))*sin(q6(i)+q7(i))-cos(q4(i))*sin(q5(i))*cos(q6(i)+q7(i)))-(q6_d(i)+q7_d(i))*cos(q4(i))*cos(q5(i))*(cos(q4(i))*sin(q6(i)+q7(i))+sin(q4(i))*sin(q5(i))*cos(q6(i)+q7(i))));
        qx_d(i)=q5_d(i)*R*sin(q4(i))+(q6_d(i)+q7_d(i))*R*cos(q4(i))*cos(q5(i));
        qy_d(i)=-q5_d(i)*R*cos(q4(i))+(q6_d(i)+q7_d(i))*R*sin(q4(i))*cos(q5(i));
    else
        %Udpate of the position and velocity of point Dstar when point Fhat is the
        %contact point
%         q1_d(i)=(q4_d(i)*cos(q4(i))*cos(q5(i))-q5_d(i)*sin(q5(i))*sin(q4(i)))*d3+Lp*((q6_d(i)+q8_d(i))*cos(q4(i))*(cos(q5(i))^2)*cos(q6(i)+q8(i))-(q6_d(i)+q8_d(i))*sin(q5(i))*(sin(q4(i))*sin(q6(i)+q8(i))-cos(q4(i))*sin(q5(i))*cos(q6(i)+q8(i))));
%         q2_d(i)=(q4_d(i)*sin(q4(i))*cos(q5(i))+q5_d(i)*sin(q5(i))*cos(q4(i)))*d3+Lp*(-(q6_d(i)+q8_d(i))*sin(q4(i))*(cos(q5(i))^2)*cos(q6(i)+q8(i))+(q6_d(i)+q8_d(i))*sin(q5(i))*(cos(q4(i))*sin(q6(i)+q8(i))+sin(q4(i))*sin(q5(i))*cos(q6(i)+q8(i))));
%         q3_d(i)=-q5_d(i)*cos(q5(i))*d3+Lp*(-(q6_d(i)+q8_d(i))*sin(q4(i))*cos(q5(i))*(sin(q4(i))*sin(q6(i)+q8(i))-cos(q4(i))*sin(q5(i))*cos(q6(i)+q8(i)))-(q6_d(i)+q8_d(i))*cos(q4(i))*cos(q5(i))*(cos(q4(i))*sin(q6(i)+q8(i))+sin(q4(i))*sin(q5(i))*cos(q6(i)+q8(i))));
        qx_d(i)=q5_d(i)*R*sin(q4(i))+(q6_d(i)+q8_d(i))*R*cos(q4(i))*cos(q5(i));
        qy_d(i)=-q5_d(i)*R*cos(q4(i))+(q6_d(i)+q8_d(i))*R*sin(q4(i))*cos(q5(i));
    end

    %Rotation matrixes
    aRb=[cos(q4(i)),-sin(q4(i)),0;sin(q4(i)),cos(q4(i)),0;0,0,1];
    bRc=[1,0,0;0,cos(q5(i)),-sin(q5(i));0,sin(q5(i)),cos(q5(i))];
    cRd=[cos(q6(i)),0,sin(q6(i));0,1,0;-sin(q6(i)),0,cos(q6(i))];
    dRe=[cos(q7(i)),0,sin(q7(i));0,1,0;-sin(q7(i)),0,cos(q7(i))];
    dRf=[cos(q8(i)),0,sin(q8(i));0,1,0;-sin(q8(i)),0,cos(q8(i))];
    %Rotation matrixes in terms of the frame A
    aRc=aRb*bRc;
    aRd=aRc*cRd;
    aRe=aRd*dRe;
    aRf=aRd*dRf;

    %Angular velocities
    AwB{i}=q4_d(i)*[0;0;1];
    BwC{i}=q5_d(i)*[1;0;0];
    CwD{i}=q6_d(i)*[0;1;0];
    DwE{i}=q7_d(i)*[0;1;0];
    DwF{i}=q8_d(i)*[0;1;0];
    %Angular velocities in frame A
    AwC{i}=AwB{i}+aRb*BwC{i};
    AwD{i}=AwC{i}+aRc*CwD{i};
    AwE{i}=AwD{i}+aRd*DwE{i};
    AwF{i}=AwD{i}+aRd*DwF{i};

    %Angular Acelerations
    AalphaB{i}=q4_dd(i)*[0;0;1];
    BalphaC{i}=q5_dd(i)*[1;0;0];
    CalphaD{i}=q6_dd(i)*[0;1;0];
    DalphaE{i}=q7_dd(i)*[0;1;0];
    DalphaF{i}=q5_dd(i)*[0;1;0];
    %Angular Acelerations in frame A
    AalphaC{i}=AalphaB{i}+aRb*BalphaC{i}+cross(AwB{i},BwC{i});
    AalphaD{i}=AalphaC{i}+aRc*CalphaD{i}+cross(AwC{i},CwD{i});
    AalphaE{i}=AalphaD{i}+aRd*DalphaE{i}+cross(AwD{i},DwE{i});
    AalphaF{i}=AalphaD{i}+aRd*DalphaF{i}+cross(AwD{i},DwF{i});

    %Velocity of point Dstar
    AvDstar{i}=[qx_d(i);qy_d(i);0]+cross(AwE{i},aRe*Lp*[0;0;1])+cross(AwD{i},aRd*d3*[0;1;0]);
    q1_d(i)=AvDstar{i}(1);
    q2_d(i)=AvDstar{i}(2);
    q3_d(i)=AvDstar{i}(3);

    q1(i+1)=q1(i)+dt*q1_d(i-1);
    q2(i+1)=q2(i)+dt*q2_d(i-1);
    q3(i+1)=q3(i)+dt*q3_d(i-1);

    q1_dd(i)=(q1_d(i)-q1_d(i-1))/dt;
    q2_dd(i)=(q2_d(i)-q2_d(i-1))/dt;
    q3_dd(i)=(q3_d(i)-q3_d(i-1))/dt;


    %Center of mass of the bodies
    rDstarO{i}=[q1(i);q2(i);q3(i)];
    rEstarO{i}=rDstarO{i}-(d3+drl3)*aRe*[0;1;0]-drl1*aRe*[0;0;1];
    rFstarO{i}=rDstarO{i}+(d3+drl3)*aRf*[0;1;0]-drl1*aRf*[0;0;1];
    rSO{i}=(1/(mD+mE+mF))*(mD*rDstarO{i}+mE*rEstarO{i}+mF*rFstarO{i});

    %Inertia Scalars with respect to points fixed in the bodies
    DID1{i}=(aRd*IDD*transpose(aRd))+mD*((d3^2)*eye(3)-(d3^2)*(aRd*[0;1;0])*transpose(aRd*[0;1;0]));
    EID1{i}=(aRe*IEE*transpose(aRe))+mE*((dot([0;-drl3;-drl1],[0;-drl3;-drl1]))*eye(3)-(aRe*[0;-drl3;-drl1])*transpose(aRe*[0;-drl3;-drl1]));
    FID2{i}=(aRf*IFF*transpose(aRf))+mF*((dot([0;drl3;-drl1],[0;drl3;-drl1]))*eye(3)-(aRf*[0;drl3;-drl1])*transpose(aRf*[0;drl3;-drl1]));

    %Inertia Scalars with respect to FIXED points
    DIO{i}=(aRd*IDD*transpose(aRd))+mD*((dot(rDstarO{i},rDstarO{i}))*eye(3)-(rDstarO{i})*transpose(rDstarO{i}));
    EIO{i}=(aRe*IEE*transpose(aRe))+mE*((dot(rEstarO{i},rEstarO{i}))*eye(3)-(rEstarO{i})*transpose(rEstarO{i}));
    FIO{i}=(aRf*IFF*transpose(aRf))+mF*((dot(rFstarO{i},rFstarO{i}))*eye(3)-(rFstarO{i})*transpose(rFstarO{i}));
end


%% Figure

%Evolution of the center of mass of body D
rdstarO=cell2mat(rDstarO);

figure
plot(t(1:end-2),rdstarO(1,:),t(1:end-2),rdstarO(2,:),t(1:end-2),rdstarO(3,:))
xlabel("$t\;(s)$")
ylabel("$\underline{r}^{D^{*}/O}$")
title("\textbf{Evolution of the center of mass of body D}")
legend(["$r^{D^{*}/O}_{x}$","$r^{D^{*}/O}_{y}$","$r^{D^{*}/O}_{z}$"])
grid on
set(gcf,'color','w')

%% Evolution of the center of mass of body E
restarO=cell2mat(rEstarO);

figure
plot(t(1:end-2),restarO(1,:),t(1:end-2),restarO(2,:),t(1:end-2),restarO(3,:))
xlabel("$t\;(s)$")
ylabel("$\underline{r}^{E^{*}/O}$")
title("\textbf{Evolution of the center of mass of body E}")
legend(["$r^{E^{*}/O}_{x}$","$r^{E^{*}/O}_{y}$","$r^{E^{*}/O}_{z}$"])
grid on
set(gcf,'color','w')

%% Evolution of the center of mass of body F
rfstarO=cell2mat(rFstarO);

figure
plot(t(1:end-2),rfstarO(1,:),t(1:end-2),rfstarO(2,:),t(1:end-2),rfstarO(3,:))
xlabel("$t\;(s)$")
ylabel("$\underline{r}^{F^{*}/O}$")
title("\textbf{Evolution of the center of mass of body F}")
legend(["$r^{F^{*}/O}_{x}$","$r^{F^{*}/O}_{y}$","$r^{F^{*}/O}_{z}$"])
grid on
set(gcf,'color','w')

%% Evolution of the center of mass of the system
rsO=cell2mat(rSO);

figure
plot(t(1:end-2),rsO(1,:),t(1:end-2),rsO(2,:),t(1:end-2),rsO(3,:))
xlabel("$t\;(s)$")
ylabel("$\underline{r}^{S/O}$")
title("\textbf{Evolution of the center of mass of the system}")
legend(["$r^{S/O}_{x}$","$r^{S/O}_{y}$","$r^{S/O}_{z}$"])
grid on
set(gcf,'color','w')

%% Evolution of the inertia scalars of body D with respect to D1
id11=zeros(1,length(t));
id12=zeros(1,length(t));
id13=zeros(1,length(t));
id22=zeros(1,length(t));
id23=zeros(1,length(t));
id33=zeros(1,length(t));

for k=2:length(t)-1
    id11(k)=DID1{k}(1,1);
    id12(k)=DID1{k}(1,2);
    id13(k)=DID1{k}(1,3);
    id22(k)=DID1{k}(2,2);
    id23(k)=DID1{k}(2,3);
    id33(k)=DID1{k}(3,3);
end
figure
plot(t,id11,t,id12,t,id13,t,id22,t,id23,t,id33)
xlabel("$t\;(s)$")
ylabel("\textit{Inertia Scalars of body D}")
title("\textbf{Evolution of the Inertia Scalars of body D with respecto to point $D_{1}$}")
legend(["$I_{11}^{D/D_{1}}$","$I_{12}^{D/D_{1}}$","$I_{13}^{D/D_{1}}$","$I_{22}^{D/D_{1}}$","$I_{23}^{D/D_{1}}$","$I_{33}^{D/D_{1}}$"])
grid on
set(gcf,'color','w')

%% Evolution of the inertia scalars of body E with respect to D1
ie11=zeros(1,length(t));
ie12=zeros(1,length(t));
ie13=zeros(1,length(t));
ie22=zeros(1,length(t));
ie23=zeros(1,length(t));
ie33=zeros(1,length(t));

for k=2:length(t)-1
    ie11(k)=EID1{k}(1,1);
    ie12(k)=EID1{k}(1,2);
    ie13(k)=EID1{k}(1,3);
    ie22(k)=EID1{k}(2,2);
    ie23(k)=EID1{k}(2,3);
    ie33(k)=EID1{k}(3,3);
end
figure
plot(t,ie11,t,ie12,t,ie13,t,ie22,t,ie23,t,ie33)
xlabel("$t\;(s)$")
ylabel("\textit{Inertia Scalars of body E}")
title("\textbf{Evolution of the Inertia Scalars of body E with respecto to point $D_{1}$}")
legend(["$I_{11}^{E/D_{1}}$","$I_{12}^{E/D_{1}}$","$I_{13}^{E/D_{1}}$","$I_{22}^{E/D_{1}}$","$I_{23}^{E/D_{1}}$","$I_{33}^{E/D_{1}}$"])
grid on
set(gcf,'color','w')

%% Evolution of the inertia scalars of body F with respect to D2
if11=zeros(1,length(t));
if12=zeros(1,length(t));
if13=zeros(1,length(t));
if22=zeros(1,length(t));
if23=zeros(1,length(t));
if33=zeros(1,length(t));

for k=2:length(t)-1
    if11(k)=FID1{k}(1,1);
    if12(k)=FID1{k}(1,2);
    if13(k)=FID1{k}(1,3);
    if22(k)=FID1{k}(2,2);
    if23(k)=FID1{k}(2,3);
    if33(k)=EID1{k}(3,3);
end
figure
plot(t,if11,t,if12,t,if13,t,if22,t,if23,t,if33)
xlabel("$t\;(s)$")
ylabel("\textit{Inertia Scalars of body F}")
title("\textbf{Evolution of the Inertia Scalars of body F with respecto to point $D_{1}$}")
legend(["$I_{11}^{F/D_{1}}$","$I_{12}^{F/D_{1}}$","$I_{13}^{F/D_{1}}$","$I_{22}^{F/D_{1}}$","$I_{23}^{F/D_{1}}$","$I_{33}^{F/D_{1}}$"])
grid on
set(gcf,'color','w')

%% Evolution of the inertia scalars of body D with respect to point O
id11O=zeros(1,length(t));
id12O=zeros(1,length(t));
id13O=zeros(1,length(t));
id22O=zeros(1,length(t));
id23O=zeros(1,length(t));
id33O=zeros(1,length(t));

for k=2:length(t)-1
    id11O(k)=DIO{k}(1,1);
    id12O(k)=DIO{k}(1,2);
    id13O(k)=DIO{k}(1,3);
    id22O(k)=DIO{k}(2,2);
    id23O(k)=DIO{k}(2,3);
    id33O(k)=DIO{k}(3,3);
end
figure
plot(t,id11O,t,id12O,t,id13O,t,id22O,t,id23O,t,id33O)
xlabel("$t\;(s)$")
ylabel("\textit{Inertia Scalars of body D}")
title("\textbf{Evolution of the Inertia Scalars of body D with respecto to point $O$}")
legend(["$I_{11}^{D/O}$","$I_{12}^{D/O}}$","$I_{13}^{D/O}$","$I_{22}^{D/O}$","$I_{23}^{D/O}$","$I_{33}^{D/O}$"])
grid on
set(gcf,'color','w')

%% Evolution of the inertia scalars of body E with respect to point O
ie11O=zeros(1,length(t));
ie12O=zeros(1,length(t));
ie13O=zeros(1,length(t));
ie22O=zeros(1,length(t));
ie23O=zeros(1,length(t));
ie33O=zeros(1,length(t));

for k=2:length(t)-1
    ie11O(k)=EIO{k}(1,1);
    ie12O(k)=EIO{k}(1,2);
    ie13O(k)=EIO{k}(1,3);
    ie22O(k)=EIO{k}(2,2);
    ie23O(k)=EIO{k}(2,3);
    ie33O(k)=EIO{k}(3,3);
end
figure
plot(t,ie11O,t,ie12O,t,ie13O,t,ie22O,t,ie23O,t,ie33O)
xlabel("$t\;(s)$")
ylabel("\textit{Inertia Scalars of body E}")
title("\textbf{Evolution of the Inertia Scalars of body E with respecto to point $O$}")
legend(["$I_{11}^{E/O}$","$I_{12}^{E/O}}$","$I_{13}^{E/O}$","$I_{22}^{E/O}$","$I_{23}^{E/O}$","$I_{33}^{E/O}$"])
grid on
set(gcf,'color','w')

%% Evolution of the inertia scalars of body F with respect to point O
if11O=zeros(1,length(t));
if12O=zeros(1,length(t));
if13O=zeros(1,length(t));
if22O=zeros(1,length(t));
if23O=zeros(1,length(t));
if33O=zeros(1,length(t));

for k=2:length(t)-1
    if11O(k)=FIO{k}(1,1);
    if12O(k)=FIO{k}(1,2);
    if13O(k)=FIO{k}(1,3);
    if22O(k)=FIO{k}(2,2);
    if23O(k)=FIO{k}(2,3);
    if33O(k)=FIO{k}(3,3);
end
figure
plot(t,if11O,t,if12O,t,if13O,t,if22O,t,if23O,t,if33O)
xlabel("$t\;(s)$")
ylabel("\textit{Inertia Scalars of body F}")
title("\textbf{Evolution of the Inertia Scalars of body F with respecto to point $O$}")
legend(["$I_{11}^{F/O}$","$I_{12}^{F/O}}$","$I_{13}^{F/O}$","$I_{22}^{F/O}$","$I_{23}^{F/O}$","$I_{33}^{F/O}$"])
grid on
set(gcf,'color','w')