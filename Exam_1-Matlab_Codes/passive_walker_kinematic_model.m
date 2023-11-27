% Implementation of a Kinematic Model for a Passive Walker with
% hemispherical feet

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

    %Velocity and aceleration of point D1
    AvD1{i}=[q1_d(i);q2_d(i);q3_d(i)]+cross(AwD{i},aRd*d3*[0;-1;0]);
    AaD1{i}=[q1_dd(i);q2_dd(i);q3_dd(i)]+cross(AwD{i},cross(AwD{i},aRd*d3*[0;-1;0]))+cross(AalphaD{i},aRd*d3*[0;-1;0]);

    %Velocity and aceleration of point D2
    AvD2{i}=[q1_d(i);q2_d(i);q3_d(i)]+cross(AwD{i},aRd*d3*[0;1;0]);
    AaD2{i}=[q1_dd(i);q2_dd(i);q3_dd(i)]+cross(AwD{i},cross(AwD{i},aRd*d3*[0;1;0]))+cross(AalphaD{i},aRd*d3*[0;1;0]);

    %Velocity and aceleration of point Estar
    AvEstar{i}=AvD1{i}+cross(AwE{i},aRe*[0;-drl3;-drl1]);
    AaEstar{i}=AaD1{i}+cross(AwE{i},cross(AwE{i},aRe*[0;-drl3;-drl1]))+cross(AalphaE{i},aRe*[0;-drl3;-drl1]);

    %Velocity and aceleration of point Fstar
    AvFstar{i}=AvD2{i}+cross(AwF{i},aRf*[0;drl3;-drl1]);
    AaFstar{i}=AaD2{i}+cross(AwF{i},cross(AwF{i},aRf*[0;drl3;-drl1]))+cross(AalphaF{i},aRf*[0;drl3;-drl1]);

end


%% Figures

%Evolution of the generalized coordinates
figure
plot(t,q4,t,q5,t,q6,t,q7,t,q8)
xlabel("$t\;(s)$")
ylabel("$q_{i}$")
title("\textbf{Evolution of $q_{i}$}")
legend(["$q_{4}$","$q_{5}$","$q_{6}$","$q_{7}$","$q_{8}$"])
grid on
set(gcf,'color','w')

%% Evolution of the time derivative of generalized coordinates
figure
plot(t,q4_d,t,q5_d,t,q6_d,t,q7_d,t,q8_d)
xlabel("$t\;(s)$")
ylabel("$\dot{q}_{i}$")
title("\textbf{Evolution of $\dot{q}_{i}$}")
legend(["$\dot{q}_{4}$","$\dot{q}_{5}$","$\dot{q}_{6}$","$\dot{q}_{7}$","$\dot{q}_{8}$"])
grid on
set(gcf,'color','w')

%% Evolution of the second time derivative of generalized coordinates
figure
plot(t,q4_dd,t,q5_dd,t,q6_dd,t,q7_dd,t,q8_dd)
xlabel("$t\;(s)$")
ylabel("$\ddot{q}_{i}$")
title("\textbf{Evolution of $\ddot{q}_{i}$}")
legend(["$\ddot{q}_{4}$","$\ddot{q}_{5}$","$\ddot{q}_{6}$","$\ddot{q}_{7}$","$\ddot{q}_{8}$"])
grid on
set(gcf,'color','w')

%% Evolution of the angular velocities
ab=cell2mat(AwB);
ac=cell2mat(AwC);
ad=cell2mat(AwD);
ae=cell2mat(AwE);
af=cell2mat(AwF);

figure
subplot(3,2,1)
plot(t(1:end-2),ab(1,:),t(1:end-2),ab(2,:),t(1:end-2),ab(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\omega}^{B}$")
title("\textbf{Evolution of $^{A}\underline{\omega}^{B}$}")
legend(["$^{A}\omega^{B}_{x}$","$^{A}\omega^{B}_{y}$","$^{A}\omega^{B}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,2)
plot(t(1:end-2),ac(1,:),t(1:end-2),ac(2,:),t(1:end-2),ac(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\omega}^{C}$")
title("\textbf{Evolution of $^{A}\underline{\omega}^{C}$}")
legend(["$^{A}\omega^{C}_{x}$","$^{A}\omega^{C}_{y}$","$^{A}\omega^{C}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,3)
plot(t(1:end-2),ad(1,:),t(1:end-2),ad(2,:),t(1:end-2),ad(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\omega}^{D}$")
title("\textbf{Evolution of $^{A}\underline{\omega}^{D}$}")
legend(["$^{A}\omega^{D}_{x}$","$^{A}\omega^{D}_{y}$","$^{A}\omega^{D}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,4)
plot(t(1:end-2),ae(1,:),t(1:end-2),ae(2,:),t(1:end-2),ae(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\omega}^{E}$")
title("\textbf{Evolution of $^{A}\underline{\omega}^{E}$}")
legend(["$^{A}\omega^{E}_{x}$","$^{A}\omega^{E}_{y}$","$^{A}\omega^{E}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,5)
plot(t(1:end-2),af(1,:),t(1:end-2),af(2,:),t(1:end-2),af(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\omega}^{F}$")
title("\textbf{Evolution of $^{A}\underline{\omega}^{F}$}")
legend(["$^{A}\omega^{F}_{x}$","$^{A}\omega^{F}_{y}$","$^{A}\omega^{F}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
sgtitle("\textbf{Evolution of the Angular Velocities}")

%% Evolution of the angular accelerations
aab=cell2mat(AalphaB);
aac=cell2mat(AalphaC);
aad=cell2mat(AalphaD);
aae=cell2mat(AalphaE);
aaf=cell2mat(AalphaF);

figure
subplot(3,2,1)
plot(t(1:end-2),aab(1,:),t(1:end-2),aab(2,:),t(1:end-2),aab(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\alpha}^{B}$")
title("\textbf{Evolution of $^{A}\underline{\alpha}^{B}$}")
legend(["$^{A}\alpha^{B}_{x}$","$^{A}\alpha^{B}_{y}$","$^{A}\alpha^{B}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,2)
plot(t(1:end-2),aac(1,:),t(1:end-2),aac(2,:),t(1:end-2),aac(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\alpha}^{C}$")
title("\textbf{Evolution of $^{A}\underline{\alpha}^{C}$}")
legend(["$^{A}\alpha^{C}_{x}$","$^{A}\alpha^{C}_{y}$","$^{A}\alpha^{C}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,3)
plot(t(1:end-2),aad(1,:),t(1:end-2),aad(2,:),t(1:end-2),aad(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\alpha}^{D}$")
title("\textbf{Evolution of $^{A}\underline{\alpha}^{D}$}")
legend(["$^{A}\alpha^{D}_{x}$","$^{A}\alpha^{D}_{y}$","$^{A}\alpha^{D}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,4)
plot(t(1:end-2),aae(1,:),t(1:end-2),aae(2,:),t(1:end-2),aae(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\alpha}^{E}$")
title("\textbf{Evolution of $^{A}\underline{\alpha}^{E}$}")
legend(["$^{A}\alpha^{E}_{x}$","$^{A}\alpha^{E}_{y}$","$^{A}\alpha^{E}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,5)
plot(t(1:end-2),aaf(1,:),t(1:end-2),aaf(2,:),t(1:end-2),aaf(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\alpha}^{F}$")
title("\textbf{Evolution of $^{A}\underline{\alpha}^{F}$}")
legend(["$^{A}\alpha^{F}_{x}$","$^{A}\alpha^{F}_{y}$","$^{A}\alpha^{F}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
sgtitle("\textbf{Evolution of the Angular Accelerations}")

%% Evolution of the velocity and acceleration of point Dstar
figure
subplot(2,1,1)
plot(t,q1_d,t,q2_d,t,q3_d)
xlabel("$t\;(s)$")
ylabel("$\underline{v}^{D^{*}}$")
title("\textbf{Evolution of $\underline{v}^{D^{*}}$}")
legend(["$\dot{q}_{1}$","$\dot{q}_{2}$","$\dot{q}_{3}$"],'FontSize',12)
grid on
set(gcf,'color','w')
subplot(2,1,2)
plot(t,q1_dd,t,q2_dd,t,q3_dd)
xlabel("$t\;(s)$")
ylabel("$\underline{a}^{D^{*}}$")
title("\textbf{Evolution of $\underline{a}^{D^{*}}$}")
legend(["$\ddot{q}_{1}$","$\ddot{q}_{2}$","$\ddot{q}_{3}$"],'FontSize',12)
grid on
set(gcf,'color','w')

%% Evolution of the velocity and aceleration of point D1
avd1=cell2mat(AvD1);
aad1=cell2mat(AaD1);

figure
subplot(2,1,1)
plot(t(1:end-2),avd1(1,:),t(1:end-2),avd1(2,:),t(1:end-2),avd1(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{v}^{D1}$")
title("\textbf{Evolution of $^{A}\underline{v}^{D1}$}")
legend(["$^{A}v^{D1}_{x}$","$^{A}v^{D}_{y}$","$^{A}v^{D}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(2,1,2)
plot(t(1:end-2),aad1(1,:),t(1:end-2),aad1(2,:),t(1:end-2),aad1(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{a}^{D1}$")
title("\textbf{Evolution of $^{A}\underline{a}^{D1}$}")
legend(["$^{A}a^{D1}_{x}$","$^{A}a^{D1}_{y}$","$^{A}a^{D1}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')

%% Evolution of the velocity and aceleration of point D2
avd2=cell2mat(AvD2);
aad2=cell2mat(AaD2);

figure
subplot(2,1,1)
plot(t(1:end-2),avd2(1,:),t(1:end-2),avd2(2,:),t(1:end-2),avd2(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{v}^{D_{2}}$")
title("\textbf{Evolution of $^{A}\underline{v}^{D_{2}}$}")
legend(["$^{A}v^{D_{2}}_{x}$","$^{A}v^{D_{2}}_{y}$","$^{A}v^{D_{2}}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(2,1,2)
plot(t(1:end-2),aad2(1,:),t(1:end-2),aad2(2,:),t(1:end-2),aad2(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{a}^{D_{2}}$")
title("\textbf{Evolution of $^{A}\underline{a}^{D_{2}}$}")
legend(["$^{A}a^{D_{2}}_{x}$","$^{A}a^{D_{2}}_{y}$","$^{A}a^{D_{2}}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')

%% Evolution of the velocity and aceleration of point Estar
avestar=cell2mat(AvEstar);
aaestar=cell2mat(AaEstar);

figure
subplot(2,1,1)
plot(t(1:end-2),avestar(1,:),t(1:end-2),avestar(2,:),t(1:end-2),avestar(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{v}^{E^{*}}$")
title("\textbf{Evolution of $^{A}\underline{v}^{E^{*}}$}")
legend(["$^{A}v^{E^{*}}_{x}$","$^{A}v^{E^{*}}_{y}$","$^{A}v^{E^{*}}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(2,1,2)
plot(t(1:end-2),aaestar(1,:),t(1:end-2),aaestar(2,:),t(1:end-2),aaestar(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{a}^{E^{*}}$")
title("\textbf{Evolution of $^{A}\underline{a}^{E^{*}}$}")
legend(["$^{A}a^{E^{*}}_{x}$","$^{A}a^{E^{*}}_{y}$","$^{A}a^{E^{*}}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')

%% Evolution of the velocity and aceleration of point Fstar
avfstar=cell2mat(AvFstar);
aafstar=cell2mat(AaFstar);

figure
subplot(2,1,1)
plot(t(1:end-2),avfstar(1,:),t(1:end-2),avfstar(2,:),t(1:end-2),avfstar(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{v}^{F^{*}}$")
title("\textbf{Evolution of $^{A}\underline{v}^{F^{*}}$}")
legend(["$^{A}v^{F^{*}}_{x}$","$^{A}v^{F^{*}}_{y}$","$^{A}v^{F^{*}}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(2,1,2)
plot(t(1:end-2),aafstar(1,:),t(1:end-2),aafstar(2,:),t(1:end-2),aafstar(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{a}^{F^{*}}$")
title("\textbf{Evolution of $^{A}\underline{a}^{F^{*}}$}")
legend(["$^{A}a^{F^{*}}_{x}$","$^{A}a^{F^{*}}_{y}$","$^{A}a^{F^{*}}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')

%% Representation of the orientations of bodys
figure
subplot(3,1,1)
plot(t,q4,t,q5,t,q6)
xlabel("$t\;(s)$")
ylabel("$q_{i}$")
title("\textbf{Evolution of orientation of body D with ZXY Euler Angles}")
legend(["$\phi$","$\theta$","$\psi$"])
grid on
set(gcf,'color','w')
subplot(3,1,2)
plot(t,q4,t,q5,t,q6+q7)
xlabel("$t\;(s)$")
ylabel("$q_{i}$")
title("\textbf{Evolution of orientation of body E with ZXY Euler Angles}")
legend(["$\phi$","$\theta$","$\psi$"])
grid on
set(gcf,'color','w')
subplot(3,1,3)
plot(t,q4,t,q5,t,q6+q8)
xlabel("$t\;(s)$")
ylabel("$q_{i}$")
title("\textbf{Evolution of orientation of body F with ZXY Euler Angles}")
legend(["$\phi$","$\theta$","$\psi$"])
grid on
set(gcf,'color','w')    