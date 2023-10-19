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

%% Definition of the generalized coordinates and their derivatives
dt=1;
t=0:dt:10;
q1=zeros(1,length(t));
q2=zeros(1,length(t));
q3=zeros(1,length(t));
q4=zeros(1,length(t));
q5=zeros(1,length(t));

q1_d=zeros(1,length(t));
q2_d=zeros(1,length(t));
q3_d=zeros(1,length(t));
q4_d=zeros(1,length(t));
q5_d=zeros(1,length(t));

q1_dd=zeros(1,length(t));
q2_dd=zeros(1,length(t));
q3_dd=zeros(1,length(t));
q4_dd=zeros(1,length(t));
q5_dd=zeros(1,length(t));

%Position and Velocity components of the point P
qx=zeros(1,length(t));
qy=zeros(1,length(t));
qx_d=zeros(1,length(t));
qy_d=zeros(1,length(t));


%% Initial conditions
q1(1)=0;
q2(1)=0;
q3(1)=0;
q4(1)=0;
q5(1)=0;

q1_d(1)=0;
q2_d(1)=0;
q3_d(1)=0;
q4_d(1)=0;
q5_d(1)=0;

q1_dd(1)=0;
q2_dd(1)=0;
q3_dd(1)=0;
q4_dd(1)=0;
q5_dd(1)=0;

qx(1)=0;
qy(1)=0;
qx_d(1)=0;
qy_d(1)=0;


%% Algorithm for the kinematically driven calculations

%Angular velocities
AwB=cell(1,length(t));
BwC=cell(1,length(t));
CwD=cell(1,length(t));
DwE=cell(1,length(t));
EwF=cell(1,length(t));

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
EalphaF=cell(1,length(t));

%Angular velocities in the frame A
AalphaC=cell(1,length(t));
AalphaD=cell(1,length(t));
AalphaE=cell(1,length(t));
AalphaF=cell(1,length(t));

%Velocity and aceleration of point Dstar
AvDstar=cell(1,length(t));
AaDstar=cell(1,length(t));

%Velocity and aceleration of point D1
AvD1=cell(1,length(t));
AaD1=cell(1,length(t));

%Velocity and aceleration of point Estar
AvEstar=cell(1,length(t));
AaEstar=cell(1,length(t));

for i=1:length(t)-1
    if q2(i)>=0       
        q1_d(i)=2*pi*1;
        q2_d(i)=0;
        q3_d(i)=0;
        q4_d(i)=0;
        q5_d(i)=0;
        
        q1_dd(i)=0;
        q2_dd(i)=0;
        q3_dd(i)=0;
        q4_dd(i)=0;
        q5_dd(i)=0;

        %Update of the generalized coordinates
        q1(i+1)=q1(i)+dt*q1_d(i);
        q2(i+1)=q2(i)+dt*q2_d(i);
        q3(i+1)=q3(i)+dt*q3_d(i);
        q4(i+1)=q4(i)+dt*q4_d(i);
        q5(i+1)=q5(i)+dt*q5_d(i);

        %Udpate of the position and velocity of point P
        qx_d(i)=q2_d(i)*R*sin(q1(i))+q3_d(i)*R*cos(q2(i))*cos(q1(i));
        qy_d(i)=-q2_d(i)*R*cos(q1(i))+q3_d(i)*R*cos(q2(i))*sin(q1(i));
        qx(i+1)=qx(i)+dt*qx_d(i);
        qy(i+1)=qy(i)+dt*qy_d(i);

        %Rotation matrixes
        aRb=[cos(q1(i)),-sin(q1(i)),0;sin(q1(i)),cos(q1(i)),0;0,0,1];
        bRc=[1,0,0;0,cos(q2(i)),-sin(q2(i));0,sin(q2(i)),cos(q2(i))];
        cRd=[cos(q3(i)),0,sin(q3(i));0,1,0;-sin(q3(i)),0,cos(q3(i))];
        dRe=[cos(q4(i)),0,sin(q4(i));0,1,0;-sin(q4(i)),0,cos(q4(i))];
        eRf=[cos(q5(i)),0,sin(q5(i));0,1,0;-sin(q5(i)),0,cos(q5(i))];
        %Rotation matrixes in terms of the frame A
        aRc=aRb*bRc;
        aRd=aRc*cRd;
        aRe=aRd*dRe;
        aRf=aRe*eRf;

        %Angular velocities
        AwB{i}=q1_d(i)*[0;0;1];
        BwC{i}=q2_d(i)*[1;0;0];
        CwD{i}=q3_d(i)*[0;1;0];
        DwE{i}=q4_d(i)*[0;1;0];
        EwF{i}=q5_d(i)*[0;1;0];
        %Angular velocities in frame A
        AwC{i}=AwB{i}+aRb*BwC{i};
        AwD{i}=AwC{i}+aRc*CwD{i};
        AwE{i}=AwD{i}+aRd*DwE{i};
        AwF{i}=AwE{i}+aRe*EwF{i};

        %Angular Acelerations
        AalphaB{i}=q1_dd(i)*[0;0;1];
        BalphaC{i}=q2_dd(i)*[1;0;0];
        CalphaD{i}=q3_dd(i)*[0;1;0];
        DalphaE{i}=q4_dd(i)*[0;1;0];
        EalphaF{i}=q5_dd(i)*[0;1;0];
        %Angular Acelerations in frame A
        AalphaC{i}=AalphaB{i}+aRb*BalphaC{i}+cross(AwB{i},BwC{i});
        AalphaD{i}=AalphaC{i}+aRc*CalphaD{i}+cross(AwC{i},CwD{i});
        AalphaE{i}=AalphaD{i}+aRd*DalphaE{i}+cross(AwD{i},DwE{i});
        AalphaF{i}=AalphaE{i}+aRe*EalphaF{i}+cross(AwE{i},EwF{i});

        %Velocity of point Dstar
        AvDstar{i}=[qx_d(i);qx_d(i);0]+cross(AwB{i},aRb*(-R)*[0;0;1]);        

        %Velocity and aceleration of point D1
        AvD1{i}=cross(AwD{i},aRd*Lp*[0;0;1]);
        AaD1{i}=cross(AwD{i},cross(AwD{i},aRd*Lp*[0;0;1]))+cross(AalphaD{i},aRd*Lp*[0;0;1]);

        %Velocity and aceleration of point Estar
        AvEstar{i}=AvD1{i}+cross(AwE{i},aRe*d3*[0;1;0]);
        AaEstar{i}=AaD1{i}+cross(AwE{i},cross(AwE{i},aRe*d3*[0;1;0]))+cross(AalphaE{i},aRe*d3*[0;1;0]);


    end
end


%% Figures

%Evolution of the generalized coordinates
figure
plot(t,q1,t,q2,t,q3,t,q4,t,q5)
xlabel("$t\;(s)$")
ylabel("$q_{i}$")
title("\textbf{Evolution of $q_{i}$}")
legend(["$q1$","$q2$","$q3$","$q4$","$q5$"])
grid on
set(gcf,'color','w')

%% Evolution of the time derivative of generalized coordinates
figure
plot(t,q1_d,t,q2_d,t,q3_d,t,q4_d,t,q5_d)
xlabel("$t\;(s)$")
ylabel("$\dot{q}_{i}$")
title("\textbf{Evolution of $\dot{q}_{i}$}")
legend(["$\dot{q}_{1}$","$\dot{q}_{2}$","$\dot{q}_{3}$","$\dot{q}_{4}$","$\dot{q}_{5}$"])
grid on
set(gcf,'color','w')

%% Evolution of the second time derivative of generalized coordinates
figure
plot(t,q1_dd,t,q2_dd,t,q3_dd,t,q4_dd,t,q5_dd)
xlabel("$t\;(s)$")
ylabel("$\ddot{q}_{i}$")
title("\textbf{Evolution of $\ddot{q}_{i}$}")
legend(["$\ddot{q}_{1}$","$\ddot{q}_{2}$","$\ddot{q}_{3}$","$\ddot{q}_{4}$","$\ddot{q}_{5}$"])
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
plot(t(1:end-1),ab(1,:),t(1:end-1),ab(2,:),t(1:end-1),ab(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\omega}^{B}$")
title("\textbf{Evolution of $^{A}\underline{\omega}^{B}$}")
legend(["$^{A}\omega^{B}_{x}$","$^{A}\omega^{B}_{y}$","$^{A}\omega^{B}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,2)
plot(t(1:end-1),ac(1,:),t(1:end-1),ac(2,:),t(1:end-1),ac(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\omega}^{C}$")
title("\textbf{Evolution of $^{A}\underline{\omega}^{C}$}")
legend(["$^{A}\omega^{C}_{x}$","$^{A}\omega^{C}_{y}$","$^{A}\omega^{C}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,3)
plot(t(1:end-1),ad(1,:),t(1:end-1),ad(2,:),t(1:end-1),ad(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\omega}^{D}$")
title("\textbf{Evolution of $^{A}\underline{\omega}^{D}$}")
legend(["$^{A}\omega^{D}_{x}$","$^{A}\omega^{D}_{y}$","$^{A}\omega^{D}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,4)
plot(t(1:end-1),ae(1,:),t(1:end-1),ae(2,:),t(1:end-1),ae(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\omega}^{E}$")
title("\textbf{Evolution of $^{A}\underline{\omega}^{E}$}")
legend(["$^{A}\omega^{E}_{x}$","$^{A}\omega^{E}_{y}$","$^{A}\omega^{E}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,5)
plot(t(1:end-1),af(1,:),t(1:end-1),af(2,:),t(1:end-1),af(3,:))
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
plot(t(1:end-1),aab(1,:),t(1:end-1),aab(2,:),t(1:end-1),aab(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\alpha}^{B}$")
title("\textbf{Evolution of $^{A}\underline{\alpha}^{B}$}")
legend(["$^{A}\alpha^{B}_{x}$","$^{A}\alpha^{B}_{y}$","$^{A}\alpha^{B}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,2)
plot(t(1:end-1),aac(1,:),t(1:end-1),aac(2,:),t(1:end-1),aac(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\alpha}^{C}$")
title("\textbf{Evolution of $^{A}\underline{\alpha}^{C}$}")
legend(["$^{A}\alpha^{C}_{x}$","$^{A}\alpha^{C}_{y}$","$^{A}\alpha^{C}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,3)
plot(t(1:end-1),aad(1,:),t(1:end-1),aad(2,:),t(1:end-1),aad(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\alpha}^{D}$")
title("\textbf{Evolution of $^{A}\underline{\alpha}^{D}$}")
legend(["$^{A}\alpha^{D}_{x}$","$^{A}\alpha^{D}_{y}$","$^{A}\alpha^{D}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,4)
plot(t(1:end-1),aae(1,:),t(1:end-1),aae(2,:),t(1:end-1),aae(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\alpha}^{E}$")
title("\textbf{Evolution of $^{A}\underline{\alpha}^{E}$}")
legend(["$^{A}\alpha^{E}_{x}$","$^{A}\alpha^{E}_{y}$","$^{A}\alpha^{E}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(3,2,5)
plot(t(1:end-1),aaf(1,:),t(1:end-1),aaf(2,:),t(1:end-1),aaf(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{\alpha}^{F}$")
title("\textbf{Evolution of $^{A}\underline{\alpha}^{F}$}")
legend(["$^{A}\alpha^{F}_{x}$","$^{A}\alpha^{F}_{y}$","$^{A}\alpha^{F}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
sgtitle("\textbf{Evolution of the Angular Accelerations}")

%% Evolution of the position and velocity of point P
figure
subplot(2,1,1)
plot(t,qx,t,qy)
xlabel("$t\;(s)$")
ylabel("$\underline{r}^{P}$")
title("\textbf{Evolution of $\underline{r}^{P}$}")
legend(["$q_{x}$","$q_{y}$"],'FontSize',12)
grid on
set(gcf,'color','w')
subplot(2,1,2)
plot(t,qx_d,t,qy_d)
xlabel("$t\;(s)$")
ylabel("$\underline{v}^{P}$")
title("\textbf{Evolution of $\underline{v}^{P}$}")
legend(["$\dot{q}_{x}$","$\dot{q}_{y}$"],'FontSize',12)
grid on
set(gcf,'color','w')

%% Evolution of the velocity and aceleration of point D1
avd1=cell2mat(AvD1);
aad1=cell2mat(AaD1);

figure
subplot(2,1,1)
plot(t(1:end-1),avd1(1,:),t(1:end-1),avd1(2,:),t(1:end-1),avd1(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{v}^{D1}$")
title("\textbf{Evolution of $^{A}\underline{v}^{D1}$}")
legend(["$^{A}v^{D1}_{x}$","$^{A}v^{D}_{y}$","$^{A}v^{D}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(2,1,2)
plot(t(1:end-1),aad1(1,:),t(1:end-1),aad1(2,:),t(1:end-1),aad1(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{a}^{D1}$")
title("\textbf{Evolution of $^{A}\underline{a}^{D1}$}")
legend(["$^{A}a^{D1}_{x}$","$^{A}a^{D1}_{y}$","$^{A}a^{D1}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')

%% Evolution of the velocity and aceleration of point Estar
avestar=cell2mat(AvEstar);
aaestar=cell2mat(AaEstar);

figure
subplot(2,1,1)
plot(t(1:end-1),avestar(1,:),t(1:end-1),avestar(2,:),t(1:end-1),avestar(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{v}^{E^{*}}$")
title("\textbf{Evolution of $^{A}\underline{v}^{E^{*}}$}")
legend(["$^{A}v^{E^{*}}_{x}$","$^{A}v^{E^{*}}_{y}$","$^{A}v^{E^{*}}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
subplot(2,1,2)
plot(t(1:end-1),aaestar(1,:),t(1:end-1),aaestar(2,:),t(1:end-1),aaestar(3,:))
xlabel("$t\;(s)$")
ylabel("$^{A}\underline{a}^{E^{*}}$")
title("\textbf{Evolution of $^{A}\underline{a}^{E^{*}}$}")
legend(["$^{A}a^{E^{*}}_{x}$","$^{A}a^{E^{*}}_{y}$","$^{A}a^{E^{*}}_{z}$"],'FontSize',8)
grid on
set(gcf,'color','w')
