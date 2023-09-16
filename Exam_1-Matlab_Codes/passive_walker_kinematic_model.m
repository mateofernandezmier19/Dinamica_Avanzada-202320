% Implementation of a Kinematic Model for a Passive Walker with
% hemispherical feet

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Definition of the generalized coordinates and their derivatives
dt=0.001;
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

for i=1:length(t)-1
    if q2(i)>=0       
        q1_d(i)=0;
        q2_d(i)=0;
        q3_d(i)=sin(2*pi*t(i));
        q4_d(i)=0;
        q5_d(i)=0;
        
        q1_dd(i)=0;
        q2_dd(i)=0;
        q3_dd(i)=2*pi*cos(2*pi*t(i));
        q4_dd(i)=0;
        q5_dd(i)=0;

        %Update of the generalized coordinates
        q1(i+1)=q1(i)+dt*q1_d(i);
        q2(i+1)=q2(i)+dt*q2_d(i);
        q3(i+1)=q3(i)+dt*q3_d(i);
        q4(i+1)=q1(i)+dt*q4_d(i);
        q5(i+1)=q5(i)+dt*q5_d(i);

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
        %Angular velocities in frame A
        AalphaC{i}=AalphaB{i}+aRb*BalphaC{i}+cross(AwB{i},BwC{i});
        AalphaD{i}=AalphaC{i}+aRc*CalphaD{i}+cross(AwC{i},CwD{i});
        AalphaE{i}=AalphaD{i}+aRd*DalphaE{i}+cross(AwD{i},DwE{i});
        AalphaF{i}=AalphaE{i}+aRe*EalphaF{i}+cross(AwE{i},EwF{i});      
    end
end



