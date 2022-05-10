%-------------------------------------------------%
% Written by - Prakrit Tyagi                      %
% Dated - 1-01-2022                               %
%-------------------------------------------------%
% Tracking(gimbal,obstacles) Multishooting, Target Prediction, LOS cost
% function 
close all
clear
clc

addpath('C:\Users\tyagi\Documents\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*
%% MPC Parameters
N = 50;%50
T = 0.1;
% Constants
m = 1.56;
S = 0.2589;
b_wingspan = 1.4224;
ar = b_wingspan^2/S;
cdo = 0.01631;
kdl = 0.04525;
rho = 1.225;
a_g = 9.8;
% Camera Parameters //to be changed in Graph_MPC_casadi//
VFOV_deg = 54.4; % vertical field of view in degree 94.4
HFOV_deg = 62.6; % horizontal field of view in degree 122.6

VFOV_rad = VFOV_deg*pi/180;
HFOV_rad = HFOV_deg*pi/180;
% Gimbal servo control limits
pitch_maxG = pi/30; roll_maxG = pi/30; yaw_maxG = pi/30;

% Obstacle parameters
obs(1) = struct('x',120,'y',300,'h',120,'r',50);%
obs(2) = struct('x',-70,'y',350,'h',120,'r',40);%
obs(3) = struct('x',150,'y',640,'h',100,'r',50);
obs(4) = struct('x',40,'y',520,'h',90,'r',30);%
obs(5) = struct('x',320,'y',730,'h',100,'r',60);%710 40
% obs(6) = struct('x',500,'y',700,'h',80,'r',30);%
obs(6) = struct('x',750,'y',630,'h',70,'r',30);%
obs(7) = struct('x',300,'y',570,'h',110,'r',30);
obs(8) = struct('x',900,'y',550,'h',75,'r',40);
obs(9) = struct('x',750,'y',470,'h',120,'r',30);%
obs(10) = struct('x',950,'y',300,'h',100,'r',30);%
obs(11) = struct('x',-60,'y',50,'h',90,'r',30);
obs(12) = struct('x',490,'y',710,'h',80,'r',40);
obs(13) = struct('x',60,'y',100,'h',70,'r',30);

n_obs = length(obs);
r = 30;

%% Defining model for UAV
x = SX.sym('x'); y = SX.sym('y'); z = SX.sym('z'); v = SX.sym('v'); gamma = SX.sym('gamma'); psi = SX.sym('psi');
states = [x;y;z;v;gamma;psi]; n_states = length(states);

alpha = SX.sym('alpha'); phi = SX.sym('phi'); Thrust = SX.sym('Thrust');
controls = [alpha;phi;Thrust]; n_controls = length(controls);

model = [v*cos(gamma)*cos(psi); v*cos(gamma)*sin(psi); v*sin(gamma); (Thrust - 0.5*rho*S*(cdo + kdl*alpha^2)*v^2 - m*a_g*sin(gamma))/m;...
    ((0.5*rho*S*alpha*v^2)*cos(phi) + Thrust*sin(alpha) - m*a_g*cos(gamma))/(m*v); (0.5*rho*S*alpha*v^2)*sin(phi)/(cos(gamma)*m*v)];

f = Function('f',{states,controls},{model});

P = SX.sym('P',n_states);
X = SX.sym('X',n_states,(N+1));
U = SX.sym('U',n_controls,N);

%% Defining model for gimbal
psiG = SX.sym('psiG'); phiG = SX.sym('phiG'); thetaG = SX.sym('thetaG');
statesG = [psiG;phiG;thetaG]; n_statesG = length(statesG);

pitchG = SX.sym('pitchG'); rollG = SX.sym('rollG'); yawG = SX.sym('yawG');
controlsG = [pitchG;rollG;yawG]; n_controlsG = length(controlsG);

modelG = [pitchG;rollG;yawG];
fG = Function('fG',{statesG,controlsG},{modelG});

PG = SX.sym('PG',n_statesG);
XG = SX.sym('XG',n_statesG,(N+1));
UG = SX.sym('UG',n_controlsG,N); 

%% Defining model for Target
xT = SX.sym('xT'); yT = SX.sym('yT'); ang = SX.sym('ang');
statesT = [xT;yT;ang]; n_statesT = length(statesT);

vT = SX.sym('vT'); angRate = SX.sym('angRate');
controlsT = [vT;angRate]; n_controlsT = length(controlsT);

modelT = [vT*cos(ang);vT*sin(ang);angRate];
fT = Function('fT',{statesT,controlsT},{modelT});

PT = SX.sym('PT',n_statesT*(N+1));
XT = SX.sym('XT',n_statesT,2);
UT = SX.sym('UT',n_controlsT);

% Calcuate next state
XT(:,1) = PT(1:n_statesT);
c_XT = XT(:,1); c_UT = UT(:,1);

fT_value = fT(c_XT,c_UT);
new_XT = c_XT + T*fT_value;
XT(:,2) = new_XT;

ffT = Function('ffT',{UT,PT(1:n_statesT)},{XT});


%% Defining the objective cost function and g constraint matrix

obj = 0; % Objective function
g = [];  % constraints vector

w_1 = [0.3;0.3];%0.2 Distance
Q = diag(w_1);
w_2 = 98;%0.8 ,6 LOS
w_3 = 0.2;%0.2 energy
w_4 = 0.9; %LOSE
w_5 = 10; % avoid obstacle increase height

p = [P;PG;PT]; % combined parameter matrix 
st  = [X(:,1);XG(:,1)]; % initial state
g = [g;st-p(1:(n_states+n_statesG))]; % initial condition constraints
for k = 1:N
    % relative distance cost function
   obj = obj + sqrt((X(1:2,k)-PT(3*k-2:3*k-1))'*Q*(X(1:2,k)-PT(3*k-2:3*k-1)));

   a = 0.5*( X(3,k)*tan(XG(1,k) + X(5,k) + VFOV_rad/2) - X(3,k)*tan(XG(1,k) + X(5,k) - VFOV_rad/2) );
   b = 0.5*( X(3,k)*tan(XG(2,k) + HFOV_rad/2) - X(3,k)*tan(XG(2,k) - HFOV_rad/2) );
   shifta = a + X(3,k)*tan(XG(1,k) + X(5,k) - VFOV_rad/2);
   shiftb = b + X(3,k)*tan(XG(2,k) - HFOV_rad/2);
   A = ((cos(X(6,k) + XG(3,k)))^2)/(a^2) + ((sin(X(6,k) + XG(3,k)))^2)/(b^2) ;
   B = 2*cos(X(6,k) + XG(3,k))*sin(X(6,k) + XG(3,k))*( (1/a^2) - (1/b^2) );
   C = (sin(X(6,k) + XG(3,k)))^2/(a^2) + (cos(X(6,k) + XG(3,k)))^2/(b^2);
   ellipse = A*(PT(3*k-2) - X(1,k) - shifta)^2 + B*(PT(3*k-2) - X(1,k) - shifta)*(PT(3*k-1) - X(2,k) - shiftb) + C*(PT(3*k-1) - X(2,k) - shiftb)^2 - 1; 
   obj = obj + w_2*ellipse;
   
%     energy cost function
    c_vel = X(4,k); c_alpha = U(1,k); c_pitch = X(5,k); c_thrust = U(3,k);
% 
%     Co_drag = cdo + kdl*c_alpha^2;
%     obj = obj + w_3*(m*c_vel*(c_thrust - 0.5*rho*S*Co_drag*c_vel^2 - m*a_g*sin(c_pitch))/m + 0.5*Co_drag*rho*S*c_vel^3 + ...
%           2*(m*a_g)^2*cos(c_pitch)^2/(pi*rho*S*ar*c_vel) + m*a_g*c_vel*sin(c_pitch));
%     obj = obj + c_thrust*c_vel;

   logic_TargetInVision = floor(ellipse/10000000)*-1; % 1 if inside ellipse else 0 for outside
%    
%    % Visibility cost function obstacle 1
   slope = (X(2,k)-PT(3*k-1))/(X(1,k)-PT(3*k-2)); % slope of ray(x-yplane) from target to UAV 
   center1_to_LOS = (abs(slope*obs(1).x - 1*obs(1).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(1).r; % distance between line LOS 2D and center of obstacle1
   slopeRef1 = (obs(1).h/(sqrt( (PT(3*k-2)-obs(1).x)^2 + (PT(3*k-1)-obs(1).y)^2 )-obs(1).r)) - (X(3,k)/(sqrt( (PT(3*k-2)-X(1,k))^2 + (PT(3*k-1)-X(2,k))^2 ))); % Slope of LOS between obstacle_1 Top and target in 3D
   logic_LosIntersectObstacle1 = floor(center1_to_LOS/100000000)*-1;
   costlos1 = logic_TargetInVision*logic_LosIntersectObstacle1*max(0,slopeRef1);
%    costlos1 = logic_LosIntersectObstacle1*max(0,slopeRef1);
%    obj = obj + w_4*center1_to_LOS;
   obj = obj + w_4*costlos1;
   
   % Visibility cost function obstacle 2
%    slope = (X(2,k)-PT(2))/(X(1,k)-PT(1)); % slope of ray(x-yplane) from target to UAV 
   center2_to_LOS = (abs(slope*obs(2).x - 1*obs(2).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(2).r; % distance between line LOS 2D and center of obstacle2
   slopeRef2 = (obs(2).h/(sqrt( (PT(3*k-2)-obs(2).x)^2 + (PT(3*k-1)-obs(2).y)^2 )-obs(2).r)) - (X(3,k)/(sqrt( (PT(3*k-2)-X(1,k))^2 + (PT(3*k-1)-X(2,k))^2 ))); % Slope of LOS between obstacle_2 Top and target in 3D
   logic_LosIntersectObstacle2 = floor(center2_to_LOS/100000000)*-1;
   costlos2 = logic_TargetInVision*logic_LosIntersectObstacle2*max(0,slopeRef2);
%    costlos2 = logic_LosIntersectObstacle2*max(0,slopeRef2);
%    obj = obj + w_4*center2_to_LOS;
   obj = obj + w_4*costlos2;
   
   % Visibility cost function obstacle 3
%    slope = (X(2,k)-PT(2))/(X(1,k)-PT(1)); % slope of ray(x-yplane) from target to UAV 
   center3_to_LOS = (abs(slope*obs(3).x - 1*obs(3).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(3).r; % distance between line LOS 2D and center of obstacle3
   slopeRef3 = (obs(3).h/(sqrt( (PT(3*k-2)-obs(3).x)^2 + (PT(3*k-1)-obs(3).y)^2 )-obs(3).r)) - (X(3,k)/(sqrt( (PT(3*k-2)-X(1,k))^2 + (PT(3*k-1)-X(2,k))^2 ))); % Slope of LOS between obstacle_3 Top and target in 3D
   logic_LosIntersectObstacle3 = floor(center3_to_LOS/100000000)*-1;
   costlos3 = logic_TargetInVision*logic_LosIntersectObstacle3*max(0,slopeRef3);
%    costlos3 = logic_LosIntersectObstacle3*max(0,slopeRef3);
%     obj = obj + w_4*center3_to_LOS;
   obj = obj + w_4*costlos3;
   
   % Visibility cost function obstacle 4
%    slope = (X(2,k)-PT(2))/(X(1,k)-PT(1)); % slope of ray(x-yplane) from target to UAV 
   center4_to_LOS = (abs(slope*obs(4).x - 1*obs(4).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(4).r; % distance between line LOS 2D and center of obstacle4
   slopeRef4 = (obs(4).h/(sqrt( (PT(3*k-2)-obs(4).x)^2 + (PT(3*k-1)-obs(4).y)^2 )-obs(4).r)) - (X(3,k)/(sqrt( (PT(3*k-2)-X(1,k))^2 + (PT(3*k-1)-X(2,k))^2 ))); % Slope of LOS between obstacle_4 Top and target in 3D
   logic_LosIntersectObstacle4 = floor(center4_to_LOS/100000000)*-1;
   costlos4 = logic_TargetInVision*logic_LosIntersectObstacle4*max(0,slopeRef4);
%    costlos4 = logic_LosIntersectObstacle4*max(0,slopeRef4);
%     obj = obj + w_4*center4_to_LOS;
   obj = obj + w_4*costlos4;
   
   % Visibility cost function obstacle 5
%    slope = (X(2,k)-PT(2))/(X(1,k)-PT(1)); % slope of ray(x-yplane) from target to UAV 
   center5_to_LOS = (abs(slope*obs(5).x - 1*obs(5).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(5).r; % distance between line LOS 2D and center of obstacle5
   slopeRef5 = (obs(5).h/(sqrt( (PT(3*k-2)-obs(5).x)^2 + (PT(3*k-1)-obs(5).y)^2 )-obs(5).r)) - (X(3,k)/(sqrt( (PT(3*k-2)-X(1,k))^2 + (PT(3*k-1)-X(2,k))^2 ))); % Slope of LOS between obstacle_5 Top and target in 3D
   logic_LosIntersectObstacle5 = floor(center5_to_LOS/100000000)*-1;
   costlos5 = logic_TargetInVision*logic_LosIntersectObstacle5*max(0,slopeRef5);
%    costlos5 = logic_LosIntersectObstacle5*max(0,slopeRef5);
%     obj = obj + w_4*center5_to_LOS;
   obj = obj + w_4*costlos5;
   
   % Visibility cost function obstacle 6
%    slope = (X(2,k)-PT(2))/(X(1,k)-PT(1)); % slope of ray(x-yplane) from target to UAV 
   center6_to_LOS = (abs(slope*obs(6).x - 1*obs(6).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(6).r; % distance between line LOS 2D and center of obstacle6
   slopeRef6 = (obs(6).h/(sqrt( (PT(3*k-2)-obs(6).x)^2 + (PT(3*k-1)-obs(6).y)^2 )-obs(6).r)) - (X(3,k)/(sqrt( (PT(3*k-2)-X(1,k))^2 + (PT(3*k-1)-X(2,k))^2 ))); % Slope of LOS between obstacle_6 Top and target in 3D
   logic_LosIntersectObstacle6 = floor(center6_to_LOS/100000000)*-1;
   costlos6 = logic_TargetInVision*logic_LosIntersectObstacle6*max(0,slopeRef6);
%    costlos6 = logic_LosIntersectObstacle6*max(0,slopeRef6);
%     obj = obj + w_4*center6_to_LOS;
   obj = obj + w_4*costlos6;
   
   % Visibility cost function obstacle 7
%    slope = (X(2,k)-PT(2))/(X(1,k)-PT(1)); % slope of ray(x-yplane) from target to UAV 
   center7_to_LOS = (abs(slope*obs(7).x - 1*obs(7).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(7).r; % distance between line LOS 2D and center of obstacle7
   slopeRef7 = (obs(7).h/(sqrt( (PT(3*k-2)-obs(7).x)^2 + (PT(3*k-1)-obs(7).y)^2 )-obs(7).r)) - (X(3,k)/(sqrt( (PT(3*k-2)-X(1,k))^2 + (PT(3*k-1)-X(2,k))^2 ))); % Slope of LOS between obstacle_7 Top and target in 3D
   logic_LosIntersectObstacle7 = floor(center7_to_LOS/100000000)*-1;
   costlos7 = logic_TargetInVision*logic_LosIntersectObstacle7*max(0,slopeRef7);
%    costlos7 = logic_LosIntersectObstacle7*max(0,slopeRef7);
%     obj = obj + w_4*center7_to_LOS;
   obj = obj + w_4*costlos7;
   
   % Visibility cost function obstacle 8
%    slope = (X(2,k)-PT(2))/(X(1,k)-PT(1)); % slope of ray(x-yplane) from target to UAV 
   center8_to_LOS = (abs(slope*obs(8).x - 1*obs(8).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(8).r; % distance between line LOS 2D and center of obstacle8
   slopeRef8 = (obs(8).h/(sqrt( (PT(3*k-2)-obs(8).x)^2 + (PT(3*k-1)-obs(8).y)^2 )-obs(8).r)) - (X(3,k)/(sqrt( (PT(3*k-2)-X(1,k))^2 + (PT(3*k-1)-X(2,k))^2 ))); % Slope of LOS between obstacle_8 Top and target in 3D
   logic_LosIntersectObstacle8 = floor(center8_to_LOS/100000000)*-1;
   costlos8 = logic_TargetInVision*logic_LosIntersectObstacle8*max(0,slopeRef8);
%    costlos8 = logic_LosIntersectObstacle8*max(0,slopeRef8);
%     obj = obj + w_4*center8_to_LOS;
   obj = obj + w_4*costlos8;
   
   % Visibility cost function obstacle 9
%    slope = (X(2,k)-PT(2))/(X(1,k)-PT(1)); % slope of ray(x-yplane) from target to UAV 
   center9_to_LOS = (abs(slope*obs(9).x - 1*obs(9).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(9).r; % distance between line LOS 2D and center of obstacle9
   slopeRef9 = (obs(9).h/(sqrt( (PT(3*k-2)-obs(9).x)^2 + (PT(3*k-1)-obs(9).y)^2 )-obs(9).r)) - (X(3,k)/(sqrt( (PT(3*k-2)-X(1,k))^2 + (PT(3*k-1)-X(2,k))^2 ))); % Slope of LOS between obstacle_9 Top and target in 3D
   logic_LosIntersectObstacle9 = floor(center9_to_LOS/100000000)*-1;
   costlos9 = logic_TargetInVision*logic_LosIntersectObstacle9*max(0,slopeRef9);
%    costlos9 = logic_LosIntersectObstacle9*max(0,slopeRef9);
%     obj = obj + w_4*center9_to_LOS;
   obj = obj + w_4*costlos9;
   
   % Visibility cost function obstacle 10
%    slope = (X(2,k)-PT(2))/(X(1,k)-PT(1)); % slope of ray(x-yplane) from target to UAV 
   center10_to_LOS = (abs(slope*obs(10).x - 1*obs(10).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(10).r; % distance between line LOS 2D and center of obstacle10
   slopeRef10 = (obs(10).h/(sqrt( (PT(3*k-2)-obs(10).x)^2 + (PT(3*k-1)-obs(10).y)^2 )-obs(10).r)) - (X(3,k)/(sqrt( (PT(3*k-2)-X(1,k))^2 + (PT(3*k-1)-X(2,k))^2 ))); % Slope of LOS between obstacle_10 Top and target in 3D
   logic_LosIntersectObstacle10 = floor(center10_to_LOS/100000000)*-1;
   costlos10 = logic_TargetInVision*logic_LosIntersectObstacle10*max(0,slopeRef10);
%    costlos10 = logic_LosIntersectObstacle10*max(0,slopeRef10);
%     obj = obj + w_4*center10_to_LOS;
   obj = obj + w_4*costlos10;

   % Visibility cost function obstacle 11
%    slope = (X(2,k)-PT(2))/(X(1,k)-PT(1)); % slope of ray(x-yplane) from target to UAV 
   center11_to_LOS = (abs(slope*obs(11).x - 1*obs(11).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(11).r; % distance between line LOS 2D and center of obstacle11
   slopeRef11 = (obs(11).h/(sqrt( (PT(3*k-2)-obs(11).x)^2 + (PT(3*k-1)-obs(11).y)^2 )-obs(11).r)) - (X(3,k)/(sqrt( (PT(3*k-2)-X(1,k))^2 + (PT(3*k-1)-X(2,k))^2 ))); % Slope of LOS between obstacle_11 Top and target in 3D
   logic_LosIntersectObstacle11 = floor(center11_to_LOS/100000000)*-1;
   costlos11 = logic_TargetInVision*logic_LosIntersectObstacle11*max(0,slopeRef11);
%    costlos11 = logic_LosIntersectObstacle11*max(0,slopeRef11);
%     obj = obj + w_4*center11_to_LOS;
   obj = obj + w_4*costlos11;
   
   % Visibility cost function obstacle 12
%    slope = (X(2,k)-PT(2))/(X(1,k)-PT(1)); % slope of ray(x-yplane) from target to UAV 
   center12_to_LOS = (abs(slope*obs(12).x - 1*obs(12).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(12).r; % distance between line LOS 2D and center of obstacle12
   slopeRef12 = (obs(12).h/(sqrt( (PT(3*k-2)-obs(12).x)^2 + (PT(3*k-1)-obs(12).y)^2 )-obs(12).r)) - (X(3,k)/(sqrt( (PT(3*k-2)-X(1,k))^2 + (PT(3*k-1)-X(2,k))^2 ))); % Slope of LOS between obstacle_12 Top and target in 3D
   logic_LosIntersectObstacle12 = floor(center12_to_LOS/100000000)*-1;
   costlos12 = logic_TargetInVision*logic_LosIntersectObstacle12*max(0,slopeRef12);
%    costlos12 = logic_LosIntersectObstacle12*max(0,slopeRef12); 
%     obj = obj + w_4*center12_to_LOS;
   obj = obj + w_4*costlos12;
   
   % Visibility cost function obstacle 13
%    slope = (X(2,k)-PT(2))/(X(1,k)-PT(1)); % slope of ray(x-yplane) from target to UAV 
   center13_to_LOS = (abs(slope*obs(13).x - 1*obs(13).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(13).r; % distance between line LOS 2D and center of obstacle13
   slopeRef13 = (obs(13).h/(sqrt( (PT(3*k-2)-obs(13).x)^2 + (PT(3*k-1)-obs(13).y)^2 )-obs(13).r)) - (X(3,k)/(sqrt( (PT(3*k-2)-X(1,k))^2 + (PT(3*k-1)-X(2,k))^2 ))); % Slope of LOS between obstacle_13 Top and target in 3D
   logic_LosIntersectObstacle13 = floor(center13_to_LOS/100000000)*-1;
   costlos13 = logic_TargetInVision*logic_LosIntersectObstacle13*max(0,slopeRef13);
%    costlos13 = logic_LosIntersectObstacle13*max(0,slopeRef13);
%     obj = obj + w_4*center13_to_LOS;
   obj = obj + w_4*costlos13;


   %obstacle collision avoidance
   obs1 = ((X(1,k)-obs(1).x)^2 + (X(2,k)-obs(1).y)^2 - r^2);
   checkobs1 = -1*floor(obs1/100000000);
%    costobs1 = checkobs1*(obs(1).h-X(3,k))*(obs(1).h-X(3,k))^2;
   costobs1 = checkobs1*(obs(1).h-X(3,k));   
   obj = obj + w_5*costobs1;
   
   obs2 = ((X(1,k)-obs(2).x)^2 + (X(2,k)-obs(2).y)^2 - r^2);
   checkobs2 = -1*floor(obs2/100000000);
%    costobs2 = checkobs2*(obs(2).h-X(3,k))*(obs(2).h-X(3,k))^2;
   costobs2 = checkobs2*(obs(2).h-X(3,k));
   obj = obj + w_5*costobs2;
   
   obs3 = ((X(1,k)-obs(3).x)^2 + (X(2,k)-obs(3).y)^2 - r^2);
   checkobs3 = -1*floor(obs3/100000000);
%    costobs3 = checkobs3*(obs(3).h-X(3,k))*(obs(3).h-X(3,k))^2;
   costobs3 = checkobs3*(obs(3).h-X(3,k));
   obj = obj + w_5*costobs3;
   
   obs4 = ((X(1,k)-obs(4).x)^2 + (X(2,k)-obs(4).y)^2 - r^2);
   checkobs4 = -1*floor(obs4/100000000);
%    costobs4 = checkobs4*(obs(4).h-X(3,k))*(obs(4).h-X(3,k))^2;
   costobs4 = checkobs4*(obs(4).h-X(3,k));
   obj = obj + w_5*costobs4;
   
   obs5 = ((X(1,k)-obs(5).x)^2 + (X(2,k)-obs(5).y)^2 - r^2);
   checkobs5 = -1*floor(obs5/100000000);
%    costobs5 = checkobs5*(obs(5).h-X(3,k))*(obs(5).h-X(3,k))^2;
   costobs5 = checkobs5*(obs(5).h-X(3,k));
   obj = obj + w_5*costobs5;

   obs6 = ((X(1,k)-obs(6).x)^2 + (X(2,k)-obs(6).y)^2 - r^2);
   checkobs6 = -1*floor(obs6/100000000);
%    costobs6 = checkobs6*(obs(6).h-X(3,k))*(obs(6).h-X(3,k))^2;
   costobs6 = checkobs6*(obs(6).h-X(3,k));
   obj = obj + w_5*costobs6;
   
   obs7 = ((X(1,k)-obs(7).x)^2 + (X(2,k)-obs(7).y)^2 - r^2);
   checkobs7 = -1*floor(obs7/100000000);
%    costobs7 = checkobs7*(obs(7).h-X(3,k))*(obs(7).h-X(3,k))^2;
   costobs7 = checkobs7*(obs(7).h-X(3,k));
   obj = obj + w_5*costobs7;

   obs8 = ((X(1,k)-obs(8).x)^2 + (X(2,k)-obs(8).y)^2 - r^2);
   checkobs8 = -1*floor(obs8/100000000);
%    costobs8 = checkobs8*(obs(8).h-X(3,k))*(obs(8).h-X(3,k))^2;
   costobs8 = checkobs8*(obs(8).h-X(3,k));
   obj = obj + w_5*costobs8;
   
   obs9 = ((X(1,k)-obs(9).x)^2 + (X(2,k)-obs(9).y)^2 - r^2);
   checkobs9 = -1*floor(obs9/100000000);
%    costobs9 = checkobs9*(obs(9).h-X(3,k))*(obs(9).h-X(3,k))^2;
   costobs9 = checkobs9*(obs(9).h-X(3,k));
   obj = obj + w_5*costobs9;
   
   obs10 = ((X(1,k)-obs(10).x)^2 + (X(2,k)-obs(10).y)^2 - r^2);
   checkobs10 = -1*floor(obs10/100000000);
%    costobs10 = checkobs10*(obs(10).h-X(3,k))*(obs(10).h-X(3,k))^2;
   costobs10 = checkobs10*(obs(10).h-X(3,k));
   obj = obj + w_5*costobs10;
   
   obs11 = ((X(1,k)-obs(11).x)^2 + (X(2,k)-obs(11).y)^2 - r^2);
   checkobs11 = -1*floor(obs11/100000000);
%    costobs11 = checkobs11*(obs(11).h-X(3,k))*(obs(11).h-X(3,k))^2;
   costobs11 = checkobs11*(obs(11).h-X(3,k));
   obj = obj + w_5*costobs11;
   
   obs12 = ((X(1,k)-obs(12).x)^2 + (X(2,k)-obs(12).y)^2 - r^2);
   checkobs12 = -1*floor(obs12/100000000);
%    costobs12 = checkobs12*(obs(12).h-X(3,k))*(obs(12).h-X(3,k))^2;
   costobs12 = checkobs12*(obs(12).h-X(3,k));
   obj = obj + w_5*costobs12;
   
   obs13 = ((X(1,k)-obs(13).x)^2 + (X(2,k)-obs(13).y)^2 - r^2);
   checkobs13 = -1*floor(obs13/100000000);
%    costobs13 = checkobs13*(obs(13).h-X(3,k))*(obs(13).h-X(3,k))^2;
   costobs13 = checkobs13*(obs(13).h-X(3,k));
   obj = obj + w_5*costobs13;

   st_u = X(:,k); c_u = U(:,k); st_g = XG(:,k); con_g = UG(:,k);
   st_u_next = X(:,k+1); st_g_next = XG(:,k+1);
   st_u_next_e = st_u + T*f(st_u,c_u); st_g_next_e = st_g + T*fG(st_g,con_g);
   st_next = [st_u_next;st_g_next]; st_next_e = [st_u_next_e;st_g_next_e];
   g = [g;st_next - st_next_e];
end

for k = 1:N+1   % box constraints due to the map margins
%     g = [g;(-sqrt((X(1,k)-obs(1).x)^2 + (X(2,k)-obs(1).y)^2) + obs(1).r)]; % obstacle 1
%     g = [g;(-sqrt((X(1,k)-obs(2).x)^2 + (X(2,k)-obs(2).y)^2) + obs(2).r)]; % obstacle 2
%     g = [g;(-sqrt((X(1,k)-obs(3).x)^2 + (X(2,k)-obs(3).y)^2) + obs(3).r)]; % obstacle 3
%     g = [g;(-sqrt((X(1,k)-obs(4).x)^2 + (X(2,k)-obs(4).y)^2) + obs(4).r)]; % obstacle 4
%     g = [g;(-sqrt((X(1,k)-obs(5).x)^2 + (X(2,k)-obs(5).y)^2) + obs(5).r)]; % obstacle 5
%     g = [g;(-sqrt((X(1,k)-obs(6).x)^2 + (X(2,k)-obs(6).y)^2) + obs(6).r)]; % obstacle 6
%     g = [g;(-sqrt((X(1,k)-obs(7).x)^2 + (X(2,k)-obs(7).y)^2) + obs(7).r)]; % obstacle 7
%     g = [g;(-sqrt((X(1,k)-obs(8).x)^2 + (X(2,k)-obs(8).y)^2) + obs(8).r)]; % obstacle 8
%     g = [g;(-sqrt((X(1,k)-obs(9).x)^2 + (X(2,k)-obs(9).y)^2) + obs(9).r)]; % obstacle 9
%     g = [g;(-sqrt((X(1,k)-obs(10).x)^2 + (X(2,k)-obs(10).y)^2) + obs(10).r)]; % obstacle 10
%     g = [g;(-sqrt((X(1,k)-obs(11).x)^2 + (X(2,k)-obs(11).y)^2) + obs(11).r)]; % obstacle 11
%     g = [g;(-sqrt((X(1,k)-obs(12).x)^2 + (X(2,k)-obs(12).y)^2) + obs(12).r)]; % obstacle 12
%     g = [g;(-sqrt((X(1,k)-obs(13).x)^2 + (X(2,k)-obs(13).y)^2) + obs(13).r)]; % obstacle 13
%     g = [g;(-sqrt((X(1,k)-obs(14).x)^2 + (X(2,k)-obs(14).y)^2) + obs(14).r)]; % obstacle 14
%       slope = (X(2,k)-PT(3*k-1))/(X(1,k)-PT(3*k-2));
%       g = [g; (abs(slope*obs(1).x - 1*obs(1).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(1).r ];  
%       g = [g; (abs(slope*obs(2).x - 1*obs(2).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(2).r ];
%       g = [g; (abs(slope*obs(3).x - 1*obs(3).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(3).r ];
%       g = [g; (abs(slope*obs(4).x - 1*obs(4).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(4).r ];
%       g = [g; (abs(slope*obs(5).x - 1*obs(5).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(5).r ];
%       g = [g; (abs(slope*obs(6).x - 1*obs(6).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(6).r ];
%       g = [g; (abs(slope*obs(7).x - 1*obs(7).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(7).r ];
%       g = [g; (abs(slope*obs(8).x - 1*obs(8).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(8).r ];
%       g = [g; (abs(slope*obs(9).x - 1*obs(9).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(9).r ];
%       g = [g; (abs(slope*obs(10).x - 1*obs(10).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(10).r ];
%       g = [g; (abs(slope*obs(11).x - 1*obs(11).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(11).r ];
%       g = [g; (abs(slope*obs(12).x - 1*obs(12).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(12).r ];
%       g = [g; (abs(slope*obs(13).x - 1*obs(13).y + PT(3*k-1) - PT(3*k-2)*slope)/sqrt(slope^2 + (-1)^2))-obs(13).r ];

end

X0 = [X;XG];
U0 = [U;UG];
OPT_variables = [reshape(X0,(n_states+n_statesG)*(N+1),1);reshape(U0,(n_controls+n_controlsG)*N,1)];
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', p);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

args = struct;
%% Constraints
args.ubg(1:(n_states+n_statesG)*(N+1)) = 0; 
args.lbg(1:(n_states+n_statesG)*(N+1)) = 0; % equality constraints
% args.ubg((n_states+n_statesG)*(N+1)+1 : (n_states+n_statesG)*(N+1)+ n_obs*(N+1)) = 0; 
% args.lbg((n_states+n_statesG)*(N+1)+1 : (n_states+n_statesG)*(N+1)+ n_obs*(N+1)) = -inf; % inequality constraints

args.ubx(1:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) = inf;  args.lbx(1:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) = -inf; %state x 
args.ubx(2:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) = inf;  args.lbx(2:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) = -inf; %state y
args.ubx(3:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) = 200;  args.lbx(3:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) =   75; %state z
args.ubx(4:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) = 25;   args.lbx(4:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) =    14; %state v
args.ubx(5:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) = pi/4; args.lbx(5:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) = -pi/4; %state gamma pitch
args.ubx(6:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) = inf;  args.lbx(6:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) =  -inf; %state psi yaw
args.ubx(7:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) = pi/6; args.lbx(7:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) = -pi/6; 
args.ubx(8:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) = pi/6; args.lbx(8:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) = -pi/6; 
args.ubx(9:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) = inf;  args.lbx(9:(n_states+n_statesG):(n_states+n_statesG)*(N+1),1) = -inf; 

args.ubx((n_states+n_statesG)*(N+1)+1:(n_controls+n_controlsG):(n_states+n_statesG)*(N+1)+(n_controls+n_controlsG)*N,1) = pi/4;  
args.lbx((n_states+n_statesG)*(N+1)+1:(n_controls+n_controlsG):(n_states+n_statesG)*(N+1)+(n_controls+n_controlsG)*N,1) = -pi/4; %state alpha
args.ubx((n_states+n_statesG)*(N+1)+2:(n_controls+n_controlsG):(n_states+n_statesG)*(N+1)+(n_controls+n_controlsG)*N,1) = pi/10;   
args.lbx((n_states+n_statesG)*(N+1)+2:(n_controls+n_controlsG):(n_states+n_statesG)*(N+1)+(n_controls+n_controlsG)*N,1) = -pi/10; %state phi
args.ubx((n_states+n_statesG)*(N+1)+3:(n_controls+n_controlsG):(n_states+n_statesG)*(N+1)+(n_controls+n_controlsG)*N,1) = 20;     
args.lbx((n_states+n_statesG)*(N+1)+3:(n_controls+n_controlsG):(n_states+n_statesG)*(N+1)+(n_controls+n_controlsG)*N,1) = 0; %state Thrust
args.ubx((n_states+n_statesG)*(N+1)+4:(n_controls+n_controlsG):(n_states+n_statesG)*(N+1)+(n_controls+n_controlsG)*N,1) = pitch_maxG;
args.lbx((n_states+n_statesG)*(N+1)+4:(n_controls+n_controlsG):(n_states+n_statesG)*(N+1)+(n_controls+n_controlsG)*N,1) = -pitch_maxG;
args.ubx((n_states+n_statesG)*(N+1)+5:(n_controls+n_controlsG):(n_states+n_statesG)*(N+1)+(n_controls+n_controlsG)*N,1) = roll_maxG;  
args.lbx((n_states+n_statesG)*(N+1)+5:(n_controls+n_controlsG):(n_states+n_statesG)*(N+1)+(n_controls+n_controlsG)*N,1) = -roll_maxG;
args.ubx((n_states+n_statesG)*(N+1)+6:(n_controls+n_controlsG):(n_states+n_statesG)*(N+1)+(n_controls+n_controlsG)*N,1) = yaw_maxG;   
args.lbx((n_states+n_statesG)*(N+1)+6:(n_controls+n_controlsG):(n_states+n_statesG)*(N+1)+(n_controls+n_controlsG)*N,1) = -yaw_maxG;
%% THE SIMULATION LOOP STARTS HERE%%
x0 = [10;1;140;18;0;pi/2]; % also change below
x0T = [1;1;pi/2];
x0G = [0;0;0];  % also change below

u0 = zeros((n_controls+n_controlsG),N);
x0c = repmat([x0;x0G],1,N+1);

velT = 12;
angRateT = 0;
u0T = [velT;angRateT];

sim_tim = 130;
x_prediction = zeros((n_states+n_statesG),N+1);
u_prediction = zeros((n_controls+n_controlsG),N);
x_historyT = x0T;
x_predictionT =  zeros(n_statesT*N,1);
% Starting MPC
count =0;
start_timer = tic;
while(norm((x0(1:2,1)-x0T(1:2,1)),2) < 1e-2 || count < sim_tim/T)

%     args.p   = [x0;x0G;x0T]; % set the values of the parameters vector
    args.p(1:(n_states+n_statesG+n_statesT)) = [x0;x0G;x0T]; 
    in_pT = x0T;
    for k =1:N
        vp_T = 12; wp_T =  -rand*pi/40;
        con_pT = [vp_T;wp_T];
        fTp_value = full(fT(in_pT,con_pT));
        x_predictionT(3*k-2:3*k,1) = in_pT + T*fTp_value;
        in_pT = x_predictionT(3*k-2:3*k);
    end
    args.p( (n_states+n_statesG+n_statesT+1):(n_states+n_statesG+n_statesT*(N+1))) = x_predictionT;
    % initial value of the optimization variables
    args.x0  = [reshape(x0c,(n_states+n_statesG)*(N+1),1);reshape(u0,(n_controls+n_controlsG)*N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    solution = full(sol.x);
    x_p = reshape(solution(1:(n_states+n_statesG)*(N+1)),(n_states+n_statesG),(N+1));
    u_p = reshape(solution((n_states+n_statesG)*(N+1)+1:end),(n_controls+n_controlsG),N); % predicted controls

    
    u0 = [u_p(:,2:N),u_p(:,N)];
    x0c = [x_p(:,2:N+1),x_p(:,N+1)];
    x0 = x_p(1:n_states,2);
    x0G = x_p((n_states+1):(n_states+n_statesG),2);
    
    ffT_value = ffT(u0T,x0T);
    x0T = full(ffT_value(:,2));
    velT = 12;
    angRateT = 0;
    u0T = [velT;angRateT];
    
    if (count >= 250 && count <= 349)
        velT = 12;
        angRateT = -pi/40;
        u0T = [velT;angRateT];
    elseif(count > 349)
        angRateT = 0;
        u0T = [velT;angRateT];
    end
    if (count >= 550 && count <= 649)
        velT = 12;
        angRateT = -pi/45;
        u0T = [velT;angRateT];
    elseif(count > 649)
        angRateT = 0;
        u0T = [velT;angRateT];
    end
    if (count >= 850 && count <= 949)
        velT = 12;
        angRateT = -pi/40;
        u0T = [velT;angRateT];
    elseif(count > 949)
        angRateT = 0;
        u0T = [velT;angRateT];
    end
    
    % storing values
    x_prediction(:,:,(count+1)) = x_p;
    u_prediction(:,:,(count+1)) = u_p;
    x_historyT = [x_historyT,x0T];
    
    count
    count = count+1;
end

end_timer = toc(start_timer);
average_time = end_timer/(count+1)
total_time = average_time*count

%% Storing data

% Storing data
x_history = zeros((n_states+n_statesG),(count+1));
u_history = zeros((n_controls+n_controlsG),(count+1));
x_history(:,1) = [[10;1;90;18;0;pi/2];[0;0;0]]; % uav + gimbal initial states
u_history(:,1) = [[0;0;0];[0;0;0]];      % uav + gimbal initial controls
for k = 1:count
   x_history(:,k+1) = x_prediction(:,2,k);
   u_history(:,k+1) = u_prediction(:,1,k);
end
save('test11OffTV.mat')
Graph_MPC_casadi11 (count,x_history,x_historyT,u_history,x_prediction,HFOV_rad,VFOV_rad,obs,T,w_2,w_3,n_obs)
