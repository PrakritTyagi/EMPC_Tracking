function Graph_MPC_casadi11 (count,x_history,x_historyT,u_history,x_prediction,HFOV_rad,VFOV_rad,obs,T,w_2,w_3,n_obs)
fig = figure;

for k=1:(count)
    clf

    x_uav = x_history(1,k);
    y_uav = x_history(2,k);
    z_uav = x_history(3,k);
    v_uav = x_history(4,k);
    pitch_uav = x_history(5,k);
    yaw_uav = x_history(6,k);
    pitchG = x_history(7,k);
    rollG = x_history(8,k);
    yawG = x_history(9,k);
   
    x_p = x_prediction(1,:,k);
    y_p = x_prediction(2,:,k);
    z_p = x_prediction(3,:,k);
    x_tar = x_historyT(1,k);
    y_tar = x_historyT(2,k);
    yaw_tar = x_historyT(3,k);
%     
    hold on
%   target trajectory     
    p1 = plot3(x_historyT(1,1:k),x_historyT(2,1:k),0*x_historyT(3,1:k),'--','color','blue','linewidth',2);
    
%    UAV trajectory
    p2 = plot3(x_history(1,1:k),x_history(2,1:k),x_history(3,1:k),'color','red','LineWidth',2);
    
%    UAV prediction
    p3 = plot3(x_p,y_p,z_p,':','linewidth',2);

%    FOV plotting 
    th = 0:pi/50:2*pi;
%     a =  z_uav*tan(VFOV_rad/2);  
%     b =  z_uav*tan(HFOV_rad/2);
    a =  0.5*( z_uav*tan( pitchG + pitch_uav + VFOV_rad/2 ) - z_uav*tan( pitchG + pitch_uav - VFOV_rad/2 ));  
    b =  0.5*( z_uav*tan( rollG + HFOV_rad/2 ) - z_uav*tan( rollG - HFOV_rad/2 ));
    
    x = a*sin(th);
    y = b*cos(th);
    xprime = x_uav + a + z_uav*tan( pitchG + pitch_uav - VFOV_rad/2 ) + x*cos(yaw_uav + yawG) - y*sin(yaw_uav + yawG);
    yprime = y_uav + b + z_uav*tan( rollG - HFOV_rad/2 ) + x*sin(yaw_uav + yawG) + y*cos(yaw_uav + yawG);
    p4 = plot3(xprime,yprime,0*cos(th));

%     rectangle
    xe = [a,a,-a,-a,a];
    ye = [b,-b,-b,b,b];
    xedge = x_uav + a + z_uav*tan( pitchG + pitch_uav - VFOV_rad/2 ) + xe*cos(yaw_uav + yawG) - ye*sin(yaw_uav + yawG);
    yedge = y_uav + b + z_uav*tan( rollG - HFOV_rad/2 ) + xe*sin(yaw_uav + yawG) + ye*cos(yaw_uav + yawG);
    plot(xedge,yedge);

%     ray
    rayx = [xedge(1:4);(x_uav + 0*[0,0,0,0])];
    rayy = [yedge(1:4);(y_uav + 0*[0,0,0,0])];
    rayz = [0*[0,0,0,0];(z_uav + 0*[0,0,0,0])];
    p4 = plot3(rayx,rayy,rayz,'color','black');

%       Los
%     ta = [5,5,-5,-5,5];
%     tb = [5,-5,-5,5,5];
%     
%     x_edgetar = x_tar + ta*cos(yaw_tar) - tb*sin(yaw_tar);
%     y_edgetar = y_tar + ta*sin(yaw_tar) + tb*cos(yaw_tar);
%     plot(x_edgetar,y_edgetar,'linewidth',4);
% %     plot LOS
%     xlos=[x_history(1,1:k);x_historyT(1,1:k)];
%     ylos=[x_history(2,1:k);x_historyT(2,1:k)];
%     zlos=[x_history(3,1:k);x_historyT(3,1:k)];
%     plot3(xlos,ylos,zlos);

    %plotting obstacles
    x_unit1 = obs(1).x + obs(1).r*cos(th);
    y_unit1 = obs(1).y + obs(1).r*sin(th);
    z_unit1 = 0*obs(1).r*cos(th);
%     plot3(x_unit1,y_unit1,z_unit1,'color','blue') % obstacle 1
%     plot3(obs(1).x,obs(1).y,0,'g*','LineWidth',3,'MarkerSize', 8) % obstacle 1
    [X1,Y1,Z1] = cylinder(obs(1).r);
    X1 = X1 + obs(1).x; Y1 = Y1 + obs(1).y; Z1 = Z1*obs(1).h;
    p5 = surf(X1,Y1,Z1);
    
    x_unit2 = obs(2).x + obs(2).r*cos(th);
    y_unit2 = obs(2).y + obs(2).r*sin(th);
    z_unit2 = 0*obs(2).r*cos(th);
%     plot3(x_unit2,y_unit2,z_unit2,'color','blue') % obstacle 2
%     plot3(obs(2).x,obs(2).y,0,'g*','LineWidth',3,'MarkerSize', 8) % obstacle 2
    [X2,Y2,Z2] = cylinder(obs(2).r);
    X2 = X2 + obs(2).x; Y2 = Y2 + obs(2).y; Z2 = Z2*obs(2).h;
    p5 = surf(X2,Y2,Z2);
    
    x_unit3 = obs(3).x + obs(3).r*cos(th);
    y_unit3 = obs(3).y + obs(3).r*sin(th);
    z_unit3 = 0*obs(3).r*cos(th);
%     plot3(x_unit3,y_unit3,z_unit3,'color','blue') % obstacle 3
    [X3,Y3,Z3] = cylinder(obs(3).r);
    X3 = X3 + obs(3).x; Y3 = Y3 + obs(3).y; Z3 = Z3*obs(3).h;
    p5 = surf(X3,Y3,Z3);
    
    x_unit4 = obs(4).x + obs(4).r*cos(th);
    y_unit4 = obs(4).y + obs(4).r*sin(th);
    z_unit4 = 0*obs(4).r*cos(th);
%     plot3(x_unit4,y_unit4,z_unit4,'color','blue') % obstacle 4
%     plot3(obs(4).x,obs(4).y,0,'g*','LineWidth',3,'MarkerSize', 8) % obstacle 4
    [X4,Y4,Z4] = cylinder(obs(4).r);
    X4 = X4 + obs(4).x; Y4 = Y4 + obs(4).y; Z4 = Z4*obs(4).h;
    p5 = surf(X4,Y4,Z4);
    
    x_unit5 = obs(5).x + obs(5).r*cos(th);
    y_unit5 = obs(5).y + obs(5).r*sin(th);
    z_unit5 = 0*obs(5).r*cos(th);
%     plot3(x_unit5,y_unit5,z_unit5,'color','blue') % obstacle 5
    [X5,Y5,Z5] = cylinder(obs(5).r);
    X5 = X5 + obs(5).x; Y5 = Y5 + obs(5).y; Z5 = Z5*obs(5).h;
    p5 = surf(X5,Y5,Z5);
    
    x_unit6 = obs(6).x + obs(6).r*cos(th);
    y_unit6 = obs(6).y + obs(6).r*sin(th);
    z_unit6 = 0*obs(6).r*cos(th);
%     plot3(x_unit6,y_unit6,z_unit6,'color','blue') % obstacle 6
%     plot3(obs(6).x,obs(6).y,0,'g*','LineWidth',3,'MarkerSize', 8) % obstacle 6
    [X6,Y6,Z6] = cylinder(obs(6).r);
    X6 = X6 + obs(6).x; Y6 = Y6 + obs(6).y; Z6 = Z6*obs(6).h;
    p5 = surf(X6,Y6,Z6);
    
    x_unit7 = obs(7).x + obs(7).r*cos(th);
    y_unit7 = obs(7).y + obs(7).r*sin(th);
    z_unit7 = 0*obs(7).r*cos(th);
%     plot3(x_unit7,y_unit7,z_unit7,'color','blue') % obstacle 7
%     plot3(obs(7).x,obs(7).y,0,'g*','LineWidth',3,'MarkerSize', 8) % obstacle 7
    [X7,Y7,Z7] = cylinder(obs(7).r);
    X7 = X7 + obs(7).x; Y7 = Y7 + obs(7).y; Z7 = Z7*obs(7).h;
    p5 = surf(X7,Y7,Z7);
    
    x_unit8 = obs(8).x + obs(8).r*cos(th);
    y_unit8 = obs(8).y + obs(8).r*sin(th);
    z_unit8 = 0*obs(8).r*cos(th);
%     plot3(x_unit8,y_unit8,z_unit8,'color','blue') % obstacle 8
    [X8,Y8,Z8] = cylinder(obs(8).r);
    X8 = X8 + obs(8).x; Y8 = Y8 + obs(8).y; Z8 = Z8*obs(8).h;
    p5 = surf(X8,Y8,Z8);
    
    x_unit9 = obs(9).x + obs(9).r*cos(th);
    y_unit9 = obs(9).y + obs(9).r*sin(th);
    z_unit9 = 0*obs(9).r*cos(th);
%     plot3(x_unit9,y_unit9,z_unit9,'color','blue') % obstacle 9
    [X9,Y9,Z9] = cylinder(obs(9).r);
    X9 = X9 + obs(9).x; Y9 = Y9 + obs(9).y; Z9 = Z9*obs(9).h;
    p5 = surf(X9,Y9,Z9);

    x_unit10 = obs(10).x + obs(10).r*cos(th);
    y_unit10 = obs(10).y + obs(10).r*sin(th);
    z_unit10 = 0*obs(10).r*cos(th);
%     plot3(x_unit10,y_unit10,z_unit10,'color','blue') % obstacle 10
%     plot3(obs(10).x,obs(10).y,0,'g*','LineWidth',3,'MarkerSize', 8) % obstacle 10
    [X10,Y10,Z10] = cylinder(obs(10).r);
    X10 = X10 + obs(10).x; Y10 = Y10 + obs(10).y; Z10 = Z10*obs(10).h;
    p5 = surf(X10,Y10,Z10);
    
    x_unit11 = obs(11).x + obs(11).r*cos(th);
    y_unit11 = obs(11).y + obs(11).r*sin(th);
    z_unit11 = 0*obs(11).r*cos(th);
%     plot3(x_unit11,y_unit11,z_unit11,'color','blue') % obstacle 11
%     plot3(obs(11).x,obs(11).y,0,'g*','LineWidth',3,'MarkerSize', 8) % obstacle 11
    [X11,Y11,Z11] = cylinder(obs(11).r);
    X11 = X11 + obs(11).x; Y11 = Y11 + obs(11).y; Z11 = Z11*obs(11).h;
    p5 = surf(X11,Y11,Z11);
    
    x_unit12 = obs(12).x + obs(12).r*cos(th);
    y_unit12 = obs(12).y + obs(12).r*sin(th);
    z_unit12 = 0*obs(12).r*cos(th);
%     plot3(x_unit12,y_unit12,z_unit12,'color','blue') % obstacle 12
    [X12,Y12,Z12] = cylinder(obs(12).r);
    X12 = X12 + obs(12).x; Y12 = Y12 + obs(12).y; Z12 = Z12*obs(12).h;
    p5 = surf(X12,Y12,Z12);
    
    x_unit13 = obs(13).x + obs(13).r*cos(th);
    y_unit13 = obs(13).y + obs(13).r*sin(th);
    z_unit13 = 0*obs(13).r*cos(th);
%     plot3(x_unit13,y_unit13,z_unit13,'color','blue') % obstacle 13
    [X13,Y13,Z13] = cylinder(obs(13).r);
    X13 = X13 + obs(13).x; Y13 = Y13 + obs(13).y; Z13 = Z13*obs(13).h;
    p5 = surf(X13,Y13,Z13);
    
%     x_unit14 = obs(14).x + r*cos(th);
%     y_unit14 = obs(14).y + r*sin(th);
%     z_unit14 = 0*r*cos(th);
% %     plot3(x_unit14,y_unit14,z_unit14,'color','blue') % obstacle 14
%     [X14,Y14,Z14] = cylinder(r);
%     X14 = X14 + obs(14).x; Y14 = Y14 + obs(14).y; Z14 = Z14*obs(14).h;
%     p4 = surf(X14,Y14,Z14);

    h = [p1(1);p2(1);p3(1);p4(1)];
    legend(h,'Target Trajectory','UAV Trajectory','UAV Future Prediction','FOV');
    pbaspect([1 1 0.25])
    axis([-200 1300 -100 1000 0 200])
    xlabel('x(m)','fontweight','bold') 
    ylabel('y(m)','fontweight','bold')
    zlabel('z(m)','fontweight','bold')


    grid on
    axis equal
%     drawnow
    
    movieVector(k) = getframe(fig, [40 40 470 400]); 
end

%% save the movie
myWriter = VideoWriter('UAV_MPC11off', 'Motion JPEG AVI');
myWriter.FrameRate = 13;
open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);

exportgraphics(fig,'trajectory2.jpg','Resolution',1800)

% energy cost

% Uav Parameter
m = 1.56;
S = 0.2589;
b = 1.4224;
ar = b^2/S;
cdo = 0.01631;
kdl = 0.04525;
rho = 1.225;
a_g = 9.8;
% w_3 = 0.1;
energy_cost = zeros(1,(count+1)); 
energy_cost2 = zeros(1,(count+1));
for k = 3:(count+1)
    Co_drag = cdo + kdl*u_history(1,k)^2;
    energy =  w_3*(m*x_history(4,k)*(u_history(3,k) - 0.5*rho*S*Co_drag*x_history(4,k)^2 - m*a_g*sin(x_history(5,k)))/m + 0.5*Co_drag*rho*S*x_history(4,k)^3 +...
    2*(m*a_g)^2*cos(x_history(5,k))^2/(pi*rho*S*ar*x_history(4,k)) + m*a_g*x_history(4,k)*sin(x_history(5,k)));
    energy2 = w_3*(x_history(4,k)*u_history(3,k));
    energy_cost(1,k) = energy;
    energy_cost2(1,k) = energy2;
end
TE = sum(energy_cost)
z = linspace(0,T*count,(count+1));
energy = figure;
plot(z,energy_cost(1,:),'lineWidth',2);
grid on
xlim([0 T*count])
xlabel('Time','fontweight','bold') 
ylabel('Energy cost function','fontweight','bold')

figure;
plot(z,energy_cost2(1,:),'lineWidth',2);
grid on
xlim([0 T*count])
xlabel('Time','fontweight','bold') 
ylabel('Energy cost function2','fontweight','bold')

% E(n) plot
vi = figure;
ellipse_history = zeros(1,(count+1));
intersect_history = ones(1,(count+1));
notinVision = 0;
for k = 1:(count+1)
intersect =1;
costlos1 = 0;
a = 0.5*( x_history(3,k)*tan(x_history(7,k) + x_history(5,k) + VFOV_rad/2) - x_history(3,k)*tan(x_history(7,k) + x_history(5,k) - VFOV_rad/2) );
b = 0.5*( x_history(3,k)*tan(x_history(8,k) + HFOV_rad/2) - x_history(3,k)*tan(x_history(8,k) - HFOV_rad/2) );
shifta = a + x_history(3,k)*tan(x_history(7,k) + x_history(5,k) - VFOV_rad/2);
shiftb = b + x_history(3,k)*tan(x_history(8,k) - HFOV_rad/2);
A = ((cos(x_history(6,k) + x_history(9,k)))^2)/(a^2) + ((sin(x_history(6,k) + x_history(9,k)))^2)/(b^2) ;
B = 2*cos(x_history(6,k) + x_history(9,k))*sin(x_history(6,k) + x_history(9,k))*( (1/a^2) - (1/b^2) );
C = (sin(x_history(6,k) + x_history(9,k)))^2/(a^2) + (cos(x_history(6,k) + x_history(9,k)))^2/(b^2);
ellipse = A*(x_historyT(1,k) - x_history(1,k) - shifta)^2 + B*(x_historyT(1,k) - x_history(1,k) - shifta)*(x_historyT(2,k) - x_history(2,k) - shiftb) + C*(x_historyT(2,k) - x_history(2,k) - shiftb)^2 - 1; 
% cost = w_2*ellipse;
logic_TargetInVision = floor(w_2*ellipse/10000000)*-1; % 1 if inside ellipse else 0 for outside
ellipse_history(1,k) = logic_TargetInVision;
if logic_TargetInVision ==0
    notinVision = notinVision +1 ;   
end
for f=1:n_obs
slope = (x_history(2,k)-x_historyT(2,k))/(x_history(1,k)-x_historyT(1,k)); % slope of ray(x-yplane) from target to UAV 
center1_to_LOS = (abs(slope*obs(f).x - 1*obs(f).y + x_historyT(2,k) - x_historyT(1,k)*slope)/sqrt(slope^2 + (-1)^2))-obs(f).r; % distance between line LOS 2D and center of obstacle1
slopeRef1 = (obs(f).h/(sqrt( (x_historyT(1,k)-obs(f).x)^2 + (x_historyT(2,k)-obs(f).y)^2 )-obs(f).r)) - (x_history(3,k)/(sqrt( (x_historyT(1,k)-x_history(1,k))^2 + (x_historyT(2,k)-x_history(2,k))^2 ))); % Slope of LOS between obstacle_1 Top and target in 3D
logic_LosIntersectObstacle1 = floor(center1_to_LOS/100000000)*-1; % if intersects 1 otherwise 0 
costlos1 = costlos1 + logic_TargetInVision*logic_LosIntersectObstacle1*max(0,slopeRef1);
end
if (costlos1 > 0)
intersect = 0;
end
intersect_history(1,k) =  intersect ;
end
notinVision
z = linspace(0,T*count,(count+1));
vision = ellipse_history.*intersect_history;

plot(z,vision(1,:),'lineWidth',2);
xlim([0 T*(count-1)])
xlabel('Time','fontweight','bold') 
ylabel('E(n)','fontweight','bold')
grid on

figure;
plot(z,ellipse_history(1,:),'lineWidth',2);
hold on;
plot(z,intersect_history(1,:),'lineWidth',2,'color','green');
grid on


