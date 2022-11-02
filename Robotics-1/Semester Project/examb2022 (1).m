%% Robotics I 2021-2022
%% Alkiviadis Panagiotis Michalitsis
%% AM : 03118868
%% Part B : Simulation 

%% *** Robot (kinematic) model parameters *** 
clc
clear all
close all

l0 = 10.0;  l1 = 15.0;  l2 = 30.0;  l3 = 30.0;

%% *** sampling period *** 
 
dt = 0.001; 

%% *** DESIRED MOTION PROFILE - TASK SPACE *** 
Tf=2.0; 	% 2sec Tf of motion 
t=dt:dt:2*Tf;

%% initial & final end-point position & acceleration (cm,sec)

xA = 15;	yA = -30;   zA = 40;

xB = 15;     yB = -30;  zB = 50;

gx = 0.05;   gy = 0.05; gz = 0.05;

%% *** DESIRED MOTION PROFILE - TASK SPACE *** 
xd  = double(trajectory(0,Tf,xA,xB,gx));    yd  = double(trajectory(0,Tf,yA,yB,gy));    zd  = double(trajectory(0,Tf,zA,zB,gz));

rxd = double(trajectory(Tf,2*Tf,xB,xA,-gx));ryd = double(trajectory(Tf,2*Tf,yB,yA,-gy));rzd = double(trajectory(Tf,2*Tf,zB,zA,-gz));

kmax = Tf/dt + 1;

x = zeros(length(t),1); y= zeros(length(t),1);  z= zeros(length(t),1);

vx= zeros(length(t),1); vy= zeros(length(t),1); vz= zeros(length(t),1);

for k = 1:2*kmax-1
   time  = (k-1)*dt;
   if (time < Tf*0.1)
        x(k) = polyval(xd(6:-1:1),time);
        vx(k)= polyval(polyder(xd(6:-1:1)),time);
        y(k) = polyval(yd(6:-1:1),time);
        vy(k) = polyval(polyder(yd(6:-1:1)),time);
        z(k) = polyval(zd(6:-1:1),time);
        vz(k) = polyval(polyder(zd(6:-1:1)),time);
   elseif (time < Tf*0.9)
        x(k) = polyval(xd(8:-1:7),time);
        vx(k) = polyval(polyder(xd(8:-1:7)),time);
        y(k) = polyval(yd(8:-1:7),time);
        vy(k) = polyval(polyder(yd(8:-1:7)),time);
        z(k) = polyval(zd(8:-1:7),time);
        vz(k) = polyval(polyder(zd(8:-1:7)),time);
   elseif (time < Tf)
        x(k) = polyval(xd(14:-1:9),time);
        vx(k) = polyval(polyder(xd(14:-1:9)),time);
        y(k) = polyval(yd(14:-1:9),time);
        vy(k) = polyval(polyder(yd(14:-1:9)),time);
        z(k) = polyval(zd(14:-1:9),time);
        vz(k) = polyval(polyder(zd(14:-1:9)),time);
   elseif (time < Tf*1.1)
        x(k) = polyval(rxd(6:-1:1),time);
        vx(k) = polyval(polyder(rxd(6:-1:1)),time);
        y(k) = polyval(ryd(6:-1:1),time);
        vy(k) = polyval(polyder(ryd(6:-1:1)),time);
        z(k) = polyval(rzd(6:-1:1),time);
        vz(k) = polyval(polyder(rzd(6:-1:1)),time);
   elseif (time < Tf*1.9)
        x(k) = polyval(rxd(8:-1:7),time);
        vx(k) = polyval(polyder(rxd(8:-1:7)),time);
        y(k) = polyval(ryd(8:-1:7),time);
        vy(k) = polyval(polyder(ryd(8:-1:7)),time);
        z(k) = polyval(rzd(8:-1:7),time);
        vz(k) = polyval(polyder(rzd(8:-1:7)),time);
   else
        x(k) = polyval(rxd(14:-1:9),time);
        vx(k) = polyval(polyder(rxd(14:-1:9)),time);
        y(k) = polyval(ryd(14:-1:9),time);
        vy(k) = polyval(polyder(ryd(14:-1:9)),time);
        z(k) = polyval(rzd(14:-1:9),time);
        vz(k) = polyval(polyder(rzd(14:-1:9)),time);
   end
end

x  = x(1:length(t));    y = y(1:length(t));     z = z(1:length(t));

vx = vx(1:length(t));   vy = vy(1:length(t));   vz = vz(1:length(t));

subplot(6,1,1); 
plot(t,x); 
ylabel('x(cm)'); 
xlabel('time(sec)');
grid on

subplot(6,1,2); 
plot(t,y); 
ylabel('y(cm)'); 
xlabel('time(sec)');
grid on

subplot(6,1,3); 
plot(t,z); 
ylabel('z(cm)'); 
xlabel('time(sec)');
grid on

subplot(6,1,4); 
plot(t,vx); 
ylabel('vx(cm/sec)'); 
xlabel('time(sec)');
grid on

subplot(6,1,5); 
plot(t,vy); 
ylabel('vy(cm/sec)'); 
xlabel('time(sec)');
grid on

subplot(6,1,6); 
plot(t,vz); 
ylabel('vz(cm/sec)'); 
xlabel('time(sec)');
grid on

%% INVESRE KINEMATICS: DESIRED MOTION -> JOINT SPACE

p1x = sqrt(x.^2+y.^2 - l1^2);
%%-------------------------------------------------------%%
q3 = asin((p1x.^2 + y.^2 -l2^2 -l3^2)./(2* l2 * l3));
q2 = atan2(sqrt(x.^2 + y.^2 - l1^2), z-l0) + atan2(l3 .* cos(q3), l2 + l3.*cos(q3));
q1 = atan2(y, x) - atan2(l2.*sin(q2) - l3.*cos(q2).*cos(q3), l1);
%%------------------------------------------------------%%
qd1 = zeros(length(t),1); 
qd2 = zeros(length(t),1); 
qd3 = zeros(length(t),1);
%%--------------------------------------------------------%%
for k = 1:length(t)
    j = inverse(l0,l1,l2,l3,q1(k),q2(k),q3(k));
    qd1(k) = j(1,1)*vx(k) + j(1,2)*vy(k);
    qd2(k) = j(2,1)*vx(k) + j(2,2)*vy(k);
    qd3(k) = j(3,1)*vx(k) + j(3,2)*vy(k);
end
%%--------------------------------------------------------%%

fig2 = figure;

subplot(6,1,1); 
plot(t,q1);
ylabel('q1(rad)'); 
xlabel('time(sec)');
grid on

subplot(6,1,2); 
plot(t,q2); 
ylabel('q2(rad)'); 
xlabel('time(sec)');
grid on

subplot(6,1,3); 
plot(t,q3);
ylabel('q3(rad)'); 
xlabel('time(sec)');
grid on

subplot(6,1,4); 
plot(t,qd1);
ylabel('qd1(rad)'); 
xlabel('time(sec)');
grid on

subplot(6,1,5); 
plot(t,qd2); 
ylabel('qd2(rad)'); 
xlabel('time(sec)');
grid on

subplot(6,1,6); 
plot(t,qd3);
ylabel('qd3(rad)'); 
xlabel('time(sec)');
grid on

c1 = cos(q1); c2 = cos(q2); c3 = cos(q3); c23 = cos(q2+q3);
s1 = sin(q1); s2 = sin(q2); s3 = sin(q3); s23 = sin(q2+q3);

x1 = l1.*s1;
y1 = zeros(length(t),1);
z1 = l1.*c1;


x2 = l2*c1.*c2 - l1.*s1;
y2 = l2*s2;
z2 = l2.*c2.*s1 + l1.*c1 + l0;


x3 = l3*c1.*c23 + l2*c1.*c2 - l1.*s1;
y3 = l3.*s23 + l2.*s2;
z3 = l3.*s1.*c23 + l2.*c2.*s1 + l1.*c1+l0;

fig3= figure; 
axis([-5 5 -5 5 0 10])                                              
axis on 
hold on 
grid on
title('Stick simulation of trajectory')
xlabel('x (cm)'); 
ylabel('y (cm)'); 
zlabel ('z (cm)');
dtk=100;                                                              



for tk=1:dtk:length(t) 
   pause(0.1);                                                           
  
   plot3([x1],[y1],[z1],'o');  
   plot3([x1(tk),x2(tk)],[y1(tk),y2(tk)], [z1(tk),z2(tk)]);
   plot3([x2(tk)],[y2(tk)], [z2(tk)],'o');  
   plot3([x2(tk),x3(tk)],[y2(tk),y3(tk)],[z2(tk),z3(tk)]);	
   plot3([x3(tk)],[y3(tk)],[z3(tk)],'o');  
   
end
plot3(x3,y3,z3,'rs'); 



function Traj = trajectory(t0,tf,xA,xB,g)
    
    syms a0 a1 a2 a3 a4 a5 b0 b1 c0 c1 c2 c3 c4 c5
    
    t1 = t0+(tf-t0)*0.1; t2 = t0+(tf-t0)*0.9;
    
    sx1 = a0 + a1*t0 + a2*t0^2 + a3*t0^3 + a4*t0^4 + a5*t0^5 - xA == 0;

    sx2 = a1 + 2*a2*t0 + 3*a3*t0^2 + 4*a4*t0^3 + 5*a5*t0^4 == 0;

    sx3 = 2*a2 + 6*a3*t0 + 12*a4*t0^2 + 20*a5*t0^3 - g == 0;

    sx4 = a0 + a1*t1 + a2*t1^2 + a3*t1^3 + a4*t1^4 + a5*t1^5 - (b0 + b1*t1) == 0;

    sx5 = a1 + 2*a2*t1 + 3*a3*t1^2 + 4*a4*t1^3 + 5*a5*t1^4 -b1 == 0;

    sx6 = 2*a2 + 6*a3*t1 + 12*a4*t1^2 + 20*a5*t1^3 == 0;

    sx7 = c0 + c1*t2 + c2*t2^2 + c3*t2^3 + c4*t2^4 + c5*t2^5 - (b0 + b1*t2) == 0;

    sx8 = c1 + 2*c2*t2 + 3*c3*t2^2 + 4*c4*t2^3 + 5*c5*t2^4 - b1 == 0;
    
    sx9 = 2*c2 + 6*c3*t2 + 12*c4*t2^2 + 20*c5*t2^3 == 0;

    sx10 = c0 + c1*tf + c2*tf^2 + c3*tf^3 + c4*tf^4 + c5*tf^5 - xB == 0;

    sx11 = c1 + 2*c2*tf + 3*c3*tf^2 + 4*c4*tf^3 + 5*c5*tf^4 == 0;

    sx12 = 2*c2 + 6*c3*tf + 12*c4*tf^2 + 20*c5*tf^3 +g == 0;

    sx13 = b0+b1*t1 - (a0 + a1*t0 + a2*t0^2 + a3*t0^3 + a4*t0^4 + a5*t0^5) - (0.5*(t1-t0)*b1) == 0;

    sx14 = (c0 + c1*tf + c2*tf^2 + c3*tf^3 + c4*tf^4 + c5*tf^5) - (b0+b1*t2) - (0.5*(t1-t0)*b1) == 0;
    
    [A,B] = equationsToMatrix([ sx1 sx2 sx3 sx4 sx5 sx6 sx7 sx8 sx9 sx10 sx11 sx12 sx13 sx14 ] , [ a0 a1 a2 a3 a4 a5 b0 b1 c0 c1 c2 c3 c4 c5 ]);
    Traj = vpa(linsolve(A,B));
    
    
end

function Jinv = inverse(l0,l1,l2,l3,q1,q2,q3)

    rotz = [sin(q1) cos(q1) 0 0; -cos(q1) sin(q1) 0 0; 0 0 1 0; 0 0 0 1];
    traz = [1 0 0 0; 0 1 0 0; 0 0 1 l0; 0 0 0 1];
    trax = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    rotx = [1 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 1];
    a01 = rotz*traz*trax*rotx;
%%----------------------------------------------------------%%
    rotz = [-sin(q2) -cos(q2) 0 0; cos(q2) -sin(q2) 0 0; 0 0 1 0; 0 0 0 1];
    traz = [1 0 0 0; 0 1 0 0; 0 0 1 -l1; 0 0 0 1];
    trax = [1 0 0 l2; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    rotx = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    a12 = rotz*traz*trax*rotx;
    a02 = a01*a12;
%%---------------------------------------------------------%%
    rotz = [cos(q3)  -sin(q3) 0 0; sin(q3) cos(q3) 0 0; 0 0 1 0; 0 0 0 1];
    traz = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    trax = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    rotx = [1 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 1];
    a23 = rotz*traz*trax*rotx;
%%-----------------------------------------------------------%%
    rotz = [0 1 0 0; -1 0 0 0; 0 0 1 0; 0 0 0 1];
    traz = [1 0 0 0; 0 1 0 0; 0 0 1 l3; 0 0 0 1];
    trax = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    rotx = [1 0 0 0; 0 0 1 0; 0 -1 0 0; 0 0 0 1];
    a3e = rotz*traz*trax*rotx;
%%-------------------------------------------------------------%%
    a = a01*a12*a23*a3e;
%%-------------------------------------------------------------%%
    b0=[0;0;1;];    b1=a01(1:3,3);  b2=a02(1:3,3);
    r0E=a(1:3,4);   r1E=a(1:3,4)-a01(1:3,4);    r2E=a(1:3,4)-a02(1:3,4);

    JL1=cross(b0,r0E);  JL2=cross(b1,r1E);  JL3=cross(b2,r2E);

    J= [JL1,JL2,JL3; b0,b1,b2];
    
    Jinv = inv(J(1:3,:));
end



