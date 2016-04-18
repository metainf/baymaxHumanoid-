function [stepLength,forwardVel] = pwSagittal(r,slope)

% PWSAGITTAL 
%   This function simulates the sagittal plane dynamics of the passive
%   walker on a ramp.  We use it to design the radius of curvature.
%
%   PWSAGITTAL(r) uses r as the radius of curvature of the foot.
%   If r is omitted, the default value (which we used for our
%   robot) is used.  
%  
%   The function outputs the approximate stepping frequency of the
%   robot.
%
% This code was written by Russ Tedrake <russt@ai.mit.edu>.  Please
% acknowlege that fact if you find this code useful.  Feel free to
% email me with questions or comments.

global m I R d b g gamma phi;
% q = [theta_stanceleg; theta_swingleg];  Q = [q;qdot;bLeftStance];

m = 2.36;
I = 9.63;
if (nargin > 0) R = r; else R = 14.8; end
if (nargin > 1) gamma = slope; else gamma = -0.027; end
d = -.2;
b = 16.35;
g = 386.088;
phi = 0.23;

Q(:,1) = zeros(5,1);
dt = 1e-3;
displaydt = 0.025;  % change the speed of the drawing
T = 10;
steplength = [];

draw(Q(:,1),0);
lastt = 0;
tcollision = 0;
for t=1:T/dt
  Q(:,t+1) = Q(:,t) + dt*[Q(3:4,t);dynamics(Q(1:2,t),Q(3:4,t));0];
  if (t - tcollision > 0.35/dt)  % 1/2 period at 1.4 Hz
    Omega_m = [ 2*b*d*cos(Q(2,t+1)-Q(1,t+1))-(b+d)*R*cos(Q(2,t+1)-gamma) ...
	  - 2*b*R*cos(Q(1,t+1)-gamma) + 2*R^2 + b^2 - b*d, ...
	  (b-d)*(b-R*cos(Q(2,t+1)-gamma)); ...
	  (b-d)*(b-R*cos(Q(2,t+1)-gamma)), 0 ];
    Omega_p = [ (b-d)*(d*cos(Q(1,t+1)-Q(2,t+1)) - ...
	  R*cos(Q(1,t+1)-gamma)+(b-d)), ...
	  -R*(b-d)*cos(Q(1,t+1)-gamma) - ...
	  R*(b+2*d)*cos(Q(2,t+1)-gamma) + d^2 + 2*R^2 + ...
	  R*b*cos(Q(2,t+1)+gamma) - b^2*cos(2*Q(2,t+1)) + ...
	  d*(b-d)*cos(Q(1,t+1)-Q(2,t+1)); ...
	  (b-d)^2, (b-d)*(d*cos(Q(1,t+1)-Q(2,t+1)) - ...
	  R*cos(Q(1,t+1)-gamma))];
    steplength = [steplength, (rot(gamma)*[R*(Q(1,t+1)-gamma);R] - ...
	rot(Q(1,t+1))*[0;d] + rot(Q(2,t+1))*[0;d] - ...
	rot(gamma)*[0;R])'*[cos(gamma); sin(gamma)]];
    Q(3:4,t+1) = inv(Omega_p)*Omega_m*Q(3:4,t+1);
    Q(:,t+1) = [Q(2,t+1); Q(1,t+1); Q(4,t+1); Q(3,t+1); 1-Q(5,t+1)];
    tcollision = t;        
  end
  if (t - lastt > displaydt/dt)
    draw(Q(:,t),t*dt);
    lastt = t;
  end
end
draw(Q(:,end),T);

ql = (1-Q(5,:)).*Q(1,:) + Q(5,:).*Q(2,:);
qldot = (1-Q(5,:)).*Q(3,:) + Q(5,:).*Q(4,:);
qr = Q(5,:).*Q(1,:) + (1-Q(5,:)).*Q(2,:);
qrdot = Q(5,:).*Q(3,:) + (1-Q(5,:)).*Q(4,:);

figure(1);
plot(ql,qldot); xlabel('\theta_{lPitch}'); ylabel('d \theta_{lPitch} / dt'); 
%plot(ql);

stepLength = mean(steplength(floor(end/2):end));

v = ((d*cos(Q(1,:)) - R*cos(gamma)).*Q(3,:))*cos(gamma) + ((d*sin(Q(1,:)) - R*cos(gamma)).*Q(3,:))*sin(gamma);
forwardVel = mean(v(floor(end/2):end));  % in inches per second

disp(['step length = ', num2str(stepLength),' inches, forward velocity = ', num2str(forwardVel), ' inches/second']);


function qddot = dynamics(q,qdot)

global m I R d b g gamma phi;

s1 = sin(q(1)); c1 = cos(q(1));
s2 = sin(q(2)); c2 = cos(q(2));
sg = sin(gamma); cg = cos(gamma);

H11 = I + m*b^2 + m*d^2 + 2*m*R^2 - 2*m*R*(b+d)*cos(q(1)-gamma);
H12 = m*(b-d)*(d*cos(q(1)-q(2))-R*cos(q(2)-gamma));
H22 = I + m*(b-d)^2;
H = [H11,H12;H12,H22];

C11 = m*R*(b+d)*qdot(1)*sin(q(1)-gamma) + 0.5*m*d*(b-d)*qdot(2)*sin(q(1)-q(2));
C12 = m*(b-d)*(d*sin(q(1)-q(2))*(qdot(2)-0.5*qdot(1)) + R*sin(q(2)-gamma)*qdot(2));
C21 = m*(b-d)*(d*sin(q(1)-q(2))*(qdot(1)-0.5*qdot(2)) - 0.5*R*sin(q(2)-gamma)*qdot(2));
C22 = 0.5*m*(b-d)*qdot(1)*(d*sin(q(1)-gamma)+R*sin(q(2)-gamma));
C = [C11,C12;C21,C22];

G = m*g*[(b+d)*sin(q(1)) - 2*R*sin(gamma); (b-d)*sin(q(2))];

qddot = inv(H)*(-C*qdot - G);

function draw(q,t)

%return;  % uncomment this to run lots of trials without displaying anything

global m I R d b g gamma phi;

alumc = [0.75, 0.75, 0.75];
woodc = [0.91,0.76,0.65];
  
alpha = max(min(q(1)-gamma,phi),-phi);
swalpha = max(min(q(2)-gamma,phi),-phi);

gc = rot(gamma)*[-R*alpha;0];  
cen = gc + rot(q(1)-alpha)*[0;R];
hip = cen - rot(q(1))*[0;d];
swingcen = hip + rot(q(2))*[0;d];
swinggc = swingcen + rot(q(2)-swalpha)*[0;-R];

axle = 0.25*[sin(0:0.1:2*pi);cos(0:0.1:2*pi)];
ankle = [0; -11.5];%[0; -7.625 ];
leg = [ 0.625*[-1,1,1,-1,-1];
      0.75,   0.75,  ankle(2),  ankle(2), 0.75];

xcen = -ankle(1); ycen = -ankle(2) + d;
footpoints = 20;
xfootlimits = R*[sin(phi),sin(-phi)]+xcen;
theta = [-phi:2*phi/footpoints:phi];
foot = [ xfootlimits(1), xfootlimits(2), R*sin(theta)+xcen, xfootlimits(1); ...
	0.0, 0.0, -R*cos(theta)+ycen, 0.0];

h = figure(25);
set(h,'DoubleBuffer','on');

clf;
hold on;

Axle = repmat(hip,1,size(axle,2)) + axle;
for i=1:size(leg,2)
  StanceLeg(:,i) = hip + rot(q(1))*leg(:,i);
  SwingLeg(:,i) = hip + rot(q(2))*leg(:,i);
end

StanceAnkle = hip + rot(q(1))*ankle;
SwingAnkle = hip + rot(q(2))*ankle;
for i=1:size(foot,2)
  StanceFoot(:,i) = StanceAnkle + rot(q(1))*foot(:,i);
  SwingFoot(:,i) = SwingAnkle + rot(q(2))*foot(:,i);
end

if (q(5)>0.5)
  fill(StanceLeg(1,:),StanceLeg(2,:),alumc);
  fill(StanceFoot(1,:),StanceFoot(2,:),woodc);
  fill(SwingLeg(1,:),SwingLeg(2,:),alumc);
  fill(SwingFoot(1,:),SwingFoot(2,:),woodc);
  plot(gc(1),gc(2),'r*');
else
  fill(SwingLeg(1,:),SwingLeg(2,:),alumc);
  fill(SwingFoot(1,:),SwingFoot(2,:),woodc);
  fill(StanceLeg(1,:),StanceLeg(2,:),alumc);
  fill(StanceFoot(1,:),StanceFoot(2,:),woodc);
  plot(gc(1),gc(2),'b*');
end

fill(Axle(1,:),Axle(2,:),alumc);

line([-15,15],[-15,15]*tan(gamma),'Color',[0 0 0]);

axis([-15 15 -9 21]);
%axis equal; 
%axis off;

title(['t = ', num2str(t)]);

drawnow;

function r = rot(q)

r = [cos(q) -sin(q); sin(q) cos(q)]; 
