function stepfreq = pwFront(r)

% PWFRONT
%   This function simulates the frontal plane dynamics of the
%   passive walker on flat terrain.  We use it to design the
%   radius of curvature.
%
%   PWFRONT(r) uses r as the radius of curvature of the foot.
%   If r is omitted, the default value (which we used for our
%   robot) is used.  
%  
%   The function outputs the approximate stepping frequency of the
%   robot.
%
% This code was written by Russ Tedrake <russt@ai.mit.edu>.  Please
% acknowlege that fact if you find this code useful.  Feel free to
% email me with questions or comments.

global m I R a g phi;
% q = [theta];  Q = [q;qdot];

m = 7.39; %pounds
I = 102.71; %lb*ft^2
if (nargin > 0) R = r; else R = 19.6; end
a = 11.6; %in
g = 386.088; %in/s^2
phi = 0.0872665; %rad

Q(:,1) = [asin(3.5/R);0];
dt = 1e-3;  % make this smaller to see limit cycles decay
displaydt = 0.05;  % change the speed of the drawing
T = 10;

draw(Q(:,1),0);
lastt = 0;
for t=1:T/dt
  Q(:,t+1) = Q(:,t) + dt*[Q(2,t);dynamics(Q(1,t),Q(2,t))];
  if (sign(Q(1,t)) ~= sign(Q(1,t+1))) % just crossed theta=0, model collision
    Q(2,t+1) = Q(2,t+1)*cos(2*atan2(R*sin(phi),(R*cos(phi)-a)));
  end
  if (t - lastt > displaydt/dt)
    draw(Q(:,t),t*dt);
    lastt = t;
  end
end
draw(Q(:,end),T);

% compute step frequency
[Ptheta,Ftheta] = psd(Q(1,:),length(Q(1,:)),1/dt);
[pi,i] = max(Ptheta);
stepfreq = Ftheta(i);


function qddot = dynamics(q,qdot)

global m I R a g phi;

if (abs(q) > phi)
  H = I + m*a^2 + m*R^2 - 2*m*a*R*cos(q);
  C = m*R*a*qdot*sin(q);
  G = m*g*a*sin(q);
else
  if (q>0) alpha = q-phi; else alpha = q+phi; end
  H = I + m*a^2 + m*R^2 - 2*m*R*a*cos(q-alpha);
  C = 0;
  G = m*g*(a*sin(q)-R*sin(alpha));
end

qddot = (-C-G)/H;



function draw(q,t)

global R phi;

alumc = [0.75, 0.75, 0.75];
woodc = [0.91,0.76,0.65];
  
ah = 9;

rot = inline('[cos(q) -sin(q); sin(q) cos(q)]','q');  

if (q(1)>0) 
  alpha = q(1)-phi; 
  sq = 1;
else 
  alpha = q(1)+phi; 
  sq = -1;
end
if (abs(q(1)) > phi)
  hip = [-R*alpha - sq*R*sin(phi) + (R-ah)*sin(q(1));R - (R-ah)*cos(q(1))];
else
  hip = [-sq*R*sin(phi);0] + rot(alpha)*[0;R] + rot(q(1))*[0;-R+ah];
end
persistent axle leg ankle foot r;

if (isempty(r) | r~=R) 
  axle = [ 2, 2,  -2, -2,  2;
    -0.125, 0.125, 0.125, -0.125, -0.125];
  leg = [ 0.5 + [0, 1.2, 1.2, 1, 1, 0, 0]; ...
	0.75,   0.75,  -8,  -8, -1.125,  -1.125, 0.75];
  
  ankle = [0.5 + 1.08; -7.625 ];
  
  xcen = -ankle(1); ycen = -ankle(2) + R - ah;
  footpoints = 10;
  xfootlimits = [-1,2.25];
  thetamin = asin((xfootlimits(1)-xcen)/R);
  thetamax = asin((xfootlimits(2)-xcen)/R);
  theta = [thetamax:-(thetamax-thetamin)/footpoints:thetamin];
  foot = [ xfootlimits(1), xfootlimits(2), R*sin(theta)+xcen, xfootlimits(1); ...
    0.0, 0.0, -R*cos(theta)+ycen, 0.0];
  r = R;
end

h = figure(25);
set(h,'DoubleBuffer','on');
clf;
hold on;

rot = rot(q(1));

for i=1:size(axle,2)
  Axle(:,i) = hip + rot*axle(:,i);
end
for i=1:size(leg,2)
  LLeg(:,i) = hip + rot*leg(:,i);
  RLeg(:,i) = hip + rot*[-leg(1,i);leg(2,i)];
end
fill(Axle(1,:),Axle(2,:),alumc);
fill(LLeg(1,:),LLeg(2,:),alumc);
fill(RLeg(1,:),RLeg(2,:),alumc);

Lankle = hip + rot*ankle;
Rankle = hip + rot*[-ankle(1);ankle(2)];
for i=1:size(foot,2)
  RFoot(:,i) = Rankle + rot*[-foot(1,i);foot(2,i)];
  LFoot(:,i) = Lankle + rot*foot(:,i);
end
fill(RFoot(1,:),RFoot(2,:),woodc);
fill(LFoot(1,:),LFoot(2,:),woodc);

line(q(1)+[-7,7],[0,0],'Color',[0 0 0]);

axis equal; 
axis([-7,7,-1,11])
%axis off;

title(['t = ', num2str(t)]);
drawnow;

