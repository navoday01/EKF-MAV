

syms px1 py2 pz3 % position
syms qx1 qy2 qz3 % orientation
syms pdotx pdoty pdotz % linear velocity
syms bgx bgy bgz % gyroscope bias
syms bax bay baz % accelerometer bias
syms ngx ngy ngz % gyroscope noise
syms nax nay naz % accelerometer noise
syms nbgx nbgy nbgz % bgdot
syms nbax nbay nbaz % badot
syms G R g wm am At Ft Ut n
syms wmx wmy wmz 
syms amx amy amz 

x3 = [pdotx; pdoty; pdotz]; % linear velocity
  
x2 = [qx1; qy2; qz3]; % orientation

G = [ cos(x2(3))*cos(x2(2)), -sin(x2(3)),   0    
      sin(x2(3))*cos(x2(2)),  cos(x2(3)),   0
      -sin(x2(2)),               0,         1]; %Euler rates

g = [0;  0; -9.81]; % gravity

%R = eul2rotm([x2(3) x2(2) x2(1)],'ZYX');
R = rotz(x2(3))*roty(x2(2))*rotx(x2(1));

wm = [wmx; wmy; wmz];

am = [amx; amy; amz];

na = [nax; nay; naz];

ng = [ngx; ngy; ngz];

nbg = [nbgx; nbgy; nbgz];

nba = [nbax; nbay; nbaz];

x4 = [bgx; bgy; bgz];

x5 = [bax; bay; baz];

xdot = [x3; 
       inv(G)*R*(wm - x4 -ng);
       g+R*(am-x5-na);
       nbg;
       nba];

I = eye(15);

x = [px1; py2; pz3; qx1; qy2; qz3; pdotx; pdoty; pdotz; bgx; bgy; bgz; bax; bay; baz];

for i = 1:15

    At(1:15,i) = diff(xdot,x(i,1));
end
%At
n = [ngx; ngy; ngz; nax; nay; naz; nbgx; nbgy; nbgz; nbax; nbay; nbaz];

for i = 1:12

    Ut(1:15,i) = diff(xdot,n(i,1));
end
%Ut
fprintf('At =')
disp(simplify(At))
fprintf('Ut =')
disp(simplify(Ut))
%{
px1 = uPrev(1,1);
py2 = uPrev(2,1);
pz3 = uPrev(3,1);
qx1 = uPrev(4,1);
qy2 = uPrev(5,1);
qz3 = uPrev(6,1);
pdotx = uPrev(7,1);
pdoty = uPrev(8,1);
pdotz = uPrev(9,1);
bgx = uPrev(10,1);
bgy = uPrev(11,1);
bgz = uPrev(12,1); 
bax = uPrev(13,1);
bay = uPrev(14,1); 
baz = uPrev(15,1);
wmx = angVel(1,1);
wmy = angVel(2,1);
wmz = angVel(3,1);
amx = acc(1,1);
amy = acc(2,1);
amz = acc(3,1);
%%
ngx = normrnd(0,1);
ngy = normrnd(0,1);
ngz = normrnd(0,1);
nax = normrnd(0,1);
nay = normrnd(0,1);
naz = normrnd(0,1);
nbgx = normrnd(0,1);
nbgy = normrnd(0,1);
nbgz = normrnd(0,1);
nbax = normrnd(0,1);
nbay = normrnd(0,1);
nbaz = normrnd(0,1);
%%

Ft = I + dt*At;

Q = diag([ngx ngy ngz nax nay naz nbgx nbgy nbgz nbax nbay nbaz]);

Qd = dt*Q;

xdot = subs(xdot);

Ft = subs(Ft);

Ut = subs(Ut);

uEst = double(uPrev + dt*xdot);

covarEst = double(Ft*covarPrev*Ft'+ Ut*Qd*Ut');
%}

 
