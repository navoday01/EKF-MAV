function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%covarPrev and uPrev are the previous mean and covariance respectively
%angVel is the angular velocity
%acc is the acceleration
%dt is the sampling time

x3 = [uPrev(7,1); uPrev(8,1); uPrev(9,1)]; % linear velocity
  
x2 = [uPrev(4,1); uPrev(5,1); uPrev(6,1)]; % orientation

% Euler rates

G = [ cos(x2(3))*cos(x2(2)), -sin(x2(3)),   0
      sin(x2(3))*cos(x2(2)),  cos(x2(3)),   0
      -sin(x2(2)),               0,         1];

% Rotation Matrix

R = [cos(x2(2))*cos(x2(3)), cos(x2(3))*sin(x2(1))*sin(x2(2)) - cos(x2(1))*sin(x2(3)), sin(x2(1))*sin(x2(3)) + cos(x2(1))*cos(x2(3))*sin(x2(2));
     cos(x2(2))*sin(x2(3)), cos(x2(1))*cos(x2(3)) + sin(x2(1))*sin(x2(2))*sin(x2(3)), cos(x2(1))*sin(x2(2))*sin(x2(3)) - cos(x2(3))*sin(x2(1));
             -sin(x2(2)),                              cos(x2(2))*sin(x2(1)),                              cos(x2(1))*cos(x2(2))];
 

g = [0;  0; -9.81]; % gravity

x4 = [uPrev(10,1); uPrev(11,1); uPrev(12,1)]; % gyroscope bias

x5 = [uPrev(13,1); uPrev(14,1); uPrev(15,1)]; % accelerometer bias

wm = [angVel(1,1); angVel(2,1); angVel(3,1)]; % angular velocity

am = [acc(1,1); acc(2,1); acc(3,1)];   % acceleration

ng = [0;0;0]; % gyroscope noise
na = [0;0;0]; % accelerometer noise
nbg = [0;0;0]; % derivative of gyroscope bias
nba = [0;0;0]; % derivative of accelerometer bias

% Process Model
xdot = [x3; 
       inv(G)*R*(wm - x4 -ng);
       g+R*(am-x5-na);
       nbg;
       nba];


% Partial derivative of xdot w.r.t. x

At =[0, 0, 0,                                                                                                                                         0,                                                                                                                                 0,                                                                                                                                                                             0, 1, 0, 0,  0,                             0,                             0,                  0,                                                0,                                                0;
     0, 0, 0,                                                                                                                                         0,                                                                                                                                 0,                                                                                                                                                                             0, 0, 1, 0,  0,                             0,                             0,                  0,                                                0,                                                0;
     0, 0, 0,                                                                                                                                         0,                                                                                                                                 0,                                                                                                                                                                             0, 0, 0, 1,  0,                             0,                             0,                  0,                                                0,                                                0;
     0, 0, 0,                            -(sin(uPrev(5,1))*(uPrev(11,1)*cos(uPrev(4,1)) + 0*cos(uPrev(4,1)) - uPrev(12,1)*sin(uPrev(4,1)) - angVel(2,1)*cos(uPrev(4,1)) - 0*sin(uPrev(4,1)) + angVel(3,1)*sin(uPrev(4,1))))/cos(uPrev(5,1)),                             -(uPrev(12,1)*cos(uPrev(4,1)) + 0*cos(uPrev(4,1)) + uPrev(11,1)*sin(uPrev(4,1)) - angVel(3,1)*cos(uPrev(4,1)) + 0*sin(uPrev(4,1)) - angVel(2,1)*sin(uPrev(4,1)))/cos(uPrev(5,1))^2,                                                                                                                                                                             0, 0, 0, 0, -1, -(sin(uPrev(4,1))*sin(uPrev(5,1)))/cos(uPrev(5,1)), -(cos(uPrev(4,1))*sin(uPrev(5,1)))/cos(uPrev(5,1)),                  0,                                                0,                                                0;
     0, 0, 0,                                                   uPrev(12,1)*cos(uPrev(4,1)) + 0*cos(uPrev(4,1)) + uPrev(11,1)*sin(uPrev(4,1)) - angVel(3,1)*cos(uPrev(4,1)) + 0*sin(uPrev(4,1)) - angVel(2,1)*sin(uPrev(4,1)),                                                                                                                                 0,                                                                                                                                                                             0, 0, 0, 0,  0,                     -cos(uPrev(4,1)),                      sin(uPrev(4,1)),                  0,                                                0,                                                0;
     0, 0, 0,                                       -(uPrev(11,1)*cos(uPrev(4,1)) + 0*cos(uPrev(4,1)) - uPrev(12,1)*sin(uPrev(4,1)) - angVel(2,1)*cos(uPrev(4,1)) - 0*sin(uPrev(4,1)) + angVel(3,1)*sin(uPrev(4,1)))/cos(uPrev(5,1)),                  -(sin(uPrev(5,1))*(uPrev(12,1)*cos(uPrev(4,1)) + 0*cos(uPrev(4,1)) + uPrev(11,1)*sin(uPrev(4,1)) - angVel(3,1)*cos(uPrev(4,1)) + 0*sin(uPrev(4,1)) - angVel(2,1)*sin(uPrev(4,1))))/cos(uPrev(5,1))^2,                                                                                                                                                                             0, 0, 0, 0,  0,            -sin(uPrev(4,1))/cos(uPrev(5,1)),            -cos(uPrev(4,1))/cos(uPrev(5,1)),                  0,                                                0,                                                0;
     0, 0, 0, - (sin(uPrev(4,1))*sin(uPrev(6,1)) + cos(uPrev(4,1))*cos(uPrev(6,1))*sin(uPrev(5,1)))*(uPrev(14,1) - acc(2,1) + 0) - (cos(uPrev(4,1))*sin(uPrev(6,1)) - cos(uPrev(6,1))*sin(uPrev(4,1))*sin(uPrev(5,1)))*(uPrev(15,1) - acc(3,1) + 0), cos(uPrev(6,1))*sin(uPrev(5,1))*(uPrev(13,1) - acc(1,1) + 0) - cos(uPrev(4,1))*cos(uPrev(5,1))*cos(uPrev(6,1))*(uPrev(15,1) - acc(3,1) + 0) - cos(uPrev(5,1))*cos(uPrev(6,1))*sin(uPrev(4,1))*(uPrev(14,1) - acc(2,1) + 0), (cos(uPrev(4,1))*cos(uPrev(6,1)) + sin(uPrev(4,1))*sin(uPrev(5,1))*sin(uPrev(6,1)))*(uPrev(14,1) - acc(2,1) + 0) - (cos(uPrev(6,1))*sin(uPrev(4,1)) - cos(uPrev(4,1))*sin(uPrev(5,1))*sin(uPrev(6,1)))*(uPrev(15,1) - acc(3,1) + 0) + cos(uPrev(5,1))*sin(uPrev(6,1))*(uPrev(13,1) - acc(1,1) + 0), 0, 0, 0,  0,                             0,                             0, -cos(uPrev(5,1))*cos(uPrev(6,1)),   cos(uPrev(4,1))*sin(uPrev(6,1)) - cos(uPrev(6,1))*sin(uPrev(4,1))*sin(uPrev(5,1)), - sin(uPrev(4,1))*sin(uPrev(6,1)) - cos(uPrev(4,1))*cos(uPrev(6,1))*sin(uPrev(5,1));
     0, 0, 0,   (cos(uPrev(6,1))*sin(uPrev(4,1)) - cos(uPrev(4,1))*sin(uPrev(5,1))*sin(uPrev(6,1)))*(uPrev(14,1) - acc(2,1) + 0) + (cos(uPrev(4,1))*cos(uPrev(6,1)) + sin(uPrev(4,1))*sin(uPrev(5,1))*sin(uPrev(6,1)))*(uPrev(15,1) - acc(3,1) + 0), sin(uPrev(5,1))*sin(uPrev(6,1))*(uPrev(13,1) - acc(1,1) + 0) - cos(uPrev(4,1))*cos(uPrev(5,1))*sin(uPrev(6,1))*(uPrev(15,1) - acc(3,1) + 0) - cos(uPrev(5,1))*sin(uPrev(4,1))*sin(uPrev(6,1))*(uPrev(14,1) - acc(2,1) + 0), (cos(uPrev(4,1))*sin(uPrev(6,1)) - cos(uPrev(6,1))*sin(uPrev(4,1))*sin(uPrev(5,1)))*(uPrev(14,1) - acc(2,1) + 0) - (sin(uPrev(4,1))*sin(uPrev(6,1)) + cos(uPrev(4,1))*cos(uPrev(6,1))*sin(uPrev(5,1)))*(uPrev(15,1) - acc(3,1) + 0) - cos(uPrev(5,1))*cos(uPrev(6,1))*(uPrev(13,1) - acc(1,1) + 0), 0, 0, 0,  0,                             0,                             0, -cos(uPrev(5,1))*sin(uPrev(6,1)), - cos(uPrev(4,1))*cos(uPrev(6,1)) - sin(uPrev(4,1))*sin(uPrev(5,1))*sin(uPrev(6,1)),   cos(uPrev(6,1))*sin(uPrev(4,1)) - cos(uPrev(4,1))*sin(uPrev(5,1))*sin(uPrev(6,1));
     0, 0, 0,                                                                 cos(uPrev(5,1))*sin(uPrev(4,1))*(uPrev(15,1) - acc(3,1) + 0) - cos(uPrev(4,1))*cos(uPrev(5,1))*(uPrev(14,1) - acc(2,1) + 0),                            cos(uPrev(5,1))*(uPrev(13,1) - acc(1,1) + 0) + cos(uPrev(4,1))*sin(uPrev(5,1))*(uPrev(15,1) - acc(3,1) + 0) + sin(uPrev(4,1))*sin(uPrev(5,1))*(uPrev(14,1) - acc(2,1) + 0),                                                                                                                                                                             0, 0, 0, 0,  0,                             0,                             0,           sin(uPrev(5,1)),                               -cos(uPrev(5,1))*sin(uPrev(4,1)),                               -cos(uPrev(4,1))*cos(uPrev(5,1));
     0, 0, 0,                                                                                                                                         0,                                                                                                                                 0,                                                                                                                                                                             0, 0, 0, 0,  0,                             0,                             0,                  0,                                                0,                                                0;
     0, 0, 0,                                                                                                                                         0,                                                                                                                                 0,                                                                                                                                                                             0, 0, 0, 0,  0,                             0,                             0,                  0,                                                0,                                                0;
     0, 0, 0,                                                                                                                                         0,                                                                                                                                 0,                                                                                                                                                                             0, 0, 0, 0,  0,                             0,                             0,                  0,                                                0,                                                0;
     0, 0, 0,                                                                                                                                         0,                                                                                                                                 0,                                                                                                                                                                             0, 0, 0, 0,  0,                             0,                             0,                  0,                                                0,                                                0;
     0, 0, 0,                                                                                                                                         0,                                                                                                                                 0,                                                                                                                                                                             0, 0, 0, 0,  0,                             0,                             0,                  0,                                                0,                                                0;
     0, 0, 0,                                                                                                                                         0,                                                                                                                                 0,                                                                                                                                                                             0, 0, 0, 0,  0,                             0,                             0,                  0,                                                0,                                                0];
 
% Partial derivative of xdot w.r.t. noise

Ut =[ 0,                             0,                             0,                  0,                                                0,                                                0, 0, 0, 0, 0, 0, 0;
      0,                             0,                             0,                  0,                                                0,                                                0, 0, 0, 0, 0, 0, 0;
      0,                             0,                             0,                  0,                                                0,                                                0, 0, 0, 0, 0, 0, 0;
     -1, -(sin(uPrev(4,1))*sin(uPrev(5,1)))/cos(uPrev(5,1)), -(cos(uPrev(4,1))*sin(uPrev(5,1)))/cos(uPrev(5,1)),                  0,                                                0,                                                0, 0, 0, 0, 0, 0, 0;
      0,                     -cos(uPrev(4,1)),                      sin(uPrev(4,1)),                  0,                                                0,                                                0, 0, 0, 0, 0, 0, 0;
      0,            -sin(uPrev(4,1))/cos(uPrev(5,1)),            -cos(uPrev(4,1))/cos(uPrev(5,1)),                  0,                                                0,                                                0, 0, 0, 0, 0, 0, 0;
      0,                             0,                             0, -cos(uPrev(5,1))*cos(uPrev(6,1)),   cos(uPrev(4,1))*sin(uPrev(6,1)) - cos(uPrev(6,1))*sin(uPrev(4,1))*sin(uPrev(5,1)), - sin(uPrev(4,1))*sin(uPrev(6,1)) - cos(uPrev(4,1))*cos(uPrev(6,1))*sin(uPrev(5,1)), 0, 0, 0, 0, 0, 0;
      0,                             0,                             0, -cos(uPrev(5,1))*sin(uPrev(6,1)), - cos(uPrev(4,1))*cos(uPrev(6,1)) - sin(uPrev(4,1))*sin(uPrev(5,1))*sin(uPrev(6,1)),   cos(uPrev(6,1))*sin(uPrev(4,1)) - cos(uPrev(4,1))*sin(uPrev(5,1))*sin(uPrev(6,1)), 0, 0, 0, 0, 0, 0;
      0,                             0,                             0,           sin(uPrev(5,1)),                               -cos(uPrev(5,1))*sin(uPrev(4,1)),                               -cos(uPrev(4,1))*cos(uPrev(5,1)), 0, 0, 0, 0, 0, 0;
      0,                             0,                             0,                  0,                                                0,                                                0, 1, 0, 0, 0, 0, 0;
      0,                             0,                             0,                  0,                                                0,                                                0, 0, 1, 0, 0, 0, 0;
      0,                             0,                             0,                  0,                                                0,                                                0, 0, 0, 1, 0, 0, 0;
      0,                             0,                             0,                  0,                                                0,                                                0, 0, 0, 0, 1, 0, 0;
      0,                             0,                             0,                  0,                                                0,                                                0, 0, 0, 0, 0, 1, 0;
      0,                             0,                             0,                  0,                                                0,                                                0, 0, 0, 0, 0, 0, 1];
 
Ft = eye(15) + dt*At;

Q = eye(12); % covariance of noise

Qd = dt*Q; % covariance of noise after discretization

uEst = double(uPrev + dt*xdot);

covarEst = double(Ft*covarPrev*Ft'+ Ut*Qd*Ut');

end

