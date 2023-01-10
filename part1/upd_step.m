function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%z_t is the measurement
%covarEst and uEst are the predicted covariance and mean respectively
%uCurr and covar_curr are the updated mean and covariance respectively

Ct = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];

R = eye(6)*0.0001; % covariance matrix for the noise of measurement

Kt = covarEst*Ct'/(Ct*covarEst*Ct'+R); % Kalman Gain

uCurr = uEst + Kt*(z_t-(Ct*uEst)) ;

covar_curr = covarEst - Kt*Ct*covarEst;
end 