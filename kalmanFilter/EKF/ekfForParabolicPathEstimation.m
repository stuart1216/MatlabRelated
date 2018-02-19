%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 % Date: 2016.10.12  
 %
 % Extended Kalman Filter to estimate parabolic path  
 %
 % the system coefficient need to adjustment carefully
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
kx = .01; ky = .05;     % damping coefficient 
g = 9.8;                % gravity 
t = 10;                 % simulation time 
Ts = 0.1;               % sample rates  
len = fix(t/Ts);        % simulation step 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% simulated orbit  
dax = 1.5; day = 1.5;  % system noise  
X = zeros(len,4); X(1,:) = [0, 50, 500, 0]; % initial value 
for k=2:len  
    x = X(k-1,1); vx = X(k-1,2); y = X(k-1,3); vy = X(k-1,4);   
    x = x + vx*Ts;  
    vx = vx + (-kx*vx^2+dax*randn(1,1))*Ts;  
    y = y + vy*Ts;  
    vy = vy + (ky*vy^2-g+day*randn(1))*Ts;  
    X(k,:) = [x, vx, y, vy];  
end  
figure(1), hold off, plot(X(:,1),X(:,3),'-b'), grid on  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% measurement
mrad = 0.001;  
dr = 10; dafa = 10*mrad; % noise  
for k=1:len  
    r = sqrt(X(k,1)^2+X(k,3)^2) + dr*randn(1,1);  
    a = atan(X(k,1)/X(k,3)) + dafa*randn(1,1);  
    Z(k,:) = [r, a];  
end  
figure(1), hold on, plot(Z(:,1).*sin(Z(:,2)), Z(:,1).*cos(Z(:,2)),'*')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% ekf 
Qk = diag([0; dax; 0; day])^2;  
Rk = diag([dr; dafa])^2;  
Xk = zeros(4,1);  
Pk = 100*eye(4);  
X_est = X;  
for k=1:len  
    Ft = JacobianF(X(k,:), kx, ky, g);  
    Hk = JacobianH(X(k,:));  
    fX = fff(X(k,:), kx, ky, g, Ts);  
    hfX = hhh(fX, Ts);  
    [Xk, Pk, Kk] = ekf(eye(4)+Ft*Ts, Qk, fX, Pk, Hk, Rk, Z(k,:)'-hfX);  
    X_est(k,:) = Xk';  
end  
figure(1), plot(X_est(:,1),X_est(:,3), '+r')  
xlabel('X'); ylabel('Y'); title('ekf simulation');  
legend('real', 'measurement', 'ekf estimated');  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 