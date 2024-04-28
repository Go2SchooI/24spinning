%% timeline init

dt = 0.005;
N = 4000;
t = 0:dt:dt*N-dt;

file_path = 'D:/RoboMaster2024/spinning_data/spinning2.csv';
in = readtable(file_path);

deltaT1 = table2array(in(3:N+2,1));
deltaT2 = table2array(in(4:N+3,1));
deltaT = (deltaT2 - deltaT1)*1e-9;

tgttheta = pi/2 - table2array(in(3:N+2,5));
framey = table2array(in(3:N+2,2));
framex = table2array(in(3:N+2,3));
%% target init
center_x = zeros(1,N);
center_y = zeros(1,N);
center_angle = zeros(1,N);
center_x_dot = zeros(1,N);
center_y_dot = zeros(1,N);
height0 = zeros(1,N);
height1 = zeros(1,N);
height_dot = zeros(1,N);
theta = zeros(1,N);
theta_dot = zeros(1,N);
theta_predict = zeros(1,N);
x_predict = zeros(1,N);
x_predict1 = zeros(1,N);
x_center_pre = zeros(1,N);
y_center_pre = zeros(1,N);
x_armor = zeros(1,N);
theta_measure = zeros(1,N);
r0 = 0.23;
r1 = 0.18+0.05*0;

pos(1,1) = framex(1);
pos(2,1) = framey(1);
theta_measure(1) = std_rad(tgttheta(1));

%% ekf init
F = diag([1,1,1,1,1,1,1,1]);
% F(1,2) = dt;
% F(3,4) = dt;
% F(5,6) = dt;

% Pinit = diag([1,1,1,1,1,1,1,1])*0.003;
Pinit = [0.001, 0, 0, 0, 0, 0, 0, 0;
    0, 0.001, 0, 0, 0, 0, 0, 0;
    0, 0, 0.001, 0, 0, 0, 0, 0;
    0, 0, 0, 0.001, 0, 0, 0, 0;
    0, 0, 0, 0, 0.01, 0, 0, 0;
    0, 0, 0, 0, 0, 0.01, 0, 0;
    0, 0, 0, 0, 0, 0, 0.0003, 0;
    0, 0, 0, 0, 0, 0, 0, 0.0003];
P = Pinit;

process_noise = [0.001, 0.001, 0.001, 0.000001];
Q = zeros(8,8);

sigmaSqY = 0.003;%0.005
sigmaSqTheta = 0.00025;
sigmaSqYaw = 0.00025;%0.0025

R = diag([1,1,1]);

H = zeros(3,8);
H(1,1) = 1;
H(2,3) = 1;
H(3,5) = 1;

xhat = zeros(8,N);
xhatminus = zeros(8,N);
z = zeros(3,N);
chisquare = zeros(1,N);
pos_est = zeros(3,N);
theta_est = zeros(1,N);
r_est = zeros(1,N);

r_init = 0.25;
r_buf = r_init;
xhat(:,1) = [pos(1,1),0,pos(2,1),0,theta_measure(1),0,r_init,r_init].';
pos_est(:,1) = [xhat(1,1) - xhat(7,1) * cos(xhat(5,1)),...
    xhat(3,1) - xhat(7,1) * sin(xhat(5,1)),...
    0];
theta_est(1) = std_rad(xhat(5,1));
switch_count = 0;
%%
for k = 2:N
    % 量测噪声更新
    dt = deltaT(k-1);
    F(1,2) = dt;
    F(3,4) = dt;
    F(5,6) = dt;
    
    Q(1,1) = dt * dt * dt / 3 * process_noise(1);
    Q(1,2) = dt * dt / 2  * process_noise(1);
    Q(2,1) = dt * dt / 2  * process_noise(1);
    Q(2,2) = dt * process_noise(1);
    Q(3,3) = dt * dt * dt / 3 * process_noise(2);
    Q(3,4) = dt * dt / 2  * process_noise(2);
    Q(4,3) = dt * dt / 2  * process_noise(2);
    Q(4,4) = dt * process_noise(2);
    Q(5,5) = dt * dt * dt / 3 * process_noise(3);
    Q(5,6) = dt * dt / 2  * process_noise(3);
    Q(6,5) = dt * dt / 2  * process_noise(3);
    Q(6,6) = dt * process_noise(3);
    Q(7,7) = dt * process_noise(4);
    Q(8,8) = dt * process_noise(4);
    
    pos(1,k) = framex(k);
    pos(2,k) = framey(k);
    theta_measure(k) = std_rad(tgttheta(k));

    xhatminus(:,k) = F*xhat(:,k-1);
    Pminus= F*P*F'+Q;
    xhatminus(5,k) = std_rad(xhatminus(5,k));

    chisquare(k) = std_rad(theta_measure(k) - theta_measure(k-1))^2;
    if chisquare(k) > 1.75
        switch_count = switch_count + 1;
    end
    % theta_z represent the angle in z vector
    theta_z = angle_process(theta_measure(k), xhat(5,k-1));
    z(:,k) = [pos(1,k), pos(2,k), theta_z];
    c_x_ = xhatminus(1,k);
    c_y_ = xhatminus(3,k);
    % theta_ represent the theta prior estimation
    theta_ = xhatminus(5,k);
    % theta_ represent the theta prior estimation match the armor
    theta__ = angle_process(xhatminus(5,k), theta_measure(k));

    y = pos(2,k);
    theta_temp = atan2(pos(1,k), pos(2,k));
    tantheta = pos(1,k) / pos(2,k);
    costheta = cos(theta_temp);
    Rc = zeros(2,2);
    Rc(1,1) = sigmaSqY * tantheta * tantheta + sigmaSqTheta * y * y / costheta^4;
    Rc(1,2) = sigmaSqY * tantheta;
    Rc(2,1) = Rc(1,2);
    Rc(2,2) = sigmaSqY;
    R = [Rc,zeros(2,1); zeros(1,2),sigmaSqYaw];
%     R = diag([1,1,1,1]);

    H = zeros(3,8);
    H(1,1) = 1;
    H(2,3) = 1;
    H(3,5) = 1;
    if mod(switch_count,2) == 0
        r_ = xhatminus(7,k);
        H(1,5) = r_ * sin(theta__);
        H(2,5) = -r_ * cos(theta__);
        H(1,7) = -cos(theta__);
        H(2,7) = -sin(theta__);
    else
        r_ = xhatminus(8,k);
        H(1,5) = r_ * sin(theta__);
        H(2,5) = -r_ * cos(theta__);
        H(1,8) = -cos(theta__);
        H(2,8) = -sin(theta__);
    end

    hx = [c_x_ - r_ * cos(theta__);
        c_y_ - r_ * sin(theta__);
        theta_];
    
    err = z(:,k)-hx;
    err(3) = std_rad(err(3));

    K = Pminus*H'/( H*Pminus*H'+ R);
    xhat(:,k) = xhatminus(:,k)+K*err;
    P = (eye(8)-K*H)*Pminus;

    xhat(5,k) = std_rad(xhat(5,k));

    theta_est(k) = std_rad(xhat(5,k));
    if mod(switch_count,2) == 0
        r_est(k) = xhat(7,k);
    else
        r_est(k) = xhat(8,k);
    end
    
    pos_est(:,k) = [xhat(1,k) - r_est(k) * cos(angle_process(theta_est(k), theta_measure(k))),...
            xhat(3,k) - r_est(k) * sin(angle_process(theta_est(k), theta_measure(k))),...
            0];
    
    % 预测效果 
    forwardTime = 0.4;
    theta_predict(k) = xhat(5,k) + forwardTime*xhat(6,k);
    theta_predict(k) = std_rad(theta_predict(k));
    theta_predict(k) = angle_process(theta_predict(k),theta_measure(k));
    
    x_armor(k) = xhat(1,k) - r_est(k)*cos(angle_process(theta_est(k),theta_measure(k)));
    x_center_pre(k) = xhat(1,k) + forwardTime*xhat(2,k);
    y_center_pre(k) = xhat(3,k) + forwardTime*xhat(4,k);
    x_predict(k) = xhat(1,k) + forwardTime*xhat(2,k) - r_est(k) * cos(theta_predict(k));
    x_predict1(k) = xhat(1,k) - r_est(k) * cos(theta_predict(k));
end

%%
figure();
subplot(3,3,1)
plot(t,xhat(1,:),t,x_center_pre)
title('center x')
subplot(3,3,2)
plot(t,xhat(3,:),t,y_center_pre)
title('center y')
subplot(3,3,3)
plot(t,xhat(5,:))
title('theta')
subplot(3,3,4)
plot(t,xhat(2,:))
title('center xdot')
subplot(3,3,5)
plot(t,xhat(4,:))
title('center ydot')
subplot(3,3,6)
plot(t,xhat(6,:))
title('theta dot')
subplot(3,3,9)
plot(t,r_est)
title('r')
figure();
for i = 1:2
    subplot(2,1,i)
    plot(t,pos(i,:),t,pos_est(i,:))
end
figure();
plot(t,theta_measure, t,xhat(5,:))
figure();
plot(t, chisquare)
figure();
plot(t,framex)
title('x')
figure();
plot(t,framey)
title('y')

function ang = angle_process(input, theta)
if abs(std_rad(input - theta)) < pi/4
    ang = input;
elseif abs(std_rad(input + pi/2 - theta)) < pi/4
    ang = input + pi/2;
elseif abs(std_rad(input + pi - theta)) < pi/4
    ang = input + pi;
else
    ang = input - pi/2;
end
end

% function ang = angle_process(input, theta)
% if abs(std_rad(input - theta)) < pi/3
%     ang = input;
% elseif abs(std_rad(input + pi/3*2 - theta)) < pi/3
%     ang = input + pi/3*2;
% else
%     ang = input - pi/3*2;
% end
% end

function ang = std_rad(input)
if input > pi
    input = input - 2*pi;
end
if input < -pi
    input = input + 2*pi;
end
ang = input;
end
