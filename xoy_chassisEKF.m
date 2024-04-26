%% timeline init

dt = 0.005;
N = 3000;
t = 0:dt:dt*N-dt;

deltaT1 = table2array(in(1:N,1));
deltaT2 = table2array(in(2:N+1,1));
deltaT = (deltaT2 - deltaT1)*1e-9;

tgttheta = table2array(in(1:N,4));
framex = table2array(in(1:N,3));
framey = -table2array(in(1:N,2));
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
x_armor = zeros(1,N);
r0 = 0.23;
r1 = 0.18+0.05*0;

center_x(1) = 0;
center_y(1) = 5;
height0(1) = 0.1;
height1(1) = 0.15;
theta(1) = pi/4*1.01;
theta_dot(1) = 0.5;
pos = zeros(3,N);
center_angle(1) = atan2(center_y(1), center_x(1));

height = zeros(1,N);
r = zeros(1,N);
theta_measure = zeros(1,N);

if abs(std_rad(theta(1)-center_angle(1))) < pi/4
    theta_measure(1) = theta(1);
    height(1) = height0(1);
    r(1) = r0;
    pos(:,1) = [center_x(1) - r0 * cos(theta(1)),...
        center_y(1) - r0 * sin(theta(1)),...
        height0(1)];
elseif abs(std_rad(theta(1)+pi/2-center_angle(1))) < pi/4
    theta_measure(1) = theta(1)+pi/2;
    height(1) = height1(1);
    r(1) = r1;
    pos(:,1) = [center_x(1) - r1 * cos(theta(1)+pi/2),...
        center_y(1) - r1 * sin(theta(1)+pi/2),...
        height1(1)];
elseif abs(std_rad(theta(1)+pi-center_angle(1))) < pi/4
    theta_measure(1) = theta(1)+pi;
    height(1) = height0(1);
    r(1) = r0;
    pos(:,1) = [center_x(1) - r0 * cos(theta(1)+pi),...
        center_y(1) - r0 * sin(theta(1)+pi),...
        height0(1)];
else
    height(1) = height1(1);
    r(1) = r1;
    theta_measure(1) = theta(1)-pi/2;
    pos(:,1) = [center_x(1) - r1 * cos(theta(1)-pi/2),...
        center_y(1) - r1 * sin(theta(1)-pi/2),...
        height1(1)];
end
% theta_measure(1) = theta(1);
% r(1) = r0;
% pos(:,1) = [center_x(1) - r0 * cos(theta_measure(1)),...
%     center_y(1) - r0 * sin(theta_measure(1)),...
%     height0(1)];
    pos(1,1) = framex(1);
    pos(2,1) = framey(1);
    theta_measure(1) = std_rad(tgttheta(1));

%% ekf init
F = diag([1,1,1,1,1,1,1,1]);
% F(1,2) = dt;
% F(3,4) = dt;
% F(5,6) = dt;

Pinit = diag([1,1,1,1,1,1,1,1])*0.1;
P = Pinit;

process_noise = [0.00001, 0.00001, 0.00005, 0.000001];
Q = zeros(8,8);
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
sigmaSqY = 0.005;
sigmaSqTheta = 0.0001;
sigmaSqYaw = 0.0025;

R = diag([1,1,1]);

% P = diag([1,1,1,1,1,1,1,1,1]);
% Q = diag([1e-2,5e-2,1e-2,5e-2,1e-2,1e-4,2e-2,4e-2,1e-3]);
% R = diag([0.1,0.1,0.1,0.2]);

% P = diag([1,1,1,1,1,1,1,1,1])*0.01;
% Q = diag([0.1,100,0.1,100,0.1,0.0001,0.2,200,0.001]);
% R = diag([50,50,10,100]);

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
pos_est(:,1) = [xhat(1,1) - xhat(7,1) * cos(xhat(7,1)),...
    xhat(3,1) - xhat(7,1) * sin(xhat(7,1)),...
    0];
theta_est(1) = std_rad(xhat(5,1));
switch_count = 0;
%%
for k = 2:N
    % 量测噪声更新
    F(1,2) = dt;
    F(3,4) = dt;
    F(5,6) = dt;
    
    
    theta_dot(k) = theta_dot(1);
    center_x_dot(k) = 3*cos(0.3*pi*k*dt)*0-1*0;
    center_y_dot(k) = -3*sin(0.2*pi*k*dt)*0-1*0;
    center_x(k) = center_x(k-1) + center_x_dot(k)*dt;
    center_y(k) = center_y(k-1) + center_y_dot(k)*dt;
    height0(k) = height0(k-1) + height_dot(k)*dt;
    height1(k) = height1(k-1) + height_dot(k)*dt;
    theta(k) = theta(k-1) + theta_dot(k)*dt;
    theta(k) = std_rad(theta(k));
    center_angle(k) = atan2(center_y(k), center_x(k));

    % theta_measure(k) represent the angle from pnp
    if abs(std_rad(theta(k)-center_angle(k))) < pi/4
        theta_measure(k) = std_rad(theta(k));
        height(k) = height0(k);
        r(k) = r0;
        pos(:,k) = [center_x(k) - r0 * cos(theta_measure(k)),...
            center_y(k) - r0 * sin(theta_measure(k)),...
            height0(k)];
    elseif abs(std_rad(theta(k)+pi/2-center_angle(k))) < pi/4
        theta_measure(k) = std_rad(theta(k)+pi/2);
        height(k) = height1(k);
        r(k) = r1;
        pos(:,k) = [center_x(k) - r1 * cos(theta_measure(k)),...
            center_y(k) - r1 * sin(theta_measure(k)),...
            height1(k)];
    elseif abs(std_rad(theta(k)+pi-center_angle(k))) < pi/4
        theta_measure(k) = std_rad(theta(k)+pi);
        height(k) = height0(k);
        r(k) = r0;
        pos(:,k) = [center_x(k) - r0 * cos(theta_measure(k)),...
            center_y(k) - r0 * sin(theta_measure(k)),...
            height0(k)];
    else
        height(k) = height1(k);
        r(k) = r1;
        theta_measure(k) = std_rad(theta(k)-pi/2);
        pos(:,k) = [center_x(k) - r1 * cos(theta(k)-pi/2),...
            center_y(k) - r1 * sin(theta(k)-pi/2),...
            height1(k)];
    end
    pos(:,k) = pos(:,k) + normrnd(0,0.01,[3,1]);
    theta_measure(k) = theta_measure(k) + normrnd(0,0.01,1);

%     height(k) = height0(k);
%     r(k) = r0;
%     pos(:,k) = [center_x(k) - r0 * cos(theta(k)),...
%         center_y(k) - r0 * sin(theta(k)),...
%         height0(k)];
%     theta_measure(k) = theta(k);

    pos(1,k) = framex(k);
    pos(2,k) = framey(k);
    theta_measure(k) = std_rad(tgttheta(k));

    xhatminus(:,k) = F*xhat(:,k-1);
    Pminus= F*P*F'+Q;
    xhatminus(5,k) = std_rad(xhatminus(5,k));

    chisquare(k) = std_rad(theta_measure(k) - theta_measure(k-1))^2;
    if chisquare(k) > 1.35
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
    pos_est(:,k) = [xhatminus(1,k) - xhatminus(7,k) * cos(theta__),...
        xhatminus(3,k) - xhatminus(7,k) * sin(theta__),...
        0];

%     plot(pos(1,1:k),pos(2,1:k), pos_est(1,1:k),pos_est(2,1:k))
%     pause(0.001);
    
    err = z(:,k)-hx;
    err(3) = std_rad(err(3));

    K = Pminus*H'/( H*Pminus*H'+ R);
    xhat(:,k) = xhatminus(:,k)+K*err;
    P = (eye(8)-K*H)*Pminus;

    xhat(5,k) = std_rad(xhat(5,k));

%     if xhat(9,k) < 0.1
%         xhat(9,k) = 0.1;
%     end
%     pos_est(:,k) = [xhat(1,k) - xhat(9,k) * cos(xhat(7,k)),...
%         xhat(3,k) - xhat(9,k) * sin(xhat(7,k)),...
%         xhat(5,k)];
    theta_est(k) = std_rad(xhat(5,k));
    if mod(switch_count,2) == 0
        r_est(k) = xhat(7,k);
    else
        r_est(k) = xhat(8,k);
    end

    forwardTime = 0.0;
    theta_predict(k) = xhat(5,k) + forwardTime*xhat(6,k);
    theta_predict(k) = std_rad(theta_predict(k));
    theta_predict(k) = angle_process(theta_predict(k),theta_measure(k));
    x_armor(k) = xhat(1,k) - xhat(7,k)*cos(angle_process(theta_est(k),theta_measure(k)));
    x_predict(k) = xhat(1,k) + forwardTime*xhat(2,k) - r_est(k) * cos(theta_predict(k));
    x_predict1(k) = xhat(1,k) - r_est(k) * cos(theta_predict(k));
%     plot xoy trojectory
%     figure(1)
%     subplot(2,1,1)
%     plot(pos(1,1:k),pos(2,1:k),center_x(1:k),center_y(1:k))
%     grid on;
%     xlim([-5,5])
%     ylim([-5,5])
%     subplot(2,1,2)
%     plot(t(1:k),xhat(6,1:k))
%     pause(0.001);
end

%%
figure();
subplot(3,3,1)
plot(t,xhat(1,:))
title('center x')
subplot(3,3,2)
plot(t,xhat(3,:))
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
plot(t,theta_measure)
figure();
plot(t,framex)
figure();
plot(t,framey)

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
