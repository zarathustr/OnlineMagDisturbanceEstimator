% Real-Time Magnetometer Disturbance 
%            Estimation via Online Nonlinear Programming
%
% author: Jin Wu
% e-mail: jin_wu_uestc@hotmail.com
% Reference: ﻿Wu, J. (2019). Real-time Magnetometer Disturbance 
%                            Estimation via Online Nonlinear Programming. 
%                            IEEE Sensors Journal, 19(12), 4405–4411. 
%                            https://doi.org/10.1109/JSEN.2019.2901925﻿


clear all
close all
clc

warning('off');

data = load('magnetic_distortion4.txt');
len = length(data(:, 1));
conv = diag([-1 1 -1]);
euler_true = data(1 : len, 1 : 3) * conv;
acc = data(1 : len, 4 : 6) * conv;
gyro = data(1 : len, 7 : 9) * conv; 
mag = data(1 : len, 10 : 12) * conv;


dt = 1 / 500;
time = dt * (1 : len);
Q1 = zeros(len, 4);
Q2 = zeros(len, 4);
q = angle2quat(- euler_true(1, 1) * pi / 180, ...
               - euler_true(1, 2) * pi / 180, ...
               euler_true(1, 3) * pi / 180 + pi, 'XYZ');
qw = q';

quaternion_w = zeros(len, 4);  % gyro-integrated quaternion
quaternion_new = zeros(len, 4); % quaternion after optimization
quaternion_prev = zeros(len, 4); % quaternion before optimization
quaternion_diff = zeros(len, 4); % difference between ground truth and estimates
norm_m = zeros(len, 1);
bias = zeros(len, 3); % estimated magnetometer disturbance
mag_prev = zeros(len, 3);
mag_new = zeros(len, 3);
gradient = zeros(len, 4); % gradient of optimization objective
eig_hess = zeros(len, 4); % eigenvalues of the hessian
iter = zeros(len, 1); % iteration numbers
fval = zeros(len, 1); % objective function value
mNs = zeros(len, 1); % estimated mNs
mDs = zeros(len, 1); % estimated mDs
states = zeros(len, 4); % the states of the optimization
global ax ay az mx my mz
global q0 q1 q2 q3

last_bmx = 0;
last_bmy = 0;
last_bmz = 0;
last_mN = 0;
last_mD = 0;

for i = 1 : len
    wx = gyro(i, 1);        wy = gyro(i, 2);        wz = gyro(i, 3);
    
    if(norm(gyro(i, :)) < 0.02) % Zero Angular Rate Update (ZARU)
        wx = 0;
        wy = 0;
        wz = 0;
    end
    
    omega4 = [  0, -wx, -wy, -wz;
               wx,   0,  wz, -wy;
               wy, -wz,   0,  wx;
               wz,  wy, -wx,   0];
    Phi = 0.5 * dt * omega4;
    qw_ = qw;
    qw = qw + Phi * qw;
    qw = qw ./ norm(qw);
    quaternion_w(i, :) = qw';
    
    C = quat2dcm(qw');    
    Ar = [0; 0; 1];
    acc_ = C * Ar;
    acc_ = acc_ ./ norm(acc_);
    
    B = 0.5 * acc_ * Ar';
    Ka = B2K(B);
    
    norm_m(i) = norm(mag(i, :));
    mag_ = mag(i, :)';
    mag_ = mag_ ./ norm(mag_);
    mag_prev(i, :) = mag_';
    mD = dot(acc_, mag_);
    Mr = [sqrt(1 - mD * mD); 0; mD];
    Bm = 0.5 * mag_ * Mr';
    Km = B2K(Bm);
    Kam = Ka + Km;
    [V, D] = eig(Kam);
    quaternion_prev(i, :) = V(:, 4)';
    quaternion_prev(i, :) = quaternion_prev(i, :) ./ norm(quaternion_prev(i, :));
    
    if(i > 1)
        ax = acc_(1); ay = acc_(2); az = acc_(3);
        mx = mag_(1); my = mag_(2); mz = mag_(3);
        qww = qw;
        q0 = qww(1); q1 = qww(2); q2 = qww(3); q3 = qww(4);
        
           
        x0 = [last_bmx, last_bmy, last_bmz, 1];
        if(abs(mD_f(x0)) >= 1)
            x0 = [0, 0, 0, 1];
        end
            
        opt = optimoptions(...
                           'fmincon', 'Display', 'off', ...
                           'MaxIterations', 20000, ...
                           'MaxFunctionEvaluations', 20000, ...
                           'FiniteDifferenceType', 'central', ...
                           'FiniteDifferenceStepSize', 1e-9, ...
                           'Algorithm', 'interior-point', ...
                           'ConstraintTolerance', 1e-100, ...
                           'OptimalityTolerance', 1e-100, ...
                           'HessianApproximation', 'bfgs', ...
                           'UseParallel', false);
            
            
        [x1, f1, ~, output, lamb, grad, hess] = fmincon(@J, x0, ...
            [], [], [], [], [], [], @quat_con, opt);
        
        [V, D] = eig(hess);
        gradient(i, :) = grad';
        eig_hess(i, :) = sort(diag(D))';
        iter(i) = output.iterations;
        fval(i) = f1;
        mNs(i) = mN_f(x1);
        mDs(i) = mD_f(x1);
        last_mD = mDs(i);
        last_mN = mNs(i);
        states(i, :) = x1;
        last_bmx = x1(1);
        last_bmy = x1(2);
        last_bmz = x1(3);
        
        bias(i, :) = [x1(1), x1(2), x1(3)];
       
        x = bias(i, :);
        mag__ = mag_ + x';
        
        norm(mag__);
        mag__ = mag__ ./ norm(mag__);
        mag_new(i, :) = mag__';
        qam = SAAM(acc_, mag__);
        quaternion_new(i, :) = qam';
        qam = qam';
        qwww = qww';
        if(f1 > 1e-15)
            qam
            qwww
            f1
        end
        
        if(mod(i, 100) == 0)
            i
            qam
            qwww
        end
        
        mD = dot(mag__, acc_);
        mN = sqrt(1 - mD * mD);
        
        if(sign(qam(1)) == sign(qwww(1)))
            quaternion_diff(i, :) = qam - qwww;
        else
            quaternion_diff(i, :) = qam + qwww;
        end
    end
    
    last_Kam = Kam;
    last_Ka = Ka;
    last_Km = Km;
    last_mag = mag_;
end