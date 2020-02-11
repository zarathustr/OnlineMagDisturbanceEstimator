% Super-fast Attitude of Accelerometer and Magnetometer (SAAM)
%
% author: Jin Wu
% e-mail: jin_wu_uestc@hotmail.com; klinsmann.zhou@gmail.com
% Reference: ﻿Wu, J., Zhou, Z., Fourati, H., & Cheng, Y. (2018). 
%                   A Super Fast Attitude Determination Algorithm for 
%                   Consumer-Level Accelerometer and Magnetometer. 
%                   IEEE Transactions on Consumer Electronics, 64(3), 
%                   375–381. https://doi.org/10.1109/TCE.2018.2859625


function q = SAAM(Ab, Mb)

    ax = Ab(1);       ay = Ab(2);       az = Ab(3);
    mx = Mb(1);       my = Mb(2);       mz = Mb(3);

    mD = ax * mx + ay * my + az * mz;
    mN = sqrt(1 - mD * mD);
    

    q = [
            ax * my - ay * mx - ay * mN;
          - mx + az * mx - ax * mz + (az - 1) * mN + ax * mD;
          - my + az * my - ay * mz + ay * mD;
          - mz + az * mD - ax * mN
    ];

    q = q ./ sqrt(q(1) * q(1) + q(2) * q(2) + q(3) * q(3) + q(4) * q(4));
end