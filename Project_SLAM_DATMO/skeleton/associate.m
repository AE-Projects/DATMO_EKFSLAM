% This function calculates the odometry information.
% Inputs:
%           mu_bar(t):          (2(Nt)+3)X1
%           sigma_bar(t):       (2(Nt)+3)X(2(Nt)+3)
%           Q(t):               2X2
%           z_i(t):             2X1
%           Nt:                 1X1

% Outputs:
%           psi(t):           2X2X(Nt)
%           pi(t):            (Nt)X1
%           H:                 2X(2Nt+3)XNt
function [psi, pi, H, z]= associate(mu_bar,sigma_bar,Q,z_i,Nt)
Fx_k = zeros(5,2*Nt+3);
Fx_k(1,1) = 1;Fx_k(2,2) = 1;Fx_k(3,3) = 1;
psi = zeros(2,2,Nt);
H = zeros(2,2*Nt+3,Nt);
pi = zeros(Nt,1);
z = zeros(2,Nt);
for k = 1:Nt
    delta_k = [mu_bar(2*k+2)- mu_bar(1);
        mu_bar(2*k+3)- mu_bar(2)];
    qk = delta_k'*delta_k;
    z_k = Measurement_model(delta_k,qk,mu_bar(3));
    z(:,k) = z_k;
    Fx_k(4:5,2*k+2:2*k+3) = eye(2,2);

    H(:,:,k) = Jacobian_H(delta_k,qk)*Fx_k;
    psi(:,:,k) = H(:,:,k)*sigma_bar*H(:,:,k)' + Q;
    pi(k,1) = (z_i - z_k)'*inv(psi(:,:,k))*(z_i - z_k);
end