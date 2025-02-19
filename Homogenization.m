%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Homogenization Methods     %
%        1. Voigt average         %
%        2. Reuss average         %
%        3. Figure                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

% Parameters
E_Al = 70; % GPa
E_SiC = 450; % GPa
nu_Al = 0.3;
nu_SiC = 0.17;
rho_Al = 2.7;
rho_SiC = 3.2;

ParticleVolFrac = linspace(0, 1, 101); % Predefined volume fraction range


% Pre-allocate arrays
CV_M3333 = zeros(1, 101);
EV_M = zeros(1, 101);
CR_M3333 = zeros(1, 101);
ER_M = zeros(1, 101);

% 1. Voigt Average
% Stiffness tensor C
C_Al = CalculateMatrixC(E_Al, nu_Al);
C_SiC = CalculateMatrixC(E_SiC, nu_SiC);

for i = 1:101
    Vol_SiC = ParticleVolFrac(i);
    Vol_Al = 1 - Vol_SiC;
    CV_M = Vol_Al * C_Al + Vol_SiC * C_SiC;
    CV_M3333(i) = CV_M(3,3);
    EV_M(i) = CalculateEM(CV_M);
end

% 2. Reuss Average
% Compliance tensor S
I = eye(6);
S_Al = C_Al\I;
S_SiC = C_SiC\I;

for j = 1:101
    Vol_SiC = ParticleVolFrac(j);
    Vol_Al = 1 - Vol_SiC;
    SR_M = Vol_Al * S_Al + Vol_SiC * S_SiC;
    CR_M = SR_M\I;
    CR_M3333(j) = CR_M(3,3);
    ER_M(j) = CalculateEM(CR_M);
end

% 3. Density Calculation
rho_M = rho_Al * (1 - ParticleVolFrac) + rho_SiC * ParticleVolFrac;

% 4. Figure - Stiffness Tensor
figure;
plot(ParticleVolFrac * 100, CV_M3333, '-o', ParticleVolFrac * 100, CR_M3333, '-*')
title('Compare Stiffness between Homogenization Methods')
xlabel('Particle Volume Fraction (%)')
ylabel('C_M_3_3_3_3 (GPa)')
legend('Voigt average', 'Reuss average', 'Location', 'northwest')

% 5. Figure - Effective Modulus
figure;
plot(ParticleVolFrac * 100, EV_M ./ rho_M, '-o', ParticleVolFrac * 100, ER_M ./ rho_M, '-*')
title('Compare Effective Modulus between Homogenization Methods')
xlabel('Particle Volume Fraction (%)')
ylabel('E_M / \rho_M')
legend('Voigt average', 'Reuss average', 'Location', 'northwest')

%% Function to calculate stiffness matrix C
function C = CalculateMatrixC(E, nu)
    I2 = eye(3); % Second-order identity tensor (3x3)
    I4 = eye(6); % Fourth-order identity tensor (6x6) in Voigt notation

    % Dyadic product I ⊗ I in Voigt notation
    % II_Voigt = I2 ⊗ I2
      II_Voigt = [1 1 1 0 0 0;
                  1 1 1 0 0 0;
                  1 1 1 0 0 0;
                  0 0 0 0 0 0;
                  0 0 0 0 0 0;
                  0 0 0 0 0 0]; 

    % Compute stiffness tensor using the given formula
    C = (E * nu / ((1 + nu) * (1 - 2 * nu))) * II_Voigt ...
      + (E / (1 + nu)) * I4;
end

%% Function to calculate EM
function E_M = CalculateEM(C_M)
    A = C_M(1,1);
    B = C_M(1,2);
    E_M = ((A - B) * (A + 2 * B)) / (A + B);
end
