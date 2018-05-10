% Title: Calculate the other three variables based on two knowing
%        thermodynamic variables. And plot the T-s Diagram of Oxygen.
% Based on: MATLAB program from << Chemical, Biochemical, and Engineering
%                                  Thermodynamics >>
% Version: 1.2, Edward Xu, 18.4.10
% SubTitle: Calculate the other three variables based on two knowing
%           thermodynamic variables.

% Constants ------------------------------------------------------------------
R = 8.314;                   % J/(mol*K), Universial Gas Constant【通用气体常数】
R_G = 8.314 / 0.032;         % R = kPa m^3/mol K ???
AAA = 25.46;                 % C_p1, heat capacity calculation parameter of Oxygen
BBB = 1.519E-2;              % C_p2, heat capacity calculation parameter of Oxygen
CCC = -0.7151E-5;            % C_p3, heat capacity calculation parameter of Oxygen
DDD = 1.311E-9;              % C_p4, heat capacity calculation parameter of Oxygen
T_c = 154.6;                 % K
p_c = 504.6;                 % kPa
% OMEGA = 0.176;
% T_boil = 243.4;            % ???
T_ref = 25 + 273.15;         % K
p_ref = 101.325;             % kPa

% Part1: Peng-Robinson Constant Calculation,
% Reference: << Chemical, Biochemical, and Engineering Thermodynamics, 5e >>
%            E6.4-2, P221, Ch6
a_Tc = 0.457235529 * (R_G * T_c)^2 ./ p_c;       % Critical Point Restriction "a(T_c)"
KAPPA = 0.4069;                                  % Dependent on OMEGA(working substance), Temperature-independent parameter in PR-EOS
%(KAPPA = 0.37464 + (1.54226 - 0.26992 * OMEGA) * OMEGA;)
T_r = T ./ T_c;                                  % Reduced Temerature
ALPHASqrt = 1 + KAPPA * (1 - sqrt(T_r));
ALPHA = ALPHASqrt^2;                             % Temperature-dependent parameter in PR-EOS
a_T = a_Tc * ALPHA;                              % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
b = 0.077796074 * R_G * T_c ./ p_c;              % m^3/mol, Critical Point Restriction "b", Temperature-independent parameter in PR-EOS
DADT = - a_Tc * KAPPA * ALPHASqrt ./ sqrt(T_c*T);
A = a_T * p ./ ((R_G * T)^2);                    % Parameters for Cubic Form of Equation of State, dimensionless form of EOS parameter a_T
B = b * p ./ (R_G * T);                          % Parameters for Cubic Form of Equation of State, dimensionless form of EOS parameter b
[ZZ,D] = ZZroot(A,B);

% Part2: Solve Peng-Robinson EOS to get compressibility fa_Tctor
% Reference: << Chemical, Biochemical, and Engineering Thermodynamics, 5e >>
%            E6.4-4, P222-223, Ch6
% Root = SolveCubic(Para_TcF(1),Para_TcF(2),Para_TcF(3));
% Root,
Z(1) = max(ZZ); % Vapor Phase, most compressible;
Z(2) = min(ZZ); % Liquid Phase, least compressible;

% Part3: Solve for Peng-Robinson compressibility factor
DepartHS = SolveDepartHS(a_T,b,B,Z,DADT,T,p,R,R_G);
DH(1) = DepartHS(1);
DH(2) = DepartHS(2);
DS(1) = DepartHS(3);
DS(2) = DepartHS(4);

% Part4: Enthalpy and Entropy of Ideal Gas Change from Reference State to (T,p)
H_IG = AAA * (T - T_ref) + BBB * (T^2 - T_ref^2)./2 + ...
       CCC * (T^3 - T_ref^3)./3 + DDD * (T^4 - T_ref^4)./4;      
S_IG = AAA * log(T ./ T_ref) + BBB * (T - T_ref) + ...
       CCC * (T^2 - T_ref^2)./2 + DDD * (T^3 - T_ref^3)./3;          
S_IG = S_IG - R * log(p/p_ref);                                 
H = DH + H_IG;
S = DS + S_IG;

% Optput the result.                                                           
Compressibility = Z;
fprintf('Temperature in this ondition is %4.1f K.\n', T);                          % K
fprintf('Pressure in this condition is %4.1f kPa.\n', p);                          % kPa
fprintf('Enthalpy of saturated vapor in this condition is %f J/mol.\n', H(1));     % J/mol
fprintf('Enthalpy of saturated liquid in this condition is %f J/mol.\n', H(2));    % J/mol
fprintf('Entropy of saturated vapor in this condition is %f J/mol*K.\n', S(1));    % J/mol*K
fprintf('Entropy of saturated liquid in this condition is %f J/mol*K.\n', S(2));   % J/mol*K

SpecifyVolume = Z * 1e3 * R_G * T ./ p;
fprintf('Specify Volume of saturated vapor in this condition is %f.\n',SpecifyVolume(1)); % m^3/kmol
fprintf('Specify Volume of saturated liquid in this condition is %f.\n',SpecifyVolume(2)); % m^3/kmol

fprintf('Compressibility Factor of saturated vapor in this condition is %f.\n',Z(1));
fprintf('Compressibility Factor of saturated liquid in this condition is %f.\n',Z(2));

% Define SubFunction area -------------------------------------------------

% SubFunction1 SolveDepartHS:
function DepartHS = SolveDepartHS(a_T,b,B,Z,DADT,T,p,R,R_G)
for i=1:2,
    ParaFuga1(i) = log((Z(i) + (1+sqrt(2))*B) ./ (Z(i) + (1-sqrt(2))*B));
    DH(i) = (T * DADT - a_T) * ParaFuga1(i)/b/sqrt(8);
    DH(i) = R * T * (Z(i)-1) + DH(i) * R ./ R_G;
    % Gas Enthalpy Departure , EQN 6.4-29
    DS(i) = DADT * ParaFuga1(i) ./ b ./ sqrt(8);
    DS(i) = R * log((Z(i)-B)) + DS(i) * R ./ R_G;
    % Gas Entropy Departure , EQN 6.4-30
end
DepartHS = [DH DS];
end

% SubFunction2 ZZroot: Solve the Equation of State.
function [ZZ,D] = ZZroot(A,B)
V(1) = 1;
V(2) = - 1 + B;
V(3) = A - B * (3 * B + 2);
V(4) = B * (B * B + B - A);
ZZ = roots(V);

v(1) = 2/27 * V(2)^3 - V(2)*V(3)/3 + V(4);
v(2) = V(3) - V(2)^2/3;
D = v(1)^2/4 + v(2)^3/27;

for i = 1:3 % Get rid off the imag root -----------------------------------
    if imag(ZZ(i)) ~= 0
        ZZ(i) = 0;
    end
end

ZZ = sort(ZZ);

if abs(ZZ(1)) < 1e-8
    ZZ(1) = ZZ(3);
end
if abs(ZZ(3)) < 1e-8
    ZZ(3) = ZZ(1);
end

end
