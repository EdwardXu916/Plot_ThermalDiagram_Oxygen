% Title: Calculate the other three variables based on two knowing
%        thermodynamic variables. And plot the T-s Diagram of R123.
% Based on: MATLAB program from << Chemical, Biochemical, and Engineering
%                                  Thermodynamics >>
% Version: 1.0, Edward Xu, 18.4.10
% SubTitle: Main.

T = ;
p = ;

% p_atm = 101.325 kPa

for i = 1:10:100
    T = 273.15 + 
    for j = 1:10:100
        
    end
    
end

%{
x(1) = - 125 + 273.15;
for i = 1:100
    x(2) = 101.325 * i;
    ThermoProp_R123_EdXu
    y1(i) = H(1);
    y2(i) = H(2);
    y3(i) = S(1);
    y4(i) = S(2);
    z(i) = x(2);
end
plot(y1,z,y2,z)
% legend(y1,y2)
%}

%{
x(2) = 101.325;
for i = 1:100
    x(1) = 273.15 - 200 + 10*i;
    ThermoProp_R123_EdXu
    y1(i) = H(1);
    y2(i) = H(2);
    y3(i) = S(1);
    y4(i) = S(2);
    z(i) = x(2);
    Z
end
plot(y1,z,y2,z)
% legend(y1,y2)
%}