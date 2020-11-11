% Idetification of a magnetic levitation/suspension plant/system.
% Measurement and estimation of parameters of a mathematical model of the
% plant.
% A levitating object is a hollow thin-walled ball made by pressing|
% stamping of mild steel sheet. It is attracted by an electromagnet with an
% E-shaped laminated transformer steel core and a enamelled solid coper
% wire winding.

%%
% Gravity of Earth:
g_E = 9.81 % m/s^4

%%
% Measurements of geometric and inertial parameters of a steel thin-walled
% ball.
% Ball mass:
m_b = 58e-3 % kg
% Ball diameter:
% d = % m

%%
% Measurements of a magnetic core and an electric winding parameters.
% Dimensions of the laminated E-shaped core:

% Number of coil turns:
% n_w = % 1
% Enamelled solid coper wire diameter:
% d_w = 

%%
% Linear viscous air drag resistance coefficient:
b_v = 0 % (can be neglected)

%%
% Finding equivalent/net/resultant/overall resistance of the electrical
% circuit comprising the electromagnet winding and other components.
% Measurements of a steady-stade relationship between the duty factor of a
% PWM control signal and an open-circuit output average voltage of a power
% stage and an output average current of the power stage loaded with the
% winding.
% Power stage output voltage should be low-pass filtered (averaged) for
% measurement. The same applies to the power stage output current.
w = 1:-0.05:0; % [0:0.05:1, 1:-0.05:0];
U = [12.19, 11.5805, 10.971, 10.3615, 9.752, 9.1425, 8.533, 7.9235, ...
     7.314, 6.7045, 6.095, 5.4855, 4.876, 4.2665, 3.657, 3.0475, 2.438, ...
     1.8285, 1.219, 0.6095, 0];
I = [2.21, 2.12, 2.03, 1.93, 1.83, 1.72, 1.62, 1.51, 1.40, 1.29, 1.18, ...
    1.07, 0.95, 0.83, 0.72, 0.60, 0.48, 0.35, 0.23, 0.11, 0.00];

% To simplify model we can assume that R is constant so wa can calculate it
% as a mean value for all measurements
R = U./I;
R_n = mean(R(1:end-1))
clear w U I R
% Temperature resistance coefficient of copper:
% alpha =  0.00393

%%
% Identification of the emf of the power supply
% Measurements of a steady-stade relationship between the duty factor of a
% PWM control signal and the open-circuit output average voltage of a power
% stage.
w = 1:-0.05:0; % [0:0.05:1, 1:-0.05:0];
U = [12.19, 11.5805, 10.971, 10.3615, 9.752, 9.1425, 8.533, 7.9235, ...
    7.314, 6.7045, 6.095, 5.4855, 4.876, 4.2665, 3.657, 3.0475, 2.438, ...
    1.8285, 1.219, 0.6095, 0];

figure
plot(w, U, '.', 'MarkerSize', 12)
grid on
xlabel('PWM duty factor')
ylabel('open-circuit output average voltage of a power stage [V]')
emf_pf = polyfit(w,U,1)
hold on
plot(w,polyval(emf_pf,w),'r-')
clear w U

%%
% Measurements of steady-state characteristic (voltage versus current) of a
% current/voltage transducer integrated into an electric unit of the
% magnetic levitation plant.
w = 1:-0.05:0;
U = [1.93, 1.85, 1.75, 1.65, 1.55, 1.47, 1.39, 1.30, 1.21, 1.12, 1.03, ...
    0.93, 0.83, 0.73, 0.63, 0.53, 0.42, 0.31, 0.21, 0.10, 0.015];
I = [2.25, 2.14, 2.03, 1.93, 1.82, 1.72, 1.61, 1.51, 1.40, 1.29, 1.19, ...
    1.07, 0.95, 0.84, 0.72, 0.60, 0.48, 0.36, 0.24, 0.11, 0.00];

clear w U I

%%
% Measurements of a steady-state characteristic (voltage versus distance)
% of an optical ball position sensor (measuring the width of an air gap
% between the face of a central limb of the electromagnet laminated core
% and the surface of a thin-walled stamped mild steel sheet ball.
% An M4 threaded/tapped rod/bar/stud is used as a length standard
% (reference). The thread pitch equals 0.7 mm.
% We shall measure bright/uncovered and dark/covered output current/voltage
% of the electrooptical ball position sensor for various ambient light(ing)
% condition -- at dark night, with electical lighting, at sunny day.
z = 1e-3*(0:0.7:21);
U = [10, 10, 10, 10, 9.98, 9.945, 9.88, 9.795, 9.645, 9.45, 9.195, ...
    8.895, 8.545, 8.18, 7.81, 7.43, 7.06, 6.7, 6.355, 6.035, 5.74, ...
    5.48, 5.25, 5.065, 4.92, 4.82, 4.76, 4.72, 4.71, 4.71, 4.71];

clear z U

%%
% Mesaurements of the relationship between the ball-core distance and the
% inductance of the electromagnet winding.
% Measurement conditions (settings of the electronic RCL "bridge"
% instrument):
% f = 
% U = 
% series equivalent circuit
% We should have use decay method to measure/estimate/identify inductance
% of the electromagnet solenoid rather than an RCL bridge. But this method
% is time-consuming and difficult to automate. It would also required firm
% immobilization of the ball as a time-varying and relatively strong
% atraction force would act on it during the test causing it to move and to
% disturb measurement results.
z = 1e-3*(0:0.7:21);
L = 1e-3*[137, 130.6, 126.5, 123.6, 121.4, 119.3, 117.8, 116.5, 115.4, ...
    114.4, 113.5, 112.7, 112, 111.5, 110.8, 110.4, 109.9, 109.5, 109.2, ...
    108.8, 108.5, 108.2, 108, 107.8, 107.6, 107.5, 107.2, 107.1, 106.9, ...
    106.8, 106.7];
R = [7.06, 5.89, 5.44, 5.15, 4.98, 4.81, 4.72, 4.6, 4.55, 4.48, 4.4, ...
    4.35, 4.34, 4.31, 4.27, 4.25, 4.21, 4.2, 4.19, 4.16, 4.14, 4.13, ...
    4.14, 4.09, 4.12, 4.12, 4.09, 4.09, 4.09, 4.06, 4.07];

fitType = fittype('a*exp(-b*x)+c+d*exp(-e*x)');
fitOptions = fitoptions(fitType);
fitOptions.StartPoint = [1, 2, 3, 4, 5];
fitOptions.Lower = [0, 0, 0, 0, 500];
fitOptions.Upper = [0.1, 400, 1, 0.1, 1000];

fitObject = fit(z',L',fitType,fitOptions);
L_fitted = feval(fitObject,z');

figure
plot(z, L, 'b.', z', L_fitted, 'r-', 'MarkerSize', 12, 'LineWidth', 1.5)
grid on
xlabel('ball-core distance [m]')
ylabel('coil inductance - measured and fitted [H]')
ind.a = fitObject.a;
ind.b = fitObject.b;
ind.c = fitObject.c;
ind.d = fitObject.d;
ind.e = fitObject.e
clear z L R
clear fitObject fitOptions fitType L_fitted

%%
% set(gcf,'PaperPositionMode','auto')

%% Equilibrum point

z_ep = 0.010; % m
v_ep = 0; % m/s
L_ep = ind.a*exp(-ind.b*z_ep) + ind.c + ind.d*exp(-ind.e*z_ep);
L_d_ep = -ind.a*ind.b*exp(-ind.b*z_ep) - ind.d*ind.e*exp(-ind.e*z_ep);
L_dd_ep = ind.a*ind.b^2*exp(-ind.b*z_ep) + ind.d*ind.e^2*exp(-ind.e*z_ep);

i_ep = sqrt(-2*m_b*g_E/L_d_ep);
u_ep = R_n * i_ep;


%% Plant model linearization at the equilibrium point

A = [0, 1, 0
    1/2/m_b*L_dd_ep*i_ep^2, -b_v/m_b, 1/m_b*L_d_ep*i_ep
    (1/L_ep*L_d_ep^2 - 1/L_ep^2*L_dd_ep^2)*v_ep*i_ep, ...
    -1/L_ep*L_d_ep*i_ep, -1/L_ep*L_d_ep*v_ep - R_n/L_ep]
B = [0; 0; 1/L_ep]
C = [1 0 0; 0 0 1]
D = [0; 0]

%%



