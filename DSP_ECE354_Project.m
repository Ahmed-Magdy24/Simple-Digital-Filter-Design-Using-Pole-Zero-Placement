%% Firtst order discrete-time lowpass filter 
clear; clc; close all;
b = [1];
a1 = [1 -0.6];
a2 = [1  0.6];

w = -pi:0.01:pi;
z = exp(1j*w);

%% Pole-zero plot and the unit circle in the Z-domain

figure('Name', 'First Order LPF - Pole-Zero Plot (a=0.6)');
zplane(b,a1);
figure('Name', 'First Order LPF - Pole-Zero Plot(a=-0.6)');
zplane(b,a2);

%% Magnitude response of the filter in dB in the frequency range (-pi < w < pi ) 

figure('Name', 'First Order LPF - Magnitude Response');

H1 = 1./(1 - 0.6*(z.^(-1)));
H2 = 1./(1 + 0.6*(z.^(-1)));

plot(w, 20*log10(abs(H1)), 'b', 'LineWidth', 1.5);
hold on;
plot(w, 20*log10(abs(H2)), 'r', 'LineWidth', 1.5);
xlabel('w'); 
ylabel('Magnitude Response (dB)');
legend('a = 0.6', 'a = -0.6');
grid on;

%% Phase response of the filter in the frequency range (-π < w ≤ π)

figure('Name', 'First Order LPF - Phase Response');
plot(w, angle(H1), 'b', 'LineWidth', 1.5); 
hold on;
plot(w, angle(H2), 'r', 'LineWidth', 1.5);
xlabel('w');
ylabel('Phase (radian)');
legend('a = 0.6', 'a = -0.6');
grid on;

%% Group delay of the filter in the frequency range (-π < w ≤ π)

figure('Name', 'First Order LPF - Group Delay');
ph1 = unwrap(angle(H1));
ph2 = unwrap(angle(H2));

gd1 = -diff(ph1)./diff(w);
gd2 = -diff(ph2)./diff(w);
w_mid = (w(1:end-1) + w(2:end))/2;

plot(w_mid, gd1, 'b', 'LineWidth', 1.5); 
hold on;
plot(w_mid, gd2, 'r', 'LineWidth', 1.5);
xlabel('w');
ylabel('Group Delay');
legend('a = 0.6', 'a = -0.6');
grid on;

%% The filter’s impulse response
[h1, n1] = impz(b, a1);
[h2, n2] = impz(b, a2);

figure('Name', 'First Order LPF - Impulse Response (a=0.6)');
stem(n1, h1, 'b');
xlabel('n');
ylabel('h_1[n]');
title('Impulse Response (a = 0.6)');
grid on;

figure('Name', 'First Order LPF - Impulse Response(a=-0.6)');
stem(n2, h2, 'r');
xlabel('n');
ylabel('h_2[n]');
title('Impulse Response (a = -0.6)');
grid on;






%% Third Order Discrete-Time Lowpass Filter
clc; clear; 

%% Initial Filter Specifications
wc = 0.25*pi;       % Cutoff frequency
wpo = 0.2*pi;       % original passband edge
wso = 0.3*pi;       % original stopband edge

w = -pi:0.01:pi;

% Pole and Zero Placement
p1 = 0.6;           % Real pole (from first-order LPF)

%% Sweep on r and wp 
         
wso = 0.3*pi;
p1 = 0.6;
wpo = 0.2*pi;

r_values  = 0.1 : 0.01 : 0.99;          % Range of r

wp1_values = 0.18*pi : 0.01*pi : 0.22*pi;   % Small sweep around 0.2*pi

optimum_r = 0;
optimum_wp1 = 0;

Optimum_ripple = inf;

for wp1 = wp1_values

        for r_trial = r_values

            % Complex poles
            p2 = r_trial * exp(1j*wp1);
            p3 = r_trial * exp(-1j*wp1);

            % Complex zeros
            z2 = exp(1j*wso);
            z3 = exp(-1j*wso);

            poles_vec = [p1; p2; p3];
            zeros_vec = [z2; z3];

            % Transfer function coefficients
            b = real(poly(zeros_vec));
            a = real(poly(poles_vec));

            % Frequency response in passband
            w_pass = linspace(0, wpo, 1000);

            H_trial = freqz(b, a, w_pass);

            % Magnitude in dB
            Hdb = 20*log10(abs(H_trial));

            % Passband ripple
            ripple = max(Hdb) - min(Hdb);

            % Save best parameters
            if ripple < Optimum_ripple

                Optimum_ripple = ripple;
                optimum_r = r_trial;
                optimum_wp1 = wp1;

            end
        end
    end
% optimum Parameters

fprintf('===== Optimum Parameters =====\n');
fprintf('optimum r  = %.3f\n', optimum_r);
fprintf('optimum wp1 = %.3f*pi\n', optimum_wp1/pi);
fprintf('Passband Ripple = %.3f dB\n', Optimum_ripple);


%% The coefficients after finding the suitable r and wp

p1 = 0.6;                            % Real pole (from first-order LPF)
p2 = optimum_r * exp(1j*optimum_wp1);      % Complex pole 1 with optimum r from the sweep and optimum wp
p3 = optimum_r * exp(-1j*optimum_wp1);     % Complex pole 2 (conjugate) with optimum r from the sweep and optimum wp
z2 = exp(1j*wso);                    % Complex zero 1 on unit circle
z3 = exp(-1j*wso);                   % Complex zero 2 on unit circle

%% Create transfer function polynomials

poles = [p1, p2, p3];
zeros = [z2, z3];

B = poly(zeros);   % Numerator coefficients
A = poly(poles);   % Denominator coefficients

% Frequency Response

[H, w] = freqz(B, A, 1024, 'whole');

H_shifted = fftshift(H);

gd = grpdelay(B, A, w);


%% 1. Pole-Zero Plot and Unit Circle
figure('Name', '3rd Order LPF - Pole-Zero Plot');
zplane(B, A);
title('Pole-Zero Plot (x: poles, o: zeros)');

%% 2. Magnitude Response in dB (-π to π)
figure('Name', '3rd Order LPF - Magnitude Response');
gain = max(abs(H_shifted));
H_shifted = H_shifted / gain;
plot((w-pi)/pi, 20*log10(abs(H_shifted)));
xlabel('Normalized Frequency (×π rad/sample)');
ylabel('Magnitude (dB)');
title('Fullband Magnitude Response');
grid on;
xlim([-1, 1]);
%% 3. Magnitude Response in dB (-ωp to ωp) for Passband Ripple
figure('Name', '3rd Order LPF - Passband Detail');
plot((w-pi)/pi, 20*log10(abs(H_shifted)));
xlabel('Normalized Frequency (×π rad/sample)');
ylabel('Magnitude (dB)');
title('Passband Magnitude Response (Zoomed)');
grid on;
xlim([-wpo/pi, wpo/pi]); % Zoom into the passband

%% 4. Phase Response in dB (-π to π)
figure('Name', '3rd Order LPF - Phase Response');
plot((w-pi)/pi, unwrap(angle(H_shifted)));
xlabel('Normalized Frequency (×π rad/sample)');
ylabel('Phase (radians)');
title('Phase Response');
grid on;
xlim([-1, 1]);

%% 5. Group Delay Response
figure('Name', '3rd Order LPF - Group Delay');
plot((w-pi)/pi, fftshift(gd));
xlabel('Normalized Frequency (×π rad/sample)');
ylabel('Group Delay (samples)');
title('Group Delay Response');
grid on;
xlim([-1, 1]);

%% 6. Impulse Response
figure('Name', '3rd Order LPF - Impulse Response');
impz(B, A);
title('Impulse Response');
grid on;








%% Fifth Order Discrete-Time Lowpass Filter- First case 
clear; clc; 

wc = 0.25 * pi;         % Cutoff frequency
delta_w = 0.1 * pi;     % Transition bandwidth
wpo = wc - delta_w/2;   % 0.2*pi
wso = wc + delta_w/2;   % 0.3*pi
optimum_wp1 = 0.21*pi ;    % optimum wp for Ap less than 1 dB for third oder

r1 = 0.88;     % Optimum magnitude of poles from third order
r2 = 0.6;      % real pole from first order
r3 = 0.88;  

%% Zeros  (On the unit circle) 
z1 = exp(1j * wso);
z2 = exp(-1j * wso);
w_z_order5 = (wso + pi) / 2; 
z3 = exp(1j * w_z_order5);
z4 = exp(-1j * w_z_order5);
zeros_vec = [z1; z2; z3; z4];

%% Poles
p1 = r2;
p2 = r1 * exp(1j * optimum_wp1);
p3 = r1 * exp(-1j * optimum_wp1);
w_p_ordr5 = (0 + optimum_wp1) / 2; 
p4 = r3 * exp(1j * w_p_ordr5);
p5 = r3 * exp(-1j * w_p_ordr5);
poles_vec = [p1; p2; p3; p4; p5];

%% Filter order5 Coefficients 
b = poly(zeros_vec); % Numerator coefficients
a = poly(poles_vec); % Denominator coefficients

%%  Full Frequency Response 
w = -pi:0.005:pi; 
H_order5_LPF = freqz(b, a, w);

%% --- Plot 1: Pole-Zero Plot ---
figure('Name', '5rd Order LPF - Pole-Zero Plot (Before adjustment)');
zplane(b, a);
grid on;
title('1. Pole-Zero Plot  (Order 5)');

%% Plot 2: Full Magnitude Response 
figure('Name', '5rd Order LPF - Magnitude Response (Before adjustment)');
plot(w, 20*log10(abs(H_order5_LPF)));
xlabel('\omega (rad/sample)');
ylabel('Magnitude Response (dB)');
title('2. Full Magnitude Response (-\pi < \omega \leq \pi)');
grid on;

%% Plot 3: Passband Magnitude Response (For Ripple Calculation) 
figure('Name', '5rd Order LPF - Passband Detail');
w_passband = -wpo:0.001:wpo;
H_passband = freqz(b, a, w_passband);
plot(w_passband, 20*log10(abs(H_passband)));
xlabel('\omega (rad/sample)');
ylabel('Magnitude Response (dB)');
title('3. Passband Magnitude Response (-\omega_p \leq \omega \leq \omega_p)');
grid on;

%%  Plot 4: Phase Response 
figure('Name', '5rd Order LPF - Phase Response (Before adjustment)');
plot(w, unwrap(angle(H_order5_LPF))); 
xlabel('\omega (rad/sample)');
ylabel('Phase (radians)');
title('4. Phase Response (-\pi < \omega \leq \pi)');
grid on;

%% Plot 5: Group Delay Response 
figure('Name', '5rd Order LPF - Group Delay Response (Before adjustment)');
[gd, w_gd] = grpdelay(b, a, w);
plot(w_gd, gd);
xlabel('\omega (rad/sample)');
ylabel('Group Delay (samples)');
title('5. Group Delay Response (-\pi < \omega \leq \pi)');
grid on;

%% --- Plot 6: Impulse Response ---
figure('Name', '5rd Order LPF - Impulse Response (Before adjustment)');
impz(b, a);
title('6. Filter''s Impulse Response h[n]');
grid on;


%% Sweep on r1 , r2 and r3

r1_values  = 0.1 : 0.01 : 0.99;
r2_values  = 0.1 : 0.01 : 0.99;
r3_values  = 0.1 : 0.01 : 0.99;

Optimum_r1 = 0;
Optimum_r2 = 0;
Optimum_r3 = 0;
Optimum_ripple = inf;

    for r1_trial = r1_values

        for r2_trial = r2_values

            for r3_trial = r3_values
               
                % Complex Zeros
                z1 = exp(1j*wso);
                z2 = exp(-1j*wso);
                w_z_order5 = (wso + pi)/2;
                z3 = exp(1j*w_z_order5);
                z4 = exp(-1j*w_z_order5);

                zeros_vec = [z1; z2; z3; z4];

                % Poles
                p1 = r2_trial;
                p2 = r1_trial * exp(1j*optimum_wp1);
                p3 = r1_trial * exp(-1j*optimum_wp1);
                w_p_order5 = (0 + optimum_wp1)/2;
                p4 = r3_trial * exp(1j*w_p_order5);
                p5 = r3_trial * exp(-1j*w_p_order5);

                poles_vec = [p1; p2; p3; p4; p5];

                % Transfer Function Coefficients
                b = real(poly(zeros_vec));
                a = real(poly(poles_vec));

                % Frequency Response in Passband   
                w_passband = linspace(0, wpo, 2000);

                H_trial = freqz(b, a, w_passband);

                % Normalize Gain
                H_trial = H_trial / abs(H_trial(1));

                % Magnitude in dB
                Hdb = 20*log10(abs(H_trial));

                % Passband Ripple
                ripple = max(Hdb) - min(Hdb);

                % Save Optimum Parameters

                if ripple < Optimum_ripple

                    Optimum_ripple = ripple;

                    Optimum_r1 = r1_trial;
                    Optimum_r2 = r2_trial;
                    Optimum_r3= r3_trial;
                end

            end
        end
    end

% Optimum Parameters

fprintf('===== Optimum Parameters =====\n');

fprintf('Optimum r1  = %.3f\n', Optimum_r1);
fprintf('Optimum r2  = %.3f\n', Optimum_r2);
fprintf('Optimum r3  = %.3f\n', Optimum_r3);
fprintf('Passband Ripple = %.3f dB\n', Optimum_ripple);





%% Fifth Order Discrete-Time Lowpass Filter- Second case with optimum r1, r2 and r3 with passband ripples = 0.444 dB

r1 = 0.9;   
r2 = 0.19; 
r3 = 0.52;  
 
%% Zeros  (On the unit circle) 
z1 = exp(1j * wso);
z2 = exp(-1j * wso);
w_z_order5 = (wso + pi) / 2; 
z3 = exp(1j * w_z_order5);
z4 = exp(-1j * w_z_order5);
zeros_vec = [z1; z2; z3; z4];

%% Poles
p1 = r2;
p2 = r1 * exp(1j * optimum_wp1);
p3 = r1 * exp(-1j * optimum_wp1);
w_p_ordr5 = (0 + optimum_wp1) / 2; 
p4 = r3 * exp(1j * w_p_ordr5);
p5 = r3 * exp(-1j * w_p_ordr5);
poles_vec = [p1; p2; p3; p4; p5];

%% Filter order5 Coefficients 
b = poly(zeros_vec);  % Numerator coefficients
a = poly(poles_vec);  % Denominator coefficients

%%  Full Frequency Response 
w = -pi:0.005:pi; 
H_order5_LPF = freqz(b, a, w);

%% --- Plot 1: Pole-Zero Plot ---
figure('Name', '5rd Order LPF - Pole-Zero Plot (After adjustment)');
zplane(b, a);
grid on;
title('1. Pole-Zero Plot  (Order 5)');

%% Plot 2: Full Magnitude Response 
figure('Name', '5rd Order LPF - Magnitude Response (After adjustment)');
plot(w, 20*log10(abs(H_order5_LPF)));
xlabel('\omega (rad/sample)');
ylabel('Magnitude Response (dB)');
title('2. Full Magnitude Response (-\pi < \omega \leq \pi)');
grid on;

%% Plot 3: Passband Magnitude Response (For Ripple Calculation) 
figure('Name', '5rd Order LPF - Passband Detail');
w_passband = -wpo :0.001:wpo;
H_passband = freqz(b, a, w_passband);
plot(w_passband, 20*log10(abs(H_passband)));
xlabel('\omega (rad/sample)');
ylabel('Magnitude Response (dB)');
title('3. Passband Magnitude Response (-\omega_p \leq \omega \leq \omega_p)');
grid on;

%%  Plot 4: Phase Response 
figure('Name', '5rd Order LPF - Phase Response (After adjustment)');
plot(w, unwrap(angle(H_order5_LPF))); 
xlabel('\omega (rad/sample)');
ylabel('Phase (radians)');
title('4. Phase Response (-\pi < \omega \leq \pi)');
grid on;

%% Plot 5: Group Delay Response 
figure('Name', '5rd Order LPF - Group Delay Response (After adjustment)');
[gd, w_gd] = grpdelay(b, a, w);
plot(w_gd, gd);
xlabel('\omega (rad/sample)');
ylabel('Group Delay (samples)');
title('5. Group Delay Response (-\pi < \omega \leq \pi)');
grid on;

%% -- Plot 6: Impulse Response ---
figure('Name', '5rd Order LPF -Impulse Response (After adjustment)');
impz(b, a);
title('6. Filter''s Impulse Response h[n]');
grid on; 







%% Frequency Transformation using Pole-Zero Pattern Rotation (HPF) 

zeros_HPF = -zeros_vec;        % Rotate all zeros by 180° on the z-plane
poles_HPF = -poles_vec;        % Rotate all poles by 180° on the z-plane

b_HPF = poly(zeros_HPF);
a_HPF = poly(poles_HPF);

%%  Full Frequency Response 
w = -pi:0.005:pi;
H_order5_HPF = freqz(b_HPF,a_HPF,w);

%% --- Plot 1: Pole-Zero Plot ---
figure('Name', '5rd Order HPF - Pole-Zero Plot');
zplane(b_HPF,a_HPF);
grid on;
title('7. Pole-Zero Plot of the HPF');

%% Plot 2: Full Magnitude Response 
figure('Name', '5rd Order HPF - Magnitude Response ');
plot(w, 20*log10(abs(H_order5_HPF)));
xlabel('\omega (rad/sample)');
ylabel('Magnitude Response (dB)');
title('2. Full Magnitude Response (-\pi < \omega \leq \pi)');
grid on;

%% Plot 3: Passband Magnitude Response (For Ripple Calculation) 
figure('Name', '5rd Order HPF - Passband Magnitude Response ');
w_passband = pi-wpo:0.001:pi ;
H_passband = freqz(b_HPF, a_HPF, w_passband);
plot(w_passband, 20*log10(abs(H_passband)));
xlabel('\omega (rad/sample)');
ylabel('Magnitude Response (dB)');
title('3. Passband Magnitude Response');
grid on;

%%  Plot 4: Phase Response 
figure('Name', '5rd Order HPF - Phase Response');
plot(w, unwrap(angle(H_order5_HPF))); 
xlabel('\omega (rad/sample)');
ylabel('Phase (radians)');
title('4. Phase Response (-\pi < \omega \leq \pi)');
grid on;

%% Plot 5: Group Delay Response 
figure('Name', '5rd Order HPF - Group Delay Response');
[gd, w_gd] = grpdelay(b_HPF, a_HPF, w);
plot(w_gd, gd);
xlabel('\omega (rad/sample)');
ylabel('Group Delay (samples)');
title('5. Group Delay Response (-\pi < \omega \leq \pi)');
grid on;

%% -- Plot 6: Impulse Response ---
figure('Name', '5rd Order HPF - Impulse Response');
impz(b_HPF, a_HPF);
title('6. Filter''s Impulse Response h[n]');
grid on; 



%% Frequency Transformation using Pole-Zero Pattern Rotation (BPF) 

% Apply Frequency Transformation (BPF at pi/2) 
% To ensure real coefficients, we must include BOTH +pi/2 and -pi/2 shifts.
% This doubles the order of the filter 

shift_pos = exp(1j * pi/2);
shift_neg = exp(-1j * pi/2);

zeros_vec_bpf = [zeros_vec * shift_pos; zeros_vec * shift_neg];
poles_vec_bpf = [poles_vec * shift_pos; poles_vec * shift_neg];

% Filter BPF Coefficients 
b_bpf = real(poly(zeros_vec_bpf)); 
a_bpf = real(poly(poles_vec_bpf)); 

%% --- 4. Plotting ---
w = -pi:0.005:pi; 
H_bpf = freqz(b_bpf, a_bpf, w);

%% Plot 1: Pole-Zero Plot 
figure('Name', '5rd Order BPF - Pole-Zero Plot');
zplane(b_bpf, a_bpf);
grid on;
title('1. Pole-Zero Plot (10th Order BPF)');

%% Plot 2: Full Magnitude Response 
figure('Name', '5rd Order BPF - Magnitude Response');
plot(w, 20*log10(abs(H_bpf)));
xlabel('\omega (rad/sample)');
ylabel('Magnitude Response (dB)');
title('2. Full Magnitude Response (-\pi < \omega \leq \pi)');
grid on;

%% Plot 3: Passband Magnitude Response (For Ripple Calculation)
% The passband is now centered around pi/2 instead of 0.
figure('Name', '5rd Order BPF - Passband Magnitude Response');
w_center = pi/2;
w_pass_bpf = (w_center - wpo):0.001:(w_center + wpo);
H_pass_bpf = freqz(b_bpf, a_bpf, w_pass_bpf);

plot(w_pass_bpf, 20*log10(abs(H_pass_bpf)));
xlabel('\omega (rad/sample)');
ylabel('Magnitude Response (dB)');
title('3. Passband Magnitude Response (BPF Passband)');
grid on;

%% Plot 4: Phase Response 
figure('Name', '5rd Order BPF - Phase Response');
plot(w, unwrap(angle(H_bpf))); 
xlabel('\omega (rad/sample)');
ylabel('Phase (radians)');
title('4. Phase Response (-\pi < \omega \leq \pi)');
grid on;

%% Plot 5: Group Delay Response 
figure('Name', '5rd Order BPF - Group Delay Response');
[gd, w_gd] = grpdelay(b_bpf, a_bpf, w);
plot(w_gd, gd);
xlabel('\omega (rad/sample)');
ylabel('Group Delay (samples)');
title('5. Group Delay Response (-\pi < \omega \leq \pi)');
grid on;

%% Plot 6: Impulse Response 
figure('Name', '5rd Order BPF - Impulse Response');
impz(b_bpf, a_bpf);
title('6. Filter''s Impulse Response h[n] (BPF)');
grid on;







%% Comb Filter 
clear; clc; 

% Optimum magnitude values from fifth order 
r1 = 0.9 ;   
r2 = 0.19 ; 
r3 = 0.52 ;  


wc = 0.25 * pi;            % Cutoff frequency
delta_w = 0.1 * pi;        % Transition bandwidth
wpo = wc - delta_w/2;      % 0.2*pi
wso = wc + delta_w/2;      % 0.3*pi
optimum_wp1 = 0.21*pi ;    % optimum wp for Ap less than 1 dB for third oder


%% Zeros  (On the unit circle) 
z1 = exp(1j * wso);
z2 = exp(-1j * wso);
w_z_order5 = (wso + pi) / 2; 
z3 = exp(1j * w_z_order5);
z4 = exp(-1j * w_z_order5);
zeros_vec = [z1; z2; z3; z4];

%% Poles
p1 = r2; 
p2 = r1 * exp(1j * optimum_wp1);
p3 = r1 * exp(-1j * optimum_wp1);
w_p_ordr5 = (0 + optimum_wp1) / 2; 
p4 = r3 * exp(1j * w_p_ordr5);
p5 = r3 * exp(-1j * w_p_ordr5);
poles_vec = [p1; p2; p3; p4; p5];

%% Filter  Coefficients 
b = poly(zeros_vec); % Numerator coefficients
a = poly(poles_vec); % Denominator coefficients

L = 8;               % Comb filter factor

b_comb = zeros(1,(length(b)-1)*L + 1);             % New numerator coefficients
a_comb = zeros(1,(length(a)-1)*L + 1);             % New denominator coefficients

b_comb(1:L:end) = b;                               % Insert zeros between numerator coefficients
a_comb(1:L:end) = a;                               % Insert zeros between denominator coefficients

%% Full Frequency Response 
w = -pi:0.005:pi; 
H_comb = freqz(b_comb, a_comb, w);

%% --- Plot 1: Pole-Zero Plot ---
figure('Name', ' Comb Filter - Pole-Zero Plot');
zplane(b_comb, a_comb);
grid on;
title('1. Pole-Zero Plot  (Comb Filter)');

%% Plot 2: Full Magnitude Response 
figure('Name', ' Comb Filter - Magnitude Response');
plot(w, 20*log10(abs(H_comb)));
xlabel('\omega (rad/sample)');
ylabel('Magnitude Response (dB)');
title('2. Full Magnitude Response (-\pi < \omega \leq \pi)');
grid on;

%% Plot 3: Passband Magnitude Response 
figure('Name', ' Comb Filter - Passband Magnitude Response');
w_pass_comb = -wpo:0.001:wpo;
H_passband = freqz(b, a, w_pass_comb);
plot(w_pass_comb, 20*log10(abs(H_passband)));
xlabel('\omega (rad/sample)');
ylabel('Magnitude Response (dB)');
title('3. Passband Magnitude Response (-\omega_p \leq \omega \leq \omega_p)');
grid on;

%%  Plot 4: Phase Response 
figure('Name', ' Comb Filter - Phase Response');
plot(w, unwrap(angle(H_comb))); 
xlabel('\omega (rad/sample)');
ylabel('Phase (radians)');
title('4. Phase Response (-\pi < \omega \leq \pi)');
grid on;

%% Plot 5: Group Delay Response 
figure('Name', ' Comb Filter - Group Delay Response');
[gd, w_gd] = grpdelay(b_comb, a_comb, w);
plot(w_gd, gd);
xlabel('\omega (rad/sample)');
ylabel('Group Delay (samples)');
title('5. Group Delay Response (-\pi < \omega \leq \pi)');
grid on;

%% -- Plot 6: Impulse Response ---
figure('Name', ' Comb Filter - Impulse Response');
impz(b_comb, a_comb);
title('6. Filter''s Impulse Response h[n]');
grid on; 






