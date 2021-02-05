% Example of 'magnetic_field' function.
% Computes the magnetic field of a circular current loop
% on the X-Y plane, and computes the loop magnetic moment (m)

% Written by Prof. Yoash Levron, Technion, Israel, 2014.

clc;

meu0 = 4*pi*1e-7; % [H/m]  (Henry / meter)
cur = 1;  % [A] loop current

% resolution in the X-Y plane
dx = 0.001;  % [m]
dy = 0.001;  % [m]

% observation region in the X-Y plane
Xmin = 0.15;  % [m]
Xmax = 0.45;  % [m]
Ymin = 0.55;   % [m]
Ymax = 0.85;   % [m]

% define current loop geometry:
cnt_x = 0.3 - dx/2;  % [m]  ring center on x
cnt_y = 0.7 - dy/2;  % [m]  ring center on y
ring_radius = 0.1;  % [m] ring radius
num_vertixes = 40;  % number of verticex on the ring

% create ring vertixes
d_teta = (2*pi)/num_vertixes;
teta = (d_teta/2):d_teta:(2*pi);
px = cnt_x + ring_radius*cos(teta);
py = cnt_y + ring_radius*sin(teta);
FROM = zeros(num_vertixes,3);
TO = zeros(num_vertixes,3);
for ii=1:(num_vertixes-1)
    FROM(ii,:) = [px(ii) py(ii) 0];
    TO(ii,:) = [px(ii+1) py(ii+1) 0];
end
FROM(num_vertixes,:) = [px(end) py(end) 0];
TO(num_vertixes,:) = [px(1) py(1) 0];
CUR =  cur*ones(num_vertixes,1);

% create observation points matrix
Xvec = Xmin:dx:Xmax;  NX = length(Xvec);
Yvec = Ymin:dy:Ymax;  NY = length(Yvec);
usermem = memory;  MaxPossibleArrayDbl=usermem.MaxPossibleArrayBytes / 8;
if (NX*NY > MaxPossibleArrayDbl / 200)
    'Warning - possibly not enough memory. Program terminated. '
    return
end
[X, Y] = meshgrid(Xvec, Yvec);
R = [X(:), Y(:), 0*X(:)];

% compute the field everywhere
Hmat = magnetic_field( FROM, TO, CUR, R );
HmatZ_vec = Hmat(:,3);  % magnetic field in z direction
HmatZ = reshape(HmatZ_vec,NY,NX);

% display field:
figure(1);
H_Z_dB = 20*log10(abs(HmatZ));
imagesc(Xvec,Yvec, H_Z_dB);
xlabel('X [m]');  ylabel('Y [m]');
title('Magnetic field magnitude in dB');
colorbar;

% MAGNETIC MOMENT
R_far_field = 1e3;  % [m]
H_far_field = magnetic_field( FROM, TO, CUR, [R_far_field 0 0] );
H_far_field_z = H_far_field(1,3);
M = abs( H_far_field_z*4*pi*R_far_field^3 );
disp('The magnetic moment should be equal to the loop-area * loop-current');
disp('magnetic moment in A-m^2 :');
Magnetic_Moment_in_Am2 = M






