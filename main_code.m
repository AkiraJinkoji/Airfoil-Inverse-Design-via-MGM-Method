clc;
clear;
fclose('all');
% generate Target Cp values, just once..
% Alternatively, we may read file with x vs. Cp for the target airfoil
%
% in the function Panel, we open the file and read the airfoil coordinates
% in the function Panel1, we supply thex and y as arrayes passed on
% to panel1 as arguments
% This is the primary difference between Panel and Panel1
%
% Generate target Cp values for the NLR1 airfoil (our target shape)
%fid = fopen('targetCp.dat','r');
%[nu, nl, alpha,xut,xlt,yut,ylt,cput,cplt]=Panel(fid);
% Read CSV (skip header if present)
data = readmatrix('targetCp.csv'); 
data = data';

mach = 0.5;
alphad = 2.0;
correction_factor_compress = sqrt(1-mach^2);

% Assign columns to variables
xlt = flip(data(1,:));   % Lower surface x-coordinates
cplt = flip(data(2,:));  % Lower surface Cp values
cplt = cplt*correction_factor_compress; 
xut = data(3,:);   % Upper surface x-coordinates
cput = data(4,:);  % Upper surface Cp values
cput = cput*correction_factor_compress;

nu = length(xut);
nl = length(xlt);

% Plot target airfoil shape
figure(1); % airfoil shape

% Begin inverse design,  20 updates, starting with NACA 0012 at the same
% alpha
%
nitmax=20; % have to vary the AoA as the LE and TE coordinates would not change
for nit = 1:nitmax
    if(nit==1)
    % Open INITIAL airfoil data and calculate starting pressure distribution
    %fid = fopen('Naca0012.dat','r');
    file_name = 'Naca0012_CorrectedFormat.dat';
    [nu_ini, nl_ini, alpha_ini,xu1,xl1,yu1,yl1,cpu1,cpl1]=Xfoil(file_name, alphad);
   %Plot initial airfoil too..
   plot(xu1,yu1,'black','LineWidth',3,'DisplayName','Starting Upper');
    hold on
    plot(xl1,yl1,'black','LineWidth',3,'DisplayName','Starting Lower');
   figure(2)
   plot(xl1,-cpl1,'green', 'DisplayName','Initial lower');
   hold on
   plot(xu1,-cpu1,'green', 'DisplayName','Initial upper');
    else
   [cpu1,cpl1] = Xfoil1(nu_ini,nl_ini,xu1,xl1,yu1,yl1,alpha_ini);
    end
    
%Interpolate the CP_initial to x locations on the target airfoil
%upper surface first
dcpu=zeros(nu_ini,1);
dcpl=zeros(nl_ini,1);

for i=2:nu_ini-1
    if (xu1(i) >= xut(1)) && (xu1(i) <= xut(nu))  % to avoid extrapolation
        cpua(i)=interp1(xut,cput,xu1(i));% if cp is not available at the control point have to set delta z =0
        dcpu(i)=cpua(i)-cpu1(i);
    else 
        dcpu(i) = 0;
    end
   % fprintf('%5d %14.8f %14.8f\n',i, xu1(i),dcpu(i));
end
dcpu(1)=0.0;
dcpu(nu)=0.0;

%lower surface next
for i=2:nl_ini-2
    if (xl1(i) >= xlt(1)) && (xl1(i) <= xlt(nu)) 
        cpla(i)=interp1(xlt,cplt,xl1(i));
        dcpl(i)=cpla(i)-cpl1(i); % same here
    else
        dcpl(i)=0;
    end 
    %  fprintf('%5d %14.8f %14.8f\n',i, xl1(i),dcpl(i));
end
dcpl(nl)=0.0;
dcpl(1)=0.0;

icount=0;
dzu=zeros(nu,1);
dzl=zeros(nl,1);
%upper surface redesign
[dzu]= mgmu(nu_ini,dcpu,xu1);
[dzl]= mgml(nl_ini,dcpl,xl1);
% Add dzu and dzl to original airfoil
% And create a new input file for the new airfoil

for i = 1:nu_ini
    yu1(i) = yu1(i) - dzu(i);
    fprintf(' %14.8f %14.8f\n',xu1(i),yu1(i));
end
for i = 1:nl_ini
    yl1(i) = yl1(i) + dzl(i);
    fprintf('%14.8f %14.8f\n',xl1(i),yl1(i));
end

end % MGM iterations 

% create the .dat file for the final airfoil
x = zeros(1,nl+nu-1);
y = zeros(1,nl+nu-1);
% Transfer supplied xu,xl etc to x and y arrays
for i = nl:nl+nu-1  
   l = i+1-nl;
   x(i) = xu1(l);
   y(i) = yu1(l);
end
    for i =1:nl 
      x(nl+1-i)=xl1(i);
      y(nl+1-i)=yl1(i);
    end

% Create .dat file for the new airfoil
dat_file_name = 'airfoil_final.dat';
% Open file for writing
fid_dat = fopen(dat_file_name, 'w');
if fid_dat == -1
    error('Failed to create XFOIL input file');
end

% Write the dat.file
fprintf(fid_dat, 'airfoil iteration \n');
% Write coordinates
if length(x)==length(y)
    for i = 1:length(x)
        fprintf(fid_dat, '%.6f %.6f\n', x(i), y(i));% 6 digits
    end
else
    error('x and y dont have the same length')
end
fclose(fid_dat);

figure(1)
plot(xu1,yu1,'--or','DisplayName','Final Upper');
hold on
plot(xl1,yl1,'--or','DisplayName','Final Lower');
title('Comparison of the airfoil geometry between the target and final version')
xlabel('x')
ylabel('y')
legend

% Compared target Cp with Cp from our design

figure (2);
plot(xlt,-cplt,'blue', 'DisplayName','Target Upper');
hold on
plot(xut,-cput,'blue', 'DisplayName','Target Upper');
plot(xl1,-cpl1,'--or', 'DisplayName','Final Upper');
plot(xu1,-cpu1,'--or', 'DisplayName','Final Upper');
ylim([-1, 1]);
xlim([0,0.99])
title('Comparison of the Cp distribution between the target and final version')
xlabel('x')
ylabel('-Cp')
legend


%% -------------------------------------------------------------------------

function [nu, nl, alphad, xu,xl,yu,yl,cpu,cpl] = Xfoil(file_name, alphad)

%Read airfoil coordinates
%alphad = fscanf(fid,'%f',1);    % Read Angle of Attack
alpha = alphad * pi /180;       % Convert to radians
%nu = fscanf(fid,'%d',1);        % Read number of points on the upper side of airfoil
%nl = fscanf(fid, '%d',1);       % Read number of points on the lower side of airfoil
%isym = fscanf(fid,'%d',1);      % Read flag that states if this airfoil is symmetric
isym = 0;
%factor=fscanf(fid,'%f',1);      % Read scaling factor % <---- which scaling factor it is ?
factor = 1;
if(isym>0)
  nl = nu;
end

fid = fopen(file_name, 'r');
fgetl(fid); % Skip header lines
data = textscan(fid, '%f %f %f'); % Read numeric data
fclose(fid);
x = data{1};
y = data{2};
nu = floor(length(x)/2)+1; % the point (0,00) belongs to both upper and lower surface
nl = floor(length(x)/2)+1;

% disp('test reading of x and y of .dat file')
% disp(length(x))
% disp(y)

xu = zeros(1,nu);
yu = zeros(1,nu);
xl = zeros(1,nl);
yl = zeros(1,nl);
for i = 1:nu              % Read the points on the upper surface
   xu(i)=x(i);
   yu(i)=y(i);
end
xu = flip(xu);
yu = flip(yu);
for i = nu:nu+nl-1 % Read the points on the upper surface
    l = i - nu+1;
    xl(l)=x(i);
   yl(l)=y(i);
end

% disp('test xupper and lower')
% disp(length(xu))
% disp(xu)
% disp(yu)
% disp(yl)

for i = 1:nu+nl-1
    fprintf('%14.8f %14.8f\n', x(i),y(i));
end


% compute Cp using xfoil

% create input file and run xfoil with it

% Define the XFOIL input filename
input_file = 'xfoil_input.txt';

% Define parameters
Re = 3000000;       % Reynolds number
alfa = alphad;           % AoA
iter_limit = 200;   % max iteration
airfoil_file = file_name; % Airfoil filename
cp_output_file = 'cp_alfa_aoa2.txt'; % Output file for Cp data

% Open file for writing
fidin = fopen('xfoil_input.txt', 'w');
if fidin == -1
    error('Failed to create XFOIL input file');
end

% Write XFOIL commands to input file
fprintf(fidin, 'load %s\n', airfoil_file);  % Load airfoil coordinates
fprintf(fidin, 'oper\n');              % command
%fprintf(fid, 'visc\n');              % command visc
%fprintf(fid, '%d\n', Re);            % Reynolds number
fprintf(fidin, 'iter %d\n', iter_limit); % Set iteration limit
fprintf(fidin, 'alfa %.2f\n', alfa);   % Set angle of attack
fprintf(fidin, 'cpwr %s\n', cp_output_file); % Save pressure coefficient data
fprintf(fidin, '\n');
fprintf(fidin, 'quit\n');

% Close the file
fclose(fidin);

% Prepare the XFOIL command
%xfoil_cmd = sprintf('xfoil.exe < %s', input_file);
xfoil_cmd = 'xfoil.exe < xfoil_input.txt > xfoiloutput';


% Display the command that will be executed
disp(['Running command: ' xfoil_cmd]);

% Run XFOIL with the input file
[status, result] = system(xfoil_cmd);

% Check if XFOIL ran successfully
if status ~= 0
    error('XFOIL execution failed with status %d: %s', status, result);
else
    disp('XFOIL simulation completed successfully');
end

% Optional: Delete the input file after execution
% delete(input_file);

% Read the output file

fidout = fopen(cp_output_file, 'r');
for i = 1:3 % 3 head lines before data
    fgetl(fidout); % Skip header lines
end
data = textscan(fidout, '%f %f %f'); % Read numeric data
fclose(fidout);
x_out = data{1};
y_out = data{2};
CP = data{3};

x_out = flip(x_out)'; % format: lower TE --> LE --> upper TE
y_out = flip(y_out)';
CP = flip(CP);

% Compute Lift and Drag Coefficients
%
cy = 0.0;
cx = 0.0;
n = nu+nl-2;
for i=1:n
dx = x_out(i+1) - x_out(i);
dy = y_out(i+1) - y_out(i);
cy = cy - CP(i,1) * dx;
cx = cx + CP(i,1) * dy;
end
%
% Print Lift and Drag coefficients on the screen
%
cl = cy * cos(alpha) - cx * sin(alpha)
cd = cy * sin(alpha) + cx * cos(alpha)
%

cpl = flip(CP(1:nl,1))'; % format LE --> TE
cpu = CP(nl:nl+nu-1,1)';% format LE --> TE
% 
% disp('test pour checker format Cp')
% disp(data)
% disp(length(CP))
% disp(CP)
% disp(length(x))
% disp(x_out)
% disp('yyyyyyy out')
% disp(y_out)
% cpl
% cpu
% 
% 
% disp('test format cpu, cpl etc...')
% disp('xu')
% disp(xu) % format: xu: 0-->1
% disp('xl')
% disp(xl)
% disp('lenght comparison')
% disp(length(xu))
% disp(length(cpu))

end

%% Xfoil1
function [cpu,cpl]=Xfoil1(nu,nl,xu,xl,yu,yl,alphad)
alpha = alphad* atan(1.0)/45.0; % convert deg --> rad
x = zeros(1,nl+nu-1);
y = zeros(1,nl+nu-1);
% Transfer supplied xu,xl etc to x and y arrays
for i = nl:nl+nu-1  
   l = i+1-nl;
   x(i) = xu(l);
   y(i) = yu(l);
end
    for i =1:nl 
      x(nl+1-i)=xl(i);
      y(nl+1-i)=yl(i);
    end
disp('a partir d ici')     
for i = 1:nu+nl-1
    fprintf('%14.8f %14.8f\n', x(i),y(i)); % format: TE lower  --> LE --> TE upper
end

% Create .dat file for the new airfoil
dat_file_name = 'airfoil_iteration.dat';
% Open file for writing
fid_dat = fopen(dat_file_name, 'w');
if fid_dat == -1
    error('Failed to create XFOIL input file');
end

% Write the dat.file
fprintf(fid_dat, 'airfoil iteration \n');
% Write coordinates
if length(x)==length(y)
    for i = 1:length(x)
        fprintf(fid_dat, '%.6f %.6f\n', x(i), y(i));% 6 digits
    end
else
    error('x and y dont have the same length')
end
fclose(fid_dat);

% Create input file for the new airfoil

% Define parameters
Re = 3000000;       % Reynolds number
alfa = alphad;           % AoA
iter_limit = 200;   % max iteration
airfoil_file = dat_file_name; % Airfoil filename
cp_output_file = 'cp_alfa_aoa2.txt'; % Output file for Cp data

% Open file for writing
fidin = fopen('xfoil_input.txt', 'w');
if fidin == -1
    error('Failed to create XFOIL input file');
end

% Write XFOIL commands to input file
fprintf(fidin, 'load %s\n', airfoil_file);  % Load airfoil coordinates
fprintf(fidin, 'oper\n');              % command
%fprintf(fid, 'visc\n');              % command visc
%fprintf(fid, '%d\n', Re);            % Reynolds number
fprintf(fidin, 'iter %d\n', iter_limit); % Set iteration limit
fprintf(fidin, 'alfa %.2f\n', alfa);   % Set angle of attack
fprintf(fidin, 'cpwr %s\n', cp_output_file); % Save pressure coefficient data
fprintf(fidin, '\n');
fprintf(fidin, 'quit\n');

% Close the file
fclose(fidin);

% Prepare the XFOIL command
xfoil_cmd = 'xfoil.exe < xfoil_input.txt > xfoiloutput';


% Display the command that will be executed
disp(['Running command: ' xfoil_cmd]);

% Run XFOIL with the input file
[status, result] = system(xfoil_cmd);

% Check if XFOIL ran successfully
if status ~= 0
    error('XFOIL execution failed with status %d: %s', status, result);
else
    disp('XFOIL simulation completed successfully');
end

% Read the output file

fidout = fopen(cp_output_file, 'r');
for i = 1:3 % 3 head lines before data
    fgetl(fidout); % Skip header lines
end
data = textscan(fidout, '%f %f %f'); % Read numeric data
fclose(fidout);
x_out = data{1};
y_out = data{2};
CP = data{3};

x_out = flip(x_out)'; % format: lower TE --> LE --> upper TE
y_out = flip(y_out)';
CP = flip(CP);


% Compute Lift and Drag Coefficients
%
cy = 0.0;
cx = 0.0;
n = nu+nl-2;
for i=1:n
dx = x(i+1) - x(i);
dy = y(i+1) - y(i);
cy = cy - CP(i,1) * dx;
cx = cx + CP(i,1) * dy;
end
%
% Print Lift and Drag coefficients on the screen
%
cl = cy * cos(alpha) - cx * sin(alpha)
cd = cy * sin(alpha) + cx * cos(alpha)
%
cpl = flip(CP(1:nl,1))'; % format LE --> TE
cpu = CP(nl:nl+nu-1,1)';% format LE --> TE

end

%% function mgm upper
function [dz]= mgmu(nu,Q,xu1)
% Apply Implicit Scheme
    % apply for lower and upper separately
    A1 = 3; B1 = 3; C1 = 3;
    B(1) = 1.0;
    C(1) = 0.0;
    B(nu)=1.0;
    Q(1)=0.0;
    dz(1)=0.0;
    dz(nu)=0.0;
    AM=zeros(nu,nu);
    AM(1,1)=1.0;
    AM(nu,nu)=1.0;
    for i =1:nu
        xu(i)=sqrt(xu1(i));
    end
    for ii=2:nu-1
        C(ii) = B1/(xu(ii+1)-xu(ii));
        C(ii) = C(ii) - C1/((xu(ii+1)-xu(ii-1))*(xu(ii+1)-xu(ii)));
        A(ii) = -C1/((xu(ii+1)-xu(ii-1))*(xu(ii)-xu(ii-1)));
        B(ii) = A1 - B1/(xu(ii+1)-xu(ii)) + C1/(xu(ii+1)-xu(ii-1))*...
            (1/(xu(ii+1)-xu(ii))+(1/(xu(ii)-xu(ii-1))));
        AM(ii,ii-1)=A(ii);
        AM(ii,ii)=B(ii);
        AM(ii,ii+1)=C(ii);
    end

    dz=AM\Q;

end

%% function mgml
function [dz]= mgml(nl,Q,xl1)
% Apply Implicit Scheme
    % apply for lower and upper separately
    A1 = 3; B1 = 3; C1 = 3;
    B(1) = 1.0;
    C(1) = 0.0;
    B(nl)=1.0;
    Q(1)=0.0;
    dz(1)=0.0;
    dz(nl)=0.0;
    AM=zeros(nl,nl);
    AM(1,1)=1.0;
    AM(nl,nl)=1.0;
    for i =1:nl
        x1(i)=sqrt(xl1(i));
    end
    for ii=2:nl-1
        C(ii) = B1/(x1(ii+1)-x1(ii));
        C(ii) = C(ii) - C1/((x1(ii+1)-x1(ii-1))*(x1(ii+1)-x1(ii)));
        A(ii) = -C1/((x1(ii+1)-x1(ii-1))*(x1(ii)-x1(ii-1)));
        B(ii) = A1 - B1/(x1(ii+1)-x1(ii)) + C1/(x1(ii+1)-x1(ii-1))*...
            (1/(x1(ii+1)-x1(ii))+(1/(x1(ii)-x1(ii-1))));
        AM(ii,ii-1)=A(ii);
        AM(ii,ii)=B(ii);
        AM(ii,ii+1)=C(ii);
       % M = A(ii)/B(ii-1);
       % B(ii) = B(ii) - M * C(ii-1);
       % Q(ii) = Q(ii) - M * Q(ii-1);
    end
%    dz(nu) = 0;
    dz=AM\Q;

end
