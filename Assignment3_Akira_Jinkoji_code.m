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

% Display sample data
% disp('Sample Lower Surface Data (xlt, cplt):');
% disp([xlt(1:5), cplt(1:5)]);
% 
% disp('Sample Upper Surface Data (xut, cput):');
% disp([xut(1:5), cput(1:5)]);

% Plot target airfoil shape
figure(1); % airfoil shape
% plot(xut,-cput,'blue','LineWidth',3,'DisplayName','Target Upper');
% hold on
% plot(xlt,-cplt,'blue','LineWidth',3,'DisplayName','Target lower');
% Keep the Figure 1 open to plot x vs y for the final design 
%fclose(fid);
%
% Begin inverse design,  20 updates, starting with NACA 0012 at the same
% alpha
%
nitmax=100; % have to vary the AoA as the LE and TE coordinates would not change
for nit = 1:nitmax
    if(nit==1)
    % Open INITIAL airfoil data and calculate starting pressure distribution
    fid = fopen('NACA0012Big_AoA2.dat','r');
   [nu_ini, nl_ini, alpha_ini,xu1,xl1,yu1,yl1,cpu1,cpl1]=Panel(fid);
   %Plot initial airfoil too..
   plot(xu1,yu1,'black','LineWidth',3,'DisplayName','Starting Upper');
hold on
plot(xl1,yl1,'black','LineWidth',3,'DisplayName','Starting Lower');
   fclose(fid);
    else
   [cpu1,cpl1] =Panel1(nu_ini,nl_ini,xu1,xl1,yu1,yl1,alpha_ini);
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
% plot to see interpolation is ok
%plot(xu1,cpu1);
%hold on
%plot(xut,cpua);
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
%hold on
%plot(xl1,cpl1);
%hold on
%plot(xlt,cpla);
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
title('Comparison of the Cp distribution between the target and final version')
xlabel('x')
ylabel('-Cp')
legend


%% -------------------------------------------------------------------------

function [nu, nl, alphad, xu,xl,yu,yl,cpu,cpl] = Panel(fid)

%Read airfoil coordinates
alphad = fscanf(fid,'%f',1);    % Read Angle of Attack
alpha = alphad * pi /180;       % Convert to radians
nu = fscanf(fid,'%d',1);        % Read number of points on the upper side of airfoil
nl = fscanf(fid, '%d',1);       % Read number of points on the lower side of airfoil
isym = fscanf(fid,'%d',1);      % Read flag that states if this airfoil is symmetric
factor=fscanf(fid,'%f',1);      % Read scaling factor % <---- which scaling factor it is ?
if(isym>0)
  nl = nu;
end

x = zeros(1,nl+nu-1);
y = zeros(1,nl+nu-1);
xu = zeros(1,nu);
yu = zeros(1,nu);
xl = zeros(1,nl);
yl = zeros(1,nl);
for i = nl:nl+nu-1              % Read the points on the upper surface
   a1=fscanf(fid,'%f',1);
   b1 = fscanf(fid,'%f',1);
   x(i) = a1;
   y(i) = b1 * factor;
   l = i+1-nl;
   xu(l)=x(i);
   yu(l)=y(i);
end

if isym == 0                    % If the airfoil is not symmetric, read lower side ordinates too..
    for i = 1:nl
      a1=fscanf(fid, '%f',1);
      b1 = fscanf(fid, '%f', 1);
      x(nl+1-i) = a1;
      y(nl+1-i) = b1 * factor;
      xl(i)=x(nl+1-i);
      yl(i)=y(nl+1-i);
    end
else
    for i =1:nl 
       x(nl+1-i) = x(nl-1+i);
       y(nl+1-i) = - y(nl-1+i);
       xl(i)=x(nl+1-i);
       yl(i)=y(nl+1-i);
    end
end
for i = 1:nu+nl-1
    fprintf('%14.8f %14.8f\n', x(i),y(i));
end
CP=zeros(1,nu+nl-1);
%Plot airfoil in window #1
%plot(x,y);
% Assemble the Influence Coefficient Matrix A
n=nu+nl-2;
ds=zeros(1,n);
A=zeros(n+1,n+1);
%pi=4. * atan(1.0);  %%%I wonder if this is required in MATLAB?
 for i = 1:n
   t1= x(i+1)-x(i);
   t2 = y(i+1)-y(i);
   ds(i) = sqrt(t1*t1+t2*t2);
 end
for j = 1:n
 A(j,n+1) = 1.0;
 for i = 1:n
   if i == j
     A(i,i) = ds(i)/(2.*pi) *(log(0.5*ds(i)) - 1.0);
   else
     xm1 = 0.5 * (x(j)+x(j+1));
     ym1 = 0.5 * (y(j)+y(j+1));
     dx  = (x(i+1)-x(i))/ds(i);
     dy  = (y(i+1)-y(i))/ds(i);
     t1  = x(i) - xm1;
     t2  = y(i) - ym1;
     t3  = x(i+1) - xm1;
     t7  = y(i+1) - ym1;
     t4  = t1 * dx + t2 * dy;
     t5  = t3 * dx + t7 * dy;
     t6  = t2 * dx - t1 * dy;
     t1  = t5 * log(t5*t5+t6*t6) - t4 * log(t4*t4+t6*t6);
     t2  = atan2(t6,t4)-atan2(t6,t5);
     A(j,i) = (0.5 * t1-t5+t4+t6*t2)/(2.*pi);
   end
 end
A(n+1,1) = 1.0;
A(n+1,n) = 1.0;
end

% Assemble the Right hand Side of the Matrix system and X midpoints
rhs=zeros(n+1,1);
Xmid=zeros(n,1);
for i = 1:n
  Xmid(i,1) = 0.5 * (x(i) + x(i+1));
  ymid = 0.5 * (y(i) + y(i+1));
  rhs(i,1) = ymid * cos(alpha) - Xmid(i) * sin(alpha);
end

% Solve the system of equations
gamma = A\rhs;                  %This solves a*gamma=rhs to get gamma
CP=zeros(n+1,1);
CP_node=zeros(n+1,1);
% Calculate coefficients of pressure
for i = 1:n
CP(i,1) = 1. - gamma(i) * gamma(i);
end
% Interpolate from panel center to panel edge
for i=2:n-1
    CP_node(i)=(CP(i-1,1)*ds(i)+CP(i,1)*ds(i-1))/(ds(i)+ds(i-1));
end
CP_node(1)= 0.5*(CP(2,1)+CP(n-1,1));
CP_node(n)=CP_node(1);
for i=1:nu
    cpu(i)=CP_node(nl-1+i);
end
for i=1:nl
    cpl(i)= CP_node(nl+1-i);
end
%figure(2);
%plot(x,CP_node);
% Compute Lift and Drag Coefficients
%
cy = 0.0;
cx = 0.0;
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
end

%% panel1
function [cpu,cpl]=Panel1(nu,nl,xu,xl,yu,yl,alphad)
alpha = alphad* atan(1.0)/45.0;
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
       
for i = 1:nu+nl-1
    fprintf('%14.8f %14.8f\n', x(i),y(i));
end
CP=zeros(1,nu+nl-1);
%Plot airfoil in window #1
%plot(x,y);
% Assemble the Influence Coefficient Matrix A
n=nu+nl-2;
ds=zeros(1,n);
A=zeros(n+1,n+1);
%pi=4. * atan(1.0);  %%%I wonder if this is required in MATLAB?
 for i = 1:n
   t1= x(i+1)-x(i);
   t2 = y(i+1)-y(i);
   ds(i) = sqrt(t1*t1+t2*t2);
 end
for j = 1:n
 A(j,n+1) = 1.0;
 for i = 1:n
   if i == j
     A(i,i) = ds(i)/(2.*pi) *(log(0.5*ds(i)) - 1.0);
   else
     xm1 = 0.5 * (x(j)+x(j+1));
     ym1 = 0.5 * (y(j)+y(j+1));
     dx  = (x(i+1)-x(i))/ds(i);
     dy  = (y(i+1)-y(i))/ds(i);
     t1  = x(i) - xm1;
     t2  = y(i) - ym1;
     t3  = x(i+1) - xm1;
     t7  = y(i+1) - ym1;
     t4  = t1 * dx + t2 * dy;
     t5  = t3 * dx + t7 * dy;
     t6  = t2 * dx - t1 * dy;
     t1  = t5 * log(t5*t5+t6*t6) - t4 * log(t4*t4+t6*t6);
     t2  = atan2(t6,t4)-atan2(t6,t5);
     A(j,i) = (0.5 * t1-t5+t4+t6*t2)/(2.*pi);
   end
 end
A(n+1,1) = 1.0;
A(n+1,n) = 1.0;
end

% Assemble the Right hand Side of the Matrix system and X midpoints
rhs=zeros(n+1,1);
Xmid=zeros(n,1);
for i = 1:n
  Xmid(i,1) = 0.5 * (x(i) + x(i+1));
  ymid = 0.5 * (y(i) + y(i+1));
  rhs(i,1) = ymid * cos(alpha) - Xmid(i) * sin(alpha);
end

% Solve the system of equations
gamma = A\rhs;                  %This solves a*gamma=rhs to get gamma
CP=zeros(n+1,1);
CP_node=zeros(n+1,1);
% Calculate coefficients of pressure
for i = 1:n
CP(i,1) = 1. - gamma(i) * gamma(i);
end
% Interpolate from panel center to panel edge
for i=2:n-1
    CP_node(i)=(CP(i-1,1)*ds(i)+CP(i,1)*ds(i-1))/(ds(i)+ds(i-1));
end
CP_node(1)= 0.5*(CP(2,1)+CP(n-1,1));
CP_node(n)=CP_node(1);
for i=1:nu
    cpu(i)=CP_node(nl-1+i);
end
for i=1:nl
    cpl(i)= CP_node(nl+1-i);
end
%figure(2);
%plot(x,CP_node);
% Compute Lift and Drag Coefficients
%
cy = 0.0;
cx = 0.0;
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
end
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

%% function mgmu
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
