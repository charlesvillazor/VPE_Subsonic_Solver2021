%Charles Villazor
%April 28th, 2021
%MAE 456
%Final Project Part 2 - Airfoil Mach 0.01
clc;clear;close all

%% Constants
p_inf = 101325; %Atmospheric pressure, Freestream Pressure (Pa)
T_inf = 300; %Atmospheric Temperature, Freestream Temperature (Kelvin)
rho_inf = 1.225; %Atmospheric Air Density (kg/m^3)
gamma = 1.4; %Ratio of specific heats for air
r = 287.05; %Specific Gas Constant for air(J/kg*K)
a_inf = sqrt(gamma*r*T_inf); %Freestream speed of sound (m/s)

%% Importing Airfoil Data
%All Airfoil stuff is done before doing any airfoil stuff
file = fopen('airfoilsmall.txt'); %opening file
fgetl(file); %read and skip first line
data = cell2mat(textscan(file, '%f %f'));
Nx = 357; %number of terms in x direction
Ny = 80; %number of terms in y direction
x = zeros(Nx,Ny); %initialize x
y = zeros(Nx,Ny); %initialize y
k = 1; %initialize variable for iteration
for i = 1:Ny
    for j = 1:Nx
        x(j,i) = data(k,1);
        y(j,i) = data(k,2);
        k = k+1; %add one to every iteration (why does this work?)
    end
end
% %% Plotting Airfoil Data
% figure(1)
% plot(data(:,1),data(:,2),'.')
% axis equal

%% Airfoil Del Array
%del = zeros(N,N,2,2);
% For N1*A1
for j = 2:Ny %interior in j direction
    for i = 1:Nx %faces in i direction
        del(i,j,1,1) = (y(i,j)-y(i,j-1)); %x component
        del(i,j,1,2) = -(x(i,j)-x(i,j-1)); %y component
    end
end
% For N2*A2
for i = 2:Nx %interior in i direction
    for j = 1:Ny %faces in j direction
        del(i,j,2,1) = -(y(i,j)-y(i-1,j)); %x component
        del(i,j,2,2) = x(i,j)-x(i-1,j); %y component
    end
end

%% Airfoil Vol Array Routine / Cell Centers Routine 
for i = 2:Nx %For area of top half
    for j = 2:Ny
        %initializing Verticies
        xa = x(i,j-1);
        ya = y(i,j-1);
        xb = x(i,j);
        yb = y(i,j);
        xc = x(i-1,j);
        yc = y(i-1,j);
        xd = x(i-1,j-1);
        yd = y(i-1,j-1);
        %Finding bottom half triangle lengths
        da = sqrt(((xb-xa)^2)+((yb-ya)^2)); %right side length
        db = sqrt(((xd-xb)^2)+((yd-yb)^2)); %diagonal
        dc = sqrt(((xa-xd)^2)+((ya-yd)^2)); %bottom length
        %Finding top half triangle lengths
        daa = sqrt(((xc-xb)^2)+((yc-yb)^2)); %top length
        dbb = sqrt(((xd-xb)^2)+((yd-yb)^2)); %diagonal
        dcc = sqrt(((xc-xd)^2)+((yc-yd)^2)); %left side length
        %Finding semi perimeters
        s = (da+db+dc)/2;
        ss = (daa+dbb+dcc)/2;
        %Finding areas and inserting into vol array
        AT1 = sqrt(s*(s-da)*(s-db)*(s-dc));
        AT2 = sqrt(ss*(ss-daa)*(ss-dbb)*(ss-dcc));
        vol(i,j) = AT1+AT2;
        %Averaging x-coordinates to get x center
        xcenter(i,j) = (xa+xb+xc+xd)/4;
        %Averaging y-coordinates to get y center
        ycenter(i,j) = (ya+yb+yc+yd)/4;
    end 
end

%% Interpolating Airfoil Cell Centers
for i = 2:Nx
    xcenter(i,1) = 2.0*xcenter(i,2)-xcenter(i,3);
    xcenter(i,Ny+1) = 2.0*xcenter(i,Ny)-xcenter(i,Ny-1);
    ycenter(i,1) = 2.0*ycenter(i,2)-ycenter(i,3);
    ycenter(i,Ny+1) = 2.0*ycenter(i,Ny)-ycenter(i,Ny-1);
end
for j=1:Ny+1
    xcenter(1,j) = 2.0*xcenter(2,j)-xcenter(3,j);
    xcenter(Nx+1,j) = 2.0*xcenter(Nx,j)-xcenter(Nx-1,j);
    ycenter(1,j) = 2.0*ycenter(2,j)-ycenter(3,j);
    ycenter(Nx+1,j) = 2.0*ycenter(Nx,j)-ycenter(Nx-1,j);
end 

%% Plotting Volume Contour for Airfoils
figure(1)
contourf(xcenter(2:Nx,2:Ny),ycenter(2:Nx,2:Ny),vol(2:Nx,2:Ny),100,'edgecolor','none')
%plot(xcenter,ycenter,'.')
colorbar 
axis([-2 3 -2 2])
title('Airfoil Contour Plot of Volume')
%% Airfoil Velocity Potential (Mach 0.01)
mach = 0.01; %mach number at freestream
u = mach; %x direction velocity
v = 0; %y direction velocity
u_inf = u*a_inf; %Freestream u velocity
v_inf = v*a_inf; %Freestream v velocity
for j = 1:Ny+1
    for i = 1:Nx+1
        phi(i,j) = u_inf*xcenter(i,j)+v_inf*ycenter(i,j);
    end
end

%% Boundary Conditions for Velocity Potential (Mach 0.01)
%far field
i = 1;
for j = 1:Ny+1
    phi(i,j) = u_inf*xcenter(i,j)+v_inf*ycenter(i,j);
end

i = Nx;         %far field
for j = 1:Ny+1
    phi(i+1,j) = u_inf*xcenter(i+1,j)+v_inf*ycenter(i+1,j);
end

j = 1;
for i = 2:100     %periodic boundary conditions in wake
    iq = i-1;
    phi(i,j) = phi(Nx-iq+1,j+1);
    phi(Nx-iq+1,j) = phi(i,j+1);
end

for i=101:258  %wall boundary conditions
    phi(i,1) = phi(i,j+1);
end

j = Ny;        %far field
for i = 1:Nx+1
    phi(i,Ny+1)=u_inf*xcenter(i,j+1)+v_inf*ycenter(i,j+1);
end

%% Alternate Airfoil Grad (Mach 0.01)
%Taken from website
for j = 2:Ny
    for i = 2:Nx
        dx1 = 0.5*(x(i,j)+x(i,j-1))-xcenter(i,j);
        dy1 = 0.5*(y(i,j)+y(i,j-1))-ycenter(i,j);
        dd1 = sqrt(dx1^2+dy1^2);
        ddd = sqrt((xcenter(i+1,j)-xcenter(i,j))^2 + (ycenter(i+1,j)-ycenter(i,j))^2);
        alp1 = dd1/ddd;
        phi1 = phi(i+1,j)*alp1 + phi(i,j)*(1.-alp1);

        dx1 = 0.5*(x(i-1,j)+x(i-1,j-1))-xcenter(i,j);
        dy1 = 0.5*(y(i-1,j)+y(i-1,j-1))-ycenter(i,j);
        dd1 = sqrt(dx1^2+dy1^2);
        ddd = sqrt((xcenter(i-1,j)-xcenter(i,j))^2 + (ycenter(i-1,j)-ycenter(i,j))^2);
        alp1 = dd1/ddd;
        phi2 = phi(i-1,j)*alp1 + phi(i,j)*(1.-alp1);

        dx1 = 0.5*(x(i,j)+x(i-1,j))-xcenter(i,j);
        dy1 = 0.5*(y(i,j)+y(i-1,j))-ycenter(i,j);
        dd1 = sqrt(dx1^2+dy1^2);
        ddd = sqrt((xcenter(i,j+1)-xcenter(i,j))^2 + (ycenter(i,j+1)-ycenter(i,j))^2);
        alp1 = dd1/ddd;
        phi3 = phi(i,j+1)*alp1 + phi(i,j)*(1.-alp1);

        dx1 = 0.5*(x(i,j-1)+x(i-1,j-1))-xcenter(i,j);
        dy1 = 0.5*(y(i,j-1)+y(i-1,j-1))-ycenter(i,j);
        dd1 = sqrt(dx1^2+dy1^2);
        ddd = sqrt((xcenter(i,j-1)-xcenter(i,j))^2 + (ycenter(i,j-1)-ycenter(i,j))^2);
        alp1 = dd1/ddd;
        phi4 = phi(i,j-1)*alp1 + phi(i,j)*(1.-alp1);

        total(:) = phi1*del(i,j,1,:)-phi2*del(i-1,j,1,:)+ phi3*del(i,j,2,:)-phi4*del(i,j-1,2,:);
        grad(i,j,:) = total(:)/vol(i,j); %grad(i,j,1) is u component, grad(i,j,2) is v component
    end
end
%grad(:,:,1) = -grad(:,:,1); %Flip U variable signs

%% Boundary Conditions for Grad Array
%taken from website
i=1;
for j=1:Ny+1
    grad(i,j,1) = u_inf;
    grad(i,j,2) = v_inf;
end

i = Nx;
for j = 1:Ny+1
    grad(i+1,j,1) = u_inf;
    grad(i+1,j,2) = v_inf;
end

j = 1;
for i=2:100
    iq = i-1;
    grad(i,j,:) = grad(Nx-iq+1,j+1,:);
    grad(Nx-iq+1,j,:) = grad(i,j+1,:);
end

j = 1;
for i=101:258
    area = sqrt((del(i,j,2,1)^2) + (del(i,j,2,2)^2));
    xnx = del(i,j,2,1)/area;
    xny = del(i,j,2,2)/area;
    vdotn = grad(i,j+1,1)*xnx + grad(i,j+1,2)*xny;
    grad(i,j,1) = grad(i,j+1,1) - 2.0*vdotn*xnx;
    grad(i,j,2) = grad(i,j+1,2) - 2.0*vdotn*xny;
end

j = Ny;
for i = 1:Nx+1
    grad(i,j+1,1) = u_inf;
    grad(i,j+1,2) = v_inf;
end

%% Contour plot of u component of Airfoil velocity (Mach 0.01)
figure(2)
%contourf(xcenter(5:Nx-5,5:Ny-5),ycenter(5:Nx-5,5:Ny-5),grad(5:Nx-5,5:Ny-5,1),100,'edgecolor','none');
contourf(xcenter(2:Nx,2:Ny),ycenter(2:Nx,2:Ny),grad(2:Nx,2:Ny,1),100,'edgecolor','none');
colorbar;
axis equal
title('Airfoil Contour Plot of u component of Velocity (Mach 0.01)')
xlabel('X Distance')
ylabel('Y Distance')

%% Contour plot of v component of Airfoil velocity (Mach 0.01)
figure(3)
contourf(xcenter(2:Nx,2:Ny),ycenter(2:Nx,2:Ny),-grad(2:Nx,2:Ny,2),100,'edgecolor','none');
colorbar;
axis equal
title('Airfoil Contour Plot of v component of Velocity (Mach 0.01)')
xlabel('X Distance')
ylabel('Y Distance')

%% Calculating Airfoil Rho (Mach 0.01)
for i = 1:Nx+1
    for j = 1:Ny+1
        rho(i,j) = rho_inf*(1+((gamma-1)/(2*(a_inf^2)))*((u_inf^2)-(((grad(i,j,1)^2)+(grad(i,j,2)^2)))))^(1/(gamma-1));
        %rho(i,j) = rho_inf;
    end
end

%% Iteration code (Mach 0.01)
tStart = tic; %measure the time it takes for iteration, start time
kmax = 10000; %Maximum amount of iterations
tolerance = 1e-4; %tolerance for resnorm
mrelax = 7; %number of relaxation steps
omega = 0.7; 
dphi = zeros(Nx+1,Ny+1); %initialize delta phi
% Calculating the residual and a array
for k = 1:kmax
    %Normal Tangential Gradient Model
    %For Faces 1 and 3 (i direction)
    for j = 2:Ny
        for i = 1:Nx
            rhoave = 0.5*(rho(i,j)+rho(i+1,j)); %Average Cell Density
            volave = 0.5*(vol(max(i,2),j)+vol(min(i+1,Nx),j)); %Average Cell Volume
            area2 = (del(i,j,1,1)^2)+(del(i,j,1,2)^2); %face area squared
            
            dtx = xcenter(i+1,j)-xcenter(i,j); %define x direction between cell centers
            dty = ycenter(i+1,j)-ycenter(i,j); %define y direction between cell centers
            dtm = sqrt((dtx^2)+(dty^2)); %Magnitude of vector between cell centers
            dtx = dtx/dtm; %unit vector for dtx
            dty = dty/dtm; %unit vector for dty
            uave = 0.5*(grad(i,j,1))+(grad(i+1,j,1)); %averaging u velocity
            
            vave = 0.5*(grad(i,j,2)+grad(i+1,j,2)); %averaging v velocity
            udott = (uave*dtx)+(vave*dty); %dot product of dtx and dty
            
            ddphi = phi(i+1,j)-phi(i,j);
            
            gradf(1) = (ddphi*(dtx/dtm))+(uave-(udott*dtx));
            gradf(2) = (ddphi*(dty/dtm))+(vave-(udott*dty));
            vdotna = (gradf(1)*del(i,j,1,1))+(gradf(2)*del(i,j,1,2));
            
            flux(i) = rhoave*vdotna;
            coeff(i) = (rhoave*area2)/volave; 
            
        end
        for i = 2:Nx
            res(i,j) = flux(i)-flux(i-1); %Add up first two parts of residual
            a(i,j,1) = coeff(i-1);
            a(i,j,5) = coeff(i);
            a(i,j,3) = -(coeff(i)+coeff(i-1));
        end
    end
    % For Faces 2 and 4 (j faces)
    for i = 2:Nx
        for j = 1:Ny
            rhoave = 0.5*(rho(i,j)+rho(i,j+1)); %Average Cell Density
            volave = 0.5*(vol(i,max(2,j))+vol(i,min(j+1,Ny))); %Average Cell Volume
            area2 = (del(i,j,2,1)^2)+(del(i,j,2,2)^2); %face area squared
            
            dtx = xcenter(i,j+1)-xcenter(i,j); %define x direction between cell centers
            dty = ycenter(i,j+1)-ycenter(i,j); %define y direction between cell centers
            dtm = sqrt((dtx^2)+(dty^2)); %Magnitude of vector between cell centers
            dtx = dtx/dtm; %unit vector for dtx
            dty = dty/dtm; %unit vector for dty
            uave = 0.5*(abs(grad(i,j,1))+abs(grad(i,j+1,1))); %averaging u velocity
            vave = 0.5*(grad(i,j,2)+grad(i,j+1,2)); %averaging v velocity
            udott = (uave*dtx)+(vave*dty); %dot product of dtx and dty
            
            ddphi = phi(i,j+1)-phi(i,j);
            
            gradf(1) = (ddphi*(dtx/dtm))+(uave-(udott*dtx));
            gradf(2) = (ddphi*(dty/dtm))+(vave-(udott*dty));
            vdotna = (gradf(1)*del(i,j,2,1))+(gradf(2)*del(i,j,2,2));
            
            flux(j) = rhoave*vdotna;
            coeff(j) = (rhoave*area2)/volave;
            
            if j == 1 && i > 100 && i < 258
                flux(j) = 0;
            end
        end
        for j = 2:Ny
            res(i,j) = res(i,j)+(flux(j)-flux(j-1)); %Add up first two parts to residual
            a(i,j,2) = coeff(j-1);
            a(i,j,4) = coeff(j);
            a(i,j,3) = a(i,j,3)-(coeff(j)+coeff(j-1));
        end
    end

%     %Normal Gradient Model (Old Code)
%     %For Faces 1 and 3
%     for j = 2:Ny
%         for i = 1:Nx
%             rhoave = 0.5*(rho(i,j)+rho(i+1,j)); %Average Cell Density
%             volave = 0.5*(vol(max(i,2),j)+vol(min(i+1,Nx),j)); %Average Cell Volume
%             area2 = (del(i,j,1,1)^2)+(del(i,j,1,2)^2); %face area squared
%             coeff(i) = (rhoave*area2)/volave; %Coefficient used for a(i,j,1:5)
%             flux(i) = coeff(i)*(phi(i+1,j)-phi(i,j)); %Mass flow rate through face
%         end
%         for i = 2:Nx
%             res(i,j) = flux(i)-flux(i-1); %Add up first two parts of residual
%             a(i,j,1) = coeff(i-1);
%             a(i,j,5) = coeff(i);
%             a(i,j,3) = -(coeff(i)+coeff(i-1));
%         end
%     end
%     % For Faces 2 and 4
%     for i = 2:Nx
%         for j = 1:Ny
%             rhoave = 0.5*(rho(i,j)+rho(i,j+1)); %Average Cell Density
%             volave = 0.5*(vol(1,max(j,2))+vol(i,min(j+1,Ny))); %Average Cell Volume
%             area2 = (del(i,j,2,1)^2)+(del(i,j,2,2)^2); %face area squared
%             coeff(j) = (rhoave*area2)/volave; %Coefficient used for a(i,j,1:5)
%             flux(j) = coeff(j)*(phi(i,j+1)-phi(i,j)); %Mass flow rate through face
%         end
%         for j = 2:Ny
%             res(i,j) = res(i,j)+flux(j)-flux(j-1); %Add up first two parts to residual
%             a(i,j,2) = coeff(j-1);
%             a(i,j,4) = coeff(j);
%             a(i,j,3) = a(i,j,3)-(coeff(j)+coeff(j-1));
%         end
%     end

   % Checking for convergence of residual and tolerance
   resnorm = 0;
    for j = 1:Ny
        for i = 1:Nx
            resnorm = resnorm +res(i,j)^2;
        end
    end
    
    %For the first iteration (k = 1)
    if k == 1
        resnorm0 = resnorm; %first resnorm
    end
    test(k) = resnorm/resnorm0;
    %k = k;
    if test(k) < tolerance %if the ratio of the resnorms is below the tolerance, it equals zero
        break
    else %if not, update phi, boundary counditions for phi, grad, boundary conditions for grad, rho array
        
        %Updating Phi
        %dphi = zeros(Nx+1,Ny+1); %initialize delta phi
        for m = 1:mrelax
            for i = 2:Nx
                for j = 2:Ny
                    dphitmp = dphi(i,j); %place holder for old solution
                    dphigs = (-res(i,j)-a(i,j,1)*dphi(i-1,j)-a(i,j,2)*dphi(i,j-1)-a(i,j,4)*dphi(i,j+1)-...
                        a(i,j,5)*dphi(i+1,j))/a(i,j,3); %Gauss-Seidal Step
                    dphi(i,j) = dphitmp +omega*(dphigs-dphitmp); %SOR Step
                end
            end
        end
        for j = 2:Ny
            for i = 2:Nx
                phi(i,j) = phi(i,j)+dphi(i,j);
            end
        end
        
        %Updating Boundary Conditions for Phi
        %far field
        i = 1;
        for j = 1:Ny+1
            phi(i,j) = u_inf*xcenter(i,j)+v_inf*ycenter(i,j);
        end

        i = Nx;         %far field
        for j = 1:Ny+1
            phi(i+1,j) = u_inf*xcenter(i+1,j)+v_inf*ycenter(i+1,j);
        end

        j = 1;
        for i = 2:100     %periodic boundary conditions in wake
            iq = i-1;
            phi(i,j) = phi(Nx-iq+1,j+1);
            phi(Nx-iq+1,j) = phi(i,j+1);
        end

        for i=101:258  %wall boundary conditions
            phi(i,1) = phi(i,2);
        end

        j = Ny;        %far field
        for i = 1:Nx+1
            phi(i,j+1)=u_inf*xcenter(i,j+1)+v_inf*ycenter(i,j+1);
        end
        
        %Calculate Grad Array in Interior (New Code)
        for j = 2:Ny
            for i = 2:Nx
                dx1 = 0.5*(x(i,j)+x(i,j-1))-xcenter(i,j);
                dy1 = 0.5*(y(i,j)+y(i,j-1))-ycenter(i,j);
                dd1 = sqrt((dx1^2)+(dy1^2));
                ddd = sqrt((xcenter(i+1,j)-xcenter(i,j))^2 + (ycenter(i+1,j)-ycenter(i,j))^2);
                alp1 = dd1/ddd;
                phi1 = phi(i+1,j)*alp1 + phi(i,j)*(1.-alp1);

                dx1 = 0.5*(x(i-1,j)+x(i-1,j-1))-xcenter(i,j);
                dy1 = 0.5*(y(i-1,j)+y(i-1,j-1))-ycenter(i,j);
                dd1 = sqrt((dx1^2)+(dy1^2));
                ddd = sqrt((xcenter(i-1,j)-xcenter(i,j))^2 + (ycenter(i-1,j)-ycenter(i,j))^2);
                alp1 = dd1/ddd;
                phi2 = phi(i-1,j)*alp1 + phi(i,j)*(1.-alp1);

                dx1 = 0.5*(x(i,j)+x(i-1,j))-xcenter(i,j);
                dy1 = 0.5*(y(i,j)+y(i-1,j))-ycenter(i,j);
                dd1 = sqrt((dx1^2)+(dy1^2));
                ddd = sqrt((xcenter(i,j+1)-xcenter(i,j))^2 + (ycenter(i,j+1)-ycenter(i,j))^2);
                alp1 = dd1/ddd;
                phi3 = phi(i,j+1)*alp1 + phi(i,j)*(1.-alp1);

                dx1 = 0.5*(x(i,j-1)+x(i-1,j-1))-xcenter(i,j);
                dy1 = 0.5*(y(i,j-1)+y(i-1,j-1))-ycenter(i,j);
                dd1 = sqrt((dx1^2)+(dy1^2));
                ddd = sqrt((xcenter(i,j-1)-xcenter(i,j))^2 + (ycenter(i,j-1)-ycenter(i,j))^2);
                alp1 = dd1/ddd;
                phi4 = phi(i,j-1)*alp1 + phi(i,j)*(1.-alp1);

                total(:) = phi1*del(i,j,1,:)-phi2*del(i-1,j,1,:)+ phi3*del(i,j,2,:)-phi4*del(i,j-1,2,:);
                grad(i,j,:) = total(:)/vol(i,j); %grad(i,j,1) is u component, grad(i,j,2) is v component
            end
        end
        %grad(:,:,1) = abs(grad(:,:,1)); %Flip U variable signs

%         % Calculating Grad Array (Old Code)
%         for j = 2:Ny
%             for i = 2:Nx
%                 phi1 = 0.5*(phi(i,j)+phi(i+1,j));
%                 phi2 = 0.5*(phi(i,j)+phi(i,j+1));
%                 phi3 = 0.5*(phi(i,j)+phi(i-1,j));
%                 phi4 = 0.5*(phi(i,j)+phi(i,j-1));
%                 total(1:2) = ((phi1*del(i,j,1,1:2))-(phi3*del(i-1,j,1,1:2))+(phi2*del(i,j,2,1:2))-(phi4*del(i,j-1,2,1:2)));
%                 grad(i,j,:) = total(:)./(vol(i,j));
%             end
%         end

        %Calculating Grad Boundary Conditions
        i=1;
        for j=1:Ny+1
            grad(i,j,1) = u_inf;
            grad(i,j,2) = v_inf;
        end

        i = Nx;
        for j = 1:Ny+1
            grad(i+1,j,1) = u_inf;
            grad(i+1,j,2) = v_inf;
        end

        j = 1;
        for i=2:100
            iq = i-1;
            grad(i,j,:) = grad(Nx-iq+1,j+1,:);
            grad(Nx-iq+1,j,:) = grad(i,j+1,:);
        end

        j = 1;
        for i=101:258
            area = sqrt((del(i,j,2,1)^2) + (del(i,j,2,2)^2));
            xnx = del(i,j,2,1)/area;
            xny = del(i,j,2,2)/area;
            vdotn = grad(i,j+1,1)*xnx + grad(i,j+1,2)*xny;
            grad(i,j,1) = grad(i,j+1,1) - 2.0*vdotn*xnx;
            grad(i,j,2) = grad(i,j+1,2) - 2.0*vdotn*xny;
        end

        j = Ny;
        for i = 1:Nx+1
            grad(i,j+1,1) = u_inf;
            grad(i,j+1,2) = v_inf;
        end
        
        %Calculating Density (rho)
        for i = 1:Nx+1
            for j = 1:Ny+1
                rho(i,j) = rho_inf*((1+((gamma-1)/(2*(a_inf^2)))*((u_inf^2)-(((grad(i,j,1)^2)+(grad(i,j,2)^2)))))^(1/(gamma-1)));
                %rho(i,j) = rho_inf;
            end
        end
        
    end
end
tEnd = toc(tStart); %Measure end time of iteration
%% Resnorm/Resnorm0 vs Iteration Plot
figure(4)
semilogy(1:1:k,test)
%plot(1:1:k,log(test))
title('Resnorm/Resnorm_0 vs Iteration Plot')
ylabel('Resnorm/Resnorm_0')
xlabel('Number of Iterations')

%% Contour plot of u component of Airfoil velocity (Mach 0.01)
figure(5)
%contourf(xcenter(5:Nx-5,5:Ny-5),ycenter(2:Nx,5:Ny-5),-grad(5:Nx-5,5:Ny-5,1),100,'edgecolor','none');
contourf(xcenter(2:Nx,2:Ny),ycenter(2:Nx,2:Ny),real(grad(2:Nx,2:Ny,1)),50,'edgecolor','none');
colorbar;
axis equal
title('Airfoil Contour Plot of u component of Velocity (Mach 0.01)(After Iterating)')
xlabel('X Distance')
ylabel('Y Distance')

%% Contour plot of v component of Airfoil velocity (Mach 0.01)
figure(6)
contourf(xcenter(2:Nx,2:Ny),ycenter(2:Nx,2:Ny),real(grad(2:Nx,2:Ny,2)),50,'edgecolor','none');
colorbar;
axis equal
title('Airfoil Contour Plot of v component of Velocity (Mach 0.01)(After Iterating)')
xlabel('X Distance')
ylabel('Y Distance')

%% Calculating Pressure (Mach 0.01)
for i = 1:Nx+1
    for j = 1:Ny+1
        p(i,j) = rho(i,j)*r*T_inf; %Pressure array (from ideal gas law)
    end
end

%% Calculating Pressure Coefficient (Mach 0.01)
for i = 1:Nx+1
    for j = 1:Ny+1
        cp(i,j) = (p(i,j)-p_inf)/(0.5*rho_inf*(u_inf^2)); %coefficient of pressure
    end
end

%% Coefficient of Pressure vs x distance plot (Along Airfoil)
figure(7)
plot((xcenter(101:258,2)/max(xcenter(101:258,2))),cp(101:258,2));
title('Coefficient of Pressure vs x/c Plot - Along Airfoil (Mach 0.01)')
ylabel('C_p')
xlabel('x/c')

% %% Coefficient of Pressure vs x distance plot (Bottom Wall)
% figure(8)
% plot((xcenter(2:Nx,Ny)/max(xcenter(2:Nx,Ny))),cp(2:Nx,Ny));
% title('Coefficient of Pressure vs X Distance Plot - Bottom Wall (Mach 0.5)')
% ylabel('Coefficient of Pressure')
% xlabel('X Distance')

%% Calculating Magnitude of velocity and converting to mach number
%Calculating velocity magnitude
for i = 1:Nx+1
    for j = 1:Ny+1
        velocity(i,j) = sqrt((grad(i,j,1)^2)+(grad(i,j,2)^2));
    end
end
%Converting to mach number
for i = 1:Nx+1
    for j = 1:Ny+1
        mach(i,j) = velocity(i,j)/a_inf;
    end
end

%% Mach Number along the airfoil (j = 2)
figure(9)
plot((x(101:258,2)/max(x(101:258,2))),mach(101:258,2))
title('Mach Number Distribution along Airfoil (Mach 0.01)')
xlabel('x/c')
ylabel('Mach Number')

%% Contour of Magnitude of Velocity
figure(10)
contourf(xcenter(2:Nx,2:Ny),ycenter(2:Nx,2:Ny),real(velocity(2:Nx,2:Ny)),50,'edgecolor','none');
colorbar;
axis equal
title('Airfoil Contour Plot of Magnitude of Velocity (Mach 0.01)(After Iterating)')
xlabel('X Distance')
ylabel('Y Distance')

% %% Mach Number along the bottom slip wall (j = Ny)
% figure(10)
% plot((xcenter(2:Nx,Ny)/max(xcenter(2:Nx,Ny))),mach(2:Nx,Ny))
% title('Mach Number Distribution along bottom Slip Wall')
% xlabel('x/c distance')
% ylabel('Mach Number')