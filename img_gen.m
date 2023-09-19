format long g

global G M c r_s r_disk inclination h outer_temp;
G = 6.67430E-11; % Gravitational constant (before I got used to natural units, I shall keep it as is as a reminder of my learning journey)
c = 299792458; % Speed of light
M = 4.154E+6 * 1.98847E+30; % mass of BH (Sagittarius A* for reference)
r_s = 2*G*M/(c^2); %Schwarzschild radius
r_disk = 20*r_s; % radius of accretion disk
inclination =  17*pi/36; % inclination of accretion disk wrt observer
h = pi/1800; % angular interval (wrt centre of BH) for numerically simulating trajectory
outer_temp = 1000; %Temperature at outer edge of accretion disk

% Image settings
dist_from_BH = 200*r_s;
fov = pi/14; % field of view, how "zoomed in" the image is to the BH
vertical_pixels = 1080;
horizontal_pixels = 1920;

pixels = real(get_image(dist_from_BH, fov, horizontal_pixels, vertical_pixels));
pixels= rot90(pixels, 2);
save('output.mat', 'pixels');

%{
%-----------------visualization part-----------------%
% This plots the trajectory of a traced ray
[r, theta, dr, condition] = trajectory(100*r_s, atan(sqrt(451^2 + 41^2)/639 * tan(pi/28)));
x = r.*sin(theta)*cos(atan(21/226));
y = r.*sin(theta)*sin(atan(21/226));
z = r.*cos(theta);

phi = linspace(0, 2*pi, 50);
theta = linspace(0, pi, 50);
r2 = linspace(3*r_s, r_disk, 25);
figure(1)
[Phi, R] = meshgrid(phi, r2);
%Accretion disk
X1 = R.*cos(Phi)*cos(inclination);
Y1 = R.*sin(Phi);
Z1 = -R.*cos(Phi)*sin(inclination);

%Event horizon
[Theta, Phi] = meshgrid(theta, phi);
X2 = r_s*sin(Theta).*cos(Phi);
Y2 = r_s*sin(Theta).*sin(Phi);
Z2 = r_s*cos(Theta);


hold on
set(gca,'DataAspectRatio',[1 1 1])
line(x, y, z)
surf(X1,Y1,Z1, 'EdgeColor', 'none', 'FaceColor', '#4DBEEE');
surf(X2,Y2,Z2, 'EdgeColor', 'none', 'FaceColor', 'black');
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
hold off
view(45,10)
%--------------end of visualization part--------------%
%}



function pixels = get_image(distance, horizontal_fov, h_res, v_res)
    tan_phi = tan(horizontal_fov/2)/(h_res - 1);

    pixels = zeros(v_res, h_res, 2);
    for i=1:h_res/2
        tic
        fprintf('\n Now plotting column: %d', i)
        for s=1:min(i,v_res/2)
            firing_angle = atan(sqrt((2*i-1)^2 + (2*s-1)^2)*tan_phi);
            azimuthal_angle = atan(s/i);
            [r, ~, dr, condition] = trajectory(distance, firing_angle);
            
            %top right quadrant
            azimuthal1 = 3*pi/2 + azimuthal_angle;
            disk_angle1 = get_polar(azimuthal1);
            
            pixels(v_res/2+s, h_res/2+i, :) = get_color(r, dr, condition, azimuthal1, disk_angle1);
            
            %top left quadrant
            azimuthal2 = pi/2 - azimuthal_angle;
            disk_angle2 = get_polar(azimuthal2);
            
            pixels(v_res/2+s, h_res/2+1-i, :) = get_color(r, dr, condition, azimuthal2, disk_angle2);
            
            %bottom left quadrant
            azimuthal3 = pi/2 + azimuthal_angle;
            disk_angle3 = get_polar(azimuthal3);
            
            pixels(v_res/2+1-s, h_res/2+1-i, :) = get_color(r, dr, condition, azimuthal3, disk_angle3);
            
            %bottom right quadrant
            azimuthal4 = 3*pi/2 - azimuthal_angle;
            disk_angle4 = get_polar(azimuthal4);
            
            pixels(v_res/2+1-s, h_res/2+i, :) = get_color(r, dr, condition, azimuthal4, disk_angle4);
            
            if (i ~= s) && (i <= v_res/2)  %second set of radial symmetry points
                %top right quadrant
                azimuthal12 = 2*pi - azimuthal_angle;
                disk_angle12 = get_polar(azimuthal12);
                
                pixels(v_res/2+i, h_res/2+s, :) = get_color(r, dr, condition, azimuthal12, disk_angle12);
                
                %top left quadrant
                azimuthal22 = azimuthal_angle;
                disk_angle22 = get_polar(azimuthal22);
                
                pixels(v_res/2+i, h_res/2+1-s, :) = get_color(r, dr, condition, azimuthal22, disk_angle22);
                
                %bottom left quadrant
                azimuthal32 = pi - azimuthal_angle;
                disk_angle32 = get_polar(azimuthal32);
                
                pixels(v_res/2+1-i, h_res/2+1-s, :) = get_color(r, dr, condition, azimuthal32, disk_angle32);
                
                %bottom right quadrant
                azimuthal42 = pi + azimuthal_angle;
                disk_angle42 = get_polar(azimuthal42);
                
                pixels(v_res/2+1-i, h_res/2+s, :) = get_color(r, dr, condition, azimuthal42, disk_angle42);
            end
  
        end
        toc
    end
end

function color = get_color(r, dr, condition, azimuthal, disk_angle)
    global G M c r_s inclination r_disk h outer_temp;
    disk_index = disk_angle/h + 1;
    linear_approx = disk_index - fix(disk_index);
    lower_disk_index = fix(disk_index);
    upper_disk_index = lower_disk_index + 1;
    
    r1 = r(lower_disk_index);
    r2 = r(upper_disk_index);
    r_at_disk = (1-linear_approx)*r1 + linear_approx*r2;
    if r1==0 || r2==0
        if condition == 1      %Case when trajectory goes into event horizon
            color = [0, 0];
        elseif condition == 2      %Case when trajectory misses accretion disk and BH
            color = [0, -1];
        end
    elseif (r_at_disk > r_disk) || (r_at_disk < 3*r_s)
        if disk_angle <= pi   %check accretion disk behind BH
            lower_back_disk_index = lower_disk_index + pi/h;
            upper_back_disk_index = lower_back_disk_index + 1;
            r1 = r(lower_back_disk_index);
            r2 = r(upper_back_disk_index);
            r_back = (1-linear_approx)*r1 + linear_approx*r2;
            if r1==0 || r2==0
                if condition == 1      %Same checking for back disk
                    color = [0, 0];
                elseif condition == 2
                    color = [0, -1];
                end
            elseif r_back > r_disk
                color = [0, -1];
            elseif r_back < 3*r_s
                if length(r) > (2*pi/h + 1)
                    color = get_color(r(2*pi/h:end), dr(2*pi/h:end), condition, azimuthal, disk_angle);
                else
                    if condition == 1
                        color = [0, 0];
                    else
                        color = [0, -1];
                    end
                end
            else
                %generate color of disk at r_back
                disk_angle = disk_angle + pi ;
                dr_back = (1-linear_approx)*dr(lower_back_disk_index) + linear_approx*dr(upper_back_disk_index);
                v_back = sqrt(G*M/(r_back-r_s));
                cos_psi = ((sin(azimuthal)/sqrt(r_back^2 + dr_back^2)) * ...
                            ((r_back*cos(disk_angle) + dr_back*sin(disk_angle))*(1-cos(inclination))*cos(azimuthal) ...
                            + (-r_back*sin(disk_angle) + dr_back*cos(disk_angle))*sin(inclination)));
                
                redshift = (1+v_back*cos_psi/c)/(1-(v_back/c)^2) * 1/sqrt(1-r_s/r_back); 
                
                temp = outer_temp*((r_disk)^(3/4))*((1-sqrt(3*r_s/r_disk))^(-1/4)) ...
                                *(r_back^(-3/4))*((1-sqrt(3*r_s/r_back))^(1/4));
                
                color = [redshift, temp];
            end
        else
            if condition == 1
                color = [0, 0];
            else
                color = [0, -1];
            end
        end
    else
        %generate color of disk at r_at_disk
        
        dr_at_disk = (1-linear_approx)*dr(lower_disk_index) + linear_approx*dr(upper_disk_index);
        v_disk = sqrt(G*M/(r_at_disk-r_s));
        cos_psi = -((sin(azimuthal)/sqrt(r_at_disk^2 + dr_at_disk^2)) * ...
                   ((r_at_disk*cos(disk_angle) + dr_at_disk*sin(disk_angle))*(1-cos(inclination))*cos(azimuthal) ...
                   + (-r_at_disk*sin(disk_angle) + dr_at_disk*cos(disk_angle))*sin(inclination)));
    
        redshift = (1+v_disk*cos_psi/c)/(1-(v_disk/c)^2) * 1/sqrt(1-r_s/r_at_disk);
        
        temp = outer_temp*((r_disk)^(3/4))*((1-sqrt(3*r_s/r_disk))^(-1/4)) ...
                                *(r_at_disk^(-3/4))*((1-sqrt(3*r_s/r_at_disk))^(1/4));
                
        color = [redshift, temp];
    end
end


function [r, theta, dr, condition]=trajectory(r0, theta_in)
    global G M c r_s r_disk h;
    u_s = 1/r_s;
    u_disk = 1/r_disk;
    u0 = 1/r0;
    du0 = initial_du(r0, theta_in);
    theta0 = 0;
    condition = 0;
    
    intervals = (6*pi)/h;
    r = zeros(1, intervals);
    theta = zeros(1, intervals);
    dr = zeros(1, intervals);
    
    condition = 0;
    for i=1:intervals
        if u0 > 2/3 * u_s     %case when trajectory goes into event horizon
            condition = 1;
            break
        elseif  (theta0 >= pi/2) && (u0 < u_disk/2)  %case when trajectory misses accretion disk
            condition = 2;
            break
        end
        r(i) = 1/u0;
        theta(i) = h*i;
        du1 = partialu(du0, u0, h, G, M, c);
        u1 = get_u(du0, du1, u0, h);
        dr(i) = (-(1/u0^2)*du0 -(1/u1^2)*du1)/2;
        u0 = u1;
        du0 = du1;
        theta0 = h*i;
    end
end
    
function du0 = initial_du(r0, theta_in)
    du0 = 1/(r0 * tan(theta_in));
end

function du1 = partialu(du0, u0, h, G, M, c)
    du1 = du0 + h*(3*G*M*u0^2 / c^2 - u0);
end

function u1 = get_u(du0, du1, u0, h)
    u1 = u0 + h*(du0 + du1)/2;
end

function theta = get_polar(phi)
    global inclination
    if (phi >= pi/2) && (phi <= 3*pi/2)
        theta = abs(atan(1./(cos(phi)*tan(inclination))));
    else 
        theta = pi - atan(1./(cos(phi)*tan(inclination)));
    end
end
