%%
% Derek Woodard, 100827083

%% Part 1

clear all
close all
set(0,'DefaultFigureWindowStyle','docked');

m0 = 9.10938356e-31;          % rest mass of an electron
mn = 0.26*m0;                 % effective mass of electrons
T = 300;                      % temperature
k = 1.38064852e-23;           % Boltzmann's Constant
Tmn = 0.2e-12;                % mean time between collisions
q = 1.60217662e-19;           % The charge of an electron
vth = sqrt((2*k*T)/mn);       % calculate thermal velocity for 1.1 with two degrees of freedom
mfp = vth * Tmn;              % calculate the mean free path for 1.2


%% The nominal size of the region is 200nm x 100nm
H = 100e-9;
W = 200e-9;

pop_tot = 5000;                  % use this variable to determine how many electrons are modeled
pop_vis = 5;                     % use this variable to change how many electrons are plotted
ts = H/vth/100;                  % The time step is set to make each plot point have an electron move 1/100th of the height of the area
iter = 500;                      % how many plot points total
Vx = 0.8;                        % voltage across the x direction is set to 0.8V
Vy = 0;                          % voltage across the y dimension not provided, so assume 0
m0 = 9.10938215e-31;             % mass of an electron
density = 1e19;                  % density of electrons
ps = 1 - exp(-((ts)/(Tmn)));     % The eletrons have a probability of scattering

animate = 1;                     % set to 1 to see each iteration
% animate = 0;                     % set to 0 to only see final product
    
%%
% part 1a) determine the electric field on the electrons
Ex = Vx/W;
Ey = Vy/H;
% E = sqrt(Ex^2 + Ey^2);
% given that no voltage is applied in the y-direction, the x-direction
% provides the full electric field for the semi-conductor.
% The electric field is: 5e5 V/m

%%
% part 1b) determine the force applied to each electron
Fx = q*Ex;
Fy = q*Ey;
% Using the electric field determined previously and the charge of an
% electron, the force applied is: -8.0109e-14 N

%%
% part 1c) determine the acceleration of the electrons
ax = Fx/m0;
ay = Fy/m0;
% Using the force found previously along with the mass of an electron the
% acceleration is: -8.7941e16 m/s^2

%%
% In order to track the electrons, the following arrays will be used
% The arrays will store the positions, velocities, and temperatures of the
% elctrons
state = zeros(pop_tot, 4);
% The rows of state represent each electron
% The columns represent the position and velocity in each direction as shown
% [x y vx vy]

traj = zeros(iter, pop_vis*2);
% The trajectory of the electrons we want to plot will be stored here

temps = zeros(iter,1);
% The temperature will be tracked in this array

J = zeros(iter,2);
% The current density is made to have an x and y column

%% 
% Start by generating the initial eletrons with a constant speed of vth
for i = 1:pop_tot
    angle = rand*2*pi;
    state(i,:) = [W*rand H*rand vth*cos(angle) vth*sin(angle)];
end


%%
% Now we iterate to update the electron positions and plot the trajectories
for i = 1:iter
    state(:,1:2) = state(:,1:2) + ts.*state(:,3:4);
    % update the state of each electron using the speed
    
    state(:,3) = state(:,3) + ts*ax;
    state(:,4) = state(:,4) + ts*ay;
    % The speed of each electron is updated using the acceleration from the
    % electric field applied
    % Note that since the elctric field is only applied in the x-direction,
    % no acceleration occurs in the y-direction
    
    % We need to ensure the boundary conditions are checked
    bcx = state(:,1) > W;
    state(bcx,1) = state(bcx,1) - W;
    % Ensure the electron is transported from the right side to the left if
    % travelling past the right X limit
    
    bcx = state(:,1) < 0;
    state(bcx,1) = state(bcx,1) + W;
    % Ensure the electron is transported from the left side to the right if
    % travelling past the left X limit
    
    bcy = state(:,2) > H;
    state(bcy,2) = 2*H - state(bcy,2);
    state(bcy,4) = -state(bcy,4);
    % Ensure the electrons bounce back down if going above high the Y limit
    
    bcy = state(:,2) < 0;
    state(bcy,2) = -state(bcy,2);
    state(bcy,4) =- state(bcy,4);
    % Ensure the electrons bounce back up if going below the low Y limit
    
    temps(i) = (sum(state(:,3).^2) + sum(state(:,4).^2)) * mn/k/2/pop_tot;
    % Track the temperature using the velocity of every electron at every
    % point
    
    for j = 1:pop_vis
        traj(i, (2*j):(2*j+1)) = state(j, 1:2);
    end
    % Track the trajectory of a small number of electrons by storing the x
    % and y position for each iteration
    
    ra = rand(pop_tot, 1) < ps;
    scattersize = size(state(ra,3:4));
    distvels = zeros(scattersize(1),1) + (randn(scattersize(1),1).*vth);
    newvels = ones(scattersize);
    
    for a = 1:scattersize(1)-1
        angle = rand*2*pi;
        vel = randn*vth;
        newvels(a, :) = [vel*cos(angle) vel*sin(angle)];
    end
 
    state(ra,3:4) = newvels;
    % Check to see if the electrons will scatter based around the
    % probability calculated
    % If they do scatter, the velocity and direction will be changed to a
    % random value within the normal distribution around vth
    
    
    J(i,1) = q*density*mean(state(:,3));
    J(i,2) = q*density*mean(state(:,4));
    % Record the current density
    
    % The animation does not need to be updated at every iteration
    % So we update it every 5 iterations to use less memory
    if animate && mod(i,5) == 0
        figure(1)
        pbaspect([2,1,1]);
        subplot(2,1,1);
        hold off;
        plot(state(1:pop_vis,1)./1e-9, state(1:pop_vis,2)./1e-9, 'o');
        title('Electrons with a Fixed Velocity with Force Applied (P.1)');
        xlabel('X (nm)');
        ylabel('Y (nm)');
        hold on;
        for j=1:pop_vis
            plot(traj(:,j*2)./1e-9, traj(:,j*2+1)./1e-9, '.');
        end
        
        if i > 1   
            subplot(2,1,2);
            hold off;
            plot(ts*(0:i-1), J(1:i,1)); 
            title('Current in x-Direction');
            ylabel('Current (A)');
            xlabel('Time (s)');
            hold on;
        end
        pause(0.05)
    end
end

% Now we need to display the trajectory after the iterations have completed
figure(1);
pbaspect([2,1,1]);
subplot(2,1,1);
title('Electrons with a Fixed Velocity with Force Applied (P.1)');
xlabel('X (nm)');
ylabel('Y (nm)');
hold on;
for i=1:pop_vis
    plot(traj(:,i*2)./1e-9, traj(:,i*2+1)./1e-9, '.');
end

% If the animation variable = 0, we need to simply plot the results
if(~animate) 
    subplot(2,1,2);
    hold off;
    plot(ts*(0:pop_vis-1), J(1:pop_vis));
    title('Current in x-Direction');
    ylabel('Current (A)');
    xlabel('Time (s)');
end

% % An electron density map is created for the bottlenecked simulation
dmap = [state(:,1), state(:,2)];
figure(2)
hist3(dmap, [200 100]);
surfHandle = get(gca, 'child');
set(surfHandle, 'FaceColor', 'interp', 'EdgeColor', 'none', 'CdataMode', 'auto');
view(2)
title('Electron Density Map');
xlabel('x (m)');
ylabel('y (m)');


%%
% Similarily to the density map, a temperature map is created
% This map is created by putting the electron velocities into areas much
% smaller than the full space called bins.
x_temp_sum = zeros(ceil(W/1e-9), ceil(H/1e-9));
y_temp_sum = zeros(ceil(W/1e-9), ceil(H/1e-9));
num_temp = zeros(ceil(W/1e-9), ceil(H/1e-9));

% Start by observing the velocities of all the electrons
for i = 1:pop_tot
    % We now need to put the electrons in their bins based on location
    x = floor(state(i,1)/1e-9);
    y = floor(state(i,2)/1e-9);
    
    if(x == 0)
        x = 1;
    end
    if(y == 0)
        y = 1;
    end
    
    % The velocities are then added to the total count for each bin
    x_temp_sum(x,y) = x_temp_sum(x,y) + state(i,3)^2;
    y_temp_sum(x,y) = y_temp_sum(x,y) + state(i,4)^2;
    num_temp(x,y) = num_temp(x,y) + 1;
end

% Using the sum of the velocities, the temperature of each bin can be
% calculated
temp = (x_temp_sum + y_temp_sum).*mn./k./2./num_temp;
temp(isnan(temp)) = 0;
temp = temp';

figure(3)
x = linspace(1, W/1e-9, (W/1e-9)+1); 
y = linspace(1, H/1e-9, (H/1e-9)+1); 
pcolor(x,y,temp)
shading interp
colormap(jet)
title('Temperature Map (Pt. 3)');
xlabel('x (nm)');
ylabel('y (nm)');
