%%
% Derek Woodard, 100827083

%% Part 3

clear all
close all
set(0,'DefaultFigureWindowStyle','docked');

m0 = 9.10938356e-31;            % rest mass of an electron
mn = 0.26*m0;                   % effective mass of electrons
T = 300;                        % temperature
k = 1.38064852e-23;             % Boltzmann's Constant
Tmn = 0.2e-12;                  % mean time between collisions
q = 1.60217662e-19;             % The charge of an electron
vth = sqrt((2*k*T)/mn);         % calculate thermal velocity for 1.1 with two degrees of freedom
mfp = vth * Tmn;                % calculate the mean free path for 1.2


%% The nominal size of the region is 200nm x 100nm
H = 100e-9;
W = 200e-9;

pop_tot = 1000;                 % use this variable to determine how many electrons are modeled
pop_vis = 5;                    % use this variable to change how many electrons are plotted
ts = H/vth/100;                 % The time step is set to make each plot point have an electron move 1/100th of the height of the area
iter = 1000;                    % how many plot points total
Vx = 0.8;                       % voltage across the x direction is set to 0.8V
Vy = 0;                         % voltage across the y dimension not provided, so assume 0
m0 = 9.10938215e-31;            % mass of an electron
density = 1e19;                 % density of electrons
ps = 1 - exp(-((ts)/(Tmn)));    % The eletrons have a probability of scattering

animate = 1;                    % set to 1 to see each iteration
% animate = 0;                    % set to 0 to only see final product
    
% Set applied voltage to 0.8
Vo = 0.8;

% For scale simplicity, the nano-factor is removed
% scale = 1e-9;
% H = round(H/scale);
% W = round(W/scale);
% ts = ts/1e-9;

% Set the spacing for the mesh and the number of nodes total based on
% region size
dx = 1;
dy = 1;
nx = W/dx;
ny = H/dy;

% There are two different sigma values based on the bottleneck
% sig1 = 1;
% sig2 = 10^-2;
% sigs = zeros(nx,ny);

% We need to set the 'bottleneck' by creating two boxes in the region
% boxes = [nx*2/5 nx*3/5 ny*2/5 ny*3/5];
% % boxes = 1e-9 .* [80 120 0 40; 80 120 60 100]; from A1
% 
% % We need to set the sigma values for the full region, giving different
% % sigma values to the boxes
% for i = 1:nx
%     for j = 1:ny
%         if i > boxes(1) && i < boxes(2) && (j < boxes(3) || j > boxes(4))
%             sigs(i,j) = sig2;
%         else
%             sigs(i,j) = sig1;
%         end
%     end
% end

%%
% The following is copied from part 2, It is not needed for part 3 and
% takes a lot of memory and time to complete.
% The G and F matrices are formed using the number of points calculated above
% G = zeros(nx*ny);
% F = zeros(1,nx*ny);
% 
% for i=1:nx
%     for j=1:ny
%         % set up the mapping variables
%         n = j+(i-1)*ny;
%         nxp = j+i*ny;
%         nxm = j+(i-2)*ny;
%         nyp = j+1+(i-1)*ny;
%         nym = j-1+(i-1)*ny;
%         
%         if i == 1
%             G(n,:) = 0;
%             G(n,n) = 1;
%             F(n) = Vo;
%         elseif i == nx
% %             G(n,:) = 0;
%             G(n,n) = 1;
% %             F(n) = 0;
%         elseif j == 1
%             G(n,nxp) = (sigs(i+1,j) + sigs(i,j))/2;
%             G(n,nxm) = (sigs(i-1,j) + sigs(i,j))/2;
%             G(n,nyp) = (sigs(i,j+1) + sigs(i,j))/2;
%             G(n,n) = -(G(n,nxp) + G(n,nxm) + G(n,nyp));
%         elseif j == ny
%             G(n,nxp) = (sigs(i+1,j) + sigs(i,j))/2;
%             G(n,nxm) = (sigs(i-1,j) + sigs(i,j))/2;
%             G(n,nym) = (sigs(i,j-1) + sigs(i,j))/2;
%             G(n,n) = -(G(n,nxp) + G(n,nxm) + G(n,nym));
%         else
%             G(n,nxp) = (sigs(i+1,j) + sigs(i,j))/2;
%             G(n,nxm) = (sigs(i-1,j) + sigs(i,j))/2;
%             G(n,nyp) = (sigs(i,j+1) + sigs(i,j))/2;
%             G(n,nym) = (sigs(i,j-1) + sigs(i,j))/2;
%             G(n,n) = -(G(n,nxp) + G(n,nxm) + G(n,nyp) + G(n,nym));
%         end
%     end
% end
% 
% V = G\F';
% S = zeros(ny,nx,1);
% 
% for i = 1:nx
%     for j = 1:ny
%         n = j+(i-1)*ny;
%         S(j,i) = V(n);
%     end
% end
% 
% figure(1)
% mesh(S)
% colormap(jet)
% xlabel('x')
% ylabel('y')
% title('Bottleneck Potential')
% 
% [Exs Eys] = gradient(S);
% Eys = -Eys;
% Exs = -Exs;
% 
% figure(2)
% quiver(Exs, Eys)
% ylabel('y');
% xlabel('x');
% title('Electric Field with Bottleneck')
% 
% Exs = Exs';
% Eys = Eys';
% Fxs = zeros(W,H);
% Fys = zeros(W,H);
% Axs = zeros(W,H);
% Ays = zeros(W,H);

% Since the electric field uses numbered blocks instead of measurements, we
% can return the scale to normal and still get the same results
% H = H*scale;
% W = W*scale;
% ts = ts*scale;
% boxes = scale .* boxes;


% This was repeated for multiple different bottleneck sizes to see how much
% the gap available affects the flow of electrons
% From figure 3 that appears, it is clear that as the bottleneck size
% decreaces, the current decreaces as well.
for bns = 0.1:0.1:0.4
    %%
    % In order to track the electrons, the following arrays will be used
    % The arrays will store the positions, velocities, and temperatures of the
    % elctrons
    state = zeros(pop_tot, 4);
    traj = zeros(iter, pop_vis*2);
    temps = zeros(iter,1);
    J = zeros(iter,2);  
    boxes = [nx*2/5 nx*3/5 ny*bns ny*(1-bns)];

    %% 
    % Start by generating the initial eletrons with a constant speed of vth
    for i = 1:pop_tot
        angle = rand*2*pi;
        state(i,:) = [W*rand H*rand vth*cos(angle) vth*sin(angle)];
            while(state(i,1) > boxes(1) && state(i,1) < boxes(2) && (state(i,2) < boxes(3) || state(i,2) > boxes(4)))
                state(i,1:2) = [W*rand H*rand];
            end
    end

    %%
    % This portion is not required for part 3, only part 2
    % Using the calculated electric field, we determine the force, then
    % acceleration at every 1x1 location
    % Fxs(:,:) = q.* Exs(:,:);
    % Fys(:,:) = q.* Eys(:,:);
    % Axs(:,:) = Fxs(:,:)./m0;
    % Ays(:,:) = Fys(:,:)./m0;
    % Given that these values are derived for incredibly small areas, the
    % resulting accelerations are not noticable when plotting the trajectories.
    % For this reason, these values are not used in the plotting, however, the
    % code that uses them is commented out to show how they can be used.
    Ex = Vx/W;
    Ey = Vy/H;
    Fx = q * Ex;
    Fy = q * Ey;
    Ax = Fx/m0;
    Ay = Fy/m0;

    xloc = zeros(pop_tot,1);
    yloc = zeros(pop_tot,1);

    %%
    % Now we iterate to update the electron positions and plot the trajectories
    for i = 1:iter
        state(:,1:2) = state(:,1:2) + ts.*state(:,3:4);
        % update the state of each electron using the speed

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

        % Check for collisions with the boxes
        for b = 1:pop_tot
            while(state(b,1) > boxes(1) && state(b,1) < boxes(2) && (state(b,2) < boxes(3) || state(b,2) > boxes(4)))
                if state(b,3) > 0
                    xd = state(b,1) - boxes(1);
                    newx = boxes(1);
                else
                    xd = boxes(2) - state(b,1);
                    newx = boxes(2);
                end

                if(state(b,4) > 0)
                    yd = state(b,2) - boxes(4);
                    newy = boxes(4);
                else
                    yd = boxes(3) - state(b,2);
                    newy = boxes(3);
                end

                if(xd < yd)
                    state(b,1) = newx;
                    state(b,3) = -state(b,3);
                else
                    state(b,2) = newy;
                    state(b,4) = -state(b,4);
                end
            end
        end

            % This code uses the acceleration matrices that do not provide
            % enough to visibly affect the plot
    %     % Determine where in the area the electons are
    %     xloc(:) = round(state(:,1));
    %     xloc(xloc==0) = 1;
    %     
    %     yloc(:) = round(state(:,2));
    %     yloc(yloc==0)=1;
    %     
    %     for b = 1:pop_tot:1
    %         state(b,3) = state(b,3) + ts*Axs(xloc(b),yloc(b));
    %         state(b,4) = state(b,4) + ts*Ays(xloc(b),yloc(b));
    %     end
    %     % The speed of each electron is updated using the acceleration from the
    %     % electric field applied
        state(:,3) = state(:,3) + ts*Ax;
        state(:,4) = state(:,4) + ts*Ay;

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
        % no need to animate for part 3
        if bns == 0.4
            if animate && mod(i,5) == 0
                figure(1)
                hold off;
                plot(state(1:pop_vis,1), state(1:pop_vis,2), 'o');
                title('Electrons with a Fixed Velocity with Force Applied (P.1)');
                xlabel('X (m)');
                ylabel('Y (m)');
                hold on;
                for j=1:pop_vis
                    plot(traj(:,j*2), traj(:,j*2+1), '.');
                end

                % Plot the boxes
                rectangle('Position',  [boxes(1) boxes(4) boxes(2)-boxes(1) ny-boxes(4)]);
                rectangle('Position',  [boxes(1) 0 boxes(2)-boxes(1) ny-boxes(4)]);
                pause(0.05)
            end
        end
    end

    % Now we need to display the trajectory after the iterations of completed
    % no need to display trajectory for part 3
    % only display the density map for the first iteration
    if bns == 0.4
        figure(1);
        title('Electrons with a Fixed Velocity with Force Applied (P.1)');
        xlabel('X (m)');
        ylabel('Y (m)');
        hold on;
        for i=1:pop_vis
            plot(traj(:,i*2), traj(:,i*2+1), '.');
        end
        % Plot the boxes
        rectangle('Position',  [boxes(1) boxes(4) boxes(2)-boxes(1) ny-boxes(4)]);
        rectangle('Position',  [boxes(1) 0 boxes(2)-boxes(1) ny-boxes(4)]);
        hold off;

    % An electron density map is created for the bottlenecked simulation
    
        dmap = [state(:,1), state(:,2)];
        figure(2)
        hist3(dmap, [200 100]);
        surfHandle = get(gca, 'child');
        set(surfHandle, 'FaceColor', 'interp', 'EdgeColor', 'none', 'CdataMode', 'auto');
        view(2)
        title('Electron Density Map');
        xlabel('x (m)');
        ylabel('y (m)');
    end
        
    Js = sqrt(J(:,1).^2 + J(:,2).^2);
    
    figure(3)
    hold on
    
    if bns == 0.1
        cur = sum(Js,2);
        tot_cur = sum(cur);
        old_cur = tot_cur;
        plot([bns,bns], [old_cur,tot_cur])
    end
    if bns > 0.1
        old_cur = tot_cur;
        cur = sum(Js,2);
        tot_cur = sum(cur);
        plot([bns-0.01,bns], [old_cur,tot_cur])
        xlabel('Bottleneck Size')
        ylabel('Current Density')
    end
    title('Effect of Bottleneck size on Current Density')
end

