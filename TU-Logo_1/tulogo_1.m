%% Formation control utilizing edge tension energy with a static, undirected
% communication topology
% Paul Glotfelter updated by Sean Wilson
% 07/2019

%% Set up Robotarium object

N = 20;
initial_conditions = generate_initial_conditions(N, 'Width', 2, 'Height', 1, 'Spacing', 0.5);
r = Robotarium('NumberOfRobots', N, 'ShowFigure', true, 'InitialConditions', initial_conditions);


% Import and scale the tu logo appropriately.
tu_img = imread('tulogo.png'); % Original input image file

% Display the image with an associated spatial referencing object.
x_img = linspace(-0.95,  1.35, size(tu_img,2));
y_img = linspace( 0.9, -0.9, size(tu_img,1)); %Note the 1 to -1 here due to the (1,1) pixel being in the top left corner.
tu_img_handle = image(x_img, y_img, tu_img,'CDataMapping','scaled');


%% Set up constants for experiment

%Gains for the transformation from single-integrator to unicycle dynamics
formation_control_gain = 10;

% Select the number of iterations for the experiment.  This value is
% arbitrary
iterations = 2000;

% Communication topology for the desired formation.  We need 2 * N - 3 = 9
% edges to ensure that the formation is rigid.
L = [3 -1 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 -1; ... %1
     -1 3 -1 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0; ... %2
     0 -1 5 -1 0 0 0 -1 -1 0 0 -1 0 0 0 0 0 0 0 0; ... %3
     0 0 -1 4 -1 0 0 0 0 0 -1 0 -1 0 0 0 0 0 0 0; ... %4
     0 0 0 -1 5 -1 0 -1 0 0 0 -1 0 0 0 0 0 0 0 0; ... %5
     0 0 0 0 -1 3 -1 0 0 0 0 0 0 -1 0 0 0 0 0 0; ... %6
     0 0 0 0 0 -1 3 0 -1 0 0 0 0 0 -1 0 0 0 0 0; ... %7
     -1 0 -1 0 -1 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0; ... %8
     0 0 -1 0 0 0 -1 0 2 -1 0 0 0 0 0 0 0 0 0 0; ... %9
     0 0 0 0 0 0 0 0 -1 2 -1 0 0 0 0 0 0 0 0 0; ... %10
     0 -1 0 -1 0 0 0 0 0 -1 5 -1 0 0 0 0 0 0 -1 0; ... %11
     0 0 -1 0 -1 0 0 0 0 0 -1 6 -1 0 0 0 0 -1 0 -1; ... %12
     0 0 0 -1 0 0 0 0 0 0 0 -1 8 -1 0 -1 -1 -1 -1 0; ... %13
     0 0 0 0 0 -1 0 0 0 0 0 0 -1 3 -1 0 0 0 0 0; ... %14
     0 0 0 0 0 0 -1 0 0 0 0 0 0 -1 3 -1 0 0 0 0; ... %15
     0 0 0 0 0 0 0 0 0 0 0 0 -1 0 -1 3 -1 0 0 0; ... %16
     0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 -1 3 -1 0 0; ... %17
     0 0 0 0 0 0 0 0 0 0 0 -1 -1 0 0 0 -1 4 -1 0; ... %18
     0 0 0 0 0 0 0 0 0 0 -1 0 -1 0 0 0 0 -1 4 -1; ... %19
     -1 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 -1 3]; %20
    
% The desired inter-agent distance for the formation
d = 0.15;

D = 3*d;

x1314 = d*sind(8.9); %2,3207
x1415 = d*sind(66.94); %13,8014

c = D - 2*(x1314 + x1415); %12,7558
a = (x1314 + x1415 + c) / sind(54.37); %35,5291
b = (x1314 + 2*x1415 + c) /sind(70.85); %45.18



% Pre-compute diagonal values for the rectangular formation
e = sqrt((2*d)^2 + d^2);
f = sqrt((3*d)^2 + d^2);
g = sqrt((6*d)^2 + (2*d)^2);
o = (x1314 + x1415 - d) / sind(1.27); %[7-9]


% Weight matrix containing the desired inter-agent distances to achieve a
% rectuangular formation
weights = [0 d 0 0 0 0 0 e 0 0 0 0 0 0 0 0 0 0 0 g; ... %1
           d 0 d 0 0 0 0 0 0 0 f 0 0 0 0 0 0 0 0 0; ... %2
           0 d 0 d 0 0 0 d d 0 0 f 0 0 0 0 0 0 0 0; ... %3
           0 0 d 0 d 0 0 0 0 0 f 0 f 0 0 0 0 0 0 0; ... %4
           0 0 0 d 0 d 0 e 0 0 0 f 0 0 0 0 0 0 0 0; ... %5
           0 0 0 0 d 0 d 0 0 0 0 0 0 D 0 0 0 0 0 0; ... %6
           0 0 0 0 0 d 0 0 o 0 0 0 0 0 D 0 0 0 0 0; ... %7
           e 0 d 0 e 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ... %8
           0 0 d 0 0 0 o 0 0 d 0 0 0 0 0 0 0 0 0 0; ... %9
           0 0 0 0 0 0 0 0 d 0 d 0 0 0 0 0 0 0 0 0; ... %10
           0 f 0 f 0 0 0 0 0 d 0 d 0 0 0 0 0 0 f 0; ... %11
           0 0 f 0 f 0 0 0 0 0 d 0 d 0 0 0 0 f 0 f; ... %12
           0 0 0 f 0 0 0 0 0 0 0 d 0 d 0 a b D f 0; ... %13
           0 0 0 0 0 D 0 0 0 0 0 0 d 0 d 0 0 0 0 0; ... %14
           0 0 0 0 0 0 D 0 0 0 0 0 0 d 0 c 0 0 0 0; ... %15
           0 0 0 0 0 0 0 0 0 0 0 0 a 0 c 0 d 0 0 0; ... %16
           0 0 0 0 0 0 0 0 0 0 0 0 b 0 0 d 0 d 0 0; ... %17
           0 0 0 0 0 0 0 0 0 0 0 f D 0 0 0 d 0 d 0; ... %18
           0 0 0 0 0 0 0 0 0 0 f 0 f 0 0 0 0 d 0 d; ... %19
           g 0 0 0 0 0 0 0 0 0 0 f 0 0 0 0 0 0 d 0]; %20
    
% Initialize velocity vector for agents.  Each agent expects a 2 x 1
% velocity vector containing the linear and angular velocity, respectively.
dx = zeros(2, N);

%% Grab tools for converting to single-integrator dynamics and ensuring safety 

uni_barrier_cert = create_uni_barrier_certificate_with_boundary('SafetyRadius', 0.01);
si_to_uni_dyn = create_si_to_uni_dynamics('LinearVelocityGain', 0.5, 'AngularVelocityLimit', pi/2);

% Iterate for the previously specified number of iterations
for t = 0:iterations
    
    
    % Retrieve the most recent poses from the Robotarium.  The time delay is
    % approximately 0.033 seconds
    x = r.get_poses();
    
    %% Algorithm
    
    %This section contains the actual algorithm for formation control!
    
    %Calculate single integrator control inputs using edge-energy consensus
    for i = 1:N
        
        % Initialize velocity to zero for each agent.  This allows us to sum
        % over agent i's neighbors
        dx(:, i) = [0 ; 0];
        
        % Get the topological neighbors of agent i from the communication
        % topology
        for j = topological_neighbors(L, i)
                
            % For each neighbor, calculate appropriate formation control term and
            % add it to the total velocity

            dx(:, i) = dx(:, i) + ...
            formation_control_gain*(norm(x(1:2, i) - x(1:2, j))^2 - weights(i, j)^2) ... 
            *(x(1:2, j) - x(1:2, i));
        end 
    end
    
    %% Avoid actuator errors
    
    % To avoid errors, we need to threshold dx
    norms = arrayfun(@(x) norm(dx(:, x)), 1:N);
    threshold = 3/4*r.max_linear_velocity;
    to_thresh = norms > threshold;
    dx(:, to_thresh) = threshold*dx(:, to_thresh)./norms(to_thresh);
    
    % Transform the single-integrator dynamics to unicycle dynamics using a provided utility function
    dx = si_to_uni_dyn(dx, x);  
    dx = uni_barrier_cert(dx, x);
    
    % Set velocities of agents 1:N
    r.set_velocities(1:N, dx);
    
    % Send the previously set velocities to the agents.  This function must be called!
    r.step();   
end

% We can call this function to debug our experiment!  Fix all the errors
% before submitting to maximize the chance that your experiment runs
% successfully.
r.debug();
