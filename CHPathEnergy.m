clc; 
clear all;
close all;
% Step 1: Start
% Step 2: Setting up WSN simulation parameters
Bits = 1000; % Replace with your specific value or calculation
networkSize = 100;  % Size of the WSN field
BS_location = [0, 0];  % Base station location
n = 50;  % Number of nodes
alivenodes = n;
Eo = 0.5;  % Initial energy of each node
ETX = 50e-9;  % Energy used in transmitting
ERX = 50e-9;  % Energy used in receiving
Eamp = 0.0013e-12;  % Transmission energy for amplifier
Efs = 10e-12;  % Reception energy for amplifier
EDA = 5e-9;  % Data aggregation energy

% Initialize a variable to store plotted paths
plottedPaths = {};

% Visualize the CH selection paths

minEnergyPath = [];  % Initialize the path with minimum energy consumption
minEnergy = inf;     % Initialize the minimum energy value to infinity


% Step 3: For every node
nodes = struct('ID', 0, 'Location', [0, 0], 'Energy', Eo, ...
    'Dtch', 0, 'Dts', 0, 'Chid', 0, 'Dkmean', 0, 'ClusterID', 0, ...
    'Cond', 1, 'Role', 0);

for i = 1:n
    % a. Set Node ID.
    nodes(i).ID = i;
    
    % b. Randomly assign coordinates as the node’s location.
    nodes(i).Location = rand(1, 2) * networkSize;
    
    % c. Set node’s initial energy = Eo.
    nodes(i).Energy = Eo;
    
    % d-h. Set other parameters as specified.
    nodes(i).Dtch = 0;
    nodes(i).Dts = 0;
    nodes(i).Chid = 0;
    nodes(i).Dkmean = 0;
    nodes(i).ClusterID = 0;
    nodes(i).Cond = 1;
    nodes(i).Role = 0;
end

% Visualize the initial network field
figure; 
locations = reshape([nodes.Location], 2, []).';  % Extract locations
scatter(locations(:, 1), locations(:, 2), 'filled');

hold on;
scatter(BS_location(1), BS_location(2), 100, 'r', 'filled', 'Marker', '*');
title('WSN field - Initial Configuration');
xlabel('X-axis');
ylabel('Y-axis');
legend('Nodes', 'Base Station');
hold off;

% Step 4: Initialize report variables
total_transmitted_data = 0;
total_received_data = 0;
total_aggregated_data = 0;
total_energy_consumed = 0;

% Step 5: Simulation
max_iterations = 15;  % Set a maximum number of iterations to avoid infinite loops
iteration_counter = 0;  % Initialize iteration counter

% Initialize variables for PDR calculation
total_transmitted_packets = 0;
total_received_packets = 0;

% Initialize a variable to store CH selection paths
chSelectionPaths = cell(n, 1);

while any([nodes.Cond]) && iteration_counter < max_iterations
    iteration_counter = iteration_counter + 1;
    
    fprintf('Iteration: %d\n', iteration_counter);

    
    % a. Setup phase
    % i-iv. Implement the setup phase as described
    
    % b. K-means Clustering
    % i. Calculate the density and the density parameter (ρ)
    alive_nodes = find([nodes.Cond]);
    Density = length(alive_nodes) / networkSize;
    Density_Parameter = 4 / Density;
    
    % ii. Calculate the number of clusters (No) to be formed
    No = round((length(alive_nodes) / 100) * Density_Parameter);
    
    % iii. Run k-means with better initialization
    data = zeros(length(alive_nodes), 2);
    for j = 1:length(alive_nodes)
        data(j, :) = nodes(alive_nodes(j)).Location;
    end
    
    % Better initialization using 'kmeans++'
    try
        [idx, centroids] = kmeans(data, No, 'MaxIter', 15, 'Replicates', 5, 'Start', 'plus');
    catch
        warning('K-means did not converge. Skipping this iteration.');
        continue;  % Skip to the next iteration if k-means does not converge
    end
    
    % iv. For each cluster, find the node closest to the centroid and elect it as CH
    for k = 1:No
        cluster_nodes = alive_nodes(idx == k);
        distances_to_centroid = zeros(length(cluster_nodes), 1);
        for l = 1:length(cluster_nodes)
            distances_to_centroid(l) = norm(nodes(cluster_nodes(l)).Location - centroids(k, :));
            % Record the CH selection path for each node in the cluster
            chSelectionPaths{cluster_nodes(l)} = [chSelectionPaths{cluster_nodes(l)}, k];
        end
        [~, idx_closest] = min(distances_to_centroid);
        elected_CH = cluster_nodes(idx_closest);
        
        % Update CH information for the elected node
        nodes(elected_CH).Role = 1;
        nodes(elected_CH).ClusterID = k;
        nodes(elected_CH).Chid = elected_CH;
    end
    
    % v. Update other node information related to clustering
    
    % b. Steady Phase
    % i-iv. Implement the steady phase as described
    for i = 1:n
        % Process each node
        if nodes(i).Cond == 1
            % Node is alive
            % Transmit data to CH
            if nodes(i).Dtch <= (Efs/Eamp)
                Etransmit = Bits * (ETX + (Efs * nodes(i).Dtch^2));
            else
                Etransmit = Bits * (ETX + (Eamp * nodes(i).Dtch));
            end
            
            % Receive data at CH
            ERecept = Bits * (ERX + EDA);
            
            % Transmit aggregated data to sink node
            if nodes(i).Dts <= (Efs/Eamp)
                EtransmitSink = Bits * ((ETX + EDA) + (Efs * nodes(i).Dts^2));
            else
                EtransmitSink = Bits * ((ETX + EDA) + (Eamp * nodes(i).Dts));
            end
            
            % Update energy levels
            nodes(i).Energy = nodes(i).Energy - (Etransmit + ERecept + EtransmitSink);
            
            % Check if the node died during transmission or reception
            if nodes(i).Energy <= 0
                nodes(i).Cond = 0;
                nodes(i).Chid = 0;
                nodes(i).Role = 0;
                % Update other parameters as needed
            else
                % Update reports
                total_transmitted_data = total_transmitted_data + Etransmit;
                total_received_data = total_received_data + ERecept;
                total_aggregated_data = total_aggregated_data + EtransmitSink;
                total_energy_consumed = total_energy_consumed + (Etransmit + ERecept + EtransmitSink);
                % Simulate packet transmission and reception
                % Assuming a probability of successful transmission and reception
                if rand < 0.8  % Example: 80% chance of successful transmission
                    total_transmitted_packets = total_transmitted_packets + 1;
                    if rand < 0.9  % Example: 90% chance of successful reception
                        total_received_packets = total_received_packets + 1;
                    end
                end
            end
        end
    end
    
    % Visualize clusters and centroids after the first round
    if iteration_counter == 1
        figure;
        gscatter(data(:, 1), data(:, 2), idx, 'rgbm', 'o', 8, 'off');
        hold on;
        scatter(centroids(:, 1), centroids(:, 2), 100, 'k', 'x');
        title('Clusters and Centroids - Round 1');
        xlabel('X-axis');
        ylabel('Y-axis');
        legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Centroids');
        hold off;
    end
end
 
% Visualize network survival over rounds
figure;
plot(1:iteration_counter, sum([nodes.Cond]), 'LineWidth', 2);

title('Network Survival Over Rounds');
xlabel('Round');
ylabel('Number of Operational Nodes');
ylim([0, n]);
grid on;



% Calculate and visualize minimum and maximum energy paths
% Initialize variables
minEnergy = inf;            % Initialize minimum energy to infinity
maxEnergy = -inf;           % Initialize maximum energy to negative infinity
minEnergyPath = [];         % Initialize minimum energy path
maxEnergyPath = [];         % Initialize maximum energy path

% Iterate through each path to calculate energy consumption
for i = 1:n
    if ~isempty(chSelectionPaths{i})
        x = nodes(i).Location(1);
        y = nodes(i).Location(2);
        path = chSelectionPaths{i};
        energyConsumption = 0;  % Initialize energy consumption for this path
        
        % Calculate energy consumption for this path
        for j = 1:length(path) - 1
            x_path = [x, centroids(path(j), 1)];
            y_path = [y, centroids(path(j), 2)];
            distance = norm([x, y] - centroids(path(j), :));
            energyConsumption = energyConsumption + distance * Eamp;
            x = centroids(path(j), 1);
            y = centroids(path(j), 2);
        end
        distanceToBS = norm([x, y] - BS_location);
        energyConsumption = energyConsumption + distanceToBS * Eamp;
        
        % Update the minimum energy path and energy if needed
        if energyConsumption < minEnergy
            minEnergy = energyConsumption;
            minEnergyPath = path;
        end
        
        % Update the maximum energy path and energy if needed
        if energyConsumption > maxEnergy
            maxEnergy = energyConsumption;
            maxEnergyPath = path;
        end
    end
end


% Plot the minimum energy path
if ~isempty(minEnergyPath)
    x = nodes(minEnergyPath(1)).Location(1);
    y = nodes(minEnergyPath(1)).Location(2);
    for j = 1:length(minEnergyPath) - 1
        x_path = [x, centroids(minEnergyPath(j), 1)];
        y_path = [y, centroids(minEnergyPath(j), 2)];
        plot(x_path, y_path, 'g-', 'LineWidth', 2);
        x = centroids(minEnergyPath(j), 1);
        y = centroids(minEnergyPath(j), 2);
    end
    x_path = [x, 0]; % Base station at (0,0)
    y_path = [y, 0]; % Base station at (0,0)
    plot(x_path, y_path, 'g-', 'LineWidth', 2); % Green color for minimum energy path
    hold on;
end


% Plot the maximum energy path
if ~isempty(maxEnergyPath)
    x = nodes(maxEnergyPath(1)).Location(1);
    y = nodes(maxEnergyPath(1)).Location(2);
    for j = 1:length(maxEnergyPath) - 1
        x_path = [x, centroids(maxEnergyPath(j), 1)];
        y_path = [y, centroids(maxEnergyPath(j), 2)];
        plot(x_path, y_path, 'b-', 'LineWidth', 2);
        x = centroids(maxEnergyPath(j), 1);
        y = centroids(maxEnergyPath(j), 2);
    end
    x_path = [x, 0]; % Base station at (0,0)
    y_path = [y, 0]; % Base station at (0,0)
    plot(x_path, y_path, 'b-', 'LineWidth', 2); % Magenta color for maximum energy path
    hold on;
end


% Plot Base Station

title('CH Selection Paths (Minimum and Maximum Energy)');
xlabel('X-axis');
ylabel('Y-axis');
legend('Minimum Energy Path', 'Maximum Energy Path');
grid on;


% Visualize all possible paths
figure;
for i = 1:n
    if ~isempty(chSelectionPaths{i})
        x = nodes(i).Location(1);
        y = nodes(i).Location(2);
        path = chSelectionPaths{i};
        
        % Check if the path has been plotted already
        pathString = num2str(path);
        if any(strcmp(plottedPaths, pathString))
            continue; % Skip if the path has been plotted before
        end
        
        % Add the path to the list of plotted paths
        plottedPaths{end+1} = pathString;
        
        for j = 1:length(path) - 1
            x_path = [x, centroids(path(j), 1)];
            y_path = [y, centroids(path(j), 2)];
            plot(x_path, y_path, 'b-', 'LineWidth', 1);
            x = centroids(path(j), 1);
            y = centroids(path(j), 2);
        end
        x_path = [x, BS_location(1)];
        y_path = [y, BS_location(2)];
        plot(x_path, y_path, 'r-', 'LineWidth', 1);
        hold on;
    end
end
scatter(BS_location(1), BS_location(2), 100, 'r', 'filled', 'Marker', '*');
title('All CH Selection Paths');
xlabel('X-axis');
ylabel('Y-axis');
legend('Paths');
grid on;


% Visualize survival nodes at the end of the simulation
figure;
scatter([nodes.Location], [nodes.Location], 'filled', 'MarkerFaceColor', 'b');
hold on;
scatter([nodes(~[nodes.Cond]).Location], [nodes(~[nodes.Cond]).Location], 'filled', 'MarkerFaceColor', 'r');
title('Survival Nodes at the End of the Simulation');
xlabel('X-axis');
ylabel('Y-axis');
legend('Survival Nodes', 'Failed Nodes');
hold off;

% Step 5: Generate final reports
fprintf('--- Simulation Reports ---\n');
fprintf('Total Transmitted Data: %f J\n', total_transmitted_data);
fprintf('Total Received Data: %f J\n', total_received_data);
fprintf('Total Aggregated Data: %f J\n', total_aggregated_data);
fprintf('Total Energy Consumed: %f J\n', total_energy_consumed);

% Print minimum and maximum energy
fprintf('Minimum Energy Required for Data Transfer: %.4e J\n', minEnergy);
fprintf('Maximum Energy Required for Data Transfer: %.4e J\n', maxEnergy);
% Step 6: Transfer Data from Base Station to Target Node

% Specify the ID of the target node
targetID = 20; % Replace with the ID of the target node

% Location of the base station
sourceLocation = BS_location;

% Initialize variables for minimum energy path
minEnergyConsumption = Inf;
bestPath = [];

% Iterate through all possible paths
for i = 1:n
    % Check if there is a path from the base station to node i
    if ~isempty(chSelectionPaths{i})
        path = chSelectionPaths{i};
        
        % Calculate energy consumption for this path
        totalEnergyConsumption = 0;
        pathCoordinates = sourceLocation;  % Start from the base station location
        
        for j = 1:length(path) - 1
            clusterID = path(j);
            centroid = centroids(clusterID, :);
            distance = norm(pathCoordinates(end, :) - centroid);
            energyConsumption = distance * Eamp;
            totalEnergyConsumption = totalEnergyConsumption + energyConsumption;
            pathCoordinates = [pathCoordinates; centroid];
        end
        
        % Calculate energy consumption from the last centroid to the target node
        lastCentroid = centroids(path(end), :);
        distanceToTarget = norm(lastCentroid - nodes(targetID).Location);
        energyConsumptionToTarget = distanceToTarget * Eamp;
        totalEnergyConsumption = totalEnergyConsumption + energyConsumptionToTarget;
        
        % Check if this path has minimum energy consumption so far
        if totalEnergyConsumption < minEnergyConsumption
            minEnergyConsumption = totalEnergyConsumption;
            bestPath = path;
        end
    end
end

% Visualize the best path
if ~isempty(bestPath)
    % Plot the path
    pathCoordinates = sourceLocation;  % Start from the base station location
    for i = 1:length(bestPath) - 1
        clusterID = bestPath(i);
        centroid = centroids(clusterID, :);
        pathCoordinates = [pathCoordinates; centroid];
    end
    lastCentroid = centroids(bestPath(end), :);
    pathCoordinates = [pathCoordinates; lastCentroid; nodes(targetID).Location];
    
    figure;
    plot(pathCoordinates(:, 1), pathCoordinates(:, 2), 'b-', 'LineWidth', 2);
    hold on;
    scatter(sourceLocation(1), sourceLocation(2), 100, 'g', 'filled', 'Marker', 'o');
    scatter(nodes(targetID).Location(1), nodes(targetID).Location(2), 100, 'r', 'filled', 'Marker', 'o');
    scatter(centroids(bestPath, 1), centroids(bestPath, 2), 100, 'k', 'filled', 'Marker', 'x');
    scatter(lastCentroid(1), lastCentroid(2), 100, 'k', 'filled', 'Marker', 'x');
    title('Best Path from Base Station to Target Node (Minimum Energy)');
    xlabel('X-axis');
    ylabel('Y-axis');
    legend('Path', 'Base Station', 'Target Node', 'Cluster Centroids');
    grid on;

    % Print energy consumption for the best path
    fprintf('Minimum energy consumed for data transfer from base station to node %d: %.4e J\n', targetID, minEnergyConsumption);
else
    fprintf('No path found from base station to node %d\n', targetID);
end

% Calculate the estimated number of transmitted packets
transmitted_packets = total_transmitted_data / Etransmit;

% Calculate the estimated number of received packets
received_packets = total_received_data / ERecept;

% Calculate Packet Delivery Ratio (PDR)
pdr = total_received_packets / total_transmitted_packets;

fprintf('Packet Delivery Ratio (PDR): %.2f%%\n', pdr * 100);

