clear all
addpath(genpath('functions'))
filename = 'pipe'; %which matlab file contains the input events?
%format is:
%       TD.x =  pixel X locations, strictly positive integers only (TD.x>0)
%       TD.y =  pixel Y locations, strictly positive integers only (TD.y>0)
%   	TD.p =  event polarity. TD.p = 0 for OFF events, TD.p = 1 for ON events
%       TD.ts = event timestamps, typically in units of microseconds

load(['Data', filesep, filename])

%make sure each struct element is a row vector
[num_rows, num_cols] = size(TD.ts);
if (num_rows<num_cols)
    TD.p = TD.p';
    TD.x = TD.x';
    TD.y = TD.y';
    TD.ts = TD.ts';
end

%generate the filter speeds and directions
filter_speeds = sqrt(2)*0.02*sqrt(2).^(0:7); %pixels per millisecond
filter_directions = -135:45:180;

%% Setup the layer parameters
% FilterSize        = 5;
filter_size         = [3,3];
spike_threshold     = 50;
reset_potential     = -128;
decay_rate          = filter_speeds*50*1e-3;


%split the events into positive and negative polarity
polarities = unique(TD.p);
pos_indices = (TD.p ~= polarities(2)); %assume the first index is the positive index
TDpos.x     = TD.x(pos_indices);
TDpos.y     = TD.y(pos_indices);
TDpos.ts    = TD.ts(pos_indices);

neg_indices = (TD.p ~= polarities(1)); %assume the first index is the positive index
TDneg.x     = TD.x(neg_indices);
TDneg.y     = TD.y(neg_indices);
TDneg.ts    = TD.ts(neg_indices);
tic
layer1_output = run_layer_1(TDneg, TDpos, filter_speeds, filter_directions, filter_size, spike_threshold, decay_rate, reset_potential, filename);

if ~exist('Results', 'dir')
    mkdir('Results')
end

save(['Results', filesep, filename, '_layer1'], 'layer1_output')

%optionally display the result
show_result(layer1_output)

%%
filter_size = [5,5];
use_delays = 0; %should second layer neurons have delays?

filter_speeds = sqrt(2)*0.02*sqrt(2).^(0:7); %pixels per millisecond
filter_directions = -135:45:180;
filters = generate_filters(filter_speeds, filter_directions, filter_size);

% set all the delays to zero
[num_speeds, num_directions] = size(filters);
if use_delays == 0
    for speed_index = 1:num_speeds
        for direction_index = 1:num_directions
            filters{speed_index,direction_index}(:) = 0;
        end
    end
end

spike_threshold      = 60;
decay_rate           = ((filter_speeds*spike_threshold*1e-3)/100); %use a much slower decay in layer 2. pretty much zero
reset_potential      = -50;

fprintf('Simulating Layer 2: Processing merged stream of %i events\n', length(layer1_output.ts));
layer2_output = spike_layer_2(layer1_output, filters, spike_threshold, decay_rate, reset_potential);

save(['Results', filesep, filename, '_layer2'], 'layer2_output');

%% display the results:

