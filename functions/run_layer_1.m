function layer1_output = run_layer_1(neg_events, pos_events, speed_vector, direction_vector, filter_size, spike_threshold, decay_rate, reset_potential, filename)
% function SpikeOut1 = run_layer_1(neg_events, pos_events, speed_vector, direction_vector, filter_size, spike_threshold, decay_rate, reset_potential, filename)
%   neg_events        -> struct containing the negative polarity spikes to process  
%   pos_events        -> struct containing the positive polarity spikes to process
%   speed_vector      -> a vector of the speeds (pix/ms)
%   direction_vector  -> a vector of the directions (degrees)
%   filter_size       -> how large should the filters be?
%   spike_threshold   -> Integrate and Fire Threshold
%   decay_rate        -> Decay Rate in units per us
%   reset_potential   -> Potential immediately after reset
%   filename          -> String prefix for saving files

%% Create the filters
disp('Generating Filters')
filters = generate_filters(speed_vector, direction_vector, filter_size);

%% Run layer 1 for each channel
if ~isempty(pos_events)
    fprintf('Implementing refraction for positive events...')
    pos_events_refracted = implement_refraction_all(pos_events, filters);
    fprintf('done\n')
    clear pos_events
    
    fprintf('Simulating Layer 1: Processing positive polarity stream, %i events\n', length(pos_events_refracted{1}.ts));
    pos_layer1_output = spike_layer_1(pos_events_refracted, filters, spike_threshold, decay_rate, reset_potential);
    clear pos_events_refracted
%     save([filename, '_pos_layer1_output'], 'pos_layer1_output')
else
    pos_layer1_output  =[];
end

if ~isempty(neg_events)
    fprintf('Implementing refraction for negative events...')
    neg_events_refracted = implement_refraction_all(neg_events, filters);
    fprintf('done\n')
    clear neg_events
    
    fprintf('Simulating Layer 1: Processing negative polarity stream, %i events\n', length(neg_events_refracted{1}.ts));
    neg_layer1_output = spike_layer_1(neg_events_refracted, filters, spike_threshold, decay_rate, reset_potential);
    clear neg_events_refracted
%     save([filename, '_neg_layer1_output'], 'neg_layer1_output')
else
    neg_layer1_output  =[];
end

%% combine pos and neg event result streams
disp('Recombining positive and negative streams')
if ~isempty(pos_layer1_output) && ~isempty(neg_layer1_output)
    layer1_output = combine_streams(pos_layer1_output, neg_layer1_output);
    clear pos_layer1_output neg_layer1_output
elseif ~isempty(neg_layer1_output)
    layer1_output = neg_layer1_output;
    clear neg_layer1_output
elseif ~isempty(pos_layer1_output)
    layer1_output = pos_layer1_output;
	clear pos_layer1_output
end

% sort into chronological order
[layer1_output.ts, order_index] = sort(layer1_output.ts);
layer1_output.x = layer1_output.x(order_index);
layer1_output.y = layer1_output.y(order_index);
layer1_output.sp = layer1_output.sp(order_index);
layer1_output.dir = layer1_output.dir(order_index);
toc

%here is where intermediate results can be saved

