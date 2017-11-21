function output_spikes = spike_layer_1 (input_spikes, filters, spike_threshold, decay_rate, reset_potential)
%function output_spikes = spike_layer_1(input_spikes, filters, spike_threshold, decay_rate, reset_potential)
%   input_spikes      -> The spikes to be filtered
%   filters           -> A 2D struct of size [num_speeds, num_directions] of 2D filters [x,y] to be used
%   spike_threshold   -> Integrate and Fire Threshold
%   decay_rate        -> Decay Rate in units per us
%   reset_potential   -> Potential immediately after reset

%% Constants
[num_speeds, num_directions] = size(filters); %how many speeds and directions are we using

filter_size = size(filters{1,1}); %how large is each filter (assume they're all the same size)
assert(rem(filter_size(1),2) == 1 & rem(filter_size(2),2) == 1); % check that the filter size is odd
filter_half_size = floor(filter_size./2);


%find how many spikes are in the input, and how big the spatial dimension is
dim_x = 0;
dim_y = 0;
num_spikes_total = 0;
for speed_index = 1:num_speeds
    num_spikes_total = max(num_spikes_total, length(input_spikes{speed_index}.ts));
    dim_x = max(dim_x, max(input_spikes{speed_index}.x(~isinf(input_spikes{speed_index}.x))));
    dim_y = max(dim_y, max(input_spikes{speed_index}.y(~isinf(input_spikes{speed_index}.y))));
end

num_spikes_per_batch = 10000; %how many spikes to process per batch?
num_batches = ceil(num_spikes_total/num_spikes_per_batch); %how many batches must we process?

synaptic_weight = 100/(filter_size(1)*filter_size(2)); %what is the weight of each synapse? (make the sum of the weights 100)

%% Organize incoming spikes
%this struct will hold the chronologically (after synapse delays) ordered spikes
ordered_spikes.ts = [];
ordered_spikes.x = [];
ordered_spikes.y = [];
ordered_spikes.sp = [];
ordered_spikes.dir = [];

%this struct will hold the output spikes
output_spikes.x     = zeros(10000, 1);
output_spikes.y     = output_spikes.x;
output_spikes.ts    = output_spikes.x;
output_spikes.sp    = output_spikes.x;
output_spikes.dir   = output_spikes.x;
%and the index where the next output spike should be placed
output_spike_number = 1;

%arrays for remembering the last time a neuron was updated, and the neuron potential at that time
neuron_time = zeros(num_speeds, num_directions, dim_x, dim_y);
neuron_potential = zeros(num_speeds, num_directions, dim_x, dim_y);

%% Generate a template which we'll add to large batches of incoming spikes to calculate their location and time
fprintf('Generating templates...')
template_index = 1;
for direction_index = 1:num_directions
    for x_index = 1:filter_size(1)
        for y_index = 1:filter_size(2)
            for speed_index = 1:num_speeds
                template{speed_index}.x(template_index:template_index+num_spikes_per_batch-1, 1)     = - x_index + (filter_half_size(1) + 1);
                template{speed_index}.y(template_index:template_index+num_spikes_per_batch-1, 1)     = - y_index + (filter_half_size(2) + 1);
                template{speed_index}.ts(template_index:template_index+num_spikes_per_batch-1, 1)    = filters{speed_index, direction_index}(x_index, y_index);
                template{speed_index}.sp(template_index:template_index+num_spikes_per_batch-1, 1)    = speed_index;
                template{speed_index}.dir(template_index:template_index+num_spikes_per_batch-1, 1)   = direction_index;
            end
            template_index =  template_index+num_spikes_per_batch;
        end
    end
end
fprintf('done!\n')

%sort the template by the length of delay to ensure chronological order of spikes
for speed_index = 1:num_speeds
    [template{speed_index}.ts, order_index]     = sort(template{speed_index}.ts);
    template{speed_index}.x                     = template{speed_index}.x(order_index);
    template{speed_index}.y                     = template{speed_index}.y(order_index);
    template{speed_index}.sp                    = template{speed_index}.sp(order_index);
    template{speed_index}.dir                   = template{speed_index}.dir(order_index);
end

num_template_repeats = (num_directions*filter_size(1)*filter_size(2));
input_spike_index = 1;
fprintf('Will run for %i loops\n', num_batches);
for batch_num = 1:num_batches
    %if we're on the last batch, the batch size will be different, so we must regenerate a smaller template
    if num_spikes_total-input_spike_index < num_spikes_per_batch
        num_spikes_per_batch = num_spikes_total-input_spike_index;
        clear template
        template_index = 1;
        for direction_index = 1:num_directions
            for x_index = 1:filter_size
                for y_index = 1:filter_size
                    for speed_index = 1:num_speeds
                        template{speed_index}.x(template_index:template_index+num_spikes_per_batch-1, 1)      = - x_index + (filter_half_size(1) + 1);
                        template{speed_index}.y(template_index:template_index+num_spikes_per_batch-1, 1)      = - y_index + (filter_half_size(2) + 1);
                        template{speed_index}.ts(template_index:template_index+num_spikes_per_batch-1, 1)     = filters{speed_index,direction_index}(x_index, y_index);
                        template{speed_index}.sp(template_index:template_index+num_spikes_per_batch-1, 1)     = speed_index;
                        template{speed_index}.dir(template_index:template_index+num_spikes_per_batch-1, 1)    = direction_index;
                    end
                    template_index =  template_index + num_spikes_per_batch;
                end
            end
        end
        
        for speed_index = 1:num_speeds
            [template{speed_index}.ts, order_index] = sort(template{speed_index}.ts);
            template{speed_index}.x                 = template{speed_index}.x(order_index);
            template{speed_index}.y                 = template{speed_index}.y(order_index);
            template{speed_index}.sp                = template{speed_index}.sp(order_index);
            template{speed_index}.dir               = template{speed_index}.dir(order_index);
        end
    end
    
    %add the template to the input spikes, and append to the ordered spikes
    for speed_index = 1:num_speeds
        ordered_spikes.ts        = [ordered_spikes.ts; repmat(input_spikes{speed_index}.ts(input_spike_index:input_spike_index+num_spikes_per_batch-1), num_template_repeats, 1) + template{speed_index}.ts];
        ordered_spikes.x         = [ordered_spikes.x; repmat(input_spikes{speed_index}.x(input_spike_index:input_spike_index+num_spikes_per_batch-1), num_template_repeats, 1) + template{speed_index}.x];
        ordered_spikes.y         = [ordered_spikes.y; repmat(input_spikes{speed_index}.y(input_spike_index:input_spike_index+num_spikes_per_batch-1), num_template_repeats, 1) + template{speed_index}.y];
        ordered_spikes.sp        = [ordered_spikes.sp; template{speed_index}.sp];
        ordered_spikes.dir       = [ordered_spikes.dir; template{speed_index}.dir];
    end
    
    % Sort the spikes in chronological order
    [ordered_spikes.ts, order_index]    = sort(ordered_spikes.ts);
    ordered_spikes.x                    = ordered_spikes.x(order_index);
    ordered_spikes.y                    = ordered_spikes.y(order_index);
    ordered_spikes.sp                   = ordered_spikes.sp(order_index);
    ordered_spikes.dir                  = ordered_spikes.dir(order_index);
    
    % remove spikes which are out of range
    order_index = (ordered_spikes.x > 0) & (ordered_spikes.y >0) & (ordered_spikes.x <= dim_x) & (ordered_spikes.y <= dim_y)& (~isinf(ordered_spikes.ts));
    ordered_spikes.ts                   = ordered_spikes.ts(order_index);
    ordered_spikes.x                    = ordered_spikes.x(order_index);
    ordered_spikes.y                    = ordered_spikes.y(order_index);
    ordered_spikes.sp                   = ordered_spikes.sp(order_index);
    ordered_spikes.dir                  = ordered_spikes.dir(order_index);
    
    input_spike_index = input_spike_index + num_spikes_per_batch;
    
    % check the smallest spike time of input spikes which have yet to be ordered
    min_event_time = inf;
    for speed_index = 1:num_speeds
        min_event_time = min(input_spikes{speed_index}.ts(input_spike_index), min_event_time);
    end
    if isempty(min_event_time)
        min_event_time = inf;
    end
    % loop through the spikes to be processed
    ordered_spikes_index = 1;
    
    while (ordered_spikes.ts(ordered_spikes_index) < min_event_time) && (ordered_spikes_index <= length(ordered_spikes.ts))
        
        % implement membrane potential decay
        if neuron_potential(ordered_spikes.sp(ordered_spikes_index), ordered_spikes.dir(ordered_spikes_index), ordered_spikes.x(ordered_spikes_index), ordered_spikes.y(ordered_spikes_index)) >= 0
            neuron_potential(ordered_spikes.sp(ordered_spikes_index), ordered_spikes.dir(ordered_spikes_index), ordered_spikes.x(ordered_spikes_index), ordered_spikes.y(ordered_spikes_index)) = max((neuron_potential(ordered_spikes.sp(ordered_spikes_index), ordered_spikes.dir(ordered_spikes_index), ordered_spikes.x(ordered_spikes_index), ordered_spikes.y(ordered_spikes_index)) - (ordered_spikes.ts(ordered_spikes_index) - neuron_time(ordered_spikes.sp(ordered_spikes_index), ordered_spikes.dir(ordered_spikes_index), ordered_spikes.x(ordered_spikes_index), ordered_spikes.y(ordered_spikes_index)))*decay_rate(ordered_spikes.sp(ordered_spikes_index))), 0);
        else
            neuron_potential(ordered_spikes.sp(ordered_spikes_index), ordered_spikes.dir(ordered_spikes_index), ordered_spikes.x(ordered_spikes_index), ordered_spikes.y(ordered_spikes_index)) = min((neuron_potential(ordered_spikes.sp(ordered_spikes_index), ordered_spikes.dir(ordered_spikes_index), ordered_spikes.x(ordered_spikes_index), ordered_spikes.y(ordered_spikes_index)) + (ordered_spikes.ts(ordered_spikes_index) - neuron_time(ordered_spikes.sp(ordered_spikes_index), ordered_spikes.dir(ordered_spikes_index), ordered_spikes.x(ordered_spikes_index), ordered_spikes.y(ordered_spikes_index)))*decay_rate(ordered_spikes.sp(ordered_spikes_index))), 0);
        end
        
        % add the synaptic weight
        neuron_potential(ordered_spikes.sp(ordered_spikes_index), ordered_spikes.dir(ordered_spikes_index), ordered_spikes.x(ordered_spikes_index), ordered_spikes.y(ordered_spikes_index)) = synaptic_weight + neuron_potential(ordered_spikes.sp(ordered_spikes_index), ordered_spikes.dir(ordered_spikes_index), ordered_spikes.x(ordered_spikes_index), ordered_spikes.y(ordered_spikes_index));
        
        % check if the threshold has been crossed
        if neuron_potential(ordered_spikes.sp(ordered_spikes_index), ordered_spikes.dir(ordered_spikes_index), ordered_spikes.x(ordered_spikes_index), ordered_spikes.y(ordered_spikes_index)) > spike_threshold
            %if so, reset all neurons at the same location (mutual inhibition)
            neuron_potential(:, :, ordered_spikes.x(ordered_spikes_index), ordered_spikes.y(ordered_spikes_index)) = reset_potential; % Inhibit all other neurons at the same location
            neuron_time(:, :, ordered_spikes.x(ordered_spikes_index), ordered_spikes.y(ordered_spikes_index)) = ordered_spikes.ts(ordered_spikes_index);
            
            %generate an output event
            output_spikes.x(output_spike_number)        = ordered_spikes.x(ordered_spikes_index);
            output_spikes.y(output_spike_number)        = ordered_spikes.y(ordered_spikes_index);
            output_spikes.ts(output_spike_number)       = ordered_spikes.ts(ordered_spikes_index);
            output_spikes.sp(output_spike_number)       = ordered_spikes.sp(ordered_spikes_index);
            output_spikes.dir(output_spike_number)      = ordered_spikes.dir(ordered_spikes_index);
            
            %increment the output spike index
            output_spike_number                         = output_spike_number + 1;
            
            %allocate another block of memory if necessary
            if(length(output_spikes.x)< output_spike_number)
                output_spikes.x        = [ordered_spikes.x; zeros(10000,1)];
                output_spikes.y        = [ordered_spikes.y; zeros(10000,1)];
                output_spikes.ts       = [ordered_spikes.ts; zeros(10000,1)];
                output_spikes.sp       = [ordered_spikes.sp; zeros(10000,1)];
                output_spikes.dir      = [ordered_spikes.dir; zeros(10000,1)];
            end
        end
        
        %remember the time at which this neuron was updated
        neuron_time(ordered_spikes.sp(ordered_spikes_index), ordered_spikes.dir(ordered_spikes_index), ordered_spikes.x(ordered_spikes_index), ordered_spikes.y(ordered_spikes_index)) = ordered_spikes.ts(ordered_spikes_index);
        
        %move on to process the next spike in the list
        ordered_spikes_index = ordered_spikes_index + 1;
    end
    
    %remove all the ordered spikes that have been processed thus far
    ordered_spikes_index = ordered_spikes_index-1;
    ordered_spikes.ts(1:ordered_spikes_index)   = [];
    ordered_spikes.x(1:ordered_spikes_index)    = [];
    ordered_spikes.y(1:ordered_spikes_index)    = [];
    ordered_spikes.dir(1:ordered_spikes_index)  = [];
    ordered_spikes.sp(1:ordered_spikes_index)   = [];
    
    disp(['Loop ', num2str(batch_num), ': ', num2str(toc), ' seconds'])
end

%remove any output spikes for which memory has been allocated but never used
output_spikes.ts(output_spike_number:end)   = [];
output_spikes.x(output_spike_number:end)    = [];
output_spikes.y(output_spike_number:end)    = [];
output_spikes.sp(output_spike_number:end)   = [];
output_spikes.dir(output_spike_number:end)  = [];

%check that the spikes are in chronological order
[output_spikes.ts, order_index] = sort(output_spikes.ts);
output_spikes.x       = output_spikes.x(order_index);
output_spikes.y       = output_spikes.y(order_index);
output_spikes.sp      = output_spikes.sp(order_index);
output_spikes.dir     = output_spikes.dir(order_index);

%give an arbitrary to the spikes. Mostly for compatibility with other functions
output_spikes.p                 = ones(size(output_spikes.dir));
