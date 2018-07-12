function events = implement_refraction(events, refractory_period)
%function events = implement_refraction(events, refractory_period)
%   events              -> a struct of events
%   refractory_period   -> the refractory period to implement (same units as "event.ts")

%an array to hold the last spike time for every pixel
last_time = ones(max(events.x), max(events.y))-(refractory_period+1);

%loop through all the events
for event_num = 1:length(events.ts)
    
    %check if the event is greater than "refractory_period" in time since the last event for this pixel
    if ((events.ts(event_num) - last_time(events.x(event_num), events.y(event_num))) > refractory_period)
        last_time(events.x(event_num), events.y(event_num)) = events.ts(event_num); %if so, update the time of the latest event at this pixel
    else
        events.ts(event_num) = 0; %otherwise set the event time to zero (which is used as a marker to later remove events)
    end
end

%remove any events with the timestamp equal to 0
invalid_indices = (events.ts == 0);
events = remove_nulls(events, invalid_indices); 