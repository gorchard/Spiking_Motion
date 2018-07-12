%SortedEvents =  SortOrder(inputEvents)
% sorts the inputEvents stream by temporal order using the 'ts' field
% This function is useful for inserting new events into the event stream and
% then automatically reordering the stream to ensure time is monotonically
% increasing
function result =  sort_order(result)

if isfield(result, 'ts')
    [~, order] = sort(result.ts);
elseif isfield(result, 'T')
    [~, order] = sort(result.T);
elseif isfield(result, 'TimeStamp')
    [~, order] = sort(result.TimeStamp);
end

fieldnames = fields(result);
for i = 1:length(fieldnames)
    if ~strcmp(fieldnames{i}, 'meta')
        result.(fieldnames{i}) = result.(fieldnames{i})(order);
    end
end