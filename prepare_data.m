function data = prepare_data(file) % filename
    column = 19; %default value for steady state sampling, 19 physical properties are recorded
    header = 2; %default value for data produced by the Bird program
    raw = importdata(file);
    raw = raw((header+1):end,:);
    data = zeros(size(raw, 1), column);
    count = 0;
    tscan = @(x) textscan(x,'%.5f');%this is the default precision of Bird output files
    for i=1:length(raw)
        temp = split(raw(i));
        if length(temp)~= column + 1 %split gives all the columns plus a label
            count = count + 1;
            continue
        else
        temp(1) = [];
        %data(i-count,:) = str2double(temp)';  %technically could be used,
        %but double precision is somewhat hard to deal with
        data(i-count,:) = cell2mat(cellfun(tscan,temp))';
        end
    end
    data = data(1:size(data,1)-count, :);
end