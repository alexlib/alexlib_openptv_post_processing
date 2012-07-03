function data = ensure_t_length(data)
fields = fieldnames(data);
for i = 1:length(data)
    data(i).t = repmat(data(i).t,length(data(i).(fields{1})),1);
end

