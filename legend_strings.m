function legend_strings(l)

labels = {};

for i=1:l
    labels{end+1} = sprintf('P%d',i);
    labels{end+1} = sprintf('V%d',i);
end

legend(labels);

end