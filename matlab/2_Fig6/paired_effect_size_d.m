function dz = paired_effect_size_d(x, y)

d = x - y;
d = d(~isnan(d));

if numel(d) < 2
    dz = NaN;
    return;
end

sd = std(d);
if sd == 0
    dz = 0;
else
    dz = mean(d) / sd;
end

end