clc
clear
i0 = [5.17 5.43 5.83 6.17 6.33];
t0 = zeros(1, 5);
Q0 = zeros(1, 5);

for i = 1:5
    [t0(i), Q0(i)] = AT1and2(i0(i));
end

plot(Q0, t0, 'o');
hold on
xx = 12.8:0.001:13.8;
yy = spline(Q0, t0, xx);
plot(xx, yy, 'linewidth', 1);

for i = 1:5
    str = sprintf("i%d=%f", i, i0(i));
    text(Q0(i), t0(i), str);
    hold on
end

title('燃油经济性-加速时间');
xlabel('百公里油耗(L/100km)');
ylabel('加速时间s');
