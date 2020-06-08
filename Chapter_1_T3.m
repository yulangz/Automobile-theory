clear
clc

Loading_weight = 2000;
Curb_weight = 1800;
Total_weight = 3880;
r = 0.367;
Eta_T = 0.85;
f = 0.013;
CDA = 2.77;
i0 = 5.83;
ig = [5.56 2.769 1.644 1.00 0.793];
If = 0.218;
Iw1 = 1.798;
Iw2 = 3.598;
L = 3.2;
a = 1.974;
hg = 0.9;
nmin = 600;
nmax = 4000;




out = "";

T_use = 0;
V_now = 0;
is_find_70 = false;
for i = 1:5
    now_ig = ig(i);
    n = nmin:1:nmax;
    Ttq = -19.313 + 295.27 * (n / 1000) - 165.44 * (n / 1000).^2 + 40.874 * (n / 1000).^3 - 3.8445 * (n / 1000).^4;
    ua = 0.377 * r * n / (now_ig * i0);
    Ft = Ttq * now_ig * i0 * Eta_T / r;


    Ff = Total_weight * 9.8 * (f + 0 * ua);
    Fw = CDA / 21.15 * ua.^2;

    % 加速度曲线
    Derta = 1 + (Iw1 + Iw2) / (Total_weight * r^2) + If * now_ig^2 * i0^2 * Eta_T / (Total_weight * r^2);
    acclerate = (Ft - Ff - Fw) / (Derta * Total_weight);
    acc_1 = 1./acclerate;
    good_index =find(acc_1 >= 10); % 删去较小值
    if isempty(good_index)
        [~ ,good_index] = size(acclerate);
    else
        good_index = good_index(1);
    end
    figure(1)
    plot(ua(1:good_index) ,acc_1(1:good_index));
    hold on
    
    % 累计计算加速时间
    uv_index = find(ua>=V_now);
    uv_index = uv_index(1);
    acc_1_use = acc_1(uv_index:end);
    ua_use = ua(uv_index:end); % m/s
    if i ~=1 % 从二挡开始累积
        T_use = (T_use + cumtrapz(ua_use./3.6, acc_1_use));
        if ~is_find_70
            figure(2)
            plot(T_use, ua_use);
            hold on
            t_index = find(ua_use>=70); % 寻找累积速度是否到达70
            if t_index
                is_find_70 = true;
                t_index = t_index(1);
                plot(T_use(t_index), ua_use(t_index), 'o');
                hold on
                out = sprintf("%s在%d档时加速到70km/h，用时%.2f秒\n",out, i, T_use(t_index));
            end
        end
    T_use = T_use(end);
    V_now = ua_use(end);
    end
end

figure(1)
legend('1/a1', '1/a2', '1/a3', '1/a4', '1/a5');
title('汽车加速度倒数图');
xlabel('ua/(km/h)');
ylabel('1/a');

figure(2)
title('二挡起步加速到70km/h');
xlabel('t/s')
ylabel('ua/(km/h)')
disp(out);