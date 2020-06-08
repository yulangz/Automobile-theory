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
G = Total_weight * 9.8;

n_ = [815 1207 1614 2012 2603 3006 3403 3804];
B_ = [1326.8 1354.7 1284.4 1122.9 1141.0 1051.2 1233.9 1129.7;
    -416.46 -303.98 -189.75 -121.59 -98.893 -73.714 -84.478 -45.291;
    72.379 36.657 14.524 7.0035 4.4763 2.8593 2.9788 0.71113;
    -5.8629 -2.0553 -0.51184 -0.18517 -0.091077 -0.05138 -0.047449 -0.00075215;
    0.17768 0.043072 0.0068164 0.0018555 0.00068906 0.00035032 0.00028230 -0.000038568];

roug = 7;

out = ""; % 文字输出

% 求油耗拟合公式
for i = 1:5
    now_ig = ig(i);
    n = nmin:1:nmax;
    ua = 0.377 * r * n / (now_ig * i0);
    Ttq = -19.313 + 295.27 * (n / 1000) - 165.44 * (n / 1000).^2 + 40.874 * (n / 1000).^3 - 3.8445 * (n / 1000).^4;

    % 利用插值法求拟合公式
    ua_ = 0.377 * r * n_ / (now_ig * i0);
    Pe = (G * f * ua_ / 3600 + CDA * ua_.^3/76140) / Eta_T;
    b_ = zeros(1, 8);

    for k = 1:8
        b_(k) = B_(:, k)' * [1; Pe(k); Pe(k)^2; Pe(k)^3; Pe(k)^4];
    end

    ua2b(i, :) = polyfit(ua_, b_, 4);

end

% 计算六工况循环耗油量
% 根据功率平衡图分析，3,4,5档均符合六工况循环条件，顾分别计算
for i = 3:5
    Q = [];
    % 工况1：  25km/h等速
    s = 50; % 行程
    ua = 25; % 起始速度
    b = polyval(ua2b(i, :), ua);
    Ff = G * (f + 0 * ua);
    Fw = CDA / 21.15 * ua.^2;
    Pf_Pw = (1 / Eta_T) * (Ff * ua / 3600 + Fw * ua / 3600);
    Pe = Pf_Pw;
    Q(1) = Pe * b * s / (102 * ua * roug);

    % 工况2
    s = 150; % 行程
    ua = 25; % 起始速度
    acc = 0.25; % 加速度

    Det_t = 1 / (3.6 * acc);
    t_ = 7.2:Det_t:23.9;
    [~, tsize] = size(t_);

    Qs = [];
    Qt = [];

    for j = 1:tsize
        % 因为功率计算公已变，故需重新插值
        Derta = 1 + (Iw1 + Iw2) / (Total_weight * r^2) + If * now_ig^2 * i0^2 * Eta_T / (Total_weight * r^2);
        Pe = (1 / Eta_T) * (G * (f + 0 * ua) / 3600 + CDA * ua^3/76140 + Derta * Total_weight * ua * acc / 3600);

        now_ig = ig(i);

        ua__ = 0.377 * r * n_ / (now_ig * i0);
        Pe_ = (1 / Eta_T) * (G * (f + 0 * ua__) / 3600 + CDA * ua__.^3/76140 + Derta * Total_weight * ua__ * acc / 3600);
        b__ = zeros(1, 8);
        for k = 1:8
            b__(k) = B_(:, k)' * [1; Pe_(k); Pe_(k)^2; Pe_(k)^3; Pe_(k)^4];
        end
        ua2bb(:) = polyfit(ua__, b__, 4);

        b = polyval(ua2bb(:), ua);
        Qt(1) = Pe * b / (361.7 * roug);

        ua = ua + 1;

        Derta = 1 + (Iw1 + Iw2) / (Total_weight * r^2) + If * now_ig^2 * i0^2 * Eta_T / (Total_weight * r^2);
        Pe = (1 / Eta_T) * (G * (f + 0 * ua) / 3600 + CDA * ua^3/76140 + Derta * Total_weight * ua * acc / 3600);

        now_ig = ig(i);

        ua__ = 0.377 * r * n_ / (now_ig * i0);
        Pe_ = (1 / Eta_T) * (G * (f + 0 * ua__) / 3600 + CDA * ua__.^3/76140 + Derta * Total_weight * ua__ * acc / 3600);
        b__ = zeros(1, 8);
        for k = 1:8
            b__(k) = B_(:, k)' * [1; Pe_(k); Pe_(k)^2; Pe_(k)^3; Pe_(k)^4];
        end
        ua2bb(:) = polyfit(ua__, b__, 4);
        
        b = polyval(ua2bb(:), ua);
        Qt(2) = Pe * b / (361.7 * roug);

        Qs(j) = 1/2 * (Qt(2) + Qt(1)) * Det_t;
    end

    Q(2) = sum(Qs);

    % 工况3
    s = 250; % 行程
    ua = 40; % 起始速度
    b = polyval(ua2b(i, :), ua);
    Ff = G * (f + 0 * ua);
    Fw = CDA / 21.15 * ua.^2;
    Pf_Pw = (1 / Eta_T) * (Ff * ua / 3600 + Fw * ua / 3600);
    Pe = Pf_Pw;
    Q(3) = Pe * b * s / (102 * ua * roug);

    % 工况4
    s = 175; % 行程
    ua = 40; % 起始速度
    acc = 0.20; % 加速度

    Det_t = 1 / (3.6 * acc);
    t_ = 46.4:Det_t:60.4;
    [~, tsize] = size(t_);

    Qs = [];
    Qt = [];

    for j = 1:tsize
        % 重新插值
        Derta = 1 + (Iw1 + Iw2) / (Total_weight * r^2) + If * now_ig^2 * i0^2 * Eta_T / (Total_weight * r^2);
        Pe = (1 / Eta_T) * (G * (f + 0 * ua) / 3600 + CDA * ua^3/76140 + Derta * Total_weight * ua * acc / 3600);

        now_ig = ig(i);

        ua__ = 0.377 * r * n_ / (now_ig * i0);
        Pe_ = (1 / Eta_T) * (G * (f + 0 * ua__) / 3600 + CDA * ua__.^3/76140 + Derta * Total_weight * ua__ * acc / 3600);
        b__ = zeros(1, 8);
        for k = 1:8
            b__(k) = B_(:, k)' * [1; Pe_(k); Pe_(k)^2; Pe_(k)^3; Pe_(k)^4];
        end
        ua2bb(:) = polyfit(ua__, b__, 4);

        b = polyval(ua2bb(:), ua);
        Qt(1) = Pe * b / (361.7 * roug);

        ua = ua + 1;

        Derta = 1 + (Iw1 + Iw2) / (Total_weight * r^2) + If * now_ig^2 * i0^2 * Eta_T / (Total_weight * r^2);
        Pe = (1 / Eta_T) * (G * (f + 0 * ua) / 3600 + CDA * ua^3/76140 + Derta * Total_weight * ua * acc / 3600);

        now_ig = ig(i);

        ua__ = 0.377 * r * n_ / (now_ig * i0);
        Pe_ = (1 / Eta_T) * (G * (f + 0 * ua__) / 3600 + CDA * ua__.^3/76140 + Derta * Total_weight * ua__ * acc / 3600);
        b__ = zeros(1, 8);
        for k = 1:8
            b__(k) = B_(:, k)' * [1; Pe_(k); Pe_(k)^2; Pe_(k)^3; Pe_(k)^4];
        end
        ua2bb(:) = polyfit(ua__, b__, 4);

        b = polyval(ua2bb(:), ua);
        Qt(2) = Pe * b / (361.7 * roug);

        Qs(j) = 1/2 * (Qt(2) + Qt(1)) * Det_t;
    end

    Q(4) = sum(Qs);

    % 工况5
    s = 250; % 行程
    ua = 50; % 起始速度
    b = polyval(ua2b(i, :), ua);
    Ff = G * (f + 0 * ua);
    Fw = CDA / 21.15 * ua.^2;
    Pf_Pw = (1 / Eta_T) * (Ff * ua / 3600 + Fw * ua / 3600);
    Pe = Pf_Pw;
    Q(5) = Pe * b * s / (102 * ua * roug);

    % 工况6
    t = 19.3;
    Qi = 0.299;
    Q(6) = t*Qi;

    out = sprintf("%s第%d档耗油量为%.2fL\n", out, i, sum(Q) * 100/1075);
end

disp(out);
