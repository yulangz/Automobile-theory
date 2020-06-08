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
G = Total_weight * 9.8


for i = 1:5
    now_ig = ig(i);
    n = nmin:1:nmax;
    Ttq = -19.313 + 295.27 * (n / 1000) - 165.44 * (n / 1000).^2 + 40.874 * (n / 1000).^3 - 3.8445 * (n / 1000).^4;
    ua = 0.377 * r * n / (now_ig * i0);

    % 驱动力、阻力与功率
    Ft = Ttq * now_ig * i0 * Eta_T / r;
    Pe = Ft .* ua ./ 3600;
    figure(1)
    h1 = plot(ua, Pe, 'red');
    hold on
    Ff = Total_weight * 9.8 * (f + 0 * ua);
    Fw = CDA / 21.15 * ua.^2;
    Pf_Pw = (1 / Eta_T) * (Ff .* ua ./ 3600 + Fw .* ua ./ 3600);
    h2 = plot(ua, Pf_Pw, 'black');
    legend([h1,h2],"Pe 1-5","(Pf+Pw)/Eta");
    hold on

end

figure(1)
title('汽车功率平衡图');
xlabel('ua/(km/h)');
ylabel('Pe/kW');
disp(out);
