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


for i = 1:5
    now_ig = ig(i);
    n = nmin:1:nmax;
    Ttq = -19.313 + 295.27 * (n / 1000) - 165.44 * (n / 1000).^2 + 40.874 * (n / 1000).^3 - 3.8445 * (n / 1000).^4;
    ua = 0.377 * r * n / (now_ig * i0);
    
    % 驱动力与行驶阻力
    Ft = Ttq * now_ig * i0 * Eta_T / r;
    figure(1)
    h1 = plot(ua, Ft, 'red');
    hold on
    Ff = Total_weight * 9.8 * (f + 0 * ua);
    Fw = CDA / 21.15 * ua.^2;
    h2 = plot(ua, Ff, 'blue');
    hold on
    h3 = plot(ua, Ff + Fw, 'black');
    hold on
    legend([h1,h2,h3],"Ft 1-5","Ff","Ff+Fw");
    
    if i==5
        % 寻找最高车速
        pos = find(Ft <= Ff+Fw);
        UMAX = ua(pos(1) - 1);
        x = UMAX;
        y = Ft(pos(1) - 1);

        figure(1)
        plot(x, y ,'ro');
        hold on
        str = sprintf("最高车速:%.2f",UMAX);
        out = sprintf("%s%s\n",out, str);
        text(x-10,y + 800,str);
    end

    if i==1
        % 寻找最大爬坡度
        max_Ft = max(Ft);
        index = find(max_Ft == Ft);
        alpha = asind( (max_Ft  - Ff(index) - Fw(index))  / (Total_weight * 9.8));
        out = sprintf("%s最大爬坡度:%.2f°\n", out, alpha);
        Fai = deg2rad(alpha) * L / a;
        out = sprintf("%s对应附着率为:%.2f\n", out, Fai);
    end
end

figure(1)
title('汽车驱动力-行驶阻力平衡图');
xlabel('ua/(km/h)');
ylabel('Ft/n, (Ff+Fw)/n');
disp(out);



