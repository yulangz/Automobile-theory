clc
clear
m = 1818.2; Iz = 3885; L = 3.048; a = 1.463; b = 1.585; k1 = -62618; k2 = -110185;
i = 20; g = 9.8; R0 = 15; u1 = 30.56;
K = m * (a / k2 - b / k1) / L^2;
Uch = (1 / K)^(1/2); %特征车速
fprintf('稳定性因数(s^2/m^2)K=%f\n', K);
fprintf('特征车速(m/s)Uch=%f\n', Uch);
u = 0:0.05:30;
S = u ./ (L * (1 + K * u.^2)); %稳态横摆角速度增益
plot(u, S);
title('汽车稳态横摆角速度增益曲线');
xlabel('车速u(m/s)');
ylabel('稳态横摆角速度增益');
fprintf('u=22.35m/s时，转向灵敏度为%f\n', S(448));
SM = k2 / (k1 + k2) - a / L;
ay = 0.4 * g;
A = K * ay * L;
B = L / R0;
R = L / (B - A);
C = R / R0; %转弯半径比
fprintf('静态储备系数S.M.=%f\n', SM);
fprintf('侧向加速度为0.4g时前、后轮侧偏角绝对值之差(rad) a1-a2=%f\n', A);
fprintf('侧向加速度为0.4g时转弯半径比值R/R0=%f\n', C);
W0 = L / u1 * (k1 * k2 / (m * Iz) * (1 + K * u1^2))^(1/2); %固有（圆）频率
D = (-m * (k1 * a^2 + k2 * b^2) - Iz * (k1 + k2)) / (2 * L * (m * Iz * k1 * k2 * (1 + K * u1^2))^(1/2)); %阻尼比
t = atan((1 - D^2)^(1/2) / (-m * u1 * a * W0 / (L * k2) - D)) / (W0 * (1 - D^2)^(1/2)); %反应时间
E = atan((1 - D^2)^(1/2) / D) / (W0 * (1 - D^2)^(1/2)) + t; %峰值反应时间
fprintf('\n车速u=30.56m/s时的瞬态响应参数分别为:\n');
fprintf('横摆角速度波动的固有(圆)频率(rad)为 %f\n', W0);
fprintf('阻尼比为 %f\n', D);
fprintf('反应时间(s)为 %f\n', t);
fprintf('峰值反应时间(s)为 %f\n', E);
