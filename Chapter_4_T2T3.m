clc
clear
mk = 4080; hgk = 0.845; Lk = 3.950; ak = 2.10; betak = 0.38; bk = Lk - ak; %空载时的参数
mm = 9290; hgm = 1.170; Lm = 3.950; am = 2.950; betam = 0.38; bm = Lm - am; %满载时的参数
z = 0:0.01:1;
fai_fk = betak * z * Lk ./ (bk + z * hgk); %空载时前轴的φf
fai_fm = betam * z * Lm ./ (bm + z * hgm); %满载时前轴的φf
fai_rk = (1 - betak) * z * Lk ./ (ak - z * hgk); %空载时后轴的φr
fai_rm = (1 - betam) * z * Lm ./ (am - z * hgm); %满载时后轴的φr
Efk = z ./ fai_fk * 100; %空载时前轴的制动效率
Efm = z ./ fai_fm * 100;
Erk = z ./ fai_rk * 100;
Erm = z ./ fai_rm * 100;
t1 = 0.02; t2 = 0.02; ua0 = 30; fai = 0.80; g = 9.8;
ak1 = Erk(81) * g * fai / 100;
am1 = Erm(81) * g * fai / 100;
Sk1 = (t1 + t2 / 2) * ua0 / 3.6 + ua0^2 / (25.92 * ak1); %制动距离
Sm1 = (t1 + t2 / 2) * ua0 / 3.6 + ua0^2 / (25.92 * am1);
fprintf('空载时，汽车制动距离Sk1=%f\n', Sk1);
fprintf('满载时，汽车制动距离Sm1=%f\n', Sm1);
ak2 = fai * g * ak / (Lk + fai * hgk);
am2 = fai * g * am / (Lm + fai * hgm);
ak3 = fai * g * bk / (Lk - fai * hgk);
am3 = fai * g * bm / (Lk - fai * hgm);
Sk2 = (t1 + t2 / 2) * ua0 / 3.6 + ua0^2 / (25.92 * ak2); %制动距离
Sm2 = (t1 + t2 / 2) * ua0 / 3.6 + ua0^2 / (25.92 * am2);
Sk3 = (t1 + t2 / 2) * ua0 / 3.6 + ua0^2 / (25.92 * ak3);
Sm3 = (t1 + t2 / 2) * ua0 / 3.6 + ua0^2 / (25.92 * am3);
fprintf('空载时，前制动器损坏，汽车制动距离Sk2=%f\n', Sk2);
fprintf('满载时，前制动器损坏，汽车制动距离Sm2=%f\n', Sm2);
fprintf('空载时，后制动器损坏，汽车制动距离Sk3=%f\n', Sk3);
fprintf('满载时，后制动器损坏，汽车制动距离Sm3=%f\n', Sm3);
