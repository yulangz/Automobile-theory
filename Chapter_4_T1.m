clc
clear
k = 4080; hgk = 0.845; Lk = 3.950; ak = 2.10; betak = 0.38; bk = Lk - ak; %空载时的参数
mm = 9290; hgm = 1.170; Lm = 3.950; am = 2.950; betam = 0.38; bm = Lm - am; %满载时的参数
z = 0:0.01:1.0;
figure(1);
fai = z;
fai_fk = betak * z * Lk ./ (bk + z * hgk); %空载时前轴的φf
fai_fm = betam * z * Lm ./ (bm + z * hgm); %满载时前轴的φf
fai_rk = (1 - betak) * z * Lk ./ (ak - z * hgk); %空载时后轴的φr
fai_rm = (1 - betam) * z * Lm ./ (am - z * hgm); %满载时后轴的φr
plot(z, fai_fk, 'b--', z, fai_fm, 'r', z, fai_rk, 'b--', z, fai_rm, 'r', z, fai, 'k');
title('利用附着系数与制动强度的关系曲线');
xlabel('制动强度(z/g)');
ylabel('利用附着系数φ');
gtext('φr(空载)'), gtext('φr(满载)'), gtext('φ=z'), gtext('φf(空载)'), gtext('φf(满载)');
figure(2);
Efk = z ./ fai_fk * 100; %空载时前轴的制动效率
Efm = z ./ fai_fm * 100;
Erk = z ./ fai_rk * 100;
Erm = z ./ fai_rm * 100;
plot(fai_fk, Efk, 'b', fai_fm, Efm, 'r', fai_rk, Erk, 'b', fai_rm, Erm, 'r');
axis([0 1 0 100]);
title('前、后制动效率曲线');
xlabel('附着系数φ');
ylabel('制动效率%');
gtext('Ef'), gtext('Er'), gtext('Er'), gtext('满载'), gtext('空载');
