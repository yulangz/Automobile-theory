clc
clear
k = 4080; hgk = 0.845; Lk = 3.950; ak = 2.10; betak = 0.38; bk = Lk - ak; %����ʱ�Ĳ���
mm = 9290; hgm = 1.170; Lm = 3.950; am = 2.950; betam = 0.38; bm = Lm - am; %����ʱ�Ĳ���
z = 0:0.01:1.0;
figure(1);
fai = z;
fai_fk = betak * z * Lk ./ (bk + z * hgk); %����ʱǰ��Ħ�f
fai_fm = betam * z * Lm ./ (bm + z * hgm); %����ʱǰ��Ħ�f
fai_rk = (1 - betak) * z * Lk ./ (ak - z * hgk); %����ʱ����Ħ�r
fai_rm = (1 - betam) * z * Lm ./ (am - z * hgm); %����ʱ����Ħ�r
plot(z, fai_fk, 'b--', z, fai_fm, 'r', z, fai_rk, 'b--', z, fai_rm, 'r', z, fai, 'k');
title('���ø���ϵ�����ƶ�ǿ�ȵĹ�ϵ����');
xlabel('�ƶ�ǿ��(z/g)');
ylabel('���ø���ϵ����');
gtext('��r(����)'), gtext('��r(����)'), gtext('��=z'), gtext('��f(����)'), gtext('��f(����)');
figure(2);
Efk = z ./ fai_fk * 100; %����ʱǰ����ƶ�Ч��
Efm = z ./ fai_fm * 100;
Erk = z ./ fai_rk * 100;
Erm = z ./ fai_rm * 100;
plot(fai_fk, Efk, 'b', fai_fm, Efm, 'r', fai_rk, Erk, 'b', fai_rm, Erm, 'r');
axis([0 1 0 100]);
title('ǰ�����ƶ�Ч������');
xlabel('����ϵ����');
ylabel('�ƶ�Ч��%');
gtext('Ef'), gtext('Er'), gtext('Er'), gtext('����'), gtext('����');
