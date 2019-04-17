r1 = [800 50];
r2 = [800 50];
r3 = [800 50];

c1 = [1e-12 2.3e-11];
c2 = [1e-12 2.3e-11];
c3 = [1e-12 2.3e-11];
c4 = [1e-12 2.3e-11];
l1 = [1e-9 4.1e-8];
l2 = [1e-9 4.1e-8];

lamdag = 3/40; 
tx1=rfckt.cpw('Height', 0.635e-3, 'EpsilonR', 10.8, 'ConductorWidth', 0.5e-3, 'LineLength', lamdag/2);
tx2=rfckt.cpw('Height', 0.635e-3, 'EpsilonR', 10.8, 'ConductorWidth', 0.5e-3, 'LineLength', lamdag/4);
f = 2e9:1.0e7:6e9;
s = 2*pi*1i*f;

analyze(tx1, f);
analyze(tx2, f);
sp1 = sparameters(tx1);
sp2 = sparameters(tx2);
tp1 = s2abcd(sp1.Parameters, 50);
tp2 = s2abcd(sp2.Parameters, 50);

halfa = 1:401;
halfb = 1:401;
halfc = 1:401;
halfd = 1:401;

for i = 1:length(f)
    halfa(i) = tp1(1,1,i);
end
for i = 1:length(f)
    halfb(i) = tp1(1,2,i);
end
for i = 1:length(f)
    halfc(i) = tp1(2,1,i);
end
for i = 1:length(f)
    halfd(i) = tp1(2,2,i);
end

quara = 1:401;
quarb = 1:401;
quarc = 1:401;
quard = 1:401;
for i = 1:length(f)
    quara(i) = tp2(1,1,i); 
end
for i = 1:length(f)
    quarb(i) = tp2(1,2,i); 
end
for i = 1:length(f)
    quarc(i) = tp2(2,1,i); 
end
for i = 1:length(f)
    quard(i) = tp2(2,2,i); 
end

a = r1(1) + s .* l1(1);
b = 1.0 ./ (s * c1(1));
c = (s .* r3(1) .* c3(1) + 1) ./ (s .* s .* r3(1) .* c3(1) .* c4(1) + s .* c3(1) + s .* c4(1));

q = r2(1) + s .* l2(1);
h = 1.0 ./ (s * c2(1));

A = (a ./ b) + (a + (q ./ b + 1) .* c) ./ q + 1;
C = (c ./ b + 1) ./ q + 1 ./ b;
D = c ./ b + h .* ((c ./ b + 1) ./ q + 1 ./ b) + 1;
B = a + c .* (a ./ b + 1) + h .* (a ./ b + (a + c .* (a ./ b + 1)) ./ q + 1);

zo = 50;

para = zeros(1,length(f));
parb = zeros(1,length(f));
parc = zeros(1,length(f));
pard = zeros(1,length(f));
for i = 1 : length(f)
    Q = [quara(i) quarb(i);quarc(i) quard(i)];
    mul = [A(i) B(i);C(i) D(i)];
    mul = Q * mul * Q;
    para(i) = mul(1,1);
    parb(i) = mul(1,2);
    parc(i) = mul(2,1);
    pard(i) = mul(2,2);
end
[ya, yb, yc, yd] = abcdtoy(para, parb, parc, pard);
[yya, yyb, yyc, yyd] = abcdtoy(halfa, halfb, halfc, halfd);
ya = ya + yya;
yb = yb + yyb;
yc = yc + yyc;
yd = yd + yyd;
[A, B, C, D] = ytoabcd(ya, yb, yc, yd);

s11 = (A + B ./ zo - C .* zo - D) ./ (A + B ./ zo + C .* zo + D);
s21 = 2 ./ (A + B ./ zo + C .* zo + D);
figure;
plot(f, log(abs(s11)));
title('s11 OFF')
xlabel('f') 
ylabel('s11')
figure;
plot(f, log(abs(s21)));
title('s21 OFF')
xlabel('f') 
ylabel('s21')

a = r1(2) + s .* l1(2);
b = 1.0 ./ (s * c1(2));
c = (s .* r3(2) .* c3(2) + 1) ./ (s .* s .* r3(2) .* c3(2) .* c4(2) + s .* c3(2) + s .* c4(2));

q = r2(2) + s .* l2(2);
h = 1.0 ./ (s * c2(2));

A = (a ./ b) + (a + (q ./ b + 1) .* c) ./ q + 1;
C = (c ./ b + 1) ./ q + 1 ./ b;
D = c ./ b + h .* ((c ./ b + 1) ./ q + 1 ./ b) + 1;
B = a + c .* (a ./ b + 1) + h .* (a ./ b + (a + c .* (a ./ b + 1)) ./ q + 1);

zo = 150;

para = zeros(1,length(f));
parb = zeros(1,length(f));
parc = zeros(1,length(f));
pard = zeros(1,length(f));
for i = 1 : length(f)
    Q = [quara(i) quarb(i);quarc(i) quard(i)];
    mul = [A(i) B(i);C(i) D(i)];
    mul = Q * mul * Q;
    para(i) = mul(1,1);
    parb(i) = mul(1,2);
    parc(i) = mul(2,1);
    pard(i) = mul(2,2);
end
[ya, yb, yc, yd] = abcdtoy(para, parb, parc, pard);
[yya, yyb, yyc, yyd] = abcdtoy(halfa, halfb, halfc, halfd);
ya = ya + yya;
yb = yb + yyb;
yc = yc + yyc;
yd = yd + yyd;
[A, B, C, D] = ytoabcd(ya, yb, yc, yd);

s11 = (A + B ./ zo - C .* zo - D) ./ (A + B ./ zo + C .* zo + D);
s21 = 2 ./ (A + B ./ zo + C .* zo + D);
figure;
plot(f, log(abs(s11)));
title('s11 ON')
xlabel('f') 
ylabel('s11')
figure;
plot(f, log(abs(s21)));
title('s21 ON')
xlabel('f')
ylabel('s21')