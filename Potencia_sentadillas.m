Sent2 = load('Sentadilla2.mat') % carga de archivo con registro de sensores
x = Sent2.Acceleration.X;
y = Sent2.Acceleration.Y;
z = Sent2.Acceleration.Z;
timestamp = Sent2.Acceleration.Timestamp;
combinedAccel = sqrt(x.^2+y.^2+z.^2);
for n = 1 : length(combinedAccel)
  Acceloff(n,1) = combinedAccel(n)-9.81;
end
t = zeros(size(timestamp));
for n = 1 : length(timestamp)
  t(n) = seconds(timestamp(n) - timestamp(1));
end
figure (1)
subplot(2,1,1)
plot(t, Acceloff)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
xlim([3.5 6])
%Registro de datos
m = 519.75; %Newtons 53 k

%Calculo Fuerza
F = m.*Acceloff;

%Calculo Potencia
V = trapz(Acceloff);
P = V * abs(F);
subplot(2,1,2)
plot(t, P)
xlim([3.5 6])

Sent1 = load('Sentadilla1.mat') % carga de archivo con registro de sensores
x = Sent1.Acceleration.X;
y = Sent1.Acceleration.Y;
z = Sent1.Acceleration.Z;
timestamp = Sent1.Acceleration.Timestamp;
combinedAccel = sqrt(x.^2+y.^2+z.^2);
for n = 1 : length(combinedAccel)
  Acceloff(n,1) = combinedAccel(n)-9.81;
end
t = zeros(size(timestamp));
for n = 1 : length(timestamp)
  t(n) = seconds(timestamp(n) - timestamp(1));
end
figure (2)
subplot(2,1,1)
plot(t, Acceloff)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
xlim([3.5 6])
%Registro de datos
m = 519.75; %Newtons 53 k

%Calculo Fuerza
F = m.*Acceloff;

%Calculo Potencia
V = trapz(Acceloff);
P = V * abs(F);
subplot(2,1,2)
plot(t, P)
xlim([3.5 6])
