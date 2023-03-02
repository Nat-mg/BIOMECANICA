%Integrantes
% Ana Lucía Soria Cardona A00827565         
% Cecilia Orozco Romo A00827013
% Natalia Montemayor Galván A00828262
% Luis Orlando Santos Cruz A00827603
% Karina Esquer Vila A00828785
clear all;
clc;

Data = dataImport('Orlando_p.csv');

markers=string(Data.Properties.VariableNames).';
idx=contains(markers,'x');
markers=erase(markers(idx),'x');

clear idx
 %% Plottin all
 
 figure (1)
 hold on
for f= 1:height(Data)
    cla
   for m= 1:length(markers)
       labelx= strcat(markers(m),'x');
       labely= strcat(markers(m),'y');
       labelz= strcat(markers(m),'z');
       marker= [Data.(labelx) Data.(labely) Data.(labelz)];
       plot3(marker(f,1),marker(f,2),marker(f,3),'k*')
   end
   axis([-1200 1200 -800 1600 0 2400])
   set(gca, 'view', [-132.3000 16.8000])
   %pause(0.01)
end
hold off

%% Fill gaps

gaps= ismissing(Data);

CutData= Data(172:881,:);

figure (2)
hold on
for f= 1:height(CutData)
    cla
   for m= 1:length(markers)
       labelx= strcat(markers(m),'x');
       labely= strcat(markers(m),'y');
       labelz= strcat(markers(m),'z');
       marker= [CutData.(labelx) CutData.(labely) CutData.(labelz)];
       plot3(marker(f,1),marker(f,2),marker(f,3),'k*')
   end
   axis([-1200 1200 -800 1600 0 2400])
   set(gca, 'view', [-132.3000 16.8000])
   % pause(0.01)
end
hold off

%% ismissing 
frames= CutData.Frame;
for m= 1:length(markers)
    labelx = strcat(markers(m),'x');
    idx= ~ismissing(CutData.(labelx));
    if sum(ismissing(CutData.(labelx))) > 0
        original_data = CutData.(labelx);
        new_data= spline(frames(idx),original_data(idx),frames);
        CutData.(labelx) = new_data;
    end
    labely = strcat(markers(m),'y');
    idx = ~ismissing(CutData.(labely));
    if sum(ismissing(CutData.(labely))) > 0 
        original_data = CutData.(labely);
        new_data= spline(frames(idx),original_data(idx),frames);
        CutData.(labely) = new_data;
    end
    labelz = strcat(markers(m),'z');
    idx = ~ismissing(CutData.(labelz));
    if sum(ismissing(CutData.(labelz))) > 0
        original_data = CutData.(labelz);
        new_data= spline(frames(idx),original_data(idx),frames);
        CutData.(labelz) = new_data;
    end
end

clear labelx labely labelz original_data new_data

%% coordinate system creation
%Pierna Derecha
%MD
md1= [CutData.Placa_MD1x CutData.Placa_MD1y CutData.Placa_MD1z]; %load variables
md2= [CutData.Placa_MD2x CutData.Placa_MD2y CutData.Placa_MD2z];
md3= [CutData.Placa_MD3x CutData.Placa_MD3y CutData.Placa_MD3z];
md4= [CutData.Placa_MD4x CutData.Placa_MD4y CutData.Placa_MD4z];

%PD
pd1= [CutData.Placa_PD1x CutData.Placa_PD1y CutData.Placa_PD1z]; %load variables
pd2= [CutData.Placa_PD2x CutData.Placa_PD2y CutData.Placa_PD2z];
pd3= [CutData.Placa_PD3x CutData.Placa_PD3y CutData.Placa_PD3z];
pd4= [CutData.Placa_PD4x CutData.Placa_PD4y CutData.Placa_PD4z];

%Pierna izquierda
%MI
mi1= [CutData.Placa_MI1x CutData.Placa_MI1y CutData.Placa_MI1z]; %load variables
mi2= [CutData.Placa_MI2x CutData.Placa_MI2y CutData.Placa_MI2z];
mi3= [CutData.Placa_MI3x CutData.Placa_MI3y CutData.Placa_MI3z];
mi4= [CutData.Placa_MI4x CutData.Placa_MI4y CutData.Placa_MI4z];

%PD
pi1= [CutData.Placa_PI1x CutData.Placa_PI1y CutData.Placa_PI1z]; %load variables
pi2= [CutData.Placa_PI2x CutData.Placa_PI2y CutData.Placa_PI2z];
pi3= [CutData.Placa_PI3x CutData.Placa_PI3y CutData.Placa_PI3z];
pi4= [CutData.Placa_PI4x CutData.Placa_PI4y CutData.Placa_PI4z];

% Creating MD vector for Coordinate System
%Pierna Derecha
MDY=[md3-md4];%position vector Local MDY
MDZ=[md2-md4];%position vector Local MDZ
MDX=cross(MDY,MDZ);%crossproduct vector Local MDX
PDY=[pd3-pd4];%position vector Local PDY
PDZ=[pd2-pd4];%position vector Local PDZ
PDX=cross(PDY,PDZ);%crossproduct vector Local PDX

%Pierna Izquierda
MIY=[mi4-mi3];%position vector Local MDY
MIZ=[mi1-md3];%position vector Local MDZ
MIX=cross(MIY,MIZ);%crossproduct vector Local MDX
PIY=[pi4-pi3];%position vector Local PDY
PIZ=[pi1-pi3];%position vector Local PDZ
PIX=cross(PIY,PIZ);%crossproduct vector Local PDX

%empty variable to allocate the angle in each frame
k_ang = nan(height(Data),3);%
%Crear variable para guardar valores de pierna derecha
k_ang_D = k_ang; 
%Crear variable para guardar valores de pierna derecha
k_ang_I = k_ang; 

for f = 1:height(frames)%
    %create coordinate system matrices
    %Pierna Derecha
     MD = [MDX(f,:).', MDY(f,:).', MDZ(f,:).'];
     PD = [PDX(f,:).', PDY(f,:).', PDZ(f,:).'];
     %Pierna Izquierda
     MI = [MIX(f,:).', MIY(f,:).', MIZ(f,:).'];
     PI = [PIX(f,:).', PIY(f,:).', PIZ(f,:).'];

    %calculate the rotation matrix
    %Pierna derecha
    RD = MD*PD';
    %Pierna derecha
    RI = MI*PI';
    
    %Alpha Derecha
    k_ang_D(f,1)= atan2d(RD(3,2),RD(3,3));
    %Betha Derecha
    k_ang_D(f,2)= atan2d(RD(2,1),RD(1,1));
    %Gamma Derecha
    if k_ang_D(f,2) == 0
        k_ang_D(f,3)= atan2d(-RD(3,1),RD(2,1)/sin(k_ang_D(f,2)));
    else
        k_ang_D(f,3)= atan2d(-RD(3,1),RD(1,1)/cos(k_ang_D(f,2)));
    end
    %Alpha Izquierda
    k_ang_I(f,1)= atan2d(RI(3,2),RI(3,3));
    %Betha Izquierda
    k_ang_I(f,2)= atan2d(RI(2,1),RI(1,1));
    %Gamma Izquierda
    if k_ang_I(f,2) == 0
        k_ang_I(f,3)= atan2d(-RI(3,1),RI(2,1)/sin(k_ang_I(f,2)));
    else
        k_ang_I(f,3)= atan2d(-RI(3,1),RI(1,1)/cos(k_ang_I(f,2)));
    end
end
hold off

%visualization Pierna Derecha
figure(3)
subplot(2,1,1)
hold on
plot([MD(1,1) 0 MD(1,2)],[MD(2,1) 0 MD(2,2)],'k')
plot([PD(1,1) 0 PD(1,2)],[PD(2,1) 0 PD(2,2)],'r')
axis([-1 1 -1 1])
title(string(k_ang_D(f)))
subtitle('Right Knee angle')
pause(0.1)
hold off

%visualization Pierna Izquierda
subplot(2,1,2)
hold on
plot([MI(1,1) 0 MI(1,2)],[MI(2,1) 0 MI(2,2)],'k')
plot([PI(1,1) 0 PI(1,2)],[PI(2,1) 0 PI(2,2)],'r')
axis([-1 1 -1 1])
title(string(k_ang_I(f)))
subtitle('Left Knee angle')
pause(0.1)
hold off

%Plot knee angle vs time right
figure (5)
t= [Data.Frame];
subplot(2,1,1)
plot(t,k_ang_D,'LineWidth',1.5)
title('Right Knee angle vs time')
xlabel('Time [s]')
ylabel('Knee angles [°]')

%Plot knee angle vs time left
t= [Data.Frame];
subplot(2,1,2)
plot(t,k_ang_I,'LineWidth',1.5)
title('Left Knee angle vs time')
xlabel('Time [s]')
ylabel('Knee angles [°]')