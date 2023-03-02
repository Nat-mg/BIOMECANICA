%Data extraction
clc
clear all

Data = dataImport('AlejandroC2.csv');

markers = string(Data.Properties.VariableNames).';
idx = contains(markers,'x');
markers = erase(markers(idx),'x');

clear idx

%% Empti matriz variables
NewData.Frame = Data.Frame;
%Cadera alta
NewData.LASIx = Data.LASIx; %adelante izq
NewData.LASIy = Data.LASIy;
NewData.LASIz = Data.LASIz;

NewData.RASIx = Data.RASIx; %adelante der
NewData.RASIy = Data.RASIy;
NewData.RASIz = Data.RASIz;

NewData.LPSIx = Data.LPSIx; %atras izq
NewData.LPSIy = Data.LPSIy;
NewData.LPSIz = Data.LPSIz;

NewData.RPSIx = Data.RPSIx; %atras der
NewData.RPSIy = Data.RPSIy;
NewData.RPSIz = Data.RPSIz;

%Pierna Izquierda
NewData.Placa_MI1x = Data.Placa_MI1x;
NewData.Placa_MI1y = Data.Placa_MI1y;
NewData.Placa_MI1z = Data.Placa_MI1z;
NewData.Placa_MI2x = Data.Placa_MI2x;
NewData.Placa_MI2y = Data.Placa_MI2y;
NewData.Placa_MI2z = Data.Placa_MI2z;
NewData.Placa_MI3x = Data.Placa_MI3x;
NewData.Placa_MI3y = Data.Placa_MI3y;
NewData.Placa_MI3z = Data.Placa_MI3z;
NewData.Placa_MI4x = Data.Placa_MI4x;
NewData.Placa_MI4y = Data.Placa_MI4y;
NewData.Placa_MI4z = Data.Placa_MI4z;

%Pierna Derecha
NewData.Placa_MD1x = Data.Placa_MD1x;
NewData.Placa_MD1y = Data.Placa_MD1y;
NewData.Placa_MD1z = Data.Placa_MD1z;
NewData.Placa_MD2x = Data.Placa_MD2x;
NewData.Placa_MD2y = Data.Placa_MD2y;
NewData.Placa_MD2z = Data.Placa_MD2z;
NewData.Placa_MD3x = Data.Placa_MD3x;
NewData.Placa_MD3y = Data.Placa_MD3y;
NewData.Placa_MD3z = Data.Placa_MD3z;
NewData.Placa_MD4x = Data.Placa_MD4x;
NewData.Placa_MD4y = Data.Placa_MD4y;
NewData.Placa_MD4z = Data.Placa_MD4z;
%Covertir structura a tabla
%NewData = struct2table(NewData)

%Calculo de punto medio marcadores cadera
%Punto medio cintura
PMLARPx = (NewData.LASIx(1) + NewData.RPSIx(1))/2;
PMLARPy = (NewData.LASIy(1) + NewData.RPSIy(1))/2;
PMLARPz = (NewData.LASIz(1) + NewData.RPSIz(1))/2;
PMRALPx = (NewData.LPSIx(1) + NewData.RASIx(1))/2;
PMRALPy = (NewData.LPSIy(1) + NewData.RASIy(1))/2;
PMRALPz = (NewData.LPSIz(1) + NewData.RASIz(1))/2;

% Distancia punto medio MD1-MI1
% Distancia punto medio MD1-MI1
MDMIx = (NewData.Placa_MI1x(1) + NewData.Placa_MD1x(1))/2;
MDMIy = (NewData.Placa_MI1y(1) + NewData.Placa_MD1y(1))/2;
MDMIz = (NewData.Placa_MI1z(1) + NewData.Placa_MD1z(1))/2;
figure
hold on
plot3(MDMIx,MDMIy,MDMIz,"Marker","o","MarkerEdgeColor","r","MarkerFaceColor","g")
plot3(NewData.Placa_MD1x(1),NewData.Placa_MD1y(1), NewData.Placa_MD1z(1),'.k')
plot3(NewData.Placa_MI1x(1),NewData.Placa_MI1y(1), NewData.Placa_MI1z(1),'.k')
plot3(NewData.LASIx(1), NewData.LASIy(1),NewData.LASIz(1),'.r')
plot3(NewData.RASIx(1), NewData.RASIy(1),NewData.RASIz(1),'.r')
plot3(NewData.LPSIx(1), NewData.LPSIy(1),NewData.LPSIz(1),'.b')
plot3(NewData.RPSIx(1), NewData.RPSIy(1),NewData.RPSIz(1),'.b')
%plot3(NewData.Placa_MD2x(1), NewData.Placa_MD2y(1), NewData.Placa_MD2z(1),'.k')
%plot3(NewData.Placa_MD3x(1), NewData.Placa_MD3y(1), NewData.Placa_MD3z(1),'.k')


%% Promedio de los puntos
PMx = (PMLARPx + PMRALPx + MDMIx)/3;
PMy = (PMLARPy + PMRALPy + MDMIy)/3;
PMz = (PMLARPz + PMRALPz + MDMIz)/3;
plot3(PMx,PMy,PMz,"Marker","*","MarkerEdgeColor","b")

%Marcadores frontales
%Distancia punto medio LARA-RASI
LARAx = (NewData.LASIx(1) + NewData.RASIx(1))/2;
LARAy = (NewData.LASIy(1) + NewData.RASIy(1))/2;
LARAz = (NewData.LASIz(1) + NewData.RASIz(1))/2;
plot3(LARAx,LARAy,LARAz,"Marker","o","MarkerEdgeColor","r","MarkerFaceColor","b")

%Marcadores posteriores
%Distancia punto medio LPSI-RPSI
LPRPx = (NewData.LPSIx(1) + NewData.RPSIx(1))/2;
LPRPy = (NewData.LPSIy(1) + NewData.RPSIy(1))/2;
LPRPz = (NewData.LPSIz(1) + NewData.RPSIz(1))/2;

plot3(LPRPx,LPRPy,LPRPz,"Marker","o","MarkerEdgeColor","r","MarkerFaceColor","r")
legend("M. Cad","P. Medio","M. Post","M. Ant")
xlabel("x label")
ylabel("y label")
zlabel("z label")
hold off

%Creación del plano
%plano en medio 
% Cooz(1,1)=(LARAx+LPRPx+MDMIx+PMx)/4;
% Cooz(1,2)=(LARAy+LPRPy+MDMIy+PMy)/4;
% Cooz(1,3)=(LARAz+LPRPz+MDMIz+PMz)/4;
Cooz(1,1)=PMx; %Coordenadas del plano en x
Cooz(1,2)=PMy; %Coordenadas del plano en y
Cooz(1,3)=PMz;

% M = ones(1000);
% A2= M.*Cooz;
% figure(2)
% hold on
% meshz(A2)
% plot3(PMx,PMy,PMz,"Marker","*","MarkerEdgeColor","b")
% hold off

%Vector en el plano
v1=[LARAx-LPRPx,LARAy-LPRPy,LARAz-LPRPz];
v2=[MDMIx-LPRPx,MDMIy-LPRPy,MDMIz-LPRPz];
    
vc=cross(v1,v2); %vector unitario, sino le estaríamos agregando más o menos !!
v=sqrt((vc(1,1)^2)+(vc(1,2)^2)+(vc(1,3)^2)); %calcular magnitud
vu=vc/v; %vector unitario

for i = 1:height(NewData.Placa_MD2x)
    MD2(i,1)= NewData.Placa_MD2x(i);
    MD2(i,2)= NewData.Placa_MD2y(i);
    MD2(i,3)= NewData.Placa_MD2z(i);
    
    MD3(i,1)= NewData.Placa_MD3x(i);
    MD3(i,2)= NewData.Placa_MD3y(i);
    MD3(i,3)= NewData.Placa_MD3z(i);
    
    MD4(i,1)= NewData.Placa_MD4x(i);
    MD4(i,2)= NewData.Placa_MD4y(i);
    MD4(i,3)= NewData.Placa_MD4z(i);

    MI2(i,1)= NewData.Placa_MI2x(i);
    MI2(i,2)= NewData.Placa_MI2y(i);
    MI2(i,3)= NewData.Placa_MI2z(i);
   
    MI3(i,1)= NewData.Placa_MI3x(i);
    MI3(i,2)= NewData.Placa_MI3y(i);
    MI3(i,3)= NewData.Placa_MI3z(i);
    
    MI4(i,1)= NewData.Placa_MI4x(i);
    MI4(i,2)= NewData.Placa_MI4y(i);
    MI4(i,3)= NewData.Placa_MI4z(i);
end

%% suma de coordenada piernas con resultado punto cruz
%Pierna derecha
dis_MD2 = Cooz + MD2;
dis_MD3 = Cooz + MD3;
dis_MD4 = Cooz + MD4;

%suma de coordenada piernas con resultado punto cruz
PMD2= MD2 + (vu.*dis_MD2);
PMD3= MD3 + (vu.*dis_MD3);
PMD4= MD4 + (vu.*dis_MD4);

%Pierna izquierda
dis_MI2 = Cooz + MI2;
dis_MI3 = Cooz + MI3;
dis_MI4 = Cooz + MI4;

%suma de coordenada piernas con resultado punto cruz
PMI2 = MI2 + (vu.*dis_MI2);
PMI3 = MI3 + (vu.*dis_MI3);
PMI4 = MI4 + (vu.*dis_MI4);

%% proyecciones
figure (3)
hold on
plot3(PMD2(:,1),PMD2(:,2),PMD2(:,3),'k')
plot3(MD2(:,1),MD2(:,2),MD2(:,3),'r')
plot3(PMD3(:,1),PMD3(:,2),PMD3(:,3),'b')
plot3(PMD4(:,1),PMD4(:,2),PMD4(:,3),'g')

plot3(PMI2(:,1),PMI2(:,2),PMI2(:,3),'k')
plot3(MI2(:,1),MI2(:,2),MI2(:,3),'r')
plot3(PMI3(:,1),PMI3(:,2),PMI3(:,3),'b')
plot3(PMI4(:,1),PMI4(:,2),PMI4(:,3),'g')

ylabel('y')
xlabel('x')
zlabel('z')
%meshz(A2)
hold off

%% Calculo de angulo
%Pierna izquierda

MI2M1 = PMI2-PMI3;
MI2M3 = PMI4-PMI3;
% angulos
for i= 1:height(MI2M1)
    pdotIzq(i,1) = dot(MI2M1(i),MI2M3(i));
    apdotIzq(i,1) = (sqrt(MI2M1(i,1)^2+MI2M1(i,2)^2)*sqrt(MI2M3(i,1)^2+MI2M3(i,2)^2));
    angleIzq(i,1) = acos(pdotIzq(i)/apdotIzq(i));
    deg_angleIzq(i,1) = rad2deg(angleIzq(i));
end
%Pierna Derecha

MD2M1 = PMD2-PMD3;
MD2M3 = PMD4-PMD3;
MD2M4 = PMD4-PMD2;
for i= 1:height(MD2M1)
    pdotDer(i,1) = dot(MD2M1(i),MD2M3(i));
    apdotDer(i,1) = (sqrt(MD2M1(i,1)^2+MD2M1(i,2)^2)*sqrt(MD2M3(i,1)^2+MD2M3(i,2)^2));
    angleDer(i,1) = acos(pdotDer(i)/apdotDer(i));
    deg_angleDer(i,1) = rad2deg(angleDer(i));
end
%% Pitagoras
%Cálculo de angulo 
%Pierna Derecha
    AB=(PMD4 - PMD3);
    AC=(PMD2 - PMD3);
    AB_AC= dot(AB,AC,2);
    AB_pitagoras=sqrt(AB(:,1).^2 + AB(:,2).^2 + AB(:,3).^2);
    AC_pitagoras=sqrt(AC(:,1).^2 + AC(:,2).^2 + AC(:,3).^2);
    ang_vectores= acosd((AB_AC)./(AC_pitagoras.*AB_pitagoras));
    
%Pierna Izquierda
    ABIzq=(PMI4 - PMI3);
    ACIzq=(PMI2 - PMI3);
    AB_ACIzq= dot(ABIzq,ACIzq,2);
    AB_pitagorasIzq=sqrt(ABIzq(:,1).^2 + ABIzq(:,2).^2 + ABIzq(:,3).^2);
    AC_pitagorasIzq=sqrt(ACIzq(:,1).^2 + ACIzq(:,2).^2 + ACIzq(:,3).^2);
    ang_vectoresIzq= acosd((AB_ACIzq)./(AC_pitagorasIzq.*AB_pitagorasIzq));