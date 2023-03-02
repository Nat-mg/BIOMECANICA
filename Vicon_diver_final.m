%Código para calculo de angulos de rodilla en movimiento diver
clc;
clear all;
i = 1;

%Data extraction
Data = dataImport('Diver_GonzaloI.csv');

%% Empti matriz variables
NewData.Frame = Data.Frame;

%Cadera alta
NewData.LASIx = Data.LASIx; %anterior izquierda
NewData.LASIy = Data.LASIy;
NewData.LASIz = Data.LASIz;

NewData.RASIx = Data.RASIx; %anterior dererecha
NewData.RASIy = Data.RASIy;
NewData.RASIz = Data.RASIz;

NewData.LPSIx = Data.LPSIx; %posterior izquierda
NewData.LPSIy = Data.LPSIy;
NewData.LPSIz = Data.LPSIz;

NewData.RPSIx = Data.RPSIx; %posterior derecha
NewData.RPSIy = Data.RPSIy;
NewData.RPSIz = Data.RPSIz;

%Pierna Izquierda
NewData.Placa_MI1x = Data.Placa_MI1x; %cadera
NewData.Placa_MI1y = Data.Placa_MI1y;
NewData.Placa_MI1z = Data.Placa_MI1z;

NewData.Placa_MI2x = Data.Placa_MI2x; %femur
NewData.Placa_MI2y = Data.Placa_MI2y;
NewData.Placa_MI2z = Data.Placa_MI2z;

NewData.Placa_MI3x = Data.Placa_MI3x; %rodilla
NewData.Placa_MI3y = Data.Placa_MI3y;
NewData.Placa_MI3z = Data.Placa_MI3z;

NewData.Placa_MI4x = Data.Placa_MI4x; %tibia
NewData.Placa_MI4y = Data.Placa_MI4y;
NewData.Placa_MI4z = Data.Placa_MI4z;

%Pierna Derecha
NewData.Placa_MD1x = Data.Placa_MD1x; %cadera
NewData.Placa_MD1y = Data.Placa_MD1y;
NewData.Placa_MD1z = Data.Placa_MD1z;

NewData.Placa_MD2x = Data.Placa_MD2x; %femur
NewData.Placa_MD2y = Data.Placa_MD2y;
NewData.Placa_MD2z = Data.Placa_MD2z;

NewData.Placa_MD3x = Data.Placa_MD3x; %rodilla
NewData.Placa_MD3y = Data.Placa_MD3y;
NewData.Placa_MD3z = Data.Placa_MD3z;

NewData.Placa_MD4x = Data.Placa_MD4x; %tibia
NewData.Placa_MD4y = Data.Placa_MD4y;
NewData.Placa_MD4z = Data.Placa_MD4z;
%Covertir structura a tabla
NewData = struct2table(NewData);

markers = string(NewData.Properties.VariableNames).';
idx = contains(markers,'x');
markers = erase(markers(idx),'x');
clear idx

%% ismissing 
frames= NewData.Frame;
for m= 1:length(markers)
    labelx = strcat(markers(m),'x');
    idx= ~ismissing(NewData.(labelx));
    if sum(ismissing(NewData.(labelx))) > 0
        original_data = NewData.(labelx);
        new_data= spline(frames(idx),original_data(idx),frames);
        NewData.(labelx) = new_data;
    end
    labely = strcat(markers(m),'y');
    idx = ~ismissing(NewData.(labely));
    if sum(ismissing(NewData.(labely))) > 0 
        original_data = NewData.(labely);
        new_data= spline(frames(idx),original_data(idx),frames);
        NewData.(labely) = new_data;
    end
    labelz = strcat(markers(m),'z');
    idx = ~ismissing(NewData.(labelz));
    if sum(ismissing(NewData.(labelz))) > 0
        original_data = NewData.(labelz);
        new_data= spline(frames(idx),original_data(idx),frames);
        NewData.(labelz) = new_data;
    end
end

clear labelx labely labelz original_data new_data
%% Calculo de punto medio marcadores cadera
%Punto medio cintura
PMLARPx = (NewData.LASIx(i) + NewData.RPSIx(i))/2; % punto medio LASI y RPSI en x
PMLARPy = (NewData.LASIy(i) + NewData.RPSIy(i))/2; % punto medio LASI y RPSI en y
PMLARPz = (NewData.LASIz(i) + NewData.RPSIz(i))/2; % punto medio LASI y RPSI en z
PMRALPx = (NewData.LPSIx(i) + NewData.RASIx(i))/2; % punto medio LPSI y RASI en x
PMRALPy = (NewData.LPSIy(i) + NewData.RASIy(i))/2; % punto medio LPSI y RASI en y
PMRALPz = (NewData.LPSIz(i) + NewData.RASIz(i))/2; % punto medio LPSI y RASI en z

%Punto medio cadera (marcadores completos)
PMMDMIx = (NewData.Placa_MI1x(i) + NewData.Placa_MD1x(i))/2; % punto medio MI1 y MD1 en x
PMMDMIy = (NewData.Placa_MI1y(i) + NewData.Placa_MD1y(i))/2; % punto medio MI1 y MD1 en y
PMMDMIz = (NewData.Placa_MI1z(i) + NewData.Placa_MD1z(i))/2; % punto medio MI1 y MD1 en z
%Punto medio cadera (marcadores creados)
% NewData.Placa_MI1x(i) = (NewData.LASIx(i) - NewData.Placa_MD1x(i))+NewData.RASIx(145);
% PMMDMIx = (NewData.Placa_MI1x(i) + NewData.Placa_MD1x(i))/2; % punto medio MI1 y MD1 en x
% PMMDMIy = (NewData.Placa_MD1y(i) + NewData.Placa_MD1y(i))/2; % punto medio MI1 y MD1 en y
% PMMDMIz = (NewData.Placa_MD1z(i) + NewData.Placa_MD1z(i))/2; % punto medio MI1 y MD1 en z

%Promedio de los puntos
Prom_x = (PMLARPx + PMRALPx + PMMDMIx)/3; %promedio de los puntos medios en x
Prom_y = (PMLARPy + PMRALPy + PMMDMIy)/3; %promedio de los puntos medios en y
Prom_z = (PMLARPz + PMRALPz + PMMDMIz)/3; %promedio de los puntos medios en z

%% Calculo de punto medio marcadores frontales
%Punto medio LARA-RASI
PMLARAx = (NewData.LASIx(i) + NewData.RASIx(i))/2; % punto medio LARA y RASI en x
PMLARAy = (NewData.LASIy(i) + NewData.RASIy(i))/2; % punto medio LARA y RASI en y
PMLARAz = (NewData.LASIz(i) + NewData.RASIz(i))/2; % punto medio LARA y RASI en z

% Calculo de punto medio marcadores posteriores
%Punto medio LPSI-RPSI
PMLPRPx = (NewData.LPSIx(i) + NewData.RPSIx(i))/2; % punto medio LPSI y RPSI en x
PMLPRPy = (NewData.LPSIy(i) + NewData.RPSIy(i))/2; % punto medio LPSI y RPSI en y
PMLPRPz = (NewData.LPSIz(i) + NewData.RPSIz(i))/2; % punto medio LPSI y RPSI en z

%% Creación del plano
% Coo(1,1)=(PMLARAx+PMLPRPx+PMMDMIx+Prom_x)/4; %Coordenadas del plano en x
% Coo(1,2)=(PMLARAy+PMLPRPy+PMMDMIy+Prom_y)/4; %Coordenadas del plano en y
% Coo(1,3)=(PMLARAz+PMLPRPz+PMMDMIz+Prom_z)/4; %Coordenadas del plano en z

Coo(1,1)=Prom_x; %Coordenadas del plano en x
Coo(1,2)=Prom_y; %Coordenadas del plano en y
Coo(1,3)=Prom_z;
%Vector en el plano
v1=[PMLARAx-PMLPRPx,PMLARAy-PMLPRPy,PMLARAz-PMLPRPz]; %Vector del punto medio LARA-RASI a LPSI-RPSI
v2=[PMMDMIx-PMLPRPx,PMMDMIy-PMLPRPy,PMMDMIz-PMLPRPz]; %Vector del punto medio LPSI-RPSI a punto medio cadera
    
vc=cross(v1,v2); %vector producto cruz
v=sqrt((vc(1,1)^2)+(vc(1,2)^2)+(vc(1,3)^2)); %calcular magnitud
vu=vc/v; %vector unitario

%% Vector con coordenadas de pierna derecha
MD2(:,1)= NewData.Placa_MD2x; %femur 
MD2(:,2)= NewData.Placa_MD2y;
MD2(:,3)= NewData.Placa_MD2z;
    
MD3(:,1)= NewData.Placa_MD3x; %rodilla
MD3(:,2)= NewData.Placa_MD3y;
MD3(:,3)= NewData.Placa_MD3z;
    
MD4(:,1)= NewData.Placa_MD4x; %tibia
MD4(:,2)= NewData.Placa_MD4y;
MD4(:,3)= NewData.Placa_MD4z;

% Vector con coordenadas de pierna derecha
MI2(:,1)= NewData.Placa_MI2x; %femur 
MI2(:,2)= NewData.Placa_MI2y;
MI2(:,3)= NewData.Placa_MI2z;
   
MI3(:,1)= NewData.Placa_MI3x; %rodilla
MI3(:,2)= NewData.Placa_MI3y;
MI3(:,3)= NewData.Placa_MI3z;
    
MI4(:,1)= NewData.Placa_MI4x; %tibia
MI4(:,2)= NewData.Placa_MI4y;
MI4(:,3)= NewData.Placa_MI4z;

%% suma de coordenada piernas con resultado punto cruz
%Pierna derecha
dis_MD2 = Coo + MD2; %distancia de femur al plano
dis_MD3 = Coo + MD3; %distancia de rodilla al plano
dis_MD4 = Coo + MD4; %distancia de tibia al plano

%suma de coordenada piernas con resultado punto cruz
PMD2= MD2 + (vu.*dis_MD2); %Coordenadas femur en el plano
PMD3= MD3 + (vu.*dis_MD3); %Coordenadas rodilla en el plano
PMD4= MD4 + (vu.*dis_MD4); %Coordenadas tibia en el plano

%Pierna izquierda
dis_MI2 = Coo + MI2;
dis_MI3 = Coo + MI3;
dis_MI4 = Coo + MI4;

%suma de coordenada piernas con resultado punto cruz
PMI2 = MI2 + (vu.*dis_MI2);
PMI3 = MI3 + (vu.*dis_MI3);
PMI4 = MI4 + (vu.*dis_MI4);


%% Cálculo de angulo (pitgoras)
%Pierna Derecha
AB=(PMD4 - PMD3); %distancia entre coordenadas de tibia y rodilla
AC=(PMD2 - PMD3); %distancia entre coordenadas de femur y rodilla
AB_AC= dot(AB,AC,2); %producto punto de AB y AC
AB_pitagoras=sqrt(AB(:,1).^2 + AB(:,2).^2 + AB(:,3).^2); %calcular pitagoras de AB
AC_pitagoras=sqrt(AC(:,1).^2 + AC(:,2).^2 + AC(:,3).^2); %calcular pitagoras de AC
ang_vectores= acosd((AB_AC)./(AC_pitagoras.*AB_pitagoras)); %obtener angulo en grados
    
%Pierna Izquierda
ABIzq=(PMI4 - PMI3);
ACIzq=(PMI2 - PMI3);
AB_ACIzq= dot(ABIzq,ACIzq,2);
AB_pitagorasIzq=sqrt(ABIzq(:,1).^2 + ABIzq(:,2).^2 + ABIzq(:,3).^2);
AC_pitagorasIzq=sqrt(ACIzq(:,1).^2 + ACIzq(:,2).^2 + ACIzq(:,3).^2);
ang_vectoresIzq= acosd((AB_ACIzq)./(AC_pitagorasIzq.*AB_pitagorasIzq));