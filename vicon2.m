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
% hecho por profe CAMBIAR

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

% Creating MD vector for Coordinate System
[md_x, md_y, md_z] = axisCreation(md4,md2,md3);

% Creating PD vector for Coordinate System
[pd_x, pd_y, pd_z] = axisCreation(pd4,pd2,pd3);

% empty variable to allocate the angle in each frame
k_ang= nan(height(Data),3);%

figure (3)
hold on
for f= 1:height(frames)%
    %create coordinate system matrices
    MD= [md_x(f,:).', md_y(f,:).', md_z(f,:).'];
    PD= [pd_x(f,:).', pd_y(f,:).', pd_z(f,:).'];
    %calculate the rotation matrix
    R= MD*PD';
    
    k_ang(f,1)= atan2d(R(3,2),R(3,3));
    k_ang(f,2)= atan2d(R(2,1),R(1,1));
    if k_ang(f,2) == 0
        k_ang(f,3)= atan2d(-R(3,1),R(2,1)/sin(k_ang(f,2)));
    else
        k_ang(f,3)= atan2d(-R(3,1),R(1,1)/cos(k_ang(f,2)));
    end
end

for f= 1:height(CutData)%
    clf
    hold on
    subplot(1,2,1)
    plot(k_ang)
    xline(f)
    axis([0 1200 -180 180])
    subplot(1,2,2)
    hold on
    for m= 1:length(markers)
        labelx= strcat(markers(m),'x');
        labely= strcat(markers(m),'y');
        labelz= strcat(markers(m),'z');
        marker= [CutData.(labelx) CutData.(labely) CutData.(labelz)];%
        plot3(marker(f,1),marker(f,2),marker(f,3),'k.')
    end
    axis([-1200 1200 -800 1600 0 2400])
    set(gca, 'view', [-132.3000 16.8000])
    %pause(0.01)
    hold off
end
hold off

%% functions
% hecho por profe CAMBIAR

function [x,y,z]= axisCreation(po,pz,py)
    %Z axis
    z= pz-po;
    z= z./norm_v(z);
    
    %Y axis
    y= py-po;
    y= y./norm_v(y);
    
    %X axis
    x= cross(y,z);
    y= cross(z,x);
end

function value= norm_v(vector)
    value= sqrt(vector(:,1).^2 + vector(:,2).^2 + vector(:,3).^2);
end