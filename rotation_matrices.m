clc;
clear all;
data = dataImport('Ejemplo2D.csv');

% MD coordinate system creation
% point extraction

md1= [data.md1x data.md1y data.md1z]; %load variables
md2= [data.md2x data.md2y data.md2z];
md3= [data.md3x data.md3y data.md3z];
md4= [data.md4x data.md4y data.md4z];

% third axis md(z)
md_z= md2-md4;
norm_v= sqrt(md_z(:,2).^2 + md_z(:,3).^2);
md_z= md_z./norm_v;

% second axis md(y)
md_y= md3-md4;
norm_v= sqrt(md_y(:,2).^2 + md_y(:,3).^2);
md_y= md_y./norm_v;

% first axis md(x)
md_x = cross(md_y,md_z);
md_y = cross(md_z,md_x);
angle = acosd(dot(md_y(1,:),md_z(1,:)));
disp(strcat("Angle between z and y axis: ",string(angle)))

%plot points MD
f=1;
figure (1)
plot([md1(f,1) md2(f,1) md3(f,1) md4(f,1)],[md1(f,2) md2(f,2) md3(f,2) md4(f,2)],'ko')
title('Plot points MD')

% PD coordinate system creation

pd1= [data.pd1x data.pd1y data.pd1z]; %load variables
pd2= [data.pd2x data.pd2y data.pd2z];
pd3= [data.pd3x data.pd3y data.pd3z];
pd4= [data.pd4x data.pd4y data.pd4z];

% third axis md(z)
pd_z= pd2-pd4;
norm_v= sqrt(pd_z(:,2).^2 + pd_z(:,3).^2);
pd_z= pd_z./norm_v;

% second axis md(y)
pd_y= pd3-pd4;
norm_v= sqrt(pd_y(:,2).^2 + pd_y(:,3).^2);
pd_y= pd_y./norm_v;

% first axis md(x)
pd_x = cross(pd_y,pd_z);
pd_y = cross(pd_z,pd_x);
angle = acosd(dot(pd_y(1,:),pd_z(1,:)));
disp(strcat("Angle between z and y axis: ",string(angle)))

%plot points MD
f=1;
figure (2)
plot([pd1(f,1) pd2(f,1) pd3(f,1) pd4(f,1)],[pd1(f,2) pd2(f,2) pd3(f,2) pd4(f,2)],'ko')
title('Plot points PD')

%empty variable to allocate the angle in each frame
k_ang= nan(height(data),1);

for f=1:height(data)
    %create coordinate system matrices
    MD=[md_y(f,2:3).',md_z(f,2:3).']; % ' transpose a vector or matrix
    PD=[pd_y(f,2:3).',pd_z(f,2:3).'];
    
    %calculate the rotation matrix
    R = MD*PD';
    
    %calculate the angle using four element of the rotation matrix
    th(1,1) = acosd(R(1,1));
    th(1,2) = -asind(R(1,2));
    th(2,1) = asind(R(2,1));
    th(2,2) = acosd(R(2,2));
    
    %calculate the correct case for each case
    if (R(1,1)>0) && (R(2,1)>0)
        angle=th(1,1);
   
    elseif (R(1,1)<0) && (R(2,1)>0)
        angle=th(1,1);
    
    elseif (R(1,1)<0) && (R(2,1)<0)
        angle=180 - th(1,2);
        
    elseif (R(1,1)>0) && (R(2,1)<0)
        angle=360 - th(1,1);
        
    elseif (R(2,1)==0)
        angle=th(1,1);
        
    elseif (R(1,1)==0) && (R(2,1)>1)
        angle=th(2,1);
        
    elseif (R(1,1)==0) && (R(2,1)<1)
        angle=360 - th(2,1);
    end
    %save the correct value
    k_ang(f)= angle;
end

figure (3)
hold on

%visualization
cla
plot([MD(1,1) 0 MD(1,2)],[MD(2,1) 0 MD(2,2)],'k')
plot([PD(1,1) 0 PD(1,2)],[PD(2,1) 0 PD(2,2)],'r')
axis([-1 1 -1 1])
title(string(k_ang(f)))
pause(0.1)

%plot the knee angle vs time
figure (4)
t= [data.Time];
plot(t,k_ang,'LineWidth',1.5)
title('Knee angle vs time')
xlabel('Time [s]')
ylabel('Knee angles [°]')