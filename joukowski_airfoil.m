clc 
clear all 
close all
%% Constants
Uinf = 110;
c = 1;
b=c/4;
alpha = [2 5 10 12];
Chmax = [0.04 0.05 0.06 0.07];
thmax = [0.05  0.07 0.1 0.12];

alpha=alpha.*(pi/180);  % Converting alpha into radians

% The four different cases
for Ca=1:4
    
    e=thmax(Ca)/1.3;
    beta=2*Chmax(Ca);
    a=b*(1+e)./cos(beta);
    x0=-b*e;
    y0=a*beta;
    
    theta_up = [0.0001:0.001:pi];
    theta_down = -theta_up;
    
    r_up = b.*(1+e.*(1-cos(theta_up))+beta.*sin(theta_up));
    r_down = b.*(1+e.*(1-cos(theta_down))+beta.*sin(theta_down));
    
    X_up = r_up.*cos(theta_up);
    X_down = r_down.*cos(theta_down);
    X_up_dash = X_up-x0;
    X_down_dash = X_down-x0;
    
    Y_up = r_up.*sin(theta_up);
    Y_down = r_down.*sin(theta_down);
    Y_up_dash = Y_up-y0;
    Y_down_dash = Y_down-y0;
    
    theta_up_dash = atan2(Y_up_dash,X_up_dash);
    theta_down_dash = atan2(Y_down_dash,X_down_dash);
    
    X1_up = X_up+((b^2*X_up)./(X_up.^2+Y_up.^2));
    X1_down = X_down+((b^2*X_down)./(X_down.^2+Y_down.^2));
    
    v_theta_up_dash = (-2*Uinf).*(sin(theta_up_dash-alpha(Ca))+sin(alpha(Ca)+beta));
    v_theta_down_dash = (-2*Uinf).*(sin(theta_down_dash-alpha(Ca))+sin(alpha(Ca)+beta));
    V1_up = sqrt(((v_theta_up_dash).^2)./(1-(2*(b^2./(r_up.^2)).*cos(2*theta_up))+(b^4./(r_up.^4))));
    V1_down = sqrt(((v_theta_down_dash).^2)./(1-(2*(b^2./(r_down.^2)).*cos(2*theta_down))+(b^4./(r_down.^4))));
    cp_up = 1-(V1_up/Uinf).^2;
    cp_down = 1-(V1_down/Uinf).^2;
    
    %% Ploting
    figure
    plot(X1_up,cp_up); % the upper face
    hold on;
    plot(X1_down,cp_down); % the bottom face
    grid on;
    title('Cp on the airfiol');
    ylabel('Cp');
    xlabel('X');
    hold off
    figure
    plot(X1_up,V1_up); % the upper face
    hold on;
    plot(X1_down,V1_down); % the bottom face
    grid on;
    title('V on the airfiol');
    ylabel('V');
    xlabel('X');
   hold off

F_up(1)=0;
dV_up_mid(1)=0;
V1_up_mid(1)=0;
F_up_mid(1)=0;
beta_up(1)=0;

F_down(1)=0;
dV_down_mid(1)=0;
V1_down_mid(1)=0;
F_down_mid(1)=0;
beta_down(1)=0;

for i=1:length(V1_up)-1
    
    ds_up=sqrt(((X1_up(i+1)-X1_up(i))^2)+((Y_up(i+1)-Y_up(i))^2));
    
    F_up(i+1)=0.5*(V1_up(i)+V1_up(i+1))*(ds_up)+F_up(i);

    dV_up_mid(i+1)=(V1_up(i+1)-V1_up(i))/(ds_up);

    V1_up_mid(i+1)=0.5*(V1_up(i)+V1_up(i+1));
   
    F_up_mid(i+1)=0.5*(F_up(i)+F_up(i+1));

    beta_up(i+1)=2.*dV_up_mid(i).*F_up_mid(i)./(V1_up_mid(i)^2);

    ds_down=sqrt(((X1_down(i+1)-X1_down(i))^2)+((Y_down(i+1)-Y_down(i))^2));
  
    F_down(i+1)=0.5*(V1_down(i)+V1_down(i+1))*(ds_down)+F_down(i);

    dV_down_mid(i+1)=(V1_down(i+1)-V1_down(i))/(ds_down);

    V1_down_mid(i+1)=0.5*(V1_down(i)+V1_down(i+1));
   
    F_down_mid(i+1)=0.5*(F_down(i)+F_down(i+1));

    beta_down(i+1)=2.*dV_down_mid(i).*F_down_mid(i)./(V1_down_mid(i)^2);

end

    figure

    plot(X1_up(1:2650),beta_up(1:2650));
    hold on;
    plot(X1_down(1:2650),beta_down(1:2650)); % the bottom face
    grid on;
    title('Beta (The Pressure-Gradient Parameter) on the airfiol');
    ylabel('Beta');
    xlabel('X');
    hold off


end