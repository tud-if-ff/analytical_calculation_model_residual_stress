%% Analytisches Berechnungsmodell zur Abschätzung von Eigenspannungen nach dem Festwalzen
% -----------------------------------------------------
% Autoren:      Thomas Werner, Fabio Simon
% Lehrstuhl:    Professur für Formgebende Fertigungsverfahren, 
%               Institut für Fertigunstechnik,
%               Technische Universität Dresden
% Datum:        12.08.2024
% -----------------------------------------------------~

% Eingabe:
% - Walzkraft (P_0)
% - Geometrie der Kontaktpartner (D_1, D_2)
% - Werkstoffkennwerte der Kontaktpartner (E_1, E_2, ny_1, ny_2, sigma_s, epsilon_s, sigma_b, epsilon_b)
% - Laufkoordinate z
% - Wahl des Referenzfalls: Fall 1 - zirkulär ausgebildete Kontaktfläche, Fall 2 - elliptisch ausgebildete Kontaktfläche

% Ergebnis:
% - Resultierende Druckeigenspannungen in Abhängigkeit von der Laufkoordinate z

% Für die Erstellung verwendete Referenzen:
% [1] H. Hertz: Über die Berührung fester elastischer Körper. In: Journal für die reine und angewandte Mathematik 92:156-171, 1881
% [2] A. P. Boresi, R. J. Schmidt, O. M. Sidebottom: Advanced mechanics of materials. John Wiley & Sons, New York, 1992
% [3] K. Han, D. Zhang, C. Yao, L. Tan, Z. Zhou, Y. Zhao: Analytical modeling of through depth strain induced by deep rolling. In: The Journal of Strain Analysis for Engineering Design 57(4):279-290, 2021
% [4] J. K. Li, Y. Mei, W. Duo: Mechanical approach to the residual stress field induced by shot peening. In: Materials Science and Engineering A147:167-173, 1991
% [5] J. Zheng, Y. Shang, Y. Guo, H. Deng, L. Ji: Analytical model of residual stress in ultrasonic rolling of 7075 aluminum alloy. In: Journal of Manufacturing Processes 80:123-140, 2022
% [6] H. Y. Miao, S. Larose, C. Perron, M. Lévesque: An analytical approach to relate shot peening parameters to Almen intensity. Surface & Coatings Technology 205: 2055-2066, 2010
% [7] M. Zhang, Z. Liu, J. Deng, M. Yang, Q. Dai, T. Zhang: Optimum design of compressive residual stress field caused by ultrasonic surface rolling with a mathematical model. In: Applied Mathematical Modeling 76:800-831, 2019

clear *; close all; clc

% Mithilfe des Kerbfaktors lässt sich der Werkzeugradius modifizieren. Dies
% hat Auswirkungen auf die sich ausbildende Kontaktfläche. Mithilfe der
% Anpassung der Variable notch_factor ändert sich die Variable contact_area_elastic. 
% Der Werkzeugradius ist so zu modifizieren, dass die Kontaktfläche den Bedingungen in der Kerbe
% entspricht, um eine Abschätzung der Eigenspannung bzw. Spannungsniveaus an dieser Stelle zu ermöglichen.
% notch_factor = 1 entspricht dem ungekerbten Fall.


%% Eingabeparameter

% Prozessparameter
P_0 = 90 % [N] Eingestellte Walzkraft laut Ecosense-Kraftanzeige bei Zustellung auf glattem Wellenbereich

% Geometrieparameter Kerbe (ROI): Winkel (Ausewerteort), Kerbradius
winkel = 20 % [°]
kerbe_r = 1*10^(-3) %[m]

% Werkzeugparameter
D_R = 40*10^(-3); % [m] Durchmesser Walzrolle
alpha_FWR = 45; % [°] Anstellwinkel Walzrolle zu Wellenachse
R_F = 4500*(1/10^(-3)); % [N/m] Federkonstante des verwendeten Walzwerkzeugs

% E-Modul
E_1 =  630 * 10^9; %[N/m^2]; Werkzeug
E_2 =  210 * 10^9; %[N/m^2]; Welle

% Querkontraktionszahl
ny_1 = 0.22;
ny_2 = 0.3;

% Durchmesser
D_1 = 0.6 * 10^(-3); %[m]
D_2 = 14 * 10^(-3); %[m]

% Radien
R_1 = D_1/2;
R_2 = D_2/2;
R_1s = R_1;
R_2s = inf;

% Winkel zwischen Hauptkrümmungsebenen
alpha_w = 0;

% Materialkennwerte
sigma_s   = 1370 * 10^6; %[Pa] Streckgrenze
epsilon_s = sigma_s/E_2; %[-] Korrespondierende Dehnung zur Streckgrenze
sigma_b = 1556 *10^6; %[Pa] Zugfestigkeit
epsilon_b = 0.03; %[-] Korrespondierende Dehnung zur Zugfestigkeit, Gleichmaßdehnung

% Laufkoordinate 
inc = 0.01 *10^(-3);
z = 0:inc: 0.8 * 10^(-3);

% Fallunterscheidung zwischen zirkulär ausgebilder Kontaktfläche und
% elliptisch ausgebildeter Kontaktfläche
Fall = 2;

% Toleranz zur Abschätzung/Iteration der Fläche
tolerance_inc=0.0001;

% Berechne Festwalzkraft und Faktor in der Kerbe zur Modifikation der Kontaktfläche
[F_W_normal, contact_area_elastic_boresi, contact_area_elastic_circular, contact_area_elastic_elliptic, factor_notch] = calculate_factor_notch(D_1,D_2,E_1,ny_1,E_2,ny_2,P_0,Fall,z,winkel,kerbe_r,D_R,alpha_FWR,R_F,tolerance_inc);

% Effektive Walzkraft in der Kerbe [N]
P_0 = F_W_normal;

% Einfluss des Kerbfaktors (Anpassung der Fläche)
R_1 = R_1*factor_notch;
R_1_s = R_1;


if Fall == 1 
    
    %% Kontaktbelastung nach Hertz
    a_e  = ((3/4)*(((1-ny_1^2)/E_1)+((1-ny_2^2)/E_2))*(R_1*P_0))^(1/3);
    q_0 = (3/2)*(P_0/(pi*a_e^2));
    
    %% %% Elastische Kontaktanalyse Kugel/Platte
    sigma_e_x = -q_0.*((1+ny_2)-(1/2).*(1./(1+(z./a_e).^2))-(1+ny_2).*(z./a_e).*atan(a_e./z));
    sigma_e_y = sigma_e_x;
    sigma_e_z = -q_0./(1+(z./a_e).^2);
    
    figure;
    plot(z/10^-3,sigma_e_x/10^6,'red',z/10^-3,sigma_e_y/10^6,'--blue',z/10^-3,sigma_e_z/10^6,'g','LineWidth', 1.5);
    legend('\sigma_x^e','\sigma_y^e','\sigma_z^e');
    xlabel('Tiefe z in mm');
    ylabel('Spannung \sigma^e in MPa');
    set(gca,'FontSize',12);
    grid on;
    saveas(gcf, '1.pdf');

    sigma_e_eq = sqrt((1/2)*((sigma_e_x-sigma_e_y).^2+(sigma_e_y-sigma_e_z).^2+(sigma_e_z-sigma_e_x).^2));
    epsilon_e_eq = sigma_e_eq/E_2;
    
    %% Elastoplastische Kontaktanalyse Kugel/Platte
    a_p = sqrt((2*P_0)/(3*pi*sigma_s));
    eta = a_p/a_e;
    contact_area_elastic = pi*a_e*a_e;
    
    epsilon_p_eq = zeros(size(epsilon_e_eq));
    for i = 1:length(epsilon_e_eq)
        if epsilon_e_eq(i) <= epsilon_s
            epsilon_p_eq(i) = epsilon_e_eq(i);
        else epsilon_e_eq(i) > epsilon_s;
            epsilon_p_eq(i) = epsilon_s+(eta*(epsilon_e_eq(i)-epsilon_s));
        end
    end

    figure
    plot(z/10^-3,epsilon_p_eq,'b',z/10^-3,epsilon_e_eq,'r','LineWidth', 1.5)
    xlabel('Tiefe z in mm')
    ylabel(['Dehnung ' char(949) '_v'])
    legend([char(949) '_v^p'], [char(949) '_v^e'])
    set(gca,'FontSize',12)
    grid on
    saveas(gcf, '2.pdf')

    H = (sigma_b-sigma_s)/(epsilon_b-epsilon_s);
    
    sigma_p_eq = zeros(size(sigma_e_eq));
    for i = 1:length(sigma_e_eq)
        if epsilon_p_eq(i) < epsilon_s
            sigma_p_eq(i) = sigma_e_eq(i);
        elseif (epsilon_p_eq(i) >= epsilon_s) && (epsilon_p_eq(i) < epsilon_b)
            sigma_p_eq(i) = sigma_s+H*(epsilon_p_eq(i)-epsilon_s);
        else 
            sigma_p_eq(i) = sigma_b;
        end
    end
    
    figure
    plot(z/10^-3,sigma_p_eq/10^6,'b',z/10^-3,sigma_e_eq/10^6,'r','LineWidth', 1.5)
    xlabel('Tiefe z in mm')
    ylabel('Spannung \sigma_v in MPa')
    legend('\sigma_v^p','\sigma_v^e')
    set(gca,'FontSize',12)
    grid on
    saveas(gcf, '3.pdf')
    
    %% Eigenspannungsberechnung    
    delta_sigma_e_eq = sigma_e_eq-(2*sigma_p_eq);
    delta_epsilon_e_eq = delta_sigma_e_eq/E_2;
    delta_epsilon_p_eq =  eta*delta_epsilon_e_eq;
    delta_sigma_p_eq = H*delta_epsilon_p_eq;
    
    sigma_r_x = zeros(size(sigma_e_eq));
    for i = 1:length(epsilon_e_eq)
        if (sigma_e_eq(i) >= sigma_s) && (sigma_e_eq(i) > 2*sigma_p_eq(i))
            sigma_r_x(i) = (1/3)*(sigma_p_eq(i)-2*sigma_p_eq(i)-delta_sigma_p_eq(i));
        elseif (sigma_e_eq(i) >= sigma_s) && (sigma_e_eq(i) <= 2*sigma_p_eq(i))
            sigma_r_x(i) = (1/3)*(sigma_p_eq(i)-sigma_e_eq(i));
        else
            sigma_r_x(i)= 0;
        end
    end
    
    sigma_r_y = zeros(size(sigma_e_eq));
    for i = 1:length(epsilon_e_eq)
        if (sigma_e_eq(i) >= sigma_s) && (sigma_e_eq(i) > 2*sigma_p_eq(i))
            sigma_r_y(i) = (1/3)*(sigma_p_eq(i)-2*sigma_p_eq(i)-delta_sigma_p_eq(i));
        elseif (sigma_e_eq(i) >= sigma_s) && (sigma_e_eq(i) <= 2*sigma_p_eq(i))
            sigma_r_y(i) = (1/3)*(sigma_p_eq(i)-sigma_e_eq(i));
        else
            sigma_r_y(i)= 0;
        end
    end
    
    sigma_r_z = zeros(size(sigma_e_eq));
    for i = 1:length(epsilon_e_eq)
        if (sigma_e_eq(i) >= sigma_s) && (sigma_e_eq(i) > 2*sigma_p_eq(i))
            sigma_r_z(i) = (-2/3)*(sigma_p_eq(i)-2*sigma_p_eq(i)-delta_sigma_p_eq(i));
        elseif (sigma_e_eq(i) >= sigma_s) && (sigma_e_eq(i) <= 2*sigma_p_eq(i))
            sigma_r_z(i) = (-2/3)*(sigma_p_eq(i)-sigma_e_eq(i));
        else
            sigma_r_z(i)= 0;
        end
    end
    

    figure
    plot(z/10^-3,sigma_e_eq/10^6,'r',z/10^-3,2*sigma_p_eq/10^6,'--b',z/10^-3,sigma_p_eq/10^6,'b','LineWidth', 1.5)
    label = {'\sigma_b','\sigma_s'};
    %yline([sigma_b/10^6 sigma_s/10^6],':',label,'LineWidth', 1.5)
    legend('\sigma_v^e','2\sigma_v^p','\sigma_v^p')
    xlabel('Tiefe z in mm')
    ylabel('Spannung \sigma_v in MPa')
    set(gca,'FontSize',12)
    grid on
    saveas(gcf, '4.pdf')

    sigma_R_x = ((1+ny_2)/(1-ny_2))*sigma_r_x;
    sigma_R_y = ((1+ny_2)/(1-ny_2))*sigma_r_x;
    sigma_R_z = 0;
    

    figure
    plot(z/10^-3, sigma_R_x/10^6, 'r', z/10^-3, sigma_R_y/10^6,'--b','LineWidth', 1.5)
    grid on
    legend('\sigma_x^R','\sigma_y^R')
    xlabel('Tiefe z in mm')
    ylabel('Eigenspannung \sigma^R  in MPa')
    set(gca,'FontSize',12)
    ylim([-1500,100])
    saveas(gcf, '5.pdf')

else
    %% Kontaktbelastung nach Hertz
    theta_r = rad2deg(acos(R_1/(R_1+(2*R_2))));
    theta = [0 10 20 30 35 40 45 50 55 60 65 70 75 80 85 90];
    alpha = [inf 6.612 3.778 2.731 2.397 2.136 1.926 1.754 1.611 1.486 1.378 1.284 1.202 1.128 1.061 1.00];
    beta  = [0 0.319 0.408 0.493 0.530 0.576 0.604 0.641 0.678 0.717 0.759 0.802 0.846 0.893 0.944 1.00];
    gamma = [nan 0.851 1.220 1.453 1.550 1.637 1.709 1.772 1.828 1.875 1.912 1.944 1.967 1.985 1.996 2.00];

    alpha_when_theta = interp1(theta,alpha,theta_r); % Berechnet alpha aus der Interpolation der Hertzschen Koeffizienten und dem Wert theta
    beta_when_theta = interp1(theta,beta,theta_r);  % Berechnet beta aus der Interpolation der Hertzschen Koeffizienten und dem Wert theta
    gamma_when_theta = interp1(theta,gamma,theta_r); % Berechnet gamma aus der Interpolation der Hertzschen Koeffizienten und dem Wert theta

    % Halbachsen und Tiefe für den elastischen Kontaktfall
    a_e = alpha_when_theta*((3/4)*((2*R_1*R_2)/(R_1+2*R_2))*(((1-ny_1^2)/E_1)+((1-ny_2^2)/E_2))*P_0)^(1/3);
    b_e = beta_when_theta*((3/4)*((2*R_1*R_2)/(R_1+2*R_2))*(((1-ny_1^2)/E_1)+((1-ny_2^2)/E_2))*P_0)^(1/3);
    c_e = gamma_when_theta*((9/128)*((R_1 + (2*R_2))/(2*R_1*R_2))*(((1-(ny_1^2))/E_1)+((1-(ny_2^2))/E_2))^2*(P_0^2))^(1/3);

    q_0 = (3*P_0)/(2*pi*a_e*b_e);

    %% Elastische Kontaktanalyse
    % Spannungen nach Boresi
    a1 = (1/4)*((1/R_1)+(1/R_2)+(1/R_1s)+(1/R_2s));
    a2 = (1/4)*sqrt((((1/R_1)-(1/R_1s))+((1/R_2)-(1/R_2s)))^2-4*((1/R_1)-(1/R_1s))*((1/R_2)-(1/R_2s)).*((sin(alpha_w)).^2));
    A = a1-a2;
    b1 = (1/4)*((1/R_1)+(1/R_2)+(1/R_1s)+(1/R_2s));
    b2 = (1/4)*sqrt((((1/R_1)-(1/R_1s))+((1/R_2)-(1/R_2s)))^2-4*((1/R_1)-(1/R_1s))*((1/R_2)-(1/R_2s)).*((sin(alpha_w)).^2));
    B = b1+b2;
    e7 = A+B;
    delta = (1/e7)*(((1-ny_1^2)/E_1)+((1-ny_2^2)/E_2));
    k = b_e/a_e;
    ks = sqrt(1-(k^2));
    n1 = k^2+k^2.*((z./b_e).^2);
    n2 = 1+k^2.*((z./b_e).^2);
    n_n = sqrt(n1./n2);
    M = (2*k)/(ks^2*ellipticE(pi/2,ks^2));
    Phi = acot(k*(z./b_e));
    omega_x  = -(1-n_n)./2+(k.*z./b_e).*(ellipticF(Phi,ks^2)-ellipticE(Phi,ks^2));
    omega_xs = -n_n./k.^2+1+(k.*z./b_e).*((1/k^2)*ellipticE(Phi,ks^2)-ellipticF(Phi,ks^2));
    omega_y = 1./(2.*n_n)+(1/2)-n_n./k^2+k.*z./b_e.*(1/k^2*ellipticE(Phi,ks^2)-ellipticF(Phi,ks^2));
    omega_ys = -1+n_n+k.*(z./b_e).*(ellipticF(Phi,ks^2)-ellipticE(Phi,ks^2));

    sigma_e_x = (M.*(omega_x+ny_2.*omega_xs))*(b_e/delta);
    sigma_e_y = (M.*(omega_y+ny_2.*omega_ys))*(b_e/delta);
    sigma_e_z = -(b_e./delta).*((1/2).*M*((1./n_n)-n_n));

    % Graphen 
    figure
    plot(z/10^-3,sigma_e_x/10^6,'r',z/10^-3,sigma_e_y/10^6,'b',z/10^-3,sigma_e_z/10^6,'g','LineWidth', 1.5)
    legend('\sigma_x^e','\sigma_y^e','\sigma_z^e')
    grid on;
    xlabel('Tiefe z in mm')
    ylabel('Spannung \sigma^e in MPa')
    set(gca,'FontSize',12)
    saveas(gcf, '1.pdf')
    
    % Vergleichspannung und -dehnung  
    sigma_e_eq = sqrt((1/2)*((sigma_e_x-sigma_e_y).^2+(sigma_e_y-sigma_e_z).^2+(sigma_e_z-sigma_e_x).^2));
    epsilon_e_eq = sigma_e_eq/E_2;
    
    %% Elastoplastische Kontaktanalyse Kugel/Zylinder
    contact_area_elastic = pi*a_e*b_e; %[m^2]
    delta_p = (P_0/(3*pi*R_1*sigma_s))*sqrt((R_1+R_2)/R_2);
    a_p = sqrt(2*R_1*delta_p);
    
    eta = a_p/a_e;
    
    epsilon_p_eq = zeros(size(epsilon_e_eq));
    for i = 1:length(epsilon_e_eq)
        if epsilon_e_eq(i) <= epsilon_s
            epsilon_p_eq(i) = epsilon_e_eq(i);
        else epsilon_e_eq(i) > epsilon_s;
            epsilon_p_eq(i) = epsilon_s+(eta*(epsilon_e_eq(i)-epsilon_s));
        end
    end
    
    figure
    plot(z/10^-3,epsilon_p_eq,'b',z/10^-3,epsilon_e_eq,'r','LineWidth', 1.5)
    grid on
    ylabel(['Dehnung ' char(949) '_v'])
    legend([char(949) '_v^p'], [char(949) '_v^e'])
    xlabel('Tiefe z in mm')
    set(gca,'FontSize',12)
    saveas(gcf, '2.pdf')
    
    H = (sigma_b-sigma_s)/(epsilon_b-epsilon_s);
    
    sigma_p_eq = zeros(size(sigma_e_eq));
    for i = 1:length(sigma_e_eq)
        if epsilon_e_eq(i) < epsilon_s
            sigma_p_eq(i) = sigma_e_eq(i);
        elseif (epsilon_p_eq(i) >= epsilon_s) && (epsilon_p_eq(i) < epsilon_b)
            sigma_p_eq(i) = sigma_s+H*(epsilon_p_eq(i)-epsilon_s);
        else 
            sigma_p_eq(i) = sigma_b;
        end
    end

    figure
    plot(z/10^-3,sigma_p_eq/10^6,'b',z/10^-3,sigma_e_eq/10^6,'r','LineWidth', 1.5)
    xlabel('Tiefe z in mm')
    ylabel('Spannung \sigma_v in MPa')
    set(gca,'FontSize',12)
    grid on
    legend('\sigma_v^p','\sigma_v^e')
    saveas(gcf, '3.pdf')
    
    %% Eigenspannungsberechnung
    delta_sigma_e_eq = sigma_e_eq-(2*sigma_p_eq);
    delta_epsilon_e_eq = delta_sigma_e_eq/E_2;
    delta_epsilon_p_eq = eta.*delta_epsilon_e_eq;
    delta_sigma_p_eq = H*delta_epsilon_p_eq;
    
    n = a_e/b_e;
    %n = sigma_e_x/sigma_e_y
    
    h = n/(sqrt(3*(n^2+n+1)));
    o = 1/(sqrt(3*(n^2+n+1)));
    j = -(n+1)/(sqrt(3*(n^2+n+1)));
    
    sigma_r_x = zeros(size(sigma_e_eq));
    for i = 1:length(epsilon_e_eq)
        if (sigma_e_eq(i) >= sigma_s) && (sigma_e_eq(i) > 2*sigma_p_eq(i))
            sigma_r_x(i) = h*(sigma_p_eq(i)-2*sigma_p_eq(i)-delta_sigma_p_eq(i));
        elseif (sigma_e_eq(i) >= sigma_s) && (sigma_e_eq(i) <= 2*sigma_p_eq(i))
            sigma_r_x(i) = h*(sigma_p_eq(i)-sigma_e_eq(i));
        else
            sigma_r_x(i)= 0;
        end
    end
    
    sigma_r_y = zeros(size(sigma_e_eq));
    for i = 1:length(epsilon_e_eq)
        if (sigma_e_eq(i) >= sigma_s) && (sigma_e_eq(i) > 2*sigma_p_eq(i))
            sigma_r_y(i) = o*(sigma_p_eq(i)-2*sigma_p_eq(i)-delta_sigma_p_eq(i));
        elseif (sigma_e_eq(i) >= sigma_s) && (sigma_e_eq(i) <= 2*sigma_p_eq(i))
            sigma_r_y(i) = o*(sigma_p_eq(i)-sigma_e_eq(i));
        else
            sigma_r_y(i)= 0;
        end
    end
    
    sigma_r_z = zeros(size(sigma_e_eq));
    for i = 1:length(epsilon_e_eq)
        if (sigma_e_eq(i) >= sigma_s) && (sigma_e_eq(i) > 2*sigma_p_eq(i))
            sigma_r_z(i) = j*(sigma_p_eq(i)-2*sigma_p_eq(i)-delta_sigma_p_eq(i));
        elseif (sigma_e_eq(i) >= sigma_s) && (sigma_e_eq(i) <= 2*sigma_p_eq(i))
            sigma_r_z(i) = j*(sigma_p_eq(i)-sigma_e_eq(i));
        else
            sigma_r_z(i)= 0;
        end
    end
    
    figure
    plot(z/10^-3,sigma_e_eq/10^6,'r',z/10^-3,2*sigma_p_eq/10^6,'--b',z/10^-3,sigma_p_eq/10^6,'b','LineWidth', 1.5)
    label = {'\sigma_b','\sigma_s'};
    %yline([sigma_b/10^6 sigma_s/10^6],':',label,'LineWidth', 1.5)
    legend('\sigma_v^e','2\sigma_v^p','\sigma_v^p')
    xlabel('Tiefe z in mm')
    ylabel('Spannung \sigma_v in MPa')
    set(gca,'FontSize',12)
    grid on
    saveas(gcf, '4.pdf')

    sigma_rel_x = ny_2/(1-ny_2)*sigma_r_z;
    sigma_rel_y = sigma_rel_x;
    sigma_R_x = sigma_r_x-sigma_rel_x;
    sigma_R_y = sigma_r_y-sigma_rel_y;
    sigma_R_z = 0;
    
    figure
    plot(z/10^-3, sigma_R_x/10^6, 'r',z/10^-3, sigma_R_y/10^6, 'b','LineWidth', 1.5)
    grid on
    legend('\sigma_x^R','\sigma_y^R')
    xlabel('Tiefe z in mm')
    ylabel('Eigenspannung \sigma^R  in MPa')
    set(gca,'FontSize',12)
    ylim([-2000,100])
    saveas(gcf, '5.pdf')
    

end


% Dateiname festlegen
filename = 'Festwalzparameter.txt';

% Datei zum Schreiben öffnen
fileID = fopen(filename, 'w');

% Schreibe die Eingabeparameter in die Datei
fprintf(fileID, 'Inputdaten:\n');
fprintf(fileID, '\n');
fprintf(fileID, 'Winkel (Auswerteort): %d °\n', winkel);
fprintf(fileID, 'Kerbradius: %d m\n', kerbe_r);
fprintf(fileID, 'Durchmesser Walzrolle: %d m\n', D_R);
fprintf(fileID, 'Anstellwinkel Walzrolle: %d °\n', alpha_FWR);
fprintf(fileID, 'Federkonstante: %d N/m\n', R_F);
fprintf(fileID, 'Walzkraft: %d N\n', P_0);
fprintf(fileID, 'Toleranz Flächenabweichung: %d\n', tolerance_inc);
fprintf(fileID, 'Elastizitätsmodul 1: %d Pa\n', E_1);
fprintf(fileID, 'Elastizitätsmodul 2: %d Pa\n', E_2);
fprintf(fileID, 'Querkontraktionszahl 1: %d \n', ny_1);
fprintf(fileID, 'Querkontraktionszahl 2: %d \n', ny_2);
fprintf(fileID, 'Durchmesser 1: %d m\n', D_1);
fprintf(fileID, 'Durchmesser 2: %d m\n', D_2);
fprintf(fileID, 'Radius 1: %d m\n', R_1);
fprintf(fileID, 'Radius 2: %d m\n', R_2);
fprintf(fileID, 'Streckgrenze: %d Pa\n', sigma_s);
fprintf(fileID, 'Korrespondierende Dehnung zur Streckgrenze: %d \n', epsilon_s);
fprintf(fileID, 'Zugfestigkeit: %d Pa\n', sigma_b);
fprintf(fileID, 'Korrespondierende Dehnung zur Zugfestigkeit: %d \n', epsilon_b);
fprintf(fileID, 'Fallunterscheidung zirkulär/elliptisch: %d \n', Fall);
fprintf(fileID, 'Faktor zur Abschätzung in der Kerbe: %d \n', factor_notch);
fprintf(fileID, 'Kontaktfläche Boresi: %d m^2\n', contact_area_elastic_boresi);
fprintf(fileID, 'Kontaktfläche Analytik: %d m^2\n', contact_area_elastic);
fprintf(fileID, 'Kontaktfläche zirkulär: %d m^2\n', contact_area_elastic_circular);
fprintf(fileID, 'Kontaktfläche elliptisch: %d m^2\n', contact_area_elastic_elliptic);

fprintf(fileID, '\n');
fprintf(fileID, '-----------------\n');
fprintf(fileID, 'Outputdaten:\n');
fprintf(fileID, '\n');
fprintf(fileID, 'Koordinate z: %d m\n', z);
fprintf(fileID, 'Eigenspannung (x): %d Pa\n', sigma_R_x);



% Datei schließen
fclose(fileID);

disp('Festwalzparameter erfolgreich in die Datei geschrieben.');