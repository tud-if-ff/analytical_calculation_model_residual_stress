function [F_W_normal, contact_area_elastic_boresi, contact_area_elastic_circular, contact_area_elastic_elliptic, factor_notch] = calculate_factor_notch(D_1, D_2, E_1, ny_1, E_2, ny_2, P_0, Fall,z,winkel,kerbe_r,D_R,alpha_FWR,R_F,tolerance_inc)
    % Übertragung aus dem Berechnungsprogramm
    winkel = winkel;
    r = kerbe_r;
    R_1 = D_1/2;

    % Winkelberechnung
    beta = 90-winkel;
    
    % Vorliegender Kontaktfall: Ellipsoid (Rolle) - Ellipsoid (Welle)
    % Hauptkrümmungsradius 1 > Hauptkrümmungsradius 2 (HK 1 > HK 2)
    r_11 = R_1 + (D_R/2 - R_1)./cosd(beta-alpha_FWR); % Rolle Hauptkrümmungsradius 1
    r_12 = R_1; % Rolle Hauptkrümmungsradius 2
    r_21 = (D_2/2+r*(1-sind(beta)))./sind(beta); % Welle Hauptkrümmungsradius 1
    r_22 = -r; % Welle Hauptkrümmungsradius 2
    alpha_w = 0;
    
    % Bestimmung der lokal wirksamen Walzkraft (Kraft in Richtung der Kontaktnormalen)
    z_R = P_0 .* sind(alpha_FWR).^2./R_F; % [m] Radiale Zustellung des Walzwerkzeugs auf glatter Welle zum Erreichen der Walzkraft F_W. Wird auch in der Kerbe jeweils in radialer Richtung eingehalten (Werkzeugkorrektur an der Drehmaschine).
    delta_F = z_R ./ cosd(beta-alpha_FWR); % [m] Auslenkung der Feder im Walzwerkzeug in Druckrichtung
    F_F = delta_F .* R_F; % [N] Druckkraft der Feder im Walzwerkzeug
    F_W_normal = F_F ./cosd(beta-alpha_FWR); % [N] Resultierende Walzkraft in Normalenrichtung
    
    % Zuordnung der geometrischen Beziehungen des Walzwerkzeuges zu den
    % geometrischen Beziehungen, welche sich nach Hertz/Boresi ergeben
    R_1_Bor = r_11;
    R_2_Bor = r_21;
    R_1s_Bor = r_12;
    R_2s_Bor = r_22;
    
    % Kontaktbelastung nach Hertz
    % Elastische Kontaktanalyse
    % Spannungen nach Boresi
    a1 = (1/4)*((1/R_1_Bor)+(1/R_2_Bor)+(1/R_1s_Bor)+(1/R_2s_Bor));
    a2 = (1/4)*sqrt((((1/R_1_Bor)-(1/R_1s_Bor))+((1/R_2_Bor)-(1/R_2s_Bor)))^2-4*((1/R_1_Bor)-(1/R_1s_Bor))*((1/R_2_Bor)-(1/R_2s_Bor)).*((sin(alpha_w)).^2));
    A = a1-a2;
    b1 = (1/4)*((1/R_1_Bor)+(1/R_2_Bor)+(1/R_1s_Bor)+(1/R_2s_Bor));
    b2 = (1/4)*sqrt((((1/R_1_Bor)-(1/R_1s_Bor))+((1/R_2_Bor)-(1/R_2s_Bor)))^2-4*((1/R_1_Bor)-(1/R_1s_Bor))*((1/R_2_Bor)-(1/R_2s_Bor)).*((sin(alpha_w)).^2));
    B = b1+b2;
    e7 = A+B;
    delta = (1/e7)*(((1-ny_1^2)/E_1)+((1-ny_2^2)/E_2));
    
    % Regressionsvorschrift für k = f(B/A) und Berechnung ks
    BA_Bor = B./A;
    k(BA_Bor<=200) = 10.^(-0.6131*log10(BA_Bor(BA_Bor<=200)) - 0.02168*log10(BA_Bor(BA_Bor<=200)));
    k(BA_Bor>200) = 10.^(-0.9111*log10(BA_Bor(BA_Bor>200)) + 0.6615*log10(BA_Bor(BA_Bor>200)));
    ks = sqrt(1-(k.^2));
    
    % Relativer Randabstand
    kzb = linspace(0,1,100); % [-] Relativer Randabstand z/b 0...1 (100 Schritte werden berechnet)
    Phi = acot(kzb); % Hilfsgröße
    
    % Elliptische Hilfsintegrale
    Int_F_Bor = ellipticF(Phi,ks^2);
    Int_H_Bor = ellipticE(Phi,ks^2);
    Int_K_Bor = ellipticK(ks^2);
    Int_E_Bor = ellipticE(ks^2);
    
    % Halbachsen und Tiefe (Abplattungsgeometrie)
    b_e = (3*k.*Int_E_Bor/(2*pi).* F_W_normal.* delta).^(1/3); % [m] Kleine Halbachse b der Kontaktellipse
    a_e = b_e./k; % [m] Große Halbachse a der Kontaktellipse
    c_e = 3*k.*F_W_normal*Int_K_Bor/(2*pi).*(e7)./b_e.*delta; % [m] Annäherung der Kontaktkörper (Tiefe der Abplattung)
    
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

    % Berechnung der Spannungen zur Überprüfung
    sigma_e_x = (M.*(omega_x+ny_2.*omega_xs))*(b_e/delta);
    sigma_e_y = (M.*(omega_y+ny_2.*omega_ys))*(b_e/delta);
    sigma_e_z = -(b_e./delta).*((1/2).*M*((1./n_n)-n_n));
    
    % Berechnung der Kontaktfläche (Modell nach Boresi)
    contact_area_elastic_boresi = a_e*b_e*pi;
    
    % Fallunterscheidung zirkulär/elliptisch und iterativer Abgleich der
    % Fläche contact area_elastic_boresi mit contact_area_elastic_circular
    % oder contact_area_elastic_elliptic
 
    if Fall==1        

        for i = 0.1:tolerance_inc:100
            P_0 = F_W_normal;
            R_1 = D_1/2;
            R_2 = D_2/2;
            R_2s = inf;
            
            R_1 = R_1*i;
            R_1s = R_1;
            
            a_e  = ((3/4)*(((1-ny_1^2)/E_1)+((1-ny_2^2)/E_2))*(R_1*P_0))^(1/3);
            contact_area_elastic_circular = pi*a_e*a_e;
            contact_area_elastic_elliptic = NaN;
            if abs(contact_area_elastic_boresi)>=abs(contact_area_elastic_circular)
                factor_notch = i;
            else
                break
            end
        end
    
    else

        for i = 0.1:tolerance_inc:100
            P_0 = F_W_normal;
            R_1 = D_1/2;
            R_2 = D_2/2;
            R_2s = inf;
            
            R_1 = R_1*i;
            R_1s = R_1;

            theta_r = rad2deg(acos(R_1/(R_1+(2*R_2))));
            theta = [0 10 20 30 35 40 45 50 55 60 65 70 75 80 85 90];
            alpha = [inf 6.612 3.778 2.731 2.397 2.136 1.926 1.754 1.611 1.486 1.378 1.284 1.202 1.128 1.061 1.00];
            beta  = [0 0.319 0.408 0.493 0.530 0.576 0.604 0.641 0.678 0.717 0.759 0.802 0.846 0.893 0.944 1.00];
            gamma = [nan 0.851 1.220 1.453 1.550 1.637 1.709 1.772 1.828 1.875 1.912 1.944 1.967 1.985 1.996 2.00];
        
            alpha_when_theta = interp1(theta,alpha,theta_r); % Berechnet alpha aus der Interpolation der Hertzschen Koeffizienten und dem Wert theta
            beta_when_theta = interp1(theta,beta,theta_r);  % Berechnet beta aus der Interpolation der Hertzschen Koeffizienten und dem Wert theta
            gamma_when_theta = interp1(theta,gamma,theta_r); % Berechnet gamma aus der Interpolation der Hertzschen Koeffizienten und dem Wert theta
            
            a_e = alpha_when_theta*((3/4)*((2*R_1*R_2)/(R_1+2*R_2))*(((1-ny_1^2)/E_1)+((1-ny_2^2)/E_2))*P_0)^(1/3);
            b_e = beta_when_theta*((3/4)*((2*R_1*R_2)/(R_1+2*R_2))*(((1-ny_1^2)/E_1)+((1-ny_2^2)/E_2))*P_0)^(1/3);
            contact_area_elastic_elliptic = pi*a_e*b_e;
            contact_area_elastic_circular = NaN;
            if abs(contact_area_elastic_boresi)>=abs(contact_area_elastic_elliptic)
                factor_notch = i;
            else
                break
            end
        end
end



