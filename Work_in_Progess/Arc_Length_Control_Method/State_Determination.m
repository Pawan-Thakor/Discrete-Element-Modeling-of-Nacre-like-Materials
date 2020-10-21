%Function defination

function [Shear_Stress_n_1, Shear_Plastic_Strain_n_1, Element_G_n_1, Colour_n_1] = State_Determination(Shear_Stress, Shear_Plastic_Strain, Shear_Delta_Strain, Element_G, Colour,p)
% All the variable are for nth Load Step and we are evaluating stress for
% n+1 load step
        
        if norm(Shear_Delta_Strain)==0
            Shear_Stress_n_1 = Shear_Stress;
            Shear_Plastic_Strain_n_1 = Shear_Plastic_Strain;
            Colour_n_1 = Colour;
            Element_G_n_1 = Element_G;
            return;
        end
        
        I3 = find(Shear_Plastic_Strain(:,1)>=p.gamma_p_s);
        I8 = setdiff([1:p.Num_Elements]', I3);
        
        % Trial Stress Evaluation
        tr_sigma(:,1) = Shear_Stress(:,1) +  p.G*Shear_Delta_Strain(:,1);
        
        % Yeild Stress Evaluation
        % For Hardening/Elastic Region
        Shear_sigma_yeild(I8,1) = p.tau_y + p.H_1*Shear_Plastic_Strain(I8,1);
        % For Softening Region
        Shear_sigma_yeild(I3,1) = p.tau_y + p.H_1*p.gamma_p_s + p.H_2*(Shear_Plastic_Strain(I3,1)-p.gamma_p_s);
   
        % Evaluation of Change in Plastic Strain 
        % (Evaluating Delta Plastic Strain)
        fr = abs(tr_sigma(:,1))-Shear_sigma_yeild(:,1);
        I0 = find(fr<=0);   
        I1 = setdiff(1:p.Num_Elements, I0);
        
        Common_Elements = intersect(I0,I3);
        % I0: Elements in Elsstic/Unloading region
        I0 = setxor(I0,Common_Elements);
        Common_Elements = intersect(I1,I3);
        % I1: Elements in Hardening region
        I1 = setxor(I1,Common_Elements);
        
        Shear_Delta_plastic_strain(I8,1) = fr(I8)/(p.G + p.H_1);
        Shear_Delta_plastic_strain(I3,1) = 0;
        Shear_Delta_plastic_strain(Shear_Delta_plastic_strain<0)=0;
        
        
        % Evaluation of Plastic Strain
        Shear_Plastic_Strain_n_1(:,1) = Shear_Plastic_Strain(:,1) + Shear_Delta_plastic_strain(:,1);
        
        % Evaluation of Stress
        Shear_Stress_n_1(I0,1) = tr_sigma(I0,1);
        Shear_Stress_n_1(I1,1) = tr_sigma(I1,1) - p.G*sign(tr_sigma(I1,1)).*Shear_Delta_plastic_strain(I1,1);
        Shear_Stress_n_1(I3,1) = Shear_Stress(I3,1) + p.G_T_2*Shear_Delta_Strain(I3,1);
        
        % Updating Stiffness values According to different regions
        Element_G_n_1(I0) = p.G;
        Element_G_n_1(I1) = p.G_T_1;
        Element_G_n_1(I3) = p.G_T_2;
        
        % Storing Colour information for Stress-Strain Curve
        % I0: Elements in Elasstic and Unloading region  (As I0 is defined first below, this I6 would overwrite the values of I0)
        % I6: Elements in Ealstic Region
        % I1: Elements in Hardening Region
        % I3: Elements in Softening Region
       % Storing Colour information for Stress-Strain Curve
        I6 = I0;
        I4 = find(Shear_Plastic_Strain_n_1(:,1)>0);
        Common_Elements = intersect(I6,I4);
        I6 = setxor(I6,Common_Elements);
        Colour_n_1(I0,1) = p.C(2);
        Colour_n_1(I6,1) = p.C(1);
        Colour_n_1(I1,1) = p.C(3);
        Colour_n_1(I3,1) = p.C(4);   

    return;
end