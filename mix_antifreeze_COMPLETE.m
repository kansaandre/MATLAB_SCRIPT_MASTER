%GIT
clear all

%ANTIFREEZE MIXTURE GEOTHERMAL - CALCULATION OF THERMAL PROPERTIES FOR A
%FREEZING POINT OF -7.7_5C OF THREE DIFFERENT ANTIFREEZE FLUIDS

%liter to m3
conv = 0.001; %liter to m3

%Water properties at 0C
k_w = 0.561 %W/m.K - Thermal conductivty 
visc_dyn_w = 1.787 %cP - Dynamic viscosity 
cp_w = 4220 %J/kg.K - Specific heat (isobaric)
rho_w = 1000; %kg/m3 - Density water


%ETHANOL & WATER MIXTURE
    %Ethanol properties at 0C
    k_eth = 0.176 %W/m.K
   
    visc_dyn_eth = 1.2318 %cP
    cp_eth = 2270 %J/kg.K 
    rho_eth = 806.7; %kg/m3 - Density ethanol 

    %In 1 liter mixture, fraction = volume, fraction used directly
    %With conversion from liter to m3 added.
    y_eth = 0.17; %Volume fraction needed to lower freezing point to -7.5C
    y_w = 1-y_eth; %Water volume fraction
    m_eth = rho_eth * y_eth * conv ;%kg - mass of ethanol 
    m_w = rho_w * y_w * conv ;%kg - mass of water
    m_tot = m_eth+m_w;%kg - total amount of mass

   x_eth = m_eth/m_tot;%Mass fraction needed to lower freezing point to -7.5c
   x_w = m_w/m_tot; %Mass fraction water

    %Ethanol neg7_5c solution (negative -7_5 celsius freezing point)
    rho_eth_mix_neg7_5c = rho_w*x_w+rho_eth*x_eth; %kg/m3
    k_eth_mix_neg7_5c = k_w*x_w+k_eth*x_eth; %W/m.K
    visc_dyn_eth_mix_neg7_5c = visc_dyn_w*x_w+visc_dyn_eth*x_eth; %cP
    visc_dyn_eth_mix_neg7_5c = visc_dyn_eth_mix_neg7_5c/1000; %N.s/m2
    cp_eth_mix_neg7_5c = cp_w*x_w+cp_eth*x_eth; %J/kg.K

    eth_neg7_5c_solution = cell(4,2); %Creating an cell array for storage

    eth_neg7_5c_solution(1,1) = {rho_eth_mix_neg7_5c}; %Store in first column
    eth_neg7_5c_solution(1,2) = {'kg/m3'};
    eth_neg7_5c_solution(2,1) = {k_eth_mix_neg7_5c}; %Store in second column
    eth_neg7_5c_solution(2,2) = {'W/m.K'};
    eth_neg7_5c_solution(3,1) = {visc_dyn_eth_mix_neg7_5c}; %...
    eth_neg7_5c_solution(3,2) = {'N.s/m2'};
    eth_neg7_5c_solution(4,1) = {cp_eth_mix_neg7_5c}; %...
    eth_neg7_5c_solution(4,2) = {'J/kg.K'};

%METHANOL & WATER MIXTURE
    %Methanol properties at 0C
    k_met = 0.207 %W/m.K
    visc_dyn_met = 0.823 %cP
    cp_met = 1370 %J/kg.K 
    rho_met = 810; %kg/m3 - Density metanol 

    %In 1 liter mixture, fraction = volume, fraction used directly
    %With conversion from liter to m3 added.
    y_met = 0.14; %Volume fraction needed to lower freezing point to -7.5C
    y_w = 1-y_met;; %Water volume fraction
    m_met = rho_met * y_met * conv ;%kg - mass of metanol 
    m_w = rho_w * y_w * conv ;%kg - mass of water
    m_tot = m_met+m_w;%kg - total amount of mass

    x_met = m_met/m_tot;%Mass fraction needed to lower freezing point -7.5c
    x_w = m_w/m_tot; %Mass fraction water

    %Methanol & water neg7_5c solution (negative -7_5 celsius freezing point)
    rho_met_mix_neg7_5c = rho_w*x_w+rho_met*x_met; %kg/m3
    k_met_mix_neg7_5c = k_w*x_w+k_met*x_met; %W/m.K
    visc_dyn_met_mix_neg7_5c = visc_dyn_w*x_w+visc_dyn_met*x_met; %cP
    visc_dyn_met_mix_neg7_5c = visc_dyn_met_mix_neg7_5c/1000; %N.s/m2
    cp_met_mix_neg7_5c = cp_w*x_w+cp_met*x_met; %J/kg.K

    met_neg7_5c_solution = cell(4,2); %Creating an cell array for storage

    met_neg7_5c_solution(1,1) = {rho_met_mix_neg7_5c}; %Store in first column
    met_neg7_5c_solution(1,2) = {'kg/m3'};
    met_neg7_5c_solution(2,1) = {k_met_mix_neg7_5c}; %Store in second column
    met_neg7_5c_solution(2,2) = {'W/m.K'};
    met_neg7_5c_solution(3,1) = {visc_dyn_met_mix_neg7_5c}; %...
    met_neg7_5c_solution(3,2) = {'N.s/m2'};
    met_neg7_5c_solution(4,1) = {cp_met_mix_neg7_5c}; %...
    met_neg7_5c_solution(4,2) = {'J/kg.K'};

%ETHYLENE GYLCOL & WATER MIXTURE
    %Ethylene glycol properties at 0C
    k_Egly = 0.256 %W/m.K
    visc_dyn_Egly = 75%cP
    cp_Egly = 2281 %J/kg.K 
    rho_Egly = 870; %kg/m3 - Density Eglyanol 

    %In 1 liter mixture, fraction = volume, fraction used directly
    %With conversion from liter to m3 added.
    y_Egly = 0.19; %Volume fraction needed to lower freezing point to -7.5C
    y_w = 1-y_Egly;; %Water volume fraction
    m_Egly = rho_Egly * y_Egly * conv ;%kg - mass of Eglyanol 
    m_w = rho_w * y_w * conv ;%kg - mass of water
    m_tot = m_Egly+m_w;%kg - total amount of mass

    x_Egly = m_Egly/m_tot;%Mass fraction needed to lower freezing point -7.5c
    x_w = m_w/m_tot; %Mass fraction water

%Ethylene glycol & water neg7_5c solution (negative -7_5 celsius freezing point)
    rho_Egly_mix_neg7_5c = rho_w*x_w+rho_Egly*x_Egly; %kg/m3
    k_Egly_mix_neg7_5c = k_w*x_w+k_Egly*x_Egly; %W/m.K
    visc_dyn_Egly_mix_neg7_5c = visc_dyn_w*x_w+visc_dyn_Egly*x_Egly; %cP
    visc_dyn_Egly_mix_neg7_5c = visc_dyn_Egly_mix_neg7_5c/1000; %N.s/m2
    cp_Egly_mix_neg7_5c = cp_w*x_w+cp_Egly*x_Egly; %J/kg.K

    Egly_neg7_5c_solution = cell(4,2); %Creating an cell array for storage

    Egly_neg7_5c_solution(1,1) = {rho_Egly_mix_neg7_5c}; %Store in first column
    Egly_neg7_5c_solution(1,2) = {'kg/m3'};
    Egly_neg7_5c_solution(2,1) = {k_Egly_mix_neg7_5c}; %Store in second column
    Egly_neg7_5c_solution(2,2) = {'W/m.K'};
    Egly_neg7_5c_solution(3,1) = {visc_dyn_Egly_mix_neg7_5c}; %...
    Egly_neg7_5c_solution(3,2) = {'N.s/m2'};
    Egly_neg7_5c_solution(4,1) = {cp_Egly_mix_neg7_5c}; %...
    Egly_neg7_5c_solution(4,2) = {'J/kg.K'};

%------------------------------------------------------------------------%
