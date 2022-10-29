%GIT

clear all

%---------------------------SOIL PRESSURE CALCULATIONS------------------%

    %-------------------SOIL/EARTH PRESSURE------------------------%
%Calculating soil pressure acting on pipe
    %https://www.engineeringtoolbox.com/underground-pipe-pressure-soil-transport-d_2145.html
    %_imp = imperical units
    %Pressure caused by soil

rho_sw = 1900; %kg/m3 - moist sand density
g = 9.81; %m/s2 - Gravitational constant

    h = 9; %m - maximum depth to pipe
    d_wt = 7; %m - Lowest water table possible to utilize     
    hw = h-d_wt; %m - depth from water table to burried pipe (water)
    hb = 3; %m - Reservoir is 3 m tall
    hs = h-hb; %m - depth from surface to backfill (soil)
    hsw = hs-d_wt; %m - depth of wet soil (water table to backfill)
      
rho_sb = 1200; %kg/m3 - Below water level soil has reduced density due to   
               %buoyancy effect.
rho_gb = 1100; %kg/m3 - gravel density due to bouyancy effect

    %Pressure from the soil alone
p_soil = rho_sw*g*d_wt+rho_sb*g*hsw+rho_gb*g*hb;%Pa (N/m2) Pressure soil

    %Pressure caused by groundwater
rho_w = 1000; %kg/m3 - density of water at 0 celsius 
p_water = rho_w*g*hw; %Pa - Pressure from the water alone


pe_tot = p_soil+p_water; %Pa - Soil pressure (combined)
PE = pe_tot/10^5; %Bar - soil pressure

    %--------------------SURCHARGE PRESSURE---------------------%
%Boussinesq's formula equation
Q = 250 * 1000 * 9.81; %N - 200 ton * 1000 kg/ton * 9.81 N/kg 
                      %This is the point load of a typical 3-4 bedroom home
z = 10;% m - maximum depth of pipe
r = 8; % m - Appropriate horizontal distance from center of gravity of a 
      % home to the pipe

    PES = ((3*Q)/(2*pi*(z^2))) * (1/(1+((r/z)^2)^(5/2)));
    PES = PES / (10^5); % bar

    %----------------TOTAL PRESSURE ACTING ON PIPE--------------------%

P = PE+PES %Bar -> Total pressure = earth pressure + surcharge pressure
P_imp = P*14.5037738; %psi

%-------------------PIPE MATERIAL/DIMENSIONS/PROPERTIES-----------------%
E_imp = 28250; %lb/in2 - Modulus of elasticity values LONG TERM HDPE
SDR = 21;%SDR - Standard dimension ratio, standard used for pipes
         %SDR = Outside Diameter / Thickness of wall
PS_imp = (2*E_imp)/(3*0.149*((SDR-1)^3)); %psi - Pipe stiffness

%--------------------RING DEFLECTION-----------------------%

    %Modifed Iowa formula from https://bit.ly/3OlEHpG
    
Dl = 1.5; %Time lag factor 
K = 0.0843;%Bedding constant 
                        %- Pressure acting on the outside of pipe in psi 
Es_imp = 3000; %psi - Modulus of soil reaction (1000 for crushed rock)
               %Also known as degree of compaction of bedding
dfx_ratio = (Dl*K*P_imp)/((0.149*PS_imp)+(0.061*Es_imp));% frac- deltax/d

dfx_ratio_percentage = dfx_ratio*100; % %
df_allowed_percentage = 7.5; %Allowable deflection of HDPE for SDR>=21(USBR)

if dfx_ratio_percentage <= df_allowed_percentage 
    Criteria_deflection_x = 1;
else 
    Criteria_deflection_x = 0;
end
%==
dfy_ratio = abs(((P_imp)/Es_imp)*(((K*Es_imp)/(0.149*PS_imp))*(((0.0595*Es_imp)/(0.149*PS_imp+0.061*Es_imp))-1)));
   
dfy_ratio_percentage = dfy_ratio*100;

if dfy_ratio_percentage <= df_allowed_percentage 
    Criteria_deflection_y = 1;
else 
    Criteria_deflection_y = 0;
end

close all

Ring_Deflection = figure;
bar(1,dfx_ratio_percentage,'FaceAlpha',0.3);
hold on
bar(1,df_allowed_percentage,'FaceAlpha',0);
bar(2,dfy_ratio_percentage,'FaceAlpha',0.3);
hold on
bar(2,df_allowed_percentage,'FaceAlpha',0);
hold off
ylabel('%')
xticklabels("")
legend('Horizontal_{actual}','Horizontal_{allowed}','Vertical_{actual}','Vertical_{allowed}')

%-------------------------BUCKLING CALCULATION--------------------------%
%Buckling analysis
 
E_imp = 28250; %psi - Modulus of elasticity values LONG TERM HDPE
hw_imp = 3.2808399 * hw; %ft - Maximum distance from pipe to water table
h_imp = h * 3.2808399; %ft - Maximum depth of pipe (surface to pipe)

Rb = 1 - 0.33*(hw_imp/h_imp); %Dimensionless - Buoyancy reduction factor
B = 1/(1+(4*exp(-0.065*h_imp)));%Dimensionless - Soil elastic support factor
    %Finding the allowable external pressure for a constrained pipe:
Pca = (5.65) * sqrt(Rb*B*Es_imp*(E_imp/(12*(SDR-1)^3))); 
Pca_SI = Pca * 0.0689475729; %Bar 

if Pca_SI >= P
    Criteria_buckling = 1;
else 
    Criteria_buckling = 0;
end
Buckling = figure;
bar(1,P,'FaceAlpha',0.5);
hold on
bar(1,Pca_SI,'FaceAlpha',0);
hold off
ylabel('bar')
xticklabels("")
legend('Buckling_{actual}','Buckling_{allowed}')

%---------------------WALL COMPRESSIVE STRESS-------------------------%
%Compressive stress analysis 

S = (P*(SDR-1))/2;

S_allowed_imp = 1000;% psi at 23C for PE3408 (HDPE class) according to 
                     %AWWA_M55_2006_PE_Pipe_Design_and_Install (page 64)
S_allowed = S_allowed_imp*0.0689475729; %Bar, maximum allowed compressive
                                        %stress (Long-term) for PE3408
Wall_Compression = figure;
bar(1,S,'FaceAlpha',0.5);
hold on
bar(1,S_allowed,'FaceAlpha',0);
hold off
ylabel('bar')
xticklabels("")
legend('Compression_{actual}','Compression_{allowed}')

if S <= S_allowed %Allowable deflection of HDPE for SDR>=21(USBR)
    Criteria_wall_compression = 1;
else 
    Criteria_wall_compression = 0;
end

%-------------------CONCLUSION-----------------%

if (Criteria_deflection_x & Criteria_deflection_y & Criteria_buckling & ...
        Criteria_wall_compression) == 1
    PIPE = 'Approved'
else
    PIPE = 'Not-Approved'

end

