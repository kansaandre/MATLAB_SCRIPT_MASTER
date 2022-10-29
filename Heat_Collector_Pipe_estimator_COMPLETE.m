%GIT
clear all

%START (1)

%HX = Collector = The collector pipe underground collecting heat from 
%the groundwater

%Overall (2)
Q_avg = 13000; % W - Heat transfer 
m_flow_min = 1; %kg/s - Mass flow rate MINIMUM through pipes 
m_flow = m_flow_min; %kg/s - variable
T_res = 275.15;%K - Reservoir temperture steady-state 2c ish, 3c at start
T_h_in = T_res ; T_h_out = T_res;%K

%Material selector
pipe_selection = 1; %0 for alu, 1 for HDPE   


%Other things (3)
g = 9.81; %m/s2 
i = 1; %Program closer 
z=1; %Flow increment/counter. Each increment = increase in flow (liter/s)
stang = 0; %step by 3 used to create arrays with thermal resistances
r=1;%Rows
c=1;%Columns
n=1;%Section value - each increment is a new section of the pipe

%Pipe dimensions 
SDR_selector = 21; %21 and 26

if SDR_selector == 26
HDPE_ROWD_COLUMNT= [0.050 0.002; 0.063 0.0025; 0.075 0.0029; 0.090 0.0035;
                  0.110 0.0042;  0.125 0.0048; 0.140 0.0054; 0.160 0.0062; 
                  0.180 0.0069;  0.200 0.0077;];

elseif SDR_selector == 21
HDPE_ROWD_COLUMNT= [0.050 0.0024; 0.063 0.0030; 0.075 0.0036; 0.090 0.0043;
                  0.110 0.0053;  0.125 0.006; 0.140 0.0067; 0.160 0.0077; 
                  0.180 0.0086;  0.200 0.0096;];
end

d_o = HDPE_ROWD_COLUMNT(r,c); %m - outside diameter (START)
t = HDPE_ROWD_COLUMNT(r,c+1); %m - thickness  (START)

%Pipe HDPE
k_hdpe = 0.51; %W/m.K
k_alu = 237; %W/m.K

%Groundwater properties (0 celsius)
    rho_fg = 1000; %kg/m3 %Density
    cp_fg = 4220; %J/kg.K %Specific heat
    u_fg = 1787 * 10^-6; %N.s/m2 %Dynamic viscosity
    k_fg = 561 * 10^-3; %W/m.K %Thermal conductivity
    v_fg = u_fg/rho_fg; %m2/s %Kinematic viscosity
    alpha_fg = k_fg/(rho_fg*cp_fg); %m2/s %Thermal diffusivity
    Pr_fg = (u_fg*cp_fg)/k_fg; %Pradtl %Prandtl number 
    beta_fg = -40*10^-6; %K^-1 %Expansion coefficent 


%Creating arrays we need (Storage space)
m=1; %Apparently Matlab adjust values itself, no need to create memory spaces

%Storage arrays created
length_neededa = zeros(1,m);
length_neededb = zeros(1,m);
qi=zeros(1,215);
q_iavg_s = zeros(1,m);
Re_ds = zeros(1,m);
es = zeros(1,m);
m_flows = zeros(1,m);
dos = zeros(1,m);
Re_d_bank = zeros(1,m);
fs = zeros(1,m);
R_storage = zeros(1,m);
delta_hs = zeros(1,m);
T_diffa = zeros(1,m);
T_c_in_storage = zeros(1,m);
T_c_out_storage = zeros(1,m);
liter_fluid_storage = zeros(1,m);
    
%Here we choose what solution we use in our system
%ATM (11.08.2022) we have either pure water [0], ethanol mix [1], 
%methanol mix [2] or ethylene glycol mixture [3].

    mix_antifreeze; %Call up the mix_antifreeze script
    
    solution = 1; 
    
    mixture_decider; %Call up the mixture_decider script

%Program (4)
while i == 1 || m_flow <= 3
    %Start values
    T_avg_GSHP = 0;
    T_diff = 0.1;%K - Diff
    l=1;%(5) %m - Reset length for new pipe
    q_iavg = 1;%(5) %W - Reset average heat for new pipe

    if m_flow > 3 % (6) 
        r = r+1; %(7) %Next pipe row
        m_flow = m_flow_min; %(8) %kg/s Reset flow to minimum
        d_o = HDPE_ROWD_COLUMNT(r,c);%(9) %m - retrieve outer diameter
        c = c+1; %(9) second column
        t = HDPE_ROWD_COLUMNT(r,c);%(10) %m - retrieve thickness
        c=1;%(11) reset column value, ready for next pipe later
        z=1;% Counter for flow increments
        stang = stang+3;%Used to organize thermal resistance storage matrix
    end

    while q_iavg < (Q_avg-100) || q_iavg > (Q_avg+100) %(12)
       
        l = l+1;%(14) %m -%Guessed length, decided after what makes sense
       
        %Pipe area etc (17)
        d_i = d_o - 2*t; %m - Outside diameter of pipe
        A_c = (pi/4)*d_i^2; %m2 - crossection of pipe
        A_s_i = pi*d_i*l; %m2 - inside surface of pipe
        A_s_o = pi*d_o*l; %m2 - outside surface of pipe
        u_m = m_flow/(rho_f*A_c);%m/s - mean velocity

        %Find a T_c_out so T_avg_GSHP = -0.5 celsius
        while T_avg_GSHP>272.66 || T_avg_GSHP<272.64 %(18)
        T_c_out = T_res-T_diff;%(21)  %K
        T_c_in = T_c_out - Q_avg/(m_flow*cp_f);%(21) %K - Reduction of temperature GSHP
        T_avg_GSHP = (T_c_in+T_c_out)/2;%(21) %K - Avg delivered to GSHP
        T_diff = T_diff+0.01;%(21) %K
        end

        %Inside pipe - Incropera et al. (2007), Ch 8.5, p.514
    Re_d = (4*m_flow)/(pi*d_i*u_f);%(23) %Reynolds number for circular pipe
      
        if Re_d > 2300 

            rough_HDPE = 0.007*10^-3; %m - absolute pipe roughness for HDPE
            rough_alu = 0.002*10^-3; %m - for aluminium
            %VALUES FROM https://www.engineeringtoolbox.com/surface-roughness-ventilation-ducts-d_209.html

            if pipe_selection == 1
                relative_rough = rough_HDPE/d_i;
            elseif pipe_selection == 0
                relative_rough = rough_alu/d_i;
                
            end
            %Haaland equation - approximation of Colebrook
            f = (1/(-1.8*log10(((relative_rough/3.7)^1.11)+(6.9/Re_d))))^2;
        elseif Re_d <= 2300
            f = 64/Re_d;
        end

        if Re_d >= 10000 %PURE TURBULENT
        nu_d_i_turb = ((f/8)*(Re_d-1000)*Pr_f)/(1+12.7*((f/8)^0.5)*((Pr_f^(2/3))-1)); %Gnielinski - Eq 8.62
        h_i_turb = (nu_d_i_turb*k_f)/d_i; %W/m2.K - Convection Coefficent Inside
        h_i = h_i_turb; %(23)
        elseif Re_d<10000 && Re_d>2300 %TRANSITION REGIME
            %TURBULENT MIX
        nu_d_i_turb = ((f/8)*(Re_d-1000)*Pr_f)/(1+12.7*((f/8)^0.5)*((Pr_f^(2/3))-1)); %Gnielinski - Eq 8.62
        h_i_turb = (nu_d_i_turb*k_f)/d_i; %W/m2.K - Convection Coefficent Inside
            %LAMINAR MIX (Ts_i ≈ constant)
        nu_d_i_lam = 3.66; 
        h_i_lam = (nu_d_i_lam*k_f)/d_i; 
        
        x_turb = (Re_d/10000); %Finding fraction of laminar flow
        x_lam = 1-x_turb;
        h_i = h_i_lam*x_lam+h_i_turb*x_turb;%(23) %An average between laminar & turbulent

        elseif Re_d <= 2300 %PURE LAMINAR
            nu_d_i_lam = 3.66; 
            h_i_lam = (nu_d_i_lam*k_f)/d_i; 
            h_i = h_i_lam; %(23)
        end

        R1 = (1/(h_i*A_s_i));%(23) %K/W - Thermal resistance inside convection
        
        %Pipe - Incropera et al. (2007), Ch 3.3.1, p.117
        if pipe_selection == 1
            k = k_hdpe;
        elseif pipe_selection == 0
            k = k_alu;
        end
        
        R2 = (log(d_o/d_i))/(2*pi*k*l);%(24) %K/W -Thermal resistance pipe wall
   
        T_s_c = zeros(1,m);%K - Array of inner surface temp
        T_s_h = zeros(1,m);%K - Array of outer surface temp
        T_c = zeros(1,m); %K - Array of mean temp of fluid inside pipe
        q = zeros(1,m); %W - Array of heat transfer local
         %storage of key data
        R1s     = zeros(1,m);%K/W - Storage of local R1 for each pipe
        R2s     = zeros(1,m);%..
        R3s     = zeros(1,m);%..
        Rtots   = zeros(1,m);%..

        %Assumptions for first iteration of a new pipe dimension
        T_s_c(1,1) = T_c_in;%(25) %K
        T_s_h(1,1) = T_res;%(25) %K
        T_c(1,1) = T_c_in;%(25) %K
    
        R3 = 1/(140*A_s_o);%(26) %140 is a good guess for conv coefficent
        rax = R3;%check
    
        Rtot = R1+R2+R3;%(27) %Start total thermal resistance
        q(1,1) = (T_res-T_c(1,1))/(Rtot);%(27) %Start heat transfer 
    
        h = 5; %s - step time (2)
        lx = u_m*h; %m - Pipe section length per step
        m_i = lx*A_c*rho_f;%(28) %kg - Water amount inside pipe section
        
        n=1;%reset previous iteration

    while T_c(1,n) <= T_c_out %(29)
        %Find energy entering stream and temp diff because of it
        E_in = q(1,n)*(lx/l)*h;%(30) %J
        T_c(1,n+1) = T_c(1,n)+E_in/(m_i*cp_f);%(31) %K
        T_s_c(1,n+1) = R1*q(1,n)+T_c(1,n+1);%(31) %K
        T_s_h(1,n+1) = R2*q(1,n)+T_s_c(1,n+1);%(31) %K
        %Find R3 (32)
        Ra_d = (g*beta_fg*(T_s_h(1,n+1)-T_res)*(d_o^3))/(v_fg*alpha_fg);%Rayleight
        Nu_d_o = (0.60+(0.387*(Ra_d^(1/6)))/(1+(0.559/Pr_fg)^(9/16))^(8/27))^2;
        h_o = (k_fg/d_o)*Nu_d_o;      
        R3 = 1/(h_o*A_s_o);%(32) %K/W - find R3 for next iteration
        R3s(1,n+1) = R3;%K/W - store local R3
        Rtot = R1+R2+R3;%(33) %K/W - Total resistance
        Rtots(1,n+1) = Rtot;%K/W - store total 

        q(1,n+1) = (T_res-T_c(1,n+1))/Rtot;%(33) %W -heat input in next iteration  
                     
        n=n+1;%(34) %next iteration
    end
        q_iavg=sum(q)/(n-1);%(35) %W-This must be equal to Q_avg, a difference 
                         %indicate wrongly guessed value l
        
    end
    l=l; %(13)


    n=n-1; %Correction due to n icreasing 1 even though while loop is false
           %, which makes an iteration without any value

    %Heat exchanger analysis with different methods
        %Logaritimic Temperature Difference (15)
        Delta_T1 = T_h_in-T_c_out;%K
        Delta_T2 = T_h_out-T_c_in;%K
        LogTemp = (Delta_T1-Delta_T2)/(log(Delta_T1/Delta_T2));%(15) %K -LMTD
        UA_log = Q_avg/LogTemp;%W/m2.K
        %NTU - Efficency data (15)
        Q_max = m_flow*cp_f*(T_res-T_c_in); %W - Qmax possible in HX
        e = Q_avg/Q_max;%(15) %Efficency of our HX, this is decided from T_c_out
        NTU = -log(1-e); %Used since Cr = 0, NTU = UA/Cmin = 1/Rtot.Cmin
        UA_e = NTU*(m_flow*cp_f);%When finding U (overall transfer coefficent)
                                %We can find needed surface area A..
        Rtotmin = 1/UA_e;%Total resistance to achieve given heat power for the
                          %efficency/temperature difference specified
    
        U_c = 1/(Rtotmin*A_s_i);%W/m2K -Needed overall heat coefficent cold
        U_h = 1/(Rtotmin*A_s_o);%W/m2K -Needed overall heat coefficent hot
 %%
    R_storage(z,r+stang) = R1;
    R_storage(z,r+1+stang) = R2;
    R_storage(z,r+2+stang) = R3;
    R_storage(z,r+3+stang) = Rtot;

    length_neededa(z,r) = l;
    length_neededb(r,z) = l;
    
    q_iavg_s(r,z) = q_iavg;
    es(r,z) = e;
    m_flows(r,z) = m_flow;
    Re_d_bank(z,r)=Re_d;
    dos(r,z) = d_o;
    fs(r,z) = f;
    T_diffa(r,z) = T_diff;
    T_c_in_storage(z,r)=T_c_in;
    T_c_out_storage(z,r)=T_c_out;

    %Time fluid use
    time=h*n;%s - Timestep (h) x iteration step (n) = time fluid used

     %Pressure loss (19) (Darcy-Weisbach
    delta_p = f*(rho_f/2)*(u_m^2)*(1/d_i)*l;% %Pa -Pressure drop pipe
    delta_pb = delta_p/10^5; %Bar - Conversion to a easier unit
    delta_h = delta_p/(rho_f*g);
    delta_hs(z,r) = delta_h;
   
    %Fluid in pipe
    liter_fluid = A_c*l*1000; %liter - crossection x length x conversion
    liter_fluid_storage(z,r) = uint16(liter_fluid); %Liter fluid needed for each
                                            %pipe and its configuration


    m_flow = m_flow+0.25;%(22) kg/s
    z=z+1;
   
    if r == size(HDPE_ROWD_COLUMNT,1) %(20)
        i = 0; %(16)
    end
    
end
    

%Extraction point
z=z-1;
m_flows=m_flows(1:r,1:z);
dos=dos(1:r,1:z);
%Re_d_bank =Re_d_bank(1:r,1:z);
length_neededa=length_neededa(1:z,1:r);
length_neededb=length_neededb(1:r,1:z);
HDPE_ROWD_COLUMNT=HDPE_ROWD_COLUMNT*1000; 
%%

close all

f1 = figure; 
bar3(m_flows(1,1:z),R_storage(1:z,1:40))
colormap("jet")
ylim([1 3])
ylabel('flow rate (liter/s)')
xlabel('Outside Diameter (mm)')
zlabel('Thermal resistance (K/W)')
xticks([4 8 12 16 20 24 28 32 36 40 44 48 52 56])
xticklabels({50,63,75,90,110,125,140,160,180,200})
grid on


f2 = figure;
bar3(m_flows(1,1:z),delta_hs(1:z,1:r))
colormap("jet")
ylabel('flow rate (liter/s)')
xlabel('Outside Diameter (mm)')
zlabel('Head loss (m)')
ylim([1 3])
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14])
xticklabels({50,63,75,90,110,125,140,160,180,200})
zticks([0 0.1 0.5 1 4 10 40 110])
set(gca, 'ZScale', 'log');

f3 = figure;
bar3(m_flows(1,1:z),length_neededa(1:z,1:r))
colormap("jet")
ylim([1 3])
ylabel('flow rate (liter/s)')
xlabel('Outside Diameter (mm)')
zlabel('Length of pipe needed (m)')
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14])
xticklabels({50,63,75,90,110,125,140,160,180,200})

f4 = figure;
bar3(m_flows(1,1:z),Re_d_bank(1:z,1:r))
colormap("jet")
ylim([1 3])
ylabel('flow rate (liter/s)')
xlabel('Outside Diameter (mm)')
zlabel('Reynolds number')
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14])
xticklabels({50,63,75,90,110,125,140,160,180,200})

magna3_25_80 = [6.22;5.22;4.30;3.44;2.81;2.11;1.56;1.04;0.7];

delta_hs = delta_hs*1.5; %increase 50% to compensate for GSHP and fittings ++

f5 = figure;
plot(m_flows(1,1:z),delta_hs(1:z,4))
xlabel('flow rate (liter/s)')
ylabel('Pressure drop in head (m)')
legend('system characteristic curve')
