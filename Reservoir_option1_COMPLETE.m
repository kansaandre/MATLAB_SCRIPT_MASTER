%RESERVOIR ESTIMATION

%DETERMINING TEMPERATURE GRADIENT AT RESERVOIR WALL

%SOIL DATA - SOLID WALL
k = 2; %W/m.K. - Thermal conductivity fully saturated condition (Johansen (1977) - master document)
Cvol_soil = 3*10^6; %J/m3.K - Volumetric heat capacity of soil (Sundberg (1988) - master document)
rho_soil = 1900; %kg/m3 - density for "sand wet" (https://www.engineeringtoolbox.com/dirt-mud-densities-d_1727.html)
Cmass_soil =Cvol_soil/rho_soil; %J/kg.K - heat capacity mass of soil 

%Boundary data
T_UA = 4 + 273.15; %Kelvin - Temperature at wall to unaffected soil (A)
T_Res = 2 + 273.15; %Kelvin - Temperature at wall to reservoir (B)
n =10; %number of nodes - note, does NOT include wall boundaries (INPUT!)
%------------------------------------------------------------------------%

%NATURAL SETUP - PURE DIFFUSION 
l_Res = 4;%m - minimum length/width of reservoir square
l = 25;%m - length from UA to reservoir wall
julia = 1;
Ts = zeros(n,1);%Creating storage space
Result_storage = zeros(1,n+2);%Creating storage space
dT_dx = zeros(1,1); %Creating storage space
A_storage = zeros(1,1);
q_storage = zeros(1,1);
l_Res_storage = zeros(1,1);
    
while l_Res <= 4
%--------------------
            %Step 1 - Grid generation
      %FVM data
    dx = l/n; %m - length between nodes (and therefore control volume)
    delta = tand(45)*dx; %m - See figure, 1/2 of width increase of square area
    
    syms T [1,n+2]; %Node that last number indicate node numbers, +2 due to boundary
    T_distri = T;
    T_distri(1,1) = T_UA; %Setting boundary A - Unaffected "wall"
    T_distri(1,n+2) = T_Res; %Setting boundary B - Reservoir wall
    
            %Area determination
        i = n+2; %Number of nodes + the two boundary wall "nodes"
        A_distri = zeros(1,i); %Areas avaliable including the boundary areas
        A_distri(1,i) = l_Res*l_Res;%m2 - area for one of the reservoir sides
        
        length = zeros(1,i); %m - length of area perpendicular side
        length(1,i) = l_Res; %m - length of side perpendicular to x direction
        
        while i ~= 1
            length(1,i-1) = length(1,i)+2*delta;%m - right to left
            A_distri(1,i-1) = length(1,i-1)*length(1,i-1);%m2 - right to left
            i = i-1;
        end
        
    %--------------------
            %Step 2 - Discretation 
    
                % NOTE that n = real nodal points while T_distri include 
                % temperature found at boundary walls "nodes" as well
    
    syms T [n 3];%Will self adjust to whatever we put inside it
    MATRIX = T;
    mxr = 1; %Row data variable 
    mxc = 1; %Column data variable
    
                %Class A - Boundary with unaffected soil
    for ma = 2 %Note that 1 is the boundary wall (Therefore 2 will be the first "real" node)
    MATRIX(mxr,mxc) = ((k/dx)*A_distri(1,ma)+(2*(k/dx)*A_distri(1,ma-1)))*T_distri(1,ma);
    MATRIX(mxr,mxc+1) = -(k/dx)*(A_distri(1,ma))*T_distri(1,ma+1);
    MATRIX(mxr,mxc+2) = (2*(k/dx)*(A_distri(1,ma-1)))*T_distri(1,ma-1);
    end
    
                %Class B - Boundary with reservoir wall
    for mb = n+1 %Note that n+2 is the boundary wall (n+1 will be the last real node)
    MATRIX(n,mxc) = (2*(k/dx)*A_distri(1,mb)+(k/dx)*A_distri(1,mb-1))*T_distri(1,mb);
    MATRIX(n,mxc+1) = -(k/dx)*(A_distri(1,mb-1))*T_distri(1,mb-1);
    MATRIX(n,mxc+2) = (2*(k/dx)*(A_distri(1,mb)))*T_distri(1,mb+1);
    end
    
                %Class non-boundary
    for mn = n:-1:3 %Nodes between boundary nodes (right to left)
       MATRIX(n-mxr,mxc) = ((k/dx)*A_distri(1,mn)+(k/dx)*A_distri(1,mn-1))*T_distri(1,mn);
       MATRIX(n-mxr,mxc+1) = -((k/dx)*A_distri(1,mn-1))*T_distri(1,mn-1);
       MATRIX(n-mxr,mxc+2) = -((k/dx)*(A_distri(1,mn))*T_distri(1,mn+1));
    mn
       mxr = mxr + 1;
    end
    MATRIX;
    
    %--------------------
    
    %Step 3 - Solving our system matrix
    
        %Creating a compare matrix
        %- This matrix can be used to see which temperature nodes/variables
        %is active in the system matrix 'MATRIX'. 
    
        %All we do below this point is to get our system matrix in the correct
        %form to be able to be solved by MATLABS linear solver
    
    syms T [n n] 
    def = T;
    for r = 1:n
    syms T [1, n+2]
    T = T(1,2:n+1);
    def(r,:) = T;
    end
    
    def = sym2cell(def);
    def = string(def);
    
        %Splitting our system matrix 'MATRIX' to be able to compare it to our 
        %compare matrix
        
    MATRIX_s = MATRIX;
    MATRIX = sym2cell(MATRIX);
    MATRIX = string(MATRIX);
    MATRIX;
    
    
        %Save our boundary condition values (not containing variables)
    BOUNDARY = zeros(n,1); %Creating storage, need these values later. Column vector.
    BOUNDARY(1,1) = double(MATRIX_s(1,3));
    BOUNDARY(n,1) = double(MATRIX_s(n,3));
    BOUNDARY;
    
    MATRIX(1,3) = "T0"; 
    MATRIX(n,3) = "T0";
    
    
        %Extract sybols (T1, T2, etc.) - Use to compare with def matrix
    pat = "T" + digitsPattern;
    MATRIX_SYM = extract(MATRIX,pat); %Symbolic part of MATRIX extracted
                                     %We need this to later compare with
                                     %default matrix to see which T component
                                     %is present/active in the system
         %Extract the numbers, needed to be put in again when we are going to
         %linsolve it --> see code at the end
    
        %We only want the actual number -> Filter, filter, filter...
    temp = replace(MATRIX,'(','');
    temp = replace(temp,')','');
    temp = replace(temp,'*','');
    temp = replace(temp,pat,'');
    MATRIX_NUM = temp;
    
        %Perform comparison to compare matrix (default)
    truth=strings(n);
    i = 1;
    complete = 0;
    while complete == 0
        for r = 1:3
            for j = 1:n
                if  def(i,j) == MATRIX_SYM(i,r)
                   truth(i,j) = MATRIX_SYM(i,r); %Matrix made to validate code
                   truth(i,j) = MATRIX_NUM(i,r);
                end
            end
        end
    
    
        if i == n
            complete = 1;
        end
        i = i+1;
    end
        
    for i = 1:n
        for j = 1:n
            if truth(i,j) == ""
                 truth(i,j) = 0;
            end
        end
    end


    for i = 1:n
        for j = 1:n
            truth_num(i,j) = str2num(truth(i,j));
        end
    end

    
        %SOLVING THE SYSTEM USING THE 'truth' MATRIX AND 'BOUNDARY' matrix 
    T = linsolve(truth_num,BOUNDARY);

    
    Ts(:,julia) = T;

            %RESULTING SYSTEM TEMPERATURE DISTRIBUITON
    Result = zeros(1,n+2);
    Result(1,1) = T_UA;
    Result(1,2:(n+1)) = T;
    Result(1,n+2) = T_Res;
    Result_storage(julia,:) = Result;

    
    %Heat transfer
    %Finding temperature gradients (dT/dx)
    l_Res_storage(julia,:) = l_Res; 

    dT_dx(julia,1) = (Result_storage(julia,n+1)-Result_storage(julia,n))/dx
    q =-k*A_distri(1,n+1)* dT_dx(julia,1) * 6; %W - Total heat transfer
    q_storage(julia,:) = q; %W - storage
    
    julia = julia + 1 %Without semicolon, able to track iteration progress
    l_Res = l_Res + 1;
    n

 
end
%%
close all
    x = 0:l/(n+1):l;
    
    f1 = figure;
    plot(x,Result_storage)
    title(["Temperture distribution from the unaffect soil to the ", "reservoir wall with different wall lengths (reservoir area)"]);
    xlabel('Distance to the reservoir from unaffected soil [m]');
    ylabel('Temperature in soil [K]');
    legend('4m','5m','6m','7m','8m');
    
    f2 = figure;
    plot(l_Res_storage,q_storage)
    title(["Total heat transfer from the unaffect soil to the reservoir","with different wall lengths (reservoir area) [pure conduction]"]);
    xlabel('Wall lengths [m]');
    ylabel('Heat transfer [W]'); 



%-----------------------------------------------
%}

