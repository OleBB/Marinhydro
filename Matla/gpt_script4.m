clear all
close all

% ---- % (frå oppgåva)
% INTEGRAL EQUATION IN THE DIFFRACTION PROBLEM
% ...A numerical solution of this equation is obtained by this script (117)
% ---- %

% Define the list of input files and corresponding titles
input_files = {'boxes_1.dat', 'boxes_2.dat', 'boxes_3.dat', 'boxes_4.dat'};
titles = {'L/D = 10', 'L/D = 2', 'L/D = 1', 'L/D = 0.1'};
Lengden = [10;2;1;0.1];

% Specify subfolder name
folderName = 'gpt_script4_output';
% Create subfolder if it doesn't exist
if ~exist(folderName, 'dir')
    mkdir(folderName);
end 

% Specify range and increment for nu
%nu_values = 0:0.05:2; % Range of nu values from 0 to 2 with increments of 0.1

% Loop over each file and process it
for idx = 1:4
    % Load the data from the current file
    file_name = input_files{idx};
    box = load('-ascii', file_name);
    
    L = Lengden(idx);

    %henter AP2 fra script 3.
    outdatscript3 = load(sprintf('outdatscript3_%d.dat', idx));

    AP2 = outdatscript3(1:end, 4);
    AM2 = outdatscript3(1:end, 5);
    
    %For at script5 skal funke, sett likt antall nu_values som diskrete
    %punkter.
    nu_values = linspace(0,2,length(box));
    
    % Extract data
    xm = box(:, 1);
    ym = box(:, 2);
    xp = box(:, 3);
    yp = box(:, 4);
    NN = size(box, 1); % Determine the number of rows dynamically
    
    % Characteristics of geometry
    dx = xp - xm;
    dy = yp - ym;
    ds = ((dx).^2 + (dy).^2).^(1/2);
    bx = 0.5 * (xm + xp);
    by = 0.5 * (ym + yp);
    n1 = -(yp - ym) ./ ds;
    n2 = (xp - xm) ./ ds;
    
    % Initialize array to hold results over all nu values
    SSX = zeros(length(nu_values),1); % Preallocate storage
    omega = zeros(length(nu_values),1);
    H1  = zeros(length(nu_values),1);
    H2 = zeros(length(nu_values),1);
    FK = zeros(length(nu_values),1);
    
    % Loop over nu values and process
    for nu_idx = 1:length(nu_values)
        nu = nu_values(nu_idx); % Current nu value
        % Incoming wave potential
        phi0 = exp(nu * (by - complex(0, 1) * bx));
        % Contributions to the integral equation
        ss = zeros(NN, NN); % Initialize the matrix
        for i = 1:NN
            for j = 1:NN
                xa = bx(j) - bx(i);
                yb = by(j) + by(i);
                zz = nu * (yb - complex(0, 1) * xa);
                f1 = -2 * exp(zz) * (expint(zz) + log(zz) - log(-zz));
                f2 = 2 * pi * exp(zz);
                % LHS contributions
                arg0 = imag(log((xm(j) - bx(i) + complex(0, 1) * (ym(j) - by(i))) /(xp(j) - bx(i) + complex(0, 1) * (yp(j) - by(i)))));
                if j - i == 0
                    arg0 = -pi;
                end
                arg1 = imag(log((xm(j) - bx(i) + complex(0, 1) * (ym(j) + by(i))) /(xp(j) - bx(i) + complex(0, 1) * (yp(j) + by(i)))));
                help1 = (n1(j) * (imag(f1) + complex(0, 1) * imag(f2)) + n2(j) * (real(f1) + complex(0, 1) * real(f2))) * nu * ds(j);
                ss(i, j) = (arg0 + arg1 + help1);
            end
        end
        rhsD = -2 * pi * phi0;
        phiD = ss \ rhsD;
        
        g = 9.81;
        ien = complex(0,1);
        w = sqrt(nu*g);
        omega(nu_idx) = w;
        een1 = exp(-1.*nu); %D=1
        sinusen = sin(0.5*nu*L)/(0.5*nu*L);%sin(0.5*K*L)}/{0.5*K*L}
        
        % Evaluating the exciting force
        XX2 = phiD .* n2 .* ds;
        sXX2 = sum(XX2); % skal være sXX2*-iwp, og så er jo omega=sqrt(gk)g
        SSX(nu_idx) = sXX2;%*(-1)*ien*w; %lese side 37, jon sier gange med -iwp
        
        %Haskind1 (158): X/g = iw*integral(phi0*dphi2dn - phi2*dphi2dn)ds
        for ii = 1
           %(phi2*dphi -phi0*n2)
           %H1ut(i,j) = arg1
           1;
        end
        H1ut = 0;
        H1svar = ien*H1ut;
        H1(nu_idx) = H1svar;        
        %Haskind2 - skal være, (162) X_2 /rhog  = i*AM2
        H2regning = ien*AM2(nu_idx)/nu; %AmplitudeMinusUendelig kommer fra script3, og (162)
        H2(nu_idx) = H2regning;
        
        %Froude-Krylov approx (152) rho*g*L*e^{-KD} * {sin(0.5*K*L)}/{0.5*K*L}
        FK(nu_idx) = L*een1*sinusen;
               
    end
    %%
    
    % Plot aggregated results for all nu values
    figure;
    plot(nu_values, (abs(SSX)));% 
    hold on;
    plot(nu_values, abs(H1),'bo');
    plot(nu_values, (abs(H2)),'rx');% 
    plot(nu_values, abs(FK),'g-.');
    legend('Direct pressure integration', 'Haskind 1', 'Haskind 2','Froude-Krylov', 'Location', 'NorthEast');
    title(['Exciting force |X_2|/gD\rho  in heave for ' titles{idx}]);
    hold off;
    
    % Save the figure
    saveas(gcf, fullfile(folderName, ['aggregated_diffractionplot_' num2str(idx) '_LD_' num2str(idx) '.png']));
    
    %eksportere verdier
    nuverdi = nu_values.'; %
    X2r40 = real(SSX); % 
    X2i40 = imag(SSX); % 
    % Combine them into a matrix
    data4 = [nuverdi, X2r40, X2i40];
    % Create a dynamic filename using sprintf
    outdatname = sprintf('outdatscript4_%d.dat', idx);
    % Save data to the dynamically named file
    %writematrix(data4, outdatname, 'Delimiter','tab');
    
    


end

%SÅ GJENSTÅR SPM: skal jeg plotte noe i 7.10.1?
% eller er det bare å lagre tallan?
%så: 7.10.2
% [ ] calc the exciting force X2/rhog for the boxes, using(155), and 
% [ ] Haskind v1 (158) and Haskind v2 (162) %ser ut som jeg må kåda disse
% selv.....

% 7.10 Plotting results
%%
    %figure;
    %hold on;
    %plot(linspace(1, 10, NN), real(sXX2), 'b');
    %plot(linspace(1, 10, NN), imag(sXX2), 'r');
    %legend('Real \phiD', 'Imag \phiD', 'Location', 'NorthEast');
    %title(['Relle og imaginære deler for ' titles{idx}]);
    %hold off;
     % Save the figure
    %saveas(gcf, fullfile(folderName, ['diffractionplot_' num2str(idx) '_LD_' num2str(idx) '.png']));