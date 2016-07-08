%------------------------------------------------
% -file = DAIS_IceFlux_model.R
% -Antarctic Ice Sheet (AIS) model
%------------------------------------------------
% -Function of the Shaffer DAIS model 2014
% -This function/DAIS model estimates the sea-level equivalence of Antarctic ice sheet melt
% -Model can be found in Shaffer_GMDD_2014
% -Matlab codes and forcings can be found at:
% -www.dcess.dk under "DCESS models"
%
% -Author: Kelsey Ruckert (klr324@psu.edu)
%------------------------------------------------
% -June 17, 2014 #Updates June 10 2015 # Coded from R to Matlab March 7th 2016
%------------------------------------------------

function llik = iceflux(theta)

global data;

% theta = parameters
% data = forcings to run the model

gamma = theta(1) ;           %sensitivity of ice flow to sea level
alpha = theta(2) ;           %sensitivity of ice flow to ocean subsurface temperature
mu = theta(3) ;              %Profile parameter related to ice stress [m^(1/2)]
eta = theta(4) ;             %Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
Po = theta(5);               %Precipitation at 0C [m of ice/yr]
kappa = theta(6) ;           %Relates precipitation to temperature [K^-1]
fo = theta(7) ;              %Constant of proportionality for ice speed [m/yr]
ho = theta(8) ;              %Initial value for runoff line calculation [m]
co = theta(9) ;              %Second value for runoff line calculation [m]
bo =  theta(10) ;            %Height of bed at the center of the continent [m]
s =  theta(11) ;             %Slope of the bed

sigma_y = theta(12);         %Standard error

TA = data(:,1) ;            %Air temperature
TO = data(:,2) ;            %Ocean temperature
GSL = data(:,3) ;           %Global sea-level rate of change
SL = data(:,4) ;            %Global sea-level rise

%standards                      %standards = [Tice, eps1, del, eps2, TOo, Volo, Roa, R]
Dsw = 1.03;                     %Density of sea water [g/cm^3]
Dice = 0.917;                   %Density of ice water [g/cm^3]
Drock = 4.0;                    %Density of rock [g/cm^3]
Tice = -1.8;                    %Freezing temperature of sea water
eps1 = Dice/(Drock-Dice);
del = Dsw/Dice;
eps2 = 1/(Drock-Dice);
TOo = 0.72;                     %Present day high latitude ocean subsurface temperature [K]
Volo = 2.4789e16;               %Steady state AIS volume for present day Ta and SL [m^3]
Roa = 1.8636e6;                 %Steady state AIS radius for present day Ta and SL [m]
R = Roa;                        %Initial condition for integration


%end information for hindcast period
endsat = 240000;
enddate = 240010;

%Simple time stepping integration with one year time step; solves for ice sheet
% radius R and from R calculates ice sheet volume Vol and sea level equivalent SLE

Rad(1:endsat)=NaN;            %Radius of ice sheet
Vol(1:endsat)=NaN;            %Ice volume
SLE(1:endsat)=NaN;            %Sea-level equivalent [m]

%Setup the model
for i=1:enddate
    %Calculate total ice flux across the grounding line
    f = fo*((1-alpha)+alpha*((TO(i)-Tice)/(TOo-Tice))^2)/((s*Roa-bo)^(gamma-1));
    hr= ho+co*TA(i); %Calculate ice sheet surface height
    %Calculate the distance from the continent center to the grounding line
    rcon = (bo-SL(i))/s;
    P = Po*exp(kappa*TA(i)); %Calculate precipitation
    beta = eta*P^1/2; %calculate the mass balance gradient
    %Calculate the distance from the continent center to where the runoff line intersects the ice sheet surface
    rR = R-(abs(hr-bo+s*R)*(hr-bo+s*R))/mu;

    %Apply whether or not the ice sheet has components of a marine ice sheet
    if R<=rR && R<=rcon
        Btot = pi*P*R*R; %Total mass accumulation on ice sheet
        F = 0;
        ISO = 0;
        fac = pi*(1+eps1)*(4/3*mu^0.5*R^1.5-s*R*R);
    elseif R>rR && R<=rcon
        Btot = pi*P*R*R-pi*beta*(hr-bo+s*R)*(R*R-rR*rR)-...
        (4*pi*beta*mu^0.5)/5*(R-rR)^2.5+...
        (4*pi*beta*mu^0.5)/3*R*(R-rR)^1.5;
        F = 0;
        ISO = 0;
        fac = pi*(1+eps1)*(4/3*mu^0.5*R^1.5-s*R*R);
    elseif R<=rR && R>=rcon
        Btot = pi*P*R*R;
        F = 2*pi*R*f*del*(s*R-bo+SL(i))^(gamma+1);
        ISO = 2*pi*eps2*(s*rcon*rcon-bo/s*rcon)*GSL(i);
        fac = pi*(1+eps1)*(4/3*mu^0.5*R^1.5-s*R*R)-2*pi*eps2*(s*R*R-bo*R);
    else
        Btot = pi*P*R*R-pi*beta*(hr-bo+s*R)*(R*R-rR*rR)-...
        (4*pi*beta*mu^0.5)/5*((R-rR)^2.5)+...
        (4*pi*beta*mu^0.5)/3*(R*(R-rR)^1.5);
        F = 2*pi*R*f*del*(s*R-bo+SL(i))^(gamma+1);
        ISO = 2*pi*eps2*(s*rcon*rcon-bo/s*rcon)*GSL(i);
        fac = pi*(1+eps1)*(4/3*mu^0.5*R^1.5-s*R*R)-2*pi*eps2*(s*R*R-bo*R);
    end

    dR = (Btot-F+ISO)/fac;
    R = R+dR; %calculate radius
    V = 8/15*pi*mu^0.5*R^2.5-1/3*pi*s*R^3; %Calculate volume
    %calculate sea volume
    Vsea = pi*(2/3*s*(R^3-rcon^3)-bo*(R^2-rcon^2));
    %calulate the volume change over time
    if R<=rcon
        Volt = (1+eps1)*V;
    else
        Volt = (1+eps1)*V-eps2*Vsea;
    end

    Rad(i) = R;
    Vol(i) = Volt;
    SLE(i) = 57*(1-Volt/Volo); %Takes steady state present day volume to correspond to 57m SLE
end

resid_y(1) = 3.9 - (SLE(120000) - mean(SLE(239961:239990)));
resid_y(2) = -11.35 - (SLE(220000) - mean(SLE(239961:239990)));
resid_y(3) = -2.625 - (SLE(234000) - mean(SLE(239961:239990)));
resid_y(4) = 0.002 - (SLE(240002) - mean(SLE(239961:239990)));

sterr_y = [2.1, 4.45, 1.375, 0.0003];

% sigma_y = theta(12);         % Variance

llik = sum(log(normpdf(resid_y, repmat(0,1,4), sqrt(sigma_y + sterr_y.^2))));






