% 240 kyr integration of time-dependent DAIS equations for calculation
% of Antarctic Ice Sheet (AIS) radius (Rad),volume (Vol)
% and sea level equivalent (SLE); see Read me file for more details

clear all
close all

% loads 240 kyr time series at one year intervals of the forcings:
% TA (Antarctic temperature reduced to sea level), SL (sea level),
% TO(high latitude subsurface ocean temperature and GSL (time rate of change
% of SL)

load ('DaisForce')

% Model parameters and their standard values

Tice=-1.8;          %freezing temperature of sea water
dsw=1.03;           % g cm^3, sea water density
dice=0.917;         % g cm^3, ice water density
drock=4.0;          % g cm^3, rock density

bo=775;             %m, bed height at the center of the continent
s=0.0006;           % slope of the bed
mu=8.7;             %m^{1/2), profile parameter
fo=1.2;             %m yr^{-1), constant of proportionality for ice speed
eta=0.012;          %m^(-1/2)yr^(-1/2), relates balance gradient to precipitation
Po=0.35;            % m ice yr^{-1) precipitation at 0 C.
kappa=0.04;         % K^{-1), relates precipitation to temperature
ho=1471;            %m,first constant for runoff line calculation
c=95;               %m K^{-1), second constant for runoff line calculation

del=dsw/dice;
eps1=dice/(drock-dice); 
eps2=1/(drock-dice);        

Roa=1.8636e6;           %m, steady state AIS Rad for present day TA and SL
Volo=2.4789e16;         % m^3, steady state AIS Vol for present day TA and SL
TOo=0.72;               %K, present day high latitude ocean subsurface temperature

% Choice of parameter values; uncomment to make choice 

gamma=1;alpha=0;          %Case 1, Original Oerlemans model (with corrections) 
%gamma=2;alpha=0;         %Case 2, Increased sensitivity of ice flow to sea level
%gamma=1;alpha=0.35;      %Case 3, Increased sensitivity of ice flow to ocean subsurface temperature 
%gamma=2;alpha=0.35;      %Case 4, Increased sensitivity of ice flow to sea level and ocean subsurface temperature

% Initial condition for integration
R=Roa;  

%Simple time stepping integration with one year time step; solves for ice sheet 
% radius R and from R calculates ice sheet volume Vol and sea level equivalent SLE

Rad(1:240000)=NaN;
Vol(1:240000)=NaN;
SLE(1:240000)=NaN;  

for i=1:240010 
  f=fo*((1-alpha)+alpha*((TO(i)-Tice)/(TOo-Tice))^2)/((s*Roa-bo)^(gamma-1));
  hr=ho+c*TA(i);
  rcon=(bo-SL(i))/s;
  P=Po*exp(kappa*TA(i));
  beta=eta*P^0.5;
  rR=R-(abs(hr-bo+s*R)*(hr-bo+s*R))/mu;
    if R<=rR && R<=rcon
       Btot=pi*P*R*R;
       F=0;
       ISO=0;
       fac=pi*(1+eps1)*(4/3*mu^.5*R^1.5-s*R*R);  
    elseif R>rR && R<=rcon
       Btot=pi*P*R*R-pi*beta*(hr-bo+s*R)*(R*R-rR*rR)-...
        (4*pi*beta*mu^.5)/5*(R-rR)^2.5+...
        (4*pi*beta*mu^.5)/3*R*(R-rR)^1.5; 
       F=0;
       ISO=0;
       fac=pi*(1+eps1)*(4/3*mu^.5*R^1.5-s*R*R); 
   elseif R<=rR && R>=rcon
       Btot=pi*P*R*R;
       F=2*pi*R*f*del*(s*R-bo+SL(i))^(gamma+1);
       ISO=2*pi*eps2*(s*rcon*rcon-bo/s*rcon)*GSL(i);
       fac=pi*(1+eps1)*(4/3*mu^.5*R^1.5-s*R*R)-...
       2*pi*eps2*(s*R*R-bo*R); 
    else
        Btot=pi*P*R*R-pi*beta*(hr-bo+s*R)*(R*R-rR*rR)-...
        (4*pi*beta*mu^.5)/5*((R-rR)^2.5)+...
        (4*pi*beta*mu^.5)/3*(R*(R-rR)^1.5); 
         F=2*pi*R*f*del*(s*R-bo+SL(i))^(gamma+1);
        ISO=2*pi*eps2*(s*rcon*rcon-bo/s*rcon)*GSL(i);
        fac=pi*(1+eps1)*(4/3*mu^.5*R^1.5-s*R*R)-...
          2*pi*eps2*(s*R*R-bo*R);  
    end
    
  dR=(Btot-F+ISO)/fac;
  R=R+dR;
  V=8/15*pi*mu^.5*R^2.5-1/3*pi*s*R^3;
  Vsea=pi*(2/3*s*(R^3-rcon^3)-bo*(R^2-rcon^2));
    if R<=rcon
      Volt=(1+eps1)*V;
    else
      Volt=(1+eps1)*V-eps2*Vsea;
    end
  Rad(i)=R;
  Vol(i)=Volt;
  SLE(i)=57*(1-Volt/Volo);           %takes Volo to correspond to 57 m SLE
end
  
save OutDais Rad Vol SLE SL
