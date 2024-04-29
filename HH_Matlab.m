%Saksham Shrestha
%Dr. Shijie Zhou
%CPB517
%Fall22
%Modeling of an action potential based on the Hodgkin Huxley model

%equilibrium potentials(mV)
ENa = 55;  
EK = -75;   
EL = -60;   

%maximum conductances(mS/cm^2)
gNamax = 120; 
gKmax = 30; 
gLmax = 0; 

%area of cell and membrane capacitance
Area = 1; %cm^2
Cap = 1.0; %uF/cm^2    

%time steps(ms)
dt = 0.1;     
tdur = 100;    
tsteps = ceil(tdur/dt);

%gating variables
m = zeros(1,tsteps);
n = zeros(1,tsteps);
h = zeros(1,tsteps);

%membrane potential and time arrays
Vm = zeros(1,tsteps);
t = (0:tsteps-1)*dt;

%initial condition
v0 = -70;
Vm(1) = v0;
inpI = zeros(1, tsteps); %injected current 

%first run
inpI(t>= 20 & t<= 30) = 10; 

%second run part a
% inpI(t>= 10 & t<= 15) = 3; 
% inpI(t>= 20 & t<= 25) = 3; 

%second run part b
% inpI(t>= 10 & t<= 15) = 3; 
% inpI(t>= 25 & t<= 30) = 3; 

%third run 
% inpI(t>= 10 & t<= 80) = 10; 

%fourth run 
% inpI(t>= 10 & t<= 80) = 3; 



m(1) = alpham(Vm(1))/(alpham(Vm(1))+betam(Vm(1)));
n(1) = alphan(Vm(1))/(alphan(Vm(1))+betan(Vm(1)));
h(1) = alphah(Vm(1))/(alphah(Vm(1))+betah(Vm(1)));

% the loop as described in the report
for i = 1: tsteps-1

  % time constants
    tau_m = 1 /(alpham(Vm(i))+betam(Vm(i)));
    tau_h = 1 /(alphah(Vm(i))+betah(Vm(i)));
    tau_n = 1 /(alphan(Vm(i))+betan(Vm(i)));
  
  % steady state values
    m_inf = alpham(Vm(i))/(alpham(Vm(i))+betam(Vm(i)));
    n_inf = alphan(Vm(i))/(alphan(Vm(i))+betan(Vm(i))); 
    h_inf = alphah(Vm(i))/(alphah(Vm(i))+betah(Vm(i)));
  
  % solving for m, h, and n using the Exponential-Euler method
    m(i+1) = m_inf+(m(i)-m_inf)*exp(-dt/tau_m);
    n(i+1) = n_inf+(n(i)-n_inf)*exp(-dt/tau_n);
    h(i+1) = h_inf+(h(i)-h_inf)*exp(-dt/tau_h);
  
  % updating the conductances
    gNa = gNamax*(m(i+1)^3)*h(i+1);   
    gK = gKmax*(n(i+1)^4);                  
    gtot = gNa+gK+gLmax;                            
    gE = gNa*ENa+gK*EK+gLmax*EL;                 
    
  % updating steady state Vm
    inf_Vm = (gE + inpI(i)/Area) / gtot;

  % updating membrane time constant
    tau_Vm = Cap/gtot;

  % solving for Vm 
    Vm(i+1) = inf_Vm + (Vm(i)-inf_Vm)*exp(-dt/tau_Vm);
end

%Plots
figure (1)
subplot(3,1,1)
plot(t, Vm); ylim([-110 40]);legend('Vm');
ylabel('V_m (mV)', 'fontsize', 20)
subplot(3,1,2)
plot(t,inpI)
ylabel('Input Current (muA/cm^2)', 'fontsize', 20)
gating = [m;n;h];
subplot(3,1,3)
plot(t, gating);legend('m', 'n', 'h');
xlabel('time(ms)', 'fontsize', 20)


%method that represents expression 1 of report
function rate=alphah(v)
const = (v+70)/20;
rate=0.07*exp(-const);
end

%method that represents expression 2 of report
function rate=betah(v)
const = (v+40)/10;
rate= 1.0 ./ (1+exp(-const));
end

%method that represents expression 3 of report
function rate=alpham(v)
const=(v+45)/10;
rate=1.0*const ./ (1-exp(-const));
rate(isnan(rate))=1;
end

%method that represents expression 4 of report
function rate=betam(v)
const = (v+70)/18;
rate=4.0*exp(-const);
end

%method that represents expression 5 of report
function rate=alphan(v)
const=(v+60)/10;
rate=0.1*const ./ (1-exp(-const));
rate(isnan(rate))=0.1;
end

%method that represents expression 6 of report
function rate=betan(v)
const = (v+70)/80;
rate=0.125*exp(-const);
end