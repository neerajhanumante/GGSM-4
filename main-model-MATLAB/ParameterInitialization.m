% || Sarvam Krishna Arpanam ||

function [ic,ecolparams,econparams] = ParameterInitialization

i=1;


ecolparams=[2.861325308717139e-002,
     1.734351784507601e-002,
     3.265224224356716e-002,
     5.732733145333690e-003,
     3.131235080180339e-001,
     9.174906638699569e-001,
     3.327275637705977e-001,
     2.350696927025789e-001,
     5.081656269272850e-002,
     8.688268022462951e-001,
     4.532e-001,%4.932828648894071e-001,
     4.658138102587053e-001,
     1.000000000000000e-003,
     5.030915711147533e-001,
     2.302639355221487e-001,
     4.286472149788781e-001,
     0.0,
     1.233437625619434e+000,
     1.000000000000000e+001,
     0.26340461,%x(1),%0.1  gP1H2 0.082553050860415                   
     0.2];
 
	ecolparams(10)=0.9;
	ecolparams(7)=1.3127275637705977e-001;
 
	ecolparams(14)=4.9030915711147533e-001;
	ecolparams(3)=3.7e-2;
	ecolparams(9)=1.081656269272850e-002;
	ecolparams(2)=2.934351784507601e-002;	
	ecolparams(3)=4.2e-2; 





	sol=[4.832248668616482e-001,
     1.357181037387526e-001,
     1.000000000000000e+000, %ap1 x(36), %
     1.000000000000000e+000, %bp1 x(37), %
     7.737088487049635e-001, %cp1 x(38), %    
     0.207876237,%x(10),%1.910770291618713e-001, ap1p, dp1h1,dh1hh,dishh0.0559187825742761        
     0.349388112,%x(11),%4.991249741897157e-002, bp1p, ep1h1,  0.165478172346092        
     0.258916193,%x(12),%3.380538111199238e-002, cp1p      0.0375155728082948       
     7.524098562805767e-001,
     1.000000000000000e-003,
     2.527165128459011e-001,
     9.926234817539428e-001,
     6.081040091094971e-001,
     2.972103070158917e-001,
     1.000000000000000e-003,
     5.646817673172643e-001,
     0.825044346,%x(14),%4.240782789383344e-001, fp1h1, kh1hh, kp1hh, kishh   0.851716686922311 
     0.219180582, %x(2), %19.000000000000000e-001, x(2) % gP1H1   
     1.286507253686182e-001,
     1.474467159083802e-002,  
     3.995283985198757e-001,
     1.555235114499388e-001, %x(17), zp1hh, zh1hh, zishh, 
     1.019919611353713e-001, %x(18),%mh1hh, mp1hh, mishh
     6.766772330024032e-001, %x(19), %np1hh,nh1hh, nishh, theta
     0.000118394,%x(3),%mP1 %1.018295456166137e-003, mP1
     9.838862467656756e-003, %0.00108434660226991 mH1
     9.000000000000011e-001];

	econparams(1:5)=sol(1:5);
	econparams(6)=3*sol(6);
    econparams(7)=3*sol(7);
	econparams(8:11)=sol(8:11);
	econparams(12)=sol(6);
	econparams(13)=sol(7);
	econparams(14:17)=sol(12:15);
	econparams(18)=sol(6);
	econparams(19)=sol(7);
	econparams(20)=sol(16);
	dP1H1=0.269204722;%x(13);%sol(6);
	econparams(21)=dP1H1/1000; %dP1H1
	econparams(25)=dP1H1/1000;  
	econparams(29)=dP1H1/1000; 
	econparams(34)=dP1H1/1000;
	econparams(39)=dP1H1/1000;
	eP1H1=sol(7);
	econparams(22)=eP1H1; %ep1h1
	econparams(26)=eP1H1; %ep1is
	fP1H1=sol(17);
	econparams(23)=fP1H1/1;
	econparams(27)=fP1H1/1;
	econparams(31)=fP1H1/1000; %kp1hh
	econparams(36)=0.5*fP1H1/1000; %kh1hh
	econparams(41)=0.5*fP1H1/1000; %kishh
	econparams(24)=sol(18);
	econparams(28)=sol(19);
	econparams(30)=sol(20)/1000*10000;
	econparams(35)=sol(20)/1000*10000;
	econparams(40)=sol(20)/1000*10000;
	mH1HH=sol(21);
	econparams(32)=0.5*mH1HH/1000;
	econparams(37)=mH1HH/1000; %mh1hh
	econparams(42)=0.5*mH1HH/1000;
	nISHH=sol(22);
	econparams(33)=0.5*nISHH/1000;
	econparams(38)=0.5*nISHH/1000;
	econparams(43)=nISHH/1000;
    
	econparams(45)=sol(23);     %theta(mass of P1 necessary to produce one unite of IS product) 
    
	econparams(46)=sol(24);
	econparams(48)=sol(25);
	econparams(49)=sol(26);
	econparams(50)=0.22;%x(4); %mHH 
 
	econparams(44)=0.3;

	econparams(47)=0.000331678;%x(5); %gRPP1 .09; 0.00453370905224315
 
	econparams(51)=0; %x(39); %p1bar
	econparams(52)=0.4;
	econparams(53)=0;
	econparams(54)=0.01*4.507354330528614e-004;%0.01*4.507354330528614e-004;
 
	econparams(55)= 0.988429965; %x(6); %1.5; etaa  0.948266028552691
	econparams(56)=1.128803224;%x(7); %0.8333;  %etab 1.09153879779457
	econparams(57)=0.01666; %etac
	econparams(58)=10.00717122;%x(8);% 1E1;	%phi 10.0080145597677
 
 
%
% initialize economic state [P1;H1;IS;HH]
%
	ic_econ=[0.80213,  %1.550317582123858e+000,%
     0.25,%1.400236056583227e-001, %
     1.008970080362582e-001,
     4.507354330528614e-001];
                   
	numHH= 10; %16; %10 %initial number of humans
 
	econparams(59)=0.20645585; %x(9); %4.055586606327234e-003; 0.00451092487952919 
%
%	Calculate ecological equilibrium, without any of the domesticated or
%	industrial agents
	ecolparams(8)=ecolparams(8)+0.25*ecolparams(8);
    
 

ic_ecol = eq21apr2004V2(ecolparams); %ic_eco=[P2;P3;H2;H3;C1;C2;RP;IRP];

% 
%	ecolparams(18)=0  %set RPIRP to zero, without changing original equilibrium.  
%	Comment this out to run ecological model itself
%                   
%     econparams=zeros(59,1)			%uncomment these two lines to run the
%     econparams(45:46)=1			    %ecological model with no economics
 
%	P1, P2, P3, H1, H2, H3, C1, C2, HH, IS, RP, IRP, P1H1massdeficit, 
%	P1ISmassdeficit,P1HHmassdeficit, H1massdeficit, ISmassdeficit, numHH, percapmass
%	
	ic(1:12)=[ic_econ(1),ic_ecol(1),ic_ecol(2),ic_econ(2),ic_ecol(3),ic_ecol(4),ic_ecol(5),ic_ecol(6),ic_econ(4),ic_econ(3),ic_ecol(7),ic_ecol(8)];

    ic(13)=0.0;
	ic(14)=0;
	ic(15)=0;
	ic(16)=0;
	ic(17)=0;
	ic(18)=numHH;
	ic(19)=4.055586606327234e-003;
    ic(20)=800;