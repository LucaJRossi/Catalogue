clear all
close all
whitebg('w')


% Plot the value of the tidal radius of a clustuster with mass M_cl along
% the orbit. Note that no mass loss is included.

M_cl= 6.4e+03; % Select the mass of th cluster

coeff_x = load('ax_coeff.dat'); 
coeff_y = load('ay_coeff.dat'); 
coeff_z = load('az_coeff.dat'); 


for j = 1:length(coeff_x)
a200 = coeff_x(j,2);
a110 = coeff_x(j,3);
a101 = coeff_x(j,4);
a100 = coeff_x(j,5);
a020 = coeff_x(j,6);
a011 = coeff_x(j,7);
a010 = coeff_x(j,8);
a002 = coeff_x(j,9);
a001 = coeff_x(j,10);
a000 = coeff_x(j,11);

b200 = coeff_y(j,2);
b110 = coeff_y(j,3);
b101 = coeff_y(j,4);
b100 = coeff_y(j,5);
b020 = coeff_y(j,6);
b011 = coeff_y(j,7);
b010 = coeff_y(j,8);
b002 = coeff_y(j,9);
b001 = coeff_y(j,10);
b000 = coeff_y(j,11);

c200 = coeff_z(j,2);
c110 = coeff_z(j,3);
c101 = coeff_z(j,4);
c100 = coeff_z(j,5);
c020 = coeff_z(j,6);
c011 = coeff_z(j,7);
c010 = coeff_z(j,8);
c002 = coeff_z(j,9);
c001 = coeff_z(j,10);
c000 = coeff_z(j,11);


t = coeff_x(j,1);


X = coeff_x(j,12);
Y = coeff_x(j,13);
Z = coeff_x(j,14);

R = sqrt(X^2 + Y^2 + Z^2);               
theta = atan2(sqrt(X^2 + Y^2),Z);
phi = atan2(Y,X);

a_x = a200*X.^2 + a110*X.*Y + a101*X.*Z + a100*X + ...
      a020*Y.^2 + a011*Y.*Z + a010*Y + ...
      a002*Z.^2 + a001*Z +a000;
  
a_y = b200*X.^2 + b110*X.*Y + b101*X.*Z + b100*X + ...
      b020*Y.^2 + b011*Y.*Z + b010*Y + ...
      b002*Z.^2 + b001*Z +b000;  

a_z = c200*X.^2 + c110*X.*Y + c101*X.*Z + c100*X + ...
      c020*Y.^2 + c011*Y.*Z + c010*Y + ...
      c002*Z.^2 + c001*Z +c000;  
  
PHI = atan2(Y,X) + pi;
GAMMA = atan2(a_y,a_x);
ALPHA = GAMMA - PHI;
BETA = Z/R;

a_gal_cl = sqrt(a_x.^2 + a_y.^2 + a_z.^2)*cos(ALPHA)*cos(BETA);


d = linspace(0,300e-3,1e4);
G = 4.302e-6;   % kpc/Mo*(km/s)^2 

a_cl = G*M_cl./d.^2;

%semilogy(d,a_cl,'k')
%axis square
    k = 1;
    r = R-d(k);
    X = r*sin(theta)*cos(phi);
    Y = r*sin(theta)*sin(phi);
    Z = r*cos(theta);
    
    a_x = a200*X.^2 + a110*X.*Y + a101*X.*Z + a100*X + ...
          a020*Y.^2 + a011*Y.*Z + a010*Y + ...
          a002*Z.^2 + a001*Z +a000;
  
    a_y = b200*X.^2 + b110*X.*Y + b101*X.*Z + b100*X + ...
          b020*Y.^2 + b011*Y.*Z + b010*Y + ...
          b002*Z.^2 + b001*Z +b000;  
  
    a_z = c200*X.^2 + c110*X.*Y + c101*X.*Z + c100*X + ...
          c020*Y.^2 + c011*Y.*Z + c010*Y + ...
          c002*Z.^2 + c001*Z +c000;  
    
    PHI = atan2(Y,X) + pi;
    GAMMA = atan2(a_y,a_x);
    ALPHA = GAMMA - PHI;
    BETA = Z/r;
    
    a_gal_star = (sqrt(a_x.^2 + a_y.^2 + a_z.^2))*cos(ALPHA)*cos(BETA);% - G*M_cl/d(i)^2 - a_gal_cl);
   
    r = R-d(k+1);
    X = r*sin(theta)*cos(phi);
    Y = r*sin(theta)*sin(phi);
    Z = r*cos(theta);

    a_x = a200*X.^2 + a110*X.*Y + a101*X.*Z + a100*X + ...
          a020*Y.^2 + a011*Y.*Z + a010*Y + ...
          a002*Z.^2 + a001*Z +a000;
  
    a_y = b200*X.^2 + b110*X.*Y + b101*X.*Z + b100*X + ...
          b020*Y.^2 + b011*Y.*Z + b010*Y + ...
          b002*Z.^2 + b001*Z +b000;  
  
    a_z = c200*X.^2 + c110*X.*Y + c101*X.*Z + c100*X + ...
          c020*Y.^2 + c011*Y.*Z + c010*Y + ...
          c002*Z.^2 + c001*Z +c000;  
      
    PHI = atan2(Y,X) + pi;
    GAMMA = atan2(a_y,a_x);
    ALPHA = GAMMA - PHI;
    BETA = Z/r;
    
    a_gal_star_plus = (sqrt(a_x.^2 + a_y.^2 + a_z.^2))*cos(ALPHA)*cos(BETA);
        
    while  abs(a_gal_star - a_gal_cl) < G*M_cl/d(k)^2 && k<length(d)-2 %&&  abs(a_gal_star_plus - a_gal_cl) > G*M_cl/d(k+1)^2 
            k = k+1;
            r = R-d(k);
            X = r*sin(theta)*cos(phi);
            Y = r*sin(theta)*sin(phi);
            Z = r*cos(theta);
    
            a_x = a200*X.^2 + a110*X.*Y + a101*X.*Z + a100*X + ...
                a020*Y.^2 + a011*Y.*Z + a010*Y + ...
                a002*Z.^2 + a001*Z +a000;
  
            a_y = b200*X.^2 + b110*X.*Y + b101*X.*Z + b100*X + ...
                b020*Y.^2 + b011*Y.*Z + b010*Y + ...
                b002*Z.^2 + b001*Z +b000;  
  
            a_z = c200*X.^2 + c110*X.*Y + c101*X.*Z + c100*X + ...
                c020*Y.^2 + c011*Y.*Z + c010*Y + ...
                c002*Z.^2 + c001*Z +c000;  
    
            PHI = atan2(Y,X) + pi;
            GAMMA = atan2(a_y,a_x);
            ALPHA = GAMMA - PHI;
            BETA = Z/r;
    
            a_gal_star = (sqrt(a_x.^2 + a_y.^2 + a_z.^2))*cos(ALPHA)*cos(BETA);% - G*M_cl/d(i)^2 - a_gal_cl);
   
            r = R-d(k+1);
            X = r*sin(theta)*cos(phi);
            Y = r*sin(theta)*sin(phi);
            Z = r*cos(theta);

            a_x = a200*X.^2 + a110*X.*Y + a101*X.*Z + a100*X + ...
                a020*Y.^2 + a011*Y.*Z + a010*Y + ...
                a002*Z.^2 + a001*Z +a000;
  
            a_y = b200*X.^2 + b110*X.*Y + b101*X.*Z + b100*X + ...
                b020*Y.^2 + b011*Y.*Z + b010*Y + ...
                b002*Z.^2 + b001*Z +b000;  
  
            a_z = c200*X.^2 + c110*X.*Y + c101*X.*Z + c100*X + ...
                c020*Y.^2 + c011*Y.*Z + c010*Y + ...
                c002*Z.^2 + c001*Z +c000;  
      
            PHI = atan2(Y,X) + pi;
            GAMMA = atan2(a_y,a_x);
            ALPHA = GAMMA - PHI;
            BETA = Z/r;
    
            a_gal_star_plus = (sqrt(a_x.^2 + a_y.^2 + a_z.^2))*cos(ALPHA)*cos(BETA);
            
    end     
            
            r_t(j) = d(k)*1000;
            RR(j) = R;

end

plot(RR,r_t,'.k')


axis square
xlabel('d (kpc)','fontsize',12)
ylabel('r_t (pc)','fontsize',12)
set(gca,'fontsize',10)
axis([0 25 0 120])