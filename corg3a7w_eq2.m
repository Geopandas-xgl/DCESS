
function [s,x,Rm,C1,D3,D2,D1,C2,C3,Res3]=...
    corg3a_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zo)

  s(1) = 0.5*(w(1)/Db + (w(1)^2/(Db^2)+4*jo2/Db)^0.5);
  s(2) = 0.5*(w(1)/Db - (w(1)^2/(Db^2)+4*jo2/Db)^0.5);
  s(3) = 0.5*(w(2)/Db + (w(2)^2/(Db^2)+4*jo2/Db)^0.5);
  s(4) = 0.5*(w(2)/Db - (w(2)^2/(Db^2)+4*jo2/Db)^0.5);
  s(5) = 0.5*(w(3)/Db + (w(3)^2/(Db^2)+4*jo2/Db)^0.5);
  s(6) = 0.5*(w(3)/Db - (w(3)^2/(Db^2)+4*jo2/Db)^0.5);
  s(7) = 0.5*(w(4)/Db + (w(4)^2/(Db^2)+4*jo2/Db)^0.5);
  s(8) = 0.5*(w(4)/Db - (w(4)^2/(Db^2)+4*jo2/Db)^0.5);
  s(9) = 0.5*(w(5)/Db + (w(5)^2/(Db^2)+4*jo2/Db)^0.5);
  s(10) = 0.5*(w(5)/Db - (w(5)^2/(Db^2)+4*jo2/Db)^0.5);

  s(11) = 0.5*(w(1)/Db + (w(1)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(12) = 0.5*(w(1)/Db - (w(1)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(13) = 0.5*(w(2)/Db + (w(2)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(14) = 0.5*(w(2)/Db - (w(2)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(15) = 0.5*(w(3)/Db + (w(3)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(16) = 0.5*(w(3)/Db - (w(3)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(17) = 0.5*(w(4)/Db + (w(4)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(18) = 0.5*(w(4)/Db - (w(4)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(19) = 0.5*(w(5)/Db + (w(5)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(20) = 0.5*(w(5)/Db - (w(5)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(21) = 0.5*(w(6)/Db + (w(6)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(22) = 0.5*(w(6)/Db - (w(6)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(23) = 0.5*(w(7)/Db + (w(7)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(24) = 0.5*(w(7)/Db - (w(7)^2/(Db^2)+4*beta*jo2/Db)^0.5); 
  
  c1=(w(1)-Db*s(1))*(1-fi(1));
  c2=(w(1)-Db*s(2))*(1-fi(1));
  c3=exp(s(1)*zb(2));
  c4=exp(s(2)*zb(2));
  c5=-exp(s(3)*zb(2));
  c6=-exp(s(4)*zb(2));
  c7=(w(1)-Db*s(1))*exp(s(1)*zb(2))*(1-fi(1));
  c8=(w(1)-Db*s(2))*exp(s(2)*zb(2))*(1-fi(1));
  c9=-(w(2)-Db*s(3))*exp(s(3)*zb(2))*(1-fi(2));
  c10=-(w(2)-Db*s(4))*exp(s(4)*zb(2))*(1-fi(2));
  c11=exp(s(3)*zb(3));
  c12=exp(s(4)*zb(3));
  c13=-exp(s(5)*zb(3));
  c14=-exp(s(6)*zb(3));
  c15=(w(2)-Db*s(3))*exp(s(3)*zb(3))*(1-fi(2));
  c16=(w(2)-Db*s(4))*exp(s(4)*zb(3))*(1-fi(2));
  c17=-(w(3)-Db*s(5))*exp(s(5)*zb(3))*(1-fi(3));
  c18=-(w(3)-Db*s(6))*exp(s(6)*zb(3))*(1-fi(3));
  c19=exp(s(5)*zo);
  c20=exp(s(6)*zo);
  c21=-exp(s(15)*zo);
  c22=-exp(s(16)*zo);
  c23=(w(3)-Db*s(5))*exp(s(5)*zo);
  c24=(w(3)-Db*s(6))*exp(s(6)*zo);
  c25=-(w(3)-Db*s(15))*exp(s(15)*zo);
  c26=-(w(3)-Db*s(16))*exp(s(16)*zo);
  c27=exp(s(15)*zb(4));
  c28=exp(s(16)*zb(4));
  c29=-exp(s(17)*zb(4));
  c30=-exp(s(18)*zb(4));
  c31=(w(3)-Db*s(15))*exp(s(15)*zb(4))*(1-fi(3));
  c32=(w(3)-Db*s(16))*exp(s(16)*zb(4))*(1-fi(3));
  c33=-(w(4)-Db*s(17))*exp(s(17)*zb(4))*(1-fi(4));
  c34=-(w(4)-Db*s(18))*exp(s(18)*zb(4))*(1-fi(4));
  c35=exp(s(17)*zb(5));
  c36=exp(s(18)*zb(5));
  c37=-exp(s(19)*zb(5));
  c38=-exp(s(20)*zb(5));
  c39=(w(4)-Db*s(17))*exp(s(17)*zb(5))*(1-fi(4));
  c40=(w(4)-Db*s(18))*exp(s(18)*zb(5))*(1-fi(4));
  c41=-(w(5)-Db*s(19))*exp(s(19)*zb(5))*(1-fi(5));
  c42=-(w(5)-Db*s(20))*exp(s(20)*zb(5))*(1-fi(5));
  c43=exp(s(19)*zb(6));
  c44=exp(s(20)*zb(6));
  c45=-exp(s(21)*zb(6));
  c46=-exp(s(22)*zb(6));
  c47=-(w(5)-Db*s(19))*exp(s(19)*zb(6))*(1-fi(5));
  c48=-(w(5)-Db*s(20))*exp(s(20)*zb(6))*(1-fi(5));
  c49=(w(6)-Db*s(21))*exp(s(21)*zb(6))*(1-fi(6));
  c50=(w(6)-Db*s(20))*exp(s(20)*zb(6))*(1-fi(6));
  c51=exp(s(21)*zb(7));
  c52=exp(s(22)*zb(7));
  c53=-exp(s(23)*zb(7));
  c54=-exp(s(24)*zb(7));
  c55=-(w(6)-Db*s(21))*exp(s(21)*zb(7))*(1-fi(6));
  c56=-(w(6)-Db*s(22))*exp(s(22)*zb(7))*(1-fi(6));
  c57=(w(7)-Db*s(23))*exp(s(23)*zb(7))*(1-fi(7));
  c58=(w(7)-Db*s(24))*exp(s(24)*zb(7))*(1-fi(7));
  c59=(Db*s(23))*exp(s(23)*zb(8));
  c60=(Db*s(24))*exp(s(24)*zb(8));
  
 
  
  A =[  c1  c2  0   0   0   0   0   0   0   0   0   0   0   0   0   0;...
        c3  c4  c5  c6  0   0   0   0   0   0   0   0   0   0   0   0;...
        c7  c8  c9  c10 0   0   0   0   0   0   0   0   0   0   0   0;...
        0   0   c11 c12 c13 c14 0   0   0   0   0   0   0   0   0   0;...
        0   0   c15 c16 c17 c18 0   0   0   0   0   0   0   0   0   0;...
        0   0   0   0   c19 c20 c21 c22 0   0   0   0   0   0   0   0;...
        0   0   0   0   c23 c24 c25 c26 0   0   0   0   0   0   0   0;...
        0   0   0   0   0   0   c27 c28 c29 c30 0   0   0   0   0   0;...
        0   0   0   0   0   0   c31 c32 c33 c34 0   0   0   0   0   0;...
        0   0   0   0   0   0   0   0   c35 c36 c37 c38 0   0   0   0;...
        0   0   0   0   0   0   0   0   c39 c40 c41 c42 0   0   0   0;...
        0   0   0   0   0   0   0   0   0   0   c43 c44 c45 c46 0   0;...
        0   0   0   0   0   0   0   0   0   0   c47 c48 c49 c50 0   0;...
        0   0   0   0   0   0   0   0   0   0   0   0   c51 c52 c53 c54;...
        0   0   0   0   0   0   0   0   0   0   0   0   c55 c56 c57 c58;...
        0   0   0   0   0   0   0   0   0   0   0   0   0   0   c59 c60];...
  
 B = [FC; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];         
           
 x = linsolve(A,B);
  
 Rm1=(1-fi(1))*jo2*(x(1)*(exp(s(1)*zb(2))-1)/s(1)+...
                    x(2)*(exp(s(2)*zb(2))-1)/s(2));
 Rm2=(1-fi(2))*jo2*(x(3)*(exp(s(3)*zb(3))-exp(s(3)*zb(2)))/s(3)+...
                    x(4)*(exp(s(4)*zb(3))-exp(s(4)*zb(2)))/s(4));
 Rm3=(1-fi(3))*jo2*(x(5)*(exp(s(5)*zo)-exp(s(5)*zb(3)))/s(5)+...
                    x(6)*(exp(s(6)*zo)-exp(s(6)*zb(3)))/s(6)+...
              beta*(x(7)*(exp(s(15)*zb(4))-exp(s(15)*zo))/s(15)+...
                    x(8)*(exp(s(16)*zb(4))-exp(s(16)*zo))/s(16)));
 Rm4=(1-fi(4))*beta*jo2*(x(9)*(exp(s(17)*zb(5))-exp(s(17)*zb(4)))/s(17)...
                        +x(10)*(exp(s(18)*zb(5))-exp(s(18)*zb(4)))/s(18));               
 Rm5=(1-fi(5))*beta*jo2*(x(11)*(exp(s(19)*zb(6))-exp(s(19)*zb(5)))/s(19)+...
                         x(12)*(exp(s(20)*zb(6))-exp(s(20)*zb(5)))/s(20));
 Rm6=(1-fi(6))*beta*jo2*(x(13)*(exp(s(21)*zb(7))-exp(s(21)*zb(6)))/s(21)+...
                         x(14)*(exp(s(22)*zb(7))-exp(s(22)*zb(6)))/s(22));
 Rm7=(1-fi(7))*beta*jo2*(x(15)*(exp(s(23)*zb(8))-exp(s(23)*zb(7)))/s(23)+...
                         x(16)*(exp(s(24)*zb(8))-exp(s(24)*zb(7)))/s(24));
 
Rm8=(1-fi(3))*jo2*(beta*(x(7)*(exp(s(15)*zb(4))-exp(s(15)*zo))/s(15)+...
                    x(8)*(exp(s(16)*zb(4))-exp(s(16)*zo))/s(16)))+Rm4+Rm5+Rm6+Rm7;
                     
 Rm = [Rm1 Rm2 Rm3 Rm4 Rm5 Rm6 Rm7 Rm8];
 
   
  C1 = IniO2 - del(1)*(x(1)/(s(1)^2)+x(2)/(s(2)^2));
  
  D3 = -del(3)*((x(5)/s(5)*exp(s(5)*zo)+x(6)/s(6)*exp(s(6)*zo)...
  + beta*(x(7)*(exp(s(15)*zb(4))-exp(s(15)*zo))/s(15)+x(8)*(exp(s(16)*zb(4))-exp(s(16)*zo))/s(16)...
    +del(4)/del(3)*(x(9)*(exp(s(17)*zb(5))-exp(s(17)*zb(4)))/s(17)+x(10)*(exp(s(18)*zb(5))-exp(s(18)*zb(4)))/s(18))...
    +del(5)/del(3)*(x(11)*(exp(s(19)*zb(6))-exp(s(19)*zb(5)))/s(19)+x(12)*(exp(s(20)*zb(6))-exp(s(20)*zb(5)))/s(20))...
    +del(6)/del(3)*(x(13)*(exp(s(21)*zb(7))-exp(s(21)*zb(6)))/s(21)+x(14)*(exp(s(22)*zb(7))-exp(s(22)*zb(6)))/s(22))...
    +del(7)/del(3)*(x(15)*(exp(s(23)*zb(8))-exp(s(23)*zb(7)))/s(23)+x(16)*(exp(s(24)*zb(8))-exp(s(24)*zb(7)))/s(24)))));
      
    %D3 = - del(3)*(x(1)/s(1)*(exp(s(1)*zo))+x(2)/s(2)*(exp(s(2)*zo)));
  
 D2 = fi(3)*F(2)/(fi(2)*F(3))*(D3 + del(3)*(x(5)*exp(s(5)*zb(3))/s(5)+x(6)*exp(s(6)*zb(3))/s(6)))-...
      del(2)*(x(3)*exp(s(3)*zb(3))/s(3)+x(4)*exp(s(4)*zb(3))/s(4));
 
 D1 = fi(2)*F(1)/(fi(1)*F(2))*(D2 + del(2)*(x(3)*exp(s(3)*zb(2))/s(3)+x(4)*exp(s(4)*zb(2))/s(4)))-...
      del(1)*(x(1)*exp(s(1)*zb(2))/s(1)+x(2)*exp(s(2)*zb(2))/s(2));     
 
  C2 = C1 - (D2-D1)*zb(2)-...
      del(2)*(x(3)*exp(s(3)*zb(2))/(s(3)^2)+x(4)*exp(s(4)*zb(2))/(s(4)^2))+...
      del(1)*(x(1)*exp(s(1)*zb(2))/(s(1)^2)+x(2)*exp(s(2)*zb(2))/(s(2)^2));    
 
  C3 = C2 - (D3-D2)*zb(3)-...
      del(3)*(x(5)*exp(s(5)*zb(3))/(s(5)^2)+x(6)*exp(s(6)*zb(3))/(s(6)^2))+...
      del(2)*(x(3)*exp(s(3)*zb(3))/(s(3)^2)+x(4)*exp(s(4)*zb(3))/(s(4)^2));

  Res3 = -(C3 + D3*zo + del(3)*(x(5)*exp(s(5)*zo)/(s(5)^2)+x(6)*exp(s(6)*zo)/(s(6)^2)));
  
 
end
