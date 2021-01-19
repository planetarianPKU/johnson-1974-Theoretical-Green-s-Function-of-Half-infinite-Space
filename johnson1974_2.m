clear all
%(18)
x1=10000;x2=0;x3=0;
x3_s=0.2*1e3;
R=sqrt(x1^2+x2^2);
phi=asind(x2/R);
r=sqrt(R^2+x3_s^2);
theta=asind(R/r);
alpha=8000;beta=4620;
rho=3300;
mu=rho*beta^2;
p=0;

Gp=0;Gs=0;
dp=1e-8;
F=[1;0;0];
Gp_time=[];
T_length=4;
dt1=(r/alpha)/1000;
for t = 0:dt1:T_length
    %P wave hasn't come
    if (t < r/alpha)
        Gp=zeros(9,1);
    %P wave has come
    else
        %P wave
        for p = dp:dp:sqrt((t/r)^2-alpha^-2)-dp
            q1=-t/r*sind(theta)+1i*sqrt(((t/r)^2-alpha^-2-p^2))*cosd(theta);           
            %(20)--(23)
            %Real part should > 0
            eta_a1=sqrt(alpha^-2+p^2-q1^2);
            eta_b1=sqrt(beta^-2+p^2-q1^2);           
            %(29)--(32)
            gamma1=eta_b1^2+p^2-q1^2;
            sigma1=gamma1^2+4*eta_a1*eta_b1*(q1^2-p^2);         
            %(33)
            M11=2*eta_b1* ((q1^2+p^2)*cosd(phi)-p^2);
            M12=2*eta_b1* (q1^2+p^2)*sind(phi)*cosd(phi);
            M13=2*q1*eta_a1*eta_b1*cosd(phi);
            M21=M12;
            M22=2*eta_b1*( (q1^2+p^2)*(sind(phi))^2 -p^2);
            M23=2*q1*eta_a1*eta_b1*sind(phi);
            M31=q1*gamma1*cosd(phi);
            M32=q1*gamma1*sind(phi);
            M33=eta_a1*gamma1;
            
            M_matrix=[M11;M12;M13;M21;M22;M23;M31;M32;M33];
            Gp=Gp+ heaviside(t-r/alpha)*real(eta_a1*sigma1^-1*((t/r)^2-alpha^-2-p^2)^-0.5*M_matrix)*dp;
            %heaviside(t-r/alpha)*
        end
    end
    Gp_time=[Gp_time Gp];
end

Gs_time=[];
%S wave
if sind(theta) < beta/alpha
    t2=r/beta;
else
    t2 = r/alpha*sind(theta)+r*(beta^-2-alpha^-2)^0.5*cosd(theta);
end
dt2=t2/1000;
for t = 0:dt2:T_length
    if sind(theta) < beta/alpha
        p2=((t/r)^2-beta^-2)^0.5;
    else
        p2 = (((t/r-(beta^-2-alpha^-2)^0.5*cosd(theta))/sind(theta))^2-alpha^-2)^0.5;
    end
    %S wave hasn't come
    if (t < t2)
        Gs=zeros(9,1);
    %S wave come after t2
    else
        for p = dp:dp:p2-dp
            q2=-t/r*sind(theta)+1i*sqrt(((t/r)^2-beta^-2-p^2))*cosd(theta);
            %(20)--(23)
            %Real part should > 0
            eta_a2=sqrt(alpha^-2+p^2-q2^2);
            eta_b2=sqrt(beta^-2+p^2-q2^2);
            %(29)--(32)
            gamma2=eta_b2^2+p^2-q2^2;
            sigma2=gamma2^2+4*eta_a2*eta_b2*(q2^2-p^2);
            
            %(34)
            N11=eta_b2^-1 * (eta_b2^2*gamma2- (gamma2-4*eta_a2*eta_b2)*((q2^2+p^2)*(sind(phi))^2 - p^2));
            N12=eta_b2^-1*(q2^2+p^2)*(gamma2-4*eta_a2*eta_b2)*sind(phi)*cosd(phi);
            N13=-q2*gamma2*cosd(phi);
            N21=N12;
            N22=eta_b2^-1*(eta_b2^2*gamma2-(gamma2-4*eta_a2*eta_b2)*((q2^2+p^2)*(cosd(phi))^2-p^2));
            N23=-q2 * gamma2 * sind(phi);
            N31=-2 * q2 * eta_a2 * eta_b2 * cosd(phi);
            N32=-2*q2*eta_a2*eta_b2*sind(phi);
            N33=2*eta_a2*(q2^2-p^2);
            
            N_matrix=[N11;N12;N13;N21;N22;N23;N31;N32;N33];
            
            Gs=Gs+heaviside(t-t2)*real(eta_b2*sigma2^-1*((t/r)^2-beta^-2-p^2)^-0.5*N_matrix)*dp;
            %heaviside(t-t2)*
        end
    end
    Gs_time=[Gs_time Gs];
end

Gp1=1/(pi^2*mu*r)*diff(Gp_time(1,:));
Gp2=1/(pi^2*mu*r)*diff(Gp_time(2,:));
Gp3=1/(pi^2*mu*r)*diff(Gp_time(3,:));
Gp4=1/(pi^2*mu*r)*diff(Gp_time(4,:));
Gp5=1/(pi^2*mu*r)*diff(Gp_time(5,:));
Gp6=1/(pi^2*mu*r)*diff(Gp_time(6,:));
Gp7=1/(pi^2*mu*r)*diff(Gp_time(7,:));
Gp8=1/(pi^2*mu*r)*diff(Gp_time(8,:));
Gp9=1/(pi^2*mu*r)*diff(Gp_time(9,:));

Gs1=1/(pi^2*mu*r)*diff(Gs_time(1,:));
Gs2=1/(pi^2*mu*r)*diff(Gs_time(2,:));
Gs3=1/(pi^2*mu*r)*diff(Gs_time(3,:));
Gs4=1/(pi^2*mu*r)*diff(Gs_time(4,:));
Gs5=1/(pi^2*mu*r)*diff(Gs_time(5,:));
Gs6=1/(pi^2*mu*r)*diff(Gs_time(6,:));
Gs7=1/(pi^2*mu*r)*diff(Gs_time(7,:));
Gs8=1/(pi^2*mu*r)*diff(Gs_time(8,:));
Gs9=1/(pi^2*mu*r)*diff(Gs_time(9,:));
%
newdt=1e-3;
Gp1=interp1(0:dt1:T_length-dt1,Gp1,0:newdt:T_length-newdt);
Gp2=interp1(0:dt1:T_length-dt1,Gp2,0:newdt:T_length-newdt);
Gp3=interp1(0:dt1:T_length-dt1,Gp3,0:newdt:T_length-newdt);
Gp4=interp1(0:dt1:T_length-dt1,Gp4,0:newdt:T_length-newdt);
Gp5=interp1(0:dt1:T_length-dt1,Gp5,0:newdt:T_length-newdt);
Gp6=interp1(0:dt1:T_length-dt1,Gp6,0:newdt:T_length-newdt);
Gp7=interp1(0:dt1:T_length-dt1,Gp7,0:newdt:T_length-newdt);
Gp8=interp1(0:dt1:T_length-dt1,Gp8,0:newdt:T_length-newdt);
Gp9=interp1(0:dt1:T_length-dt1,Gp9,0:newdt:T_length-newdt);


Gs1=interp1(0:dt2:T_length-dt2,Gs1,0:newdt:T_length-newdt);
Gs2=interp1(0:dt2:T_length-dt2,Gs2,0:newdt:T_length-newdt);
Gs3=interp1(0:dt2:T_length-dt2,Gs3,0:newdt:T_length-newdt);
Gs4=interp1(0:dt2:T_length-dt2,Gs4,0:newdt:T_length-newdt);
Gs5=interp1(0:dt2:T_length-dt2,Gs5,0:newdt:T_length-newdt);
Gs6=interp1(0:dt2:T_length-dt2,Gs6,0:newdt:T_length-newdt);
Gs7=interp1(0:dt2:T_length-dt2,Gs7,0:newdt:T_length-newdt);
Gs8=interp1(0:dt2:T_length-dt2,Gs8,0:newdt:T_length-newdt);
Gs9=interp1(0:dt2:T_length-dt2,Gs9,0:newdt:T_length-newdt);

Gp_final=[Gp1;Gp2;Gp3;Gp4;Gp5;Gp6;Gp7;Gp8;Gp9];
Gs_final=[Gs1;Gs2;Gs3;Gs4;Gs5;Gs6;Gs7;Gs8;Gs9];

G_final=Gp_final+Gs_final;
for i = 1:9
    T=0:newdt:T_length-newdt;
    if max(abs(G_final(i,:)))>0
    x=G_final(i,:)/max(abs(G_final(i,:)));
    else
        x=0*G_final(i,:);
    end
plot(T,x-1.2*i,'LineWidth',2,'color','black');hold on;
end
legend('1','2','3','4','5','6','7','8','9')
