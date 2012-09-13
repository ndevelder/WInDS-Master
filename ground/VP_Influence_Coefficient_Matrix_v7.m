function [Gamma] = VP_Influence_Coefficient_Matrix_v7(V_ext,Xp,Yp,Zp,co)

Np = size(Xp,1)-1;
N = Np^2;
P = zeros(N,3);
dx = Xp(1,2,1)-Xp(1,1,1);
dy = Yp(2,1,1)-Yp(1,1,1);

i=1;

for j=1:Np
    for k=1:Np
        P(i,1) = Xp(j,k,1)+0.5*(Xp(j,k+1,1)-Xp(j,k,1));
        P(i,2) = Yp(j,k,1)+0.5*(Yp(j+1,1,1)-Yp(j,1,1));
        P(i,3) = Zp(1,1,k)+0.5*(Zp(1,1,k+1)-Zp(1,1,k));
        i=i+1;
    end
end

r1 = zeros(N,3);
r2 = zeros(N,3);
r3 = zeros(N,3);
r4 = zeros(N,3);

n1 = zeros(N,1);
n2 = zeros(N,1);
n3 = zeros(N,1);
n4 = zeros(N,1);

t1 = zeros(N,1);
t2 = zeros(N,1);
t3 = zeros(N,1);
t4 = zeros(N,1);

c1 = zeros(N,3);
c2 = zeros(N,3);
c3 = zeros(N,3);
c4 = zeros(N,3);

A1 = zeros(N,3);
A2 = zeros(N,3);
A3 = zeros(N,3);
A4 = zeros(N,3);

A = zeros(N,Np^2);
n_vec = zeros(N,3);
B = zeros(N,Np^2);
V_n = zeros(N,1);
Gamma = zeros(N,1);

x = P(:,1);   %x-coordinate of point in space
y = P(:,2);   %y-coordinate
z = P(:,3);   %z-coordinate

co = co*(dx+dy)/2;

j = 1;

for n = 1:Np  % count through y direction
    for m = 1:Np % count through x direction

        r1 = [x - Xp(n,m,1),y - Yp(n,m,1),z - Zp(1,1,m)];
        r2 = [x - Xp(n,m+1,1),y - Yp(n,m+1,1),z - Zp(1,1,m+1)];
        r3 = [x - Xp(n+1,m+1,1),y - Yp(n+1,m+1,1),z - Zp(1,1,m+1)];
        r4 = [x - Xp(n+1,m,1),y - Yp(n+1,m,1),z - Zp(1,1,m)];

        n1 = (sum((r1.^2),2)).^.5;
        n2 = (sum((r2.^2),2)).^.5;
        n3 = (sum((r3.^2),2)).^.5;
        n4 = (sum((r4.^2),2)).^.5;

        t1 = (n1+n2) ./ (n1.*n2.*(n1.*n2+sum(r1.*r2,2))+co);
        c1 = cross(r1,r2,2);
        A1 = c1.*[t1 t1 t1];

        t2 = (n2+n3) ./ (n2.*n3.*(n2.*n3+sum(r2.*r3,2))+co);
        c2 = cross(r2,r3,2);
        A2 = c2.*[t2 t2 t2];

        t3 = (n3+n4) ./ (n3.*n4.*(n3.*n4+sum(r3.*r4,2))+co);
        c3 = cross(r3,r4,2);
        A3 = c3.*[t3 t3 t3];

        t4 = (n4+n1) ./ (n4.*n1.*(n4.*n1+sum(r4.*r1,2))+co);
        c4 = cross(r4,r1,2);
        A4 = c4.*[t4 t4 t4];
                                                
        n_vec(j,:) = c2(1,:)/norm(c2(1,:));
        V_n(j,1) = -dot(V_ext(j,:),n_vec(j,:));
        
        A(:,j) = (A1+A2+A3+A4)/4/pi*[n_vec(j,:)]';

        j=j+1;

    end

end

B = inv(A);
Gamma = B*V_n;