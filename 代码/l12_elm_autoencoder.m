function x1 = l12_elm_autoencoder(A,b,lam,itrs)
AA = (A') * A;
Lf = max(eig(AA));
Li = 1/Lf;
m = size(A,2);
n = size(b,2);
x = zeros(m,n);
value=((54^(1/3))/4)*((lam*Li)^(2/3));
x1=x;
L1 = Li * AA;
L2 = Li * A' * b;
for i=1:itrs
    g2=x1+L2-L1*x1;
    Phi=(lam/8)*((abs(g2)/3).^(-3/2));
    Phi=mapminmax(Phi',0,1)';
    Phi(Phi(:,:)==inf)=1;
    x1=((2/3)*g2).*(1+cos((2*pi/3)-((2/3)*(acos(Phi)))));
    index=sign(max(abs(g2)-value,0));
    x1=x1.*index;
end

