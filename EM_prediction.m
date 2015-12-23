Xtrain=csvread('final_data_Xtrain');
Ytrain=csvread('final_data_Ytrain');
[N,d]=size(Xtrain);
N0=length(find(Ytrain==0));
X0=Xtrain(:,1:N0);
X1=Xtrain(:,N0:N);
w=zeros(1,d)';
sig=1.5;
lam=1;
E=zeros (1,N);
xE=zeros (d,1);
xxT=zeros(d,d);
wf=zeros(d,N);
AB=zeros(1,N);
for i=1:N
    xxT=Xtrain(:,i)*Xtrain(:,i)'+xxT;
end
ww=zeros(1,N);
I=eye(d);
T=100;
lp=zeros(T,1);
for j=1:T
%E-step
a0= X0' * w ./1.5;
a1= X1' * w ./1.5;
Q0=normcdf(-a0);
M0=normpdf(-a0);
Q1=normcdf(-a1);%
M1=normpdf(-a1);%

E1=a1.*1.5+1.5.*M1./(1-Q1);
E0=a0.*1.5-1.5.*M0./Q0;
E=[E0' E1']';
xE=Xtrain*E;
w=inv(I+xxT./(1.5^2))*xE/(1.5^2);
aa0= X0' * w ./1.5;
aa1= X1' * w ./1.5;
QQ0=normcdf(aa0);
QQ1=normcdf(aa1);
lp(j)=d/2*log(1/(2*pi))-0.5*(w')*w+sum(log(QQ1))+sum(log(1-QQ0));
AB(j)=sum(log(QQ1));
wf(:,j)=w;
end


Xtest=csvread('final_data_Xtest');
Ytest=csvread('final_data_Ytest');
xtw=Xtest'*w/1.5;
p0=1-normcdf((xtw)./1.5);
p1=normcdf((xtw)./1.5);
% 
[Nt,dt]=size(Xtest);
result=ones(Nt,1);
k=1
A=0;%p0 t0
B=0;%p1 t0
C=0;%p1 t1
D=0;%p0 t1
for i=1:Nt
    if p0(i)>p1(i) 
        result(i)=0;
        if ytest(i)==0 A=A+1;
        else
            D=D+1;
           
            k=k+1;
        end
    else
        result(1)=1;
        if ytest(i)==0 
            B=B+1;
            
            k=k+1;
        else C=C+1;
        end
    end
    
end











        

