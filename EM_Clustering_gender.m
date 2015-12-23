%This is the algorithm for clustering height and weight data for male and
%female seperately
clear
close all
%read data
data=xlsread('final_data');
%distinguish male and female
X=data';
X0=X(:,X(3,:)==0);
X=X(:,X(3,:)==1);
X0=[X0(1,:)',X0(2,:)']';
X=[X(1,:)',X(2,:)']';
K=4;
[d,N]=size(X);
pi=zeros(1,K);
%male
Y=X';
phi=zeros(K,N);
[idx,mu] = kmeans(Y,K);
sigma = zeros(d,d,K);
time=300;
for k = 1:K
    sigma(:,:,k) = cov(Y(idx==k,:));
    pi(k)=sum(idx==k)/N;
end
mu=mu';
MAX=zeros(1,time);
%run iterations for 300 time
for T=1:time
for i=1:N
    for j=1:K
        phi(j,i)=pi(j)*mvnpdf(X(:,i),mu(:,j),sigma(:,:,j));
    end
end
A=sum(phi);
%likelyhood
L=0;
for l=1:N
    L=log(A(l))+L;
end
MAX(T)=L;
A=repmat(A,[K,1]);
phi=phi./A;
n=sum(phi,2);
for i=1:K
    mu(:,i)=(phi(i,:)*X')'./n(i);
    cbk=repmat(mu(:,i),[1,N]);
    sigma(:,:,i)=(X-cbk)*diag(phi(i,:))*(X-cbk)'./n(i);
end
pi=n./N;
end

%female
Y0=X0';
[d,N0]=size(X0);

pi0=zeros(1,K);
phi0=zeros(K,N0);
[idx,mu0] = kmeans(Y0,K);
sigma0 = zeros(d,d,K);
time=300;
for k = 1:K
    sigma0(:,:,k) = cov(Y0(idx==k,:));
    pi0(k)=sum(idx==k)/N0;
end
mu0=mu0';
MAX0=zeros(1,time);
for T=1:time
for i=1:N0
    for j=1:K
        phi0(j,i)=pi0(j)*mvnpdf(X0(:,i),mu0(:,j),sigma0(:,:,j));
    end
end
A0=sum(phi0);
%likelyhood
L0=0;
for l=1:N0
    L0=log(A0(l))+L0;
end
MAX0(T)=L0;
A0=repmat(A0,[K,1]);
phi0=phi0./A0;
n0=sum(phi0,2);
for i=1:K
    mu0(:,i)=(phi0(i,:)*X0')'./n0(i);
    cbk=repmat(mu0(:,i),[1,N0]);
    sigma0(:,:,i)=(X0-cbk)*diag(phi0(i,:))*(X0-cbk)'./n0(i);
end
pi0=n0./N0;
end


[a,index]=max(phi0);
col=['b','g','m','r','c','y','k']; 
shape=['s','o','d','p'];
for i=1:N0
    plot3(0,X0(1,i),X0(2,i),shape(index(i)),'Color',col(index(i)),'MarkerFaceColor',col(index(i)),'MarkerSize',12)

    hold on
end
[a,index]=max(phi);
col=['g','b','r','m','c','y','k']; 
shape=['o','s','p','d'];
for i=1:N
    plot3(1,X(1,i),X(2,i),shape(index(i)),'Color',col(index(i)),'MarkerFaceColor',col(index(i)),'MarkerSize',12)

    hold on
end
rotate3d on;
grid on;
for k=1:K
    plot3(0,mu0(1,k),mu0(2,k),'o','Color','k','MarkerSize',20,'LineWidth',3)
    hold on
end
for k=1:K
    plot3(1,mu(1,k),mu(2,k),'o','Color','k','MarkerSize',15,'LineWidth',3)
    hold on
end
