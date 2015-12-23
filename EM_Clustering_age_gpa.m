%3-dimensional height-weight-gpa/age clustering
clear
close all
data=xlsread('final_data_2');
X=data';
K=4;
[d,N]=size(X);
pi=zeros(1,K);
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

[a,index]=max(phi);
col=['g','b','r','m','c','y','k']; 
shape=['^','s','p','d'];
for i=1:N
    plot3(X(3,i),X(1,i),X(2,i),shape(index(i)),'Color',col(index(i)),'MarkerFaceColor',col(index(i)),'MarkerSize',12)
    
    hold on
end
rotate3d on;
grid on;
for k=1:K
    plot3(mu(3,k),mu(1,k),mu(2,k),'o','Color','k','MarkerSize',15,'LineWidth',3)
    hold on
end 
