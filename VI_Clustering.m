clear;
close all

data=xlsread('final_data_2');
X=data';
[d,N]=size(X);
K=4;
alpha=ones(1,K);
%m,sigma
 Y=X';
[idx,mu] = kmeans(Y,K);
sigma = zeros(d,d,K);
for k = 1:K
    sigma(:,:,k) = cov(Y(idx==k,:));
end
mu=mu';
T=500;
%a,B
a=d.*ones(1,K);
B=zeros(d,d,K);
for k=1:K
B(:,:,k)=0.2.*cov(Y);
end

t1=zeros(1,K);
t2=zeros(1,K);
t3=zeros(1,K);
t4=zeros(1,K);
phi=zeros(N,K);
L=zeros(1,T);
for t=1:T
%q(c)
for j=1:K
for i=1:N
t1(j)=psi(a(j)/2)+psi((a(j)-1)/2)-log(det(B(:,:,j)));
t2(j)=(X(:,i)-mu(:,j))'*(a(j).*pinv(B(:,:,j)))*(X(:,i)-mu(:,j));
t3(j)=trace(a(j).*pinv(B(:,:,j))*sigma(:,:,j));
t4(j)=psi(alpha(j))-psi(sum(alpha));
phi(i,j)=exp(0.5*t1(j)-0.5*t2(j)-0.5*t3(j)+t4(j));
end
end
phi=phi./repmat(sum(phi,2),[1,K]);
%(1) and (2) and (3) and(5)

n=sum(phi);
for j=1:K
    %q(pi)
    alpha(j)=1+n(j);
    %q(mu)
    sigma(:,:,j)=pinv(0.1.*eye(d)+n(j)*a(j).*pinv(B(:,:,j)));
    mu(:,j)=sigma(:,:,j)*(a(j)*(B(:,:,j)\eye(d))*(X*phi(:,j)));
    %q(sigma)
    a(j)=d+n(j);
    cbk=repmat(mu(:,j),[1,N]);
    zzb=zeros(d,d);%B
    for i=1:N
        zzb=zzb+phi(i,j).*sigma(:,:,j);
    end
    B(:,:,j)=0.2.*cov(Y)+(X-cbk)*diag(phi(:,j))*(X-cbk)'+zzb;
end
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

