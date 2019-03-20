%������ԭʼ�Ķ�ά�����ͽ��ƺ�tau_star��lambertw��ʾ��
% Pmd��gamma��ʾ�ķ�������ȫ�����ţ�������ȷ��

% ��ά�������ŵ�np��alpha
%Pfa��Pmd����PDF�ۻ���ã�tau_star��һ��forѭ���õ���ʾ
clear all
vw=5;vb=1;%ע�⣬��Щ���Ǳ�׼��
n=100;
rho=1;lamda_ab=1;
epsilon=0.1;
alpha=linspace(0.01,0.999,20);
 np=linspace(1,99,20);
% np=linspace(3,99,33);
tau=linspace(0.1,3,20);%tauҪ����rhod+vw^2�����ܱ�֤Pmd������
nd=n-np;
np_op=zeros(1,length(alpha));
capacity_pair=zeros(1,length(alpha));
Pmd_star=zeros(length(alpha),length(np));
Pfa_star=zeros(length(alpha),length(np));
xi_star=zeros(length(alpha),length(np));
xi_star_test=zeros(1,length(np));
f_md=zeros(1,length(tau));
f_fa=zeros(1,length(tau));

for i=1:length(alpha)
    i
    %����ÿ��alpha������һ��np��һ��kexi���ҳ�xi>=1-epsilon����Щnp������Щλ�����ҳ�ʹ��
    %�ŵ���������np

%��ֵ�ó�tau_star,����ÿ��np������һ��tau�õ�һ����Ӧ��xi��
%�ҵ�min xi��tau_star,�ó���Ӧ��xi_star
for j=1:length(np)
    j
    rhop1=(1-alpha(i)).*rho.*n./np(j);
    rhod1=alpha(i).*rho.*n./nd(j);
    nd(j)=n-np(j);
    for m=1:length(tau)
    
    a=n*tau(m).*(rhop1-rhod1)/((rhod1+vw^2).*(rhop1+vw^2));%�������κ�����Ĳ���
    f_md(m)=(n^n*tau(m).^(n-1)*exp(-n*tau(m)./(rhod1+vw^2))./((rhop1+vw^2).^np(j)))./((rhod1+vw^2).^nd(j)).*...
        hypergeom(np(j),n,a)/gamma(n);%T��PDF
    f_fa(m)=(n/vw^2).^n./gamma(n).*tau(m).^(n-1).*exp(-n.*tau(m)./vw^2);
    end
    %��PDF��õ�Pmd ��Pfa
P_md=cumsum( f_md)/sum( f_md);  %�����ֲ���Pmd
P_fa=1-cumsum( f_fa)/sum( f_fa);  %�����ֲ���Pfa
xi=P_md+P_fa;

    loca=(xi==min(xi));
    %�ҵ�ʹ��loca��һ�γ���1��λ��
    location=find(loca,1);
    tau_star=tau(location);
%     xi_star_test(j)=xi(location);
    xi_star(i,j)=xi(location);
    
end
%�ҵ�����������xi_star��λ��
plot(np,xi_star(i,:));
pos=(xi_star(i,:)>=1-epsilon);
xi_star_ok=xi_star(i,pos);
% pos=(xi_star(i,:)==xi_star_ok);
%��û�����������ģ�������np��capacityΪ0
if any(pos)==false
    np_op=0;
    capacity_pair(i)=0;
    continue
end
%�ѿ��е�np��λ���ҳ���
np_ok=np(pos);

% �������������rho_eff������
rho_eff=(lamda_ab^2.*rho.*n.*alpha(i).*(1-alpha(i)))./...
    ((n-np_ok).*vb^2.*(lamda_ab.*(1-alpha(i))+vb^2./rho./n+lamda_ab.*alpha(i)./(n-np_ok)));
capacity=(n-np_ok)./n.*log2(1+rho_eff); 
%�ҵ�ʹ����������np��λ��
pos2=(capacity==max(capacity));
np_op(i)=np_ok(pos2);

%��ÿ��alpha��Ӧ������np����һ�ԣ������Ӧ��capacity
capacity_pair(i)=max(capacity); %�ŵ�����
end
pos3=(capacity_pair==max(capacity_pair));
pos3=find(pos3,1);
alpha_op_whole=alpha(pos3);
np_op_whole=np_op(pos3);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��ά�������ŵ�np��alpha
%tau_star��lambertw������ʾ��Pmd��igamma��ʾ

vw=5;vb=1;n=100;
rho=1;lamda_ab=1;
epsilon=0.1;
alpha=linspace(0.001,0.999,20);
np=linspace(1,99,20);
mut=rho+vw^2;
nd=n-np;
np_op_orig=zeros(1,length(alpha));
rho_eff_pair=zeros(1,length(alpha));
capacity_pair_orig=zeros(1,length(alpha));
tau_star=zeros(1,length(np));
sigma_T_squ=zeros(length(alpha),length(np));
rho_eff_max=zeros(length(alpha),length(np));
for i=1:length(alpha)
    i
    rhop=(1-alpha(i)).*rho.*n./np;
rhod=alpha(i).*rho.*n./nd;
    %����ÿ��alpha������һ��np��һ��kexi���ҳ�xi>=1-epsilon����Щnp������Щλ�����ҳ�ʹ��
    %�ŵ���������np
       sigma_T_squ(i,:)=(1/n^2).*(np.*((1-alpha(i)).*rho.*n./np+vw^2).^2+...
    nd.*(alpha(i).*rho.*n./nd+vw^2).^2);
%����tau_star�ı��ʽ
A_num=-rhop.*exp((rhod+vw^2)./(rhop+vw^2)).*(gamma(n)./gamma(np)).^(1/n)./(n.*(rhop+vw^2));
for j=1:length(A_num)
    if A_num(j)<-0.3679
        A_num(j)=-1/exp(1);
    end
end
tau_star=-vw^2.*(rhop+vw^2)./rhop.*lambertw(-1,A_num);
%    tau_star(i,:)=(vw^2.*(n.*sigma_T_squ(i,:)-mut.*vw^2))./(n.*sigma_T_squ(i,:)-vw^4)+...
%     (sqrt(sigma_T_squ(i,:).*vw^4.*(n.*rho^2+(vw^4-n.*sigma_T_squ(i,:))...
%     .*log(vw^4./n./sigma_T_squ))))./(n.*sigma_T_squ-vw^4);  
% Pmd_satr(i,:)=0.5.*(erf(mut./sqrt(2.*sigma_T_squ(i,:)))-erf((mut-tau_star)./sqrt(2.*sigma_T_squ(i,:))));
Pmd_star=1-igamma(np,n.*(tau_star-rhod-vw^2)./(rhop+vw^2))./gamma(np);
for m=1:length(Pmd_star)
    if Pmd_star(m)<0
        Pmd_star(m)=0;
    elseif Pmd_star(m)>1
        Pmd_star(m)=1;
    end
end
Pfa_star=qfunc((tau_star-vw^2)./(vw^2./sqrt(n)));
xi_star=Pmd_star+Pfa_star;
% %ֱ�Ӽ���xi_star=1-epsilon��np������ȫ�����������ǲ��Ǹý���
% fun=@(np)xi_star+epsilon-1;
% op_op2=solve(fun,np);
% op_op2=double(op_op2);
plot(np,xi_star);
%ylim([0.99 1.01])
hold on; plot([1 10],[1-epsilon 1-epsilon]);
pos=(xi_star>=1-epsilon);
if any(pos)==false
    np_op_orig(i)=0;
    capacity_pair_orig(i)=0;
    continue
end
xi_star_ok=(xi_star(pos));
%�ѿ��е�np��λ���ҳ���
np_ok=np(pos);

% �������е�rho_eff������
rho_eff=(lamda_ab^2.*rho.*n.*alpha(i).*(1-alpha(i)))./...
    ((n-np).*vb^2.*(lamda_ab.*(1-alpha(i))+vb^2./rho./n+lamda_ab.*alpha(i)./(n-np)));
capacity=(n-np)./n.*log2(1+rho_eff); 
%�ҳ�����������rho_eff������
rho_eff_ok=rho_eff(pos);
capacity_ok=capacity(pos);
%�ҵ�ʹ����������np��λ��
pos2=(capacity==max(capacity_ok));
np_op_orig(i)=np(pos2);

%��ÿ��alpha��Ӧ������np����һ�ԣ������Ӧ��capacity
rho_eff_pair(i)=(lamda_ab^2.*rho.*n.*alpha(i).*(1-alpha(i)))./...
    ((n-np_op_orig(i)).*vb^2.*(lamda_ab.*(1-alpha(i))+vb^2./rho./n+lamda_ab.*alpha(i)./(n-np_op_orig(i))));
capacity_pair_orig(i)=(n-np_op_orig(i))./n.*log2(1+rho_eff_pair(i)); %�ŵ�����
end
pos3=(capacity_pair_orig==max(capacity_pair_orig));
alpha_op_whole_orig=alpha(pos3);
np_op_whole_orig=np_op_orig(pos3);





plot(np_op,capacity_pair);xlabel('n_p-op');ylabel('capacity C')
figure(2),plot(alpha,np_op);xlabel('\alpha');ylabel('optimal op')

