%用最最原始的二维搜索和近似后，tau_star用lambertw表示，
% Pmd用gamma表示的方法，求全局最优，看看精确度

% 二维搜索最优的np和alpha
%Pfa和Pmd都用PDF累积获得，tau_star用一层for循环得到表示
clear all
vw=5;vb=1;%注意，这些都是标准差
n=100;
rho=1;lamda_ab=1;
epsilon=0.1;
alpha=linspace(0.01,0.999,20);
 np=linspace(1,99,20);
% np=linspace(3,99,33);
tau=linspace(0.1,3,20);%tau要大于rhod+vw^2，才能保证Pmd大于零
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
    %对于每个alpha，给出一串np和一串kexi，找出xi>=1-epsilon的那些np，在这些位置中找出使得
    %信道容量最大的np

%数值得出tau_star,对于每个np，代入一串tau得到一串相应的xi，
%找到min xi的tau_star,得出相应的xi_star
for j=1:length(np)
    j
    rhop1=(1-alpha(i)).*rho.*n./np(j);
    rhod1=alpha(i).*rho.*n./nd(j);
    nd(j)=n-np(j);
    for m=1:length(tau)
    
    a=n*tau(m).*(rhop1-rhod1)/((rhod1+vw^2).*(rhop1+vw^2));%合流几何函数里的参数
    f_md(m)=(n^n*tau(m).^(n-1)*exp(-n*tau(m)./(rhod1+vw^2))./((rhop1+vw^2).^np(j)))./((rhod1+vw^2).^nd(j)).*...
        hypergeom(np(j),n,a)/gamma(n);%T的PDF
    f_fa(m)=(n/vw^2).^n./gamma(n).*tau(m).^(n-1).*exp(-n.*tau(m)./vw^2);
    end
    %用PDF获得的Pmd 和Pfa
P_md=cumsum( f_md)/sum( f_md);  %卡方分布的Pmd
P_fa=1-cumsum( f_fa)/sum( f_fa);  %卡方分布的Pfa
xi=P_md+P_fa;

    loca=(xi==min(xi));
    %找到使得loca第一次出现1的位置
    location=find(loca,1);
    tau_star=tau(location);
%     xi_star_test(j)=xi(location);
    xi_star(i,j)=xi(location);
    
end
%找到满足条件的xi_star的位置
plot(np,xi_star(i,:));
pos=(xi_star(i,:)>=1-epsilon);
xi_star_ok=xi_star(i,pos);
% pos=(xi_star(i,:)==xi_star_ok);
%若没有满足条件的，则设置np和capacity为0
if any(pos)==false
    np_op=0;
    capacity_pair(i)=0;
    continue
end
%把可行的np的位置找出来
np_ok=np(pos);

% 计算所有满足的rho_eff和容量
rho_eff=(lamda_ab^2.*rho.*n.*alpha(i).*(1-alpha(i)))./...
    ((n-np_ok).*vb^2.*(lamda_ab.*(1-alpha(i))+vb^2./rho./n+lamda_ab.*alpha(i)./(n-np_ok)));
capacity=(n-np_ok)./n.*log2(1+rho_eff); 
%找到使得容量最大的np的位置
pos2=(capacity==max(capacity));
np_op(i)=np_ok(pos2);

%把每个alpha对应的最优np当成一对，求出相应的capacity
capacity_pair(i)=max(capacity); %信道容量
end
pos3=(capacity_pair==max(capacity_pair));
pos3=find(pos3,1);
alpha_op_whole=alpha(pos3);
np_op_whole=np_op(pos3);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 二维搜索最优的np和alpha
%tau_star用lambertw函数表示，Pmd用igamma表示

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
    %对于每个alpha，给出一串np和一串kexi，找出xi>=1-epsilon的那些np，在这些位置中找出使得
    %信道容量最大的np
       sigma_T_squ(i,:)=(1/n^2).*(np.*((1-alpha(i)).*rho.*n./np+vw^2).^2+...
    nd.*(alpha(i).*rho.*n./nd+vw^2).^2);
%给出tau_star的表达式
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
% %直接计算xi_star=1-epsilon的np，看看全局搜索到的是不是该交点
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
%把可行的np的位置找出来
np_ok=np(pos);

% 计算所有的rho_eff和容量
rho_eff=(lamda_ab^2.*rho.*n.*alpha(i).*(1-alpha(i)))./...
    ((n-np).*vb^2.*(lamda_ab.*(1-alpha(i))+vb^2./rho./n+lamda_ab.*alpha(i)./(n-np)));
capacity=(n-np)./n.*log2(1+rho_eff); 
%找出满足条件的rho_eff和容量
rho_eff_ok=rho_eff(pos);
capacity_ok=capacity(pos);
%找到使得容量最大的np的位置
pos2=(capacity==max(capacity_ok));
np_op_orig(i)=np(pos2);

%把每个alpha对应的最优np当成一对，求出相应的capacity
rho_eff_pair(i)=(lamda_ab^2.*rho.*n.*alpha(i).*(1-alpha(i)))./...
    ((n-np_op_orig(i)).*vb^2.*(lamda_ab.*(1-alpha(i))+vb^2./rho./n+lamda_ab.*alpha(i)./(n-np_op_orig(i))));
capacity_pair_orig(i)=(n-np_op_orig(i))./n.*log2(1+rho_eff_pair(i)); %信道容量
end
pos3=(capacity_pair_orig==max(capacity_pair_orig));
alpha_op_whole_orig=alpha(pos3);
np_op_whole_orig=np_op_orig(pos3);





plot(np_op,capacity_pair);xlabel('n_p-op');ylabel('capacity C')
figure(2),plot(alpha,np_op);xlabel('\alpha');ylabel('optimal op')

