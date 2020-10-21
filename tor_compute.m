function [tau_theo,load_torque_est,motor_torque_est] = tor_compute(q,dq,ddq)
%#codegen
% 输入参数: 关节角度q: 6xn; 关节角速度q: 6xn; 关节角加速度q: 6xn;
% 输出参数: 理论计算力矩tau_theo: 6xn; 辨识参数估计负载力矩load_torque_est: 6xn; 辨识参数估计电机力矩motor_torque_est: 6xn;
n = size(q,2); % n为轨迹点数
m = size(q,1); % m为自由度个数
%% 名义动力学模型 理论力矩估计
%名义动力学参数设置
g = 9.81; %重力加速度g
%改进DH参数
alpha=[0;pi/2;0;0;pi/2;-pi/2];    
a=[0;0;-0.264;-0.237;0;0];         
d=[0.144;0;0;0.1065;0.114;0.09];           
theta=[0;-pi/2;0;-pi/2;0;0];
dh_list = [alpha a d theta];

%连杆质量(kg)
mass_list = [2.9196; 6.6409; 2.8613; 1.7066; 1.7066; 0.1759]; 

%连杆质心相对于坐标系{i}坐标(m)
mass_center_list = [0, -0.0031, -0.0139;
                   -0.132, 0, 0.1078;
                   -0.1643, -0.0001, 0.0095;
                   -0.0001, -0.0210, -0.0025;
                   0.0001, 0.0210, -0.0025;
                   -0.0001, -0.0004, -0.0137];
                        
%inertia_tensor_list 连杆相对于质心坐标系的惯量张量(kg*m^2)
inertia_tensor_list = zeros(3,3,6);
inertia_tensor_list(:,:,1)  = [43,0,0; 0,41,-1; 0,-1,32]*10^-4;
inertia_tensor_list(:,:,2)  = [98,0,0; 0,1084,0; 0,0,1061]*10^-4;
inertia_tensor_list(:,:,3)  = [38,0,-4; 0,291,0; -4,0,285]*10^-4;
inertia_tensor_list(:,:,4)  = [21,0,0; 0,17,-1; 0,-1,20]*10^-4;
inertia_tensor_list(:,:,5)  = [21,0,0; 0,17,1; 0,1,20]*10^-4;
inertia_tensor_list(:,:,6)  = [1,0,0; 0,1,0; 0,0,2]*10^-4; 

%f_tip: 机械臂末端施加外力和力矩
f_tip = [0 0 0; 0 0 0];

% 理论力矩计算
tau_theo = zeros(m,n);
for i = 1:n
    taulist_temp = Newtown_InverseDynamics(q(:,i), dq(:,i), ddq(:,i), g,...
                                   dh_list, mass_list, mass_center_list, inertia_tensor_list, f_tip);
    tau_theo(1:m,i)=taulist_temp;                      
%     tau_theo = [tau_theo,taulist_temp];
end
% 最终计算力矩与实际力矩正负相反
tau_theo = -tau_theo;

%% 辨识模型估计关节力矩
%摩擦力项数
nf = 5;  % 正反(库伦+粘滞)摩擦 + 电机惯量力矩补偿
% 导入辨识所得最小参数数据
loadpath1 = 'para_20s_1_0kg_5f.mat'; 
data1 = load(loadpath1);

% 计算线性回归矩阵
% J = []; K = []; 
% Kf = [];
J=zeros(6*n,6);K=zeros(6*n,60+nf*6);Kf=zeros(6*n,nf*6);
for i = 1:n
    [J_temp,K_temp,Kf_temp] = Compute_Dynmatrix(q(:,i), dq(:,i), ddq(:,i), g, dh_list, nf);
    J(6*(i-1)+1:6*i,:)=J_temp;
    K(6*(i-1)+1:6*i,1:60)=K_temp;
    Kf(6*(i-1)+1:6*i,:)=Kf_temp;
%     J = [J; J_temp]; % 与外力有关
%     K = [K;K_temp]; % 与机器人惯性参数有关
%     Kf = [Kf; Kf_temp]; % 与摩擦力有关
end
% 考虑摩擦力
K(:,61:60+nf*6)=Kf;
% K = [K,Kf];

%% 线性回归矩阵消除全为0的列, 未注释代码等同于注释代码，因为不知道find函数代码生成
% C_k=Self_Find(K);
% K_d=zeros(6*n,length(C_k));
% for i=1:length(C_k)
%     K_d(:,i)=K(:,C_k(i));
% end
% K_d = K(:,find(sum(abs(K))'>sqrt(eps)));
% disp(isequal(K_d,K_d1));      %对比结果
K_d = [K(:,10),K(:,12:13),K(:,15:end)];
%%
% QR分解，得到与最小参数集对应的线性回归矩阵
%[~,R,P] = qr(K_d);

% 使用与辨识最小参数集的同一个列转置矩阵P
P = data1.P;
%
n_d = rank(K_d);

KK = K_d*P;
% KK_m = KK(:,1:n_d);
KK_m = KK(:,1:size(data1.para_all,1));
% 代入辨识的最小参数集
para_all = data1.para_all;

% 估计关节力矩
load_torque_est = KK_m * para_all;
load_torque_est = reshape(load_torque_est,[size(load_torque_est,1)/n,n]);

%% 辨识模型估计电机力矩
%摩擦力项数
nf = 4;  % 正反(库伦+粘滞)摩擦
% 导入辨识所得最小参数数据
loadpath1 = 'para_20s_0kg_1_4f.mat'; 
data1 = load(loadpath1);

% 计算线性回归矩阵
% J = []; K = []; Kf = [];
J=zeros(6*n,6);K=zeros(6*n,60+nf*6);Kf=zeros(6*n,nf*6);
for i = 1:n
    [J_temp,K_temp,Kf_temp] = Compute_Dynmatrix(q(:,i), dq(:,i), ddq(:,i), g, dh_list, nf);
    J(6*(i-1)+1:6*i,:)=J_temp;
    K(6*(i-1)+1:6*i,1:60)=K_temp;
    Kf(6*(i-1)+1:6*i,:)=Kf_temp;
%     J = [J; J_temp]; % 与外力有关
%     K = [K;K_temp]; % 与机器人惯性参数有关
%     Kf = [Kf; Kf_temp]; % 与摩擦力有关
end
% 考虑摩擦力
%    K = [K,Kf];
K(:,61:60+nf*6)=Kf;
%% 线性回归矩阵消除全为0的列, 未注释代码等同于注释代码，因为不知道find函数代码生成C_k=Self_Find(K);
% K_d=zeros(6*n,length(C_k));
% for i=1:length(C_k)
%     K_d(:,i)=K(:,C_k(i));
% end
% K_d = K(:,find(sum(abs(K))'>sqrt(eps)));
% disp(isequal(K_d,K_d1));      %对比结果
K_d = [K(:,10),K(:,12:13),K(:,15:end)];
%%
% QR分解，得到与最小参数集对应的线性回归矩阵
%[~,R,P] = qr(K_d);

% 使用与辨识最小参数集的同一个列转置矩阵P
P = data1.P;
%
n_d = rank(K_d);

KK = K_d*P;
% KK_m = KK(:,1:n_d);
KK_m = KK(:,1:size(data1.para_all,1));
% 代入辨识的最小参数集
para_all = data1.para_all;

% 估计电机力矩
motor_torque_est = KK_m * para_all;
motor_torque_est = reshape(motor_torque_est,[size(motor_torque_est,1)/n,n]);
end

function C_k=Self_Find(K)
    SumK=sum(abs(K))';
    num=1;
    for i=1:length(SumK)
        if(SumK(i)>sqrt(eps))
            num=num+1;    %求满足SumK(i)>sqrt(eps)的数据个数
        end
    end
    C_k=zeros(1,num-1);  
    num=1;
    for i=1:length(SumK)
        if(SumK(i)>sqrt(eps))
            C_k(num)=i;  %将满足SumK(i)>sqrt(eps)的SumK下标记录
            num=num+1;
        end
    end
%     num=1;
end