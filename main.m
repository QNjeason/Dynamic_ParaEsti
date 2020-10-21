clc;clear;close all
%% 导入轨迹数据
 loadpath = 'mod_data_15s_0kg_2.mat';
data = load(loadpath);
%关节位置、角速度和角加速度
q = data.q_f';
dq = data.dq_calculate1';
ddq = data.ddq_calculate1';

%关节负载力矩(力矩传感器测得)
load_torque = data.force';
load_torque_f = data.force_f';

%关节电机力矩(电流计算得到)
motor_torque = data.motor_torque';
motor_torque_f = data.motor_torque_f';

% 计算轨迹运行时间t，降采样率为5
down_rate = 5;
dt = 0.002*down_rate;
t = 0:dt:dt*(size(q,2)-1);
n = size(q,2);

%% 计算理论力矩,辨识模型估计的负载力矩和电机力矩
%理论计算力矩tau_theo; 辨识参数估计负载力矩load_torque_est; 辨识参数估计电机力矩motor_torque_est;
[tau_theo,load_torque_est,motor_torque_est] = tor_compute(q,dq,ddq);
%% load_torque plot
%close all
 figure
for cnt = 1:6
subplot(3,2,cnt)
plot(t,load_torque(cnt,:),'g'),hold on
plot(t,tau_theo(cnt,:),'b','linewidth',1),hold on
%plot(t,efforts_f(cnt,:),'k','linewidth',1),hold on
plot(t,load_torque_est(cnt,:),'r','linewidth',1),hold on
xlabel('时间t(s)'), ylabel('负载力矩(N*m)')
legend('实际力矩','理论力矩','辨识力矩','Location','southeast')
title(['关节',num2str(cnt)])
end
%% motor_torque plot
%close all
figure
for cnt = 1:6
subplot(3,2,cnt)
plot(t,motor_torque(cnt,:),'g'),hold on
plot(t,motor_torque_est(cnt,:),'r','linewidth',1),hold on
xlabel('时间t(s)'), ylabel('电机力矩(N*m)')
legend('实际力矩','辨识力矩','Location','southeast')
title(['关节',num2str(cnt)])
end