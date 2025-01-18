
    % -------------------------
    % 初始参数设置
    % -------------------------
    % a:射线电波波束的方位角(正北方向为0度，顺时针为正方向)
    % b:发射仰角
    % 发射点经纬度分别为longtitude\latitude
    a=0;
    
    %仅计算单个仰角启用
    %b=0.8*pi;
    
    longtitude=0;
    latitude=0;

    c = 3e8; % 光速 (单位：m/s)
    r_m = 6371 + 310; % 电离层最大电子密度高度，单位：km
    r_b = 6371+ 169; % 电离层F底部高度,仅考虑夜间情况（不存在D层），单位：km
    N_m = 1.1e12;       % F2层最大电子密度 (单位：电子/m^3)
    %fre = 2.5e6;         % 工作频率 (单位：Hz)

    t0 = 0;           % 初始群路径
    t_end = 10000;     % 最大群路径
    h = 0.001;            % 初始步长
    tol = 1e-6;       % 容许误差

    r0 = 6371.1; % 起始高度 (地球半径 + 0.5km)，单位：km
    
    
    theta0 = 0.5*pi-latitude;    % 起始纬度角，单位：弧度
    phi0 = longtitude;       % 起始经度角，单位：弧度
    k_r0 = 2*pi*fre/c*sin(b);       % 初始波矢径向分量
    k_theta0 = -2*pi*fre/c*cos(b)*cos(a);  % 初始波矢纬度分量
    k_phi0 = 2*pi*fre/c*cos(b)*sin(a);  % 初始波矢经度分量
    d_k_r0=0; % 径向方向的波矢量随群路径的变化率
    d_k_theta0=0; % theta方向的波矢量随群路径的变化率
    d_k_phi0=0; % phi方向的波矢量随群路径的变化率
    y0 = [r0; theta0; phi0;  k_r0; k_theta0; k_phi0]; % 初始状态向量
    dk_dp0 = [d_k_r0; d_k_theta0; d_k_phi0]; % 波矢量随群路径变化率的初始量
    

   % -------------------------
% 使用自适应Runge-Kutta方法计算
% -------------------------
% 定义不同的仰角范围（仅计算单个仰角时禁用）
elevations = [0.1, 0.2, 0.3, 0.4, 0.5] * pi;  % 例如，仰角从 0.1π 到 0.5π 弧度
colors = ['r', 'g', 'b', 'm', 'c'];  % 每个仰角对应的颜色

% 创建一个新的图形
figure;

% 循环计算不同仰角下的射线轨迹
hold on;
for i = 1:length(elevations)
    b = elevations(i);  % 当前仰角
    k_r0 = 2*pi*fre/c*sin(b);  % 更新初始波矢径向分量
    k_theta0 = -2*pi*fre/c*cos(b)*cos(a);  % 更新初始波矢纬度分量
    k_phi0 = 2*pi*fre/c*cos(b)*sin(a);  % 更新初始波矢经度分量
    y0 = [r0; theta0; phi0; k_r0; k_theta0; k_phi0];  % 更新初始状态向量

    % 计算当前仰角的射线轨迹
    [T, Y] = adaptive_runge_kutta(@(t, y) haselgrove_eq(t, y, r_m, r_b, N_m, fre), t0, y0, t_end, h, tol, r_m, r_b, N_m, fre);

    % 计算沿地球表面传播的距离
    r = Y(:, 1);      % r (球坐标半径)
    theta = Y(:, 2);  % theta (纬度角)
    phi = Y(:, 3);    % phi (经度角)

    % 计算沿地球表面传播的距离（考虑纬度和经度）
    d_surface = 6371e3 * sqrt((theta - theta(1)).^2 + (phi - phi(1)).^2);

    % 绘制当前仰角下的射线轨迹
    plot(d_surface, r - 6371, colors(i), 'LineWidth', 1);  % 绘制射线的高度与传播距离关系
end

% 图形设置
xlabel('传播距离 (km)');  % 横坐标为传播距离
ylabel('射线高度 (km)');   % 纵坐标为射线的高度
% 设置坐标轴范围
xlim([0, 1000]);  % 设置传播距离的范围为 0 到 1000 公里
ylim([0, 500]);   % 设置射线高度的范围为 0 到 500 公里
title('不同仰角下射线高度随传播距离的变化');
grid on;
legend(arrayfun(@(x) sprintf('\\theta = %.2f\\pi', x/pi), elevations, 'UniformOutput', false));  % 图例显示仰角

hold off;

function dydt = haselgrove_eq(t, y, r_m, r_b, N_m, f)
    % Haselgrove方程实现
    c = 3e8; % 光速 (单位：m/s)
    r = y(1)*1e3;
    r_b=r_b*1e3;
    r_m=r_m*1e3;
    theta = y(2); phi = y(3);
    k_r = y(4); k_theta = y(5); k_phi = y(6);

    % 计算折射率n及X梯度
    %n = calc_refraction_index(r, r_m, r_b, N_m, f);
    grad_x = calc_gradient_x(r, r_m, r_b, N_m, f);
    %disp("n:"); disp(n); % 折射率
    %disp("grad_n:"); disp(grad_n); % 梯度

    % 防止sin(theta) = 0
    if abs(sin(theta)) < 1e-6
        sin_theta = 1e-6; % 避免除零
    else
        sin_theta = sin(theta);
    end
    % 防止tan(theta) = 0
    if abs(tan(theta)) < 1e-6
        tan_theta = 1e-6; % 避免除零
    else
        tan_theta = tan(theta);
    end

    % 构建Haselgrove方程
    dydt = zeros(6,1);
    dydt(1) = c/(2*pi*f) * k_r; % dr/dP'
    dydt(2) = c/(r*2*pi*f) * k_theta; % dtheta/dP'
    dydt(3) = c/(r*2*pi*f*sin_theta) * k_phi; % dphi/dP'
    dydt(4) = -pi*f/c * grad_x(1) + c/(r*2*pi*f) * k_theta^2 + c/(r*2*pi*f) * k_phi^2; % dk_r/dP'
    dydt(5) = 1/r * (-pi * f/c * grad_x(2)- c/(2*pi*f) * k_r * k_theta + c/(2*pi*f*tan_theta)*k_phi^2); % dk_theta/dP'
    dydt(6) = 1/(r*sin_theta) * (-pi * f/c * grad_x(3) - c/(2*pi*f)*k_r*k_phi*sin_theta - c/(2*pi*f) * k_theta*k_phi*cos(theta)); % dk_phi/dP'
end

function Ne = calc_refraction_index(r, r_m, r_b, N_m, f)
    % 计算准抛物模型下的每立方米的电子浓度
    % r距离地表的高度，r_m电子浓度最大的高度（一般为F2层最大电子浓度的高度）
    % rb电离层底部的高度，一般为D层底部（白天），E层（夜间）底部高度
    % N_m最大电子浓度（一般为F2层的最大电子浓度）
    % 电波的工作频率，用于计算折射率n
     y_m = r_m - r_b; % 电离层半厚度
    if r < r_b || r > r_m*(r_b/(r_b-y_m))
        Ne=0;
    else
        Ne = N_m * (1 - ((r - r_m) / r_m)^2 * (r_b/ r)^2);% 每立方米电子浓度
        %f_p2 = 81 * N_e; % 临界频率的平方
        %n2 = 1 - f_p2 / f^2;
        %n = sqrt(max(n2, 0)); % 确保折射率非负
    end
end

function grad_X = calc_gradient_x(r, r_m, r_b, N_m, f)
    % 计算准抛物模型下x的梯度
    %r_b=r_b*1e-3;r_m=r_m*1e-3;r=r*1e-3;f=f*1e-6;
    % 电子的电荷量 e (单位：库伦)
    e = 1.602e-19; % 库伦

    % 电子的质量 m_e (单位：千克)
    m_e = 9.11e-31; % 千克

    % 真空电容率 ε₀ (单位：法拉每米)
    epsilon_0 = 8.854e-12; % 法拉每米
    y_m = r_m - r_b; % 电离层半厚度
    if r < r_b 
        grad_X = [0; 0; 0]; % 假设对流层中的折射率始终为1，梯度为0
    else
        y_m = r_m - r_b; % 电离层半厚度
      %   %电子密度的径向导数
      %   
        dN_e_dr = 2*N_m * (-(r - r_m) / y_m^2 * (r_b / r)^2 + ((r - r_m)^2 / y_m^2) * (r_b^2 / r^3));
      %   %X梯度的径向分量
      %  dX_dr = 80.6 * 1/f^2 *dNe_dr;
      %  grad_X = [dX_dr; 0; 0]; % 在QP准抛物线模型下仅需计算X径向分量
        % 计算电子浓度 N_e(r) 对 r 的导数
    %dN_e_dr = (2 * N_m * (r - r_m) / y_m^2) * (r_b^2 / r^3 - (r_b / r)^2);

   

    % 计算等离子体频率 f_n 和它的导数
    % f_n^2 = N_e / (epsilon_0 * m_e)
    f_n2 = (N_m * e^2/ (epsilon_0 * m_e)); % f_n^2 的常数部分
    
    % 计算 f_n 对 r 的导数
    df_n2_dr = dN_e_dr * e^2/ (epsilon_0 * m_e);  % f_n^2 的导数

    % 计算 X = (f_n / f)^2
    % 因为 f 是常数，可以直接计算 df_n^2 / dr
    
    %X = f_n2 / f^2;

    % X 对 r 的导数
    dX_dr = df_n2_dr / (f^2);
    grad_X = [dX_dr; 0; 0];
    end
        
end


function [T, Y] = adaptive_runge_kutta(f, t0, y0, t_end, h, tol,r_m,r_b,N_m,fre)
    T = t0; Y = y0';
    t = t0; y = y0;
    %dk_r=dk_p0(1);
    %dk_theta=dk_p0(2);
    %dk_phi=dk_p0(3);
    r=Y(1);theta=Y(2);phi=Y(3);
    k_r=Y(4); k_theta=Y(5); k_phi=Y(6);
    h_min=1e-3;
    h_max=10;
    % 电子的电荷量 e (单位：库伦)
    e = 1.602e-19; % 库伦

    % 电子的质量 m_e (单位：千克)
    m_e = 9.11e-31; % 千克

    % 真空电容率 ε₀ (单位：法拉每米)
    epsilon_0 = 8.854e-12; % 法拉每米
    while t < t_end
        if t + h > t_end
            h = t_end - t;
        end
        % 自适应步长h
        Ne = calc_refraction_index(r*1e3, r_m*1e3, r_b*1e3, N_m, fre);% 计算该高度下的电子浓度
        Wp_2=Ne*e^2/(epsilon_0*m_e); % 计算该高度下的离子谐振频率
        grad_x = calc_gradient_x(r*1e3, r_m*1e3, r_b*1e3, N_m, fre);% 计算该高度下的电子梯度
        he=tol/(Wp_2+grad_x(1)+grad_x(2)+grad_x(3));
     
        accept = true;
        while accept
            [y4, y5] = runge_kutta_step(f, t, y, h);
            error = norm(y5 - y4, 2);
            if error <= tol
                accept = false; % 如果误差满足容忍度
            else
                h45 = h * 0.9 * (tol / error)^(1/5); % 调整步长

                h = min(h_max,max(h_min,min(h45,he)));
            end
        end
      

        % 更新状态
        t = t + h;

        % 直接使用 y4 中的值更新状态
        r = y4(1);        % 更新半径
        theta = y4(2);    % 更新纬度角
        phi = y4(3);      % 更新经度角
        k_r = y4(4);      % 更新径向波矢分量
        k_theta = y4(5);  % 更新纬度波矢分量
        k_phi = y4(6);    % 更新经度波矢分量

        % 更新状态向量
        y = y4;           % 将新的状态向量赋值给 y

        T = [T; t];
        Y = [Y; y'];
        if y(1)<=6371
            % 仅计算单跳波束，再次到达地面后停止计算
            break
        end
        % 预测下一个步长
        h = max(h_min, min(h_max, h * 0.9 * (tol / error)^(1/5)));
    end
end



function [y4, y5] = runge_kutta_step(f, t, y, h)
    % 自适应Runge-Kutta Fehlberg方法
    k1 = h * f(t, y);
    k2 = h * f(t + h/4, y + k1/4);
    %disp("k2:"); disp(k2);
    k3 = h * f(t + 3*h/8, y + 3*k1/32 + 9*k2/32);
    %disp("k3:"); disp(k3);
    k4 = h * f(t + 12*h/13, y + 1932*k1/2197 - 7200*k2/2197 + 7296*k3/2197);
    k5 = h * f(t + h, y + 439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104);
    k6 = h * f(t + h/2, y - 8*k1/27 + 2*k2 - 3544*k3/2565 + 1859*k4/4104 - 11*k5/40);

    y4 = y + (25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5);
    y5 = y + (16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55);
end
